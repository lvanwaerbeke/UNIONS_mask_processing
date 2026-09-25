#!/usr/bin/env python3
"""Build a trim-edge FITS mask from tile sky boundaries and image WCS.

For each tile, pixels outside the RA/Dec box listed in ``FullSky_tiles_cuts.txt``
receive the value 512 (2^9). Pixels inside the box receive 0.
The output FITS keeps the original header/WCS unchanged.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


TRIM_MASK_VALUE = np.int16(512)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile", help="Single tile ID, e.g. 122.316.")
    parser.add_argument("--tile-list-file", help="Text file containing one tile ID per line.")
    parser.add_argument(
        "--tile-cuts",
        default="FullSky_tiles_cuts.txt",
        help="File containing tile ID, RAmin, RAmax, DECmin, DECmax.",
    )
    parser.add_argument(
        "--source-templates",
        nargs="+",
        default=[
            "UNIONS.{tile}_ugriz.mask.fits",
            "UNIONS.{tile}_r.fits",
        ],
        help="Source FITS templates checked in order for WCS and image size.",
    )
    parser.add_argument(
        "--output-template",
        default="UNIONS.{tile}_r_masktrim.fits.gz",
        help="Output trim-mask FITS filename template. The script always writes gzip-compressed FITS.",
    )
    parser.add_argument(
        "--row-chunk-size",
        type=int,
        default=256,
        help="Number of image rows to transform per WCS chunk.",
    )
    parser.add_argument(
        "--missing-tiles-log",
        default="skipped_trim_tiles.log",
        help="Path to a text file where skipped tile IDs will be written.",
    )
    parser.add_argument(
        "--problem-log",
        default="trim_edges_mask_problems.log",
        help="Log file where processing failures other than missing tile cuts will be recorded.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of tile worker processes. Use 16 for a 16-CPU machine.",
    )
    return parser.parse_args()


def read_tile_list_file(path: Path) -> list[str]:
    tiles: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tiles.append(stripped.split()[0])
    return tiles


def collect_tiles(args: argparse.Namespace) -> list[str]:
    if args.tile_list_file:
        return read_tile_list_file(Path(args.tile_list_file))
    if args.tile:
        return [args.tile]
    raise SystemExit("Provide either --tile or --tile-list-file.")


def parse_tile_cuts(path: Path) -> dict[str, tuple[float, float, float, float]]:
    tile_cuts: dict[str, tuple[float, float, float, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 5:
            continue
        tile_cuts[parts[0]] = (
            float(parts[1]),
            float(parts[2]),
            float(parts[3]),
            float(parts[4]),
        )
    return tile_cuts


def choose_source_fits(tile: str, templates: list[str]) -> Path:
    for template in templates:
        candidate = Path(template.format(tile=tile))
        if candidate.exists():
            return candidate
    attempted = ", ".join(template.format(tile=tile) for template in templates)
    raise FileNotFoundError(f"No source FITS found for tile {tile}. Tried: {attempted}")


def ra_in_interval(ra: np.ndarray, ra_min: float, ra_max: float) -> np.ndarray:
    if ra_min <= ra_max:
        return (ra >= ra_min) & (ra <= ra_max)
    return (ra >= ra_min) | (ra <= ra_max)


def build_trim_mask(
    header: fits.Header,
    ra_min: float,
    ra_max: float,
    dec_min: float,
    dec_max: float,
    row_chunk_size: int,
) -> np.ndarray:
    width = int(header["NAXIS1"])
    height = int(header["NAXIS2"])
    wcs = WCS(header)
    mask = np.full((height, width), TRIM_MASK_VALUE, dtype=np.int16)
    x_pixels = np.arange(width, dtype=np.float64)

    for y0 in range(0, height, row_chunk_size):
        y1 = min(y0 + row_chunk_size, height)
        y_pixels = np.arange(y0, y1, dtype=np.float64)
        xx, yy = np.meshgrid(x_pixels, y_pixels)
        ra, dec = wcs.all_pix2world(xx, yy, 0)
        inside = ra_in_interval(ra, ra_min, ra_max) & (dec >= dec_min) & (dec <= dec_max)
        mask[y0:y1, :] = np.where(inside, 0, TRIM_MASK_VALUE).astype(np.int16, copy=False)

    return mask


def write_mask_fits(output_path: Path, header: fits.Header, mask: np.ndarray) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=mask.astype(np.int16, copy=False), header=header).writeto(
        output_path, overwrite=True
    )


def append_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")
        handle.flush()


def process_tile(
    tile: str,
    args: argparse.Namespace,
    cuts: tuple[float, float, float, float],
) -> dict[str, object]:
    ra_min, ra_max, dec_min, dec_max = cuts
    source_fits = choose_source_fits(tile, args.source_templates)
    output_path = Path(args.output_template.format(tile=tile))
    if output_path.suffix != ".gz":
        output_path = Path(f"{output_path}.gz")

    with fits.open(source_fits, memmap=False) as hdul:
        header = hdul[0].header.copy()

    mask = build_trim_mask(header, ra_min, ra_max, dec_min, dec_max, args.row_chunk_size)
    write_mask_fits(output_path, header, mask)

    return {
        "tile": tile,
        "ra_min": ra_min,
        "ra_max": ra_max,
        "dec_min": dec_min,
        "dec_max": dec_max,
        "source_fits": str(source_fits),
        "output_path": str(output_path),
        "non_zero_pixels": int(np.count_nonzero(mask)),
        "unique_values": sorted(int(v) for v in np.unique(mask)),
    }


def print_result(result: dict[str, object]) -> None:
    print(f"Tile {result['tile']}")
    print(
        f"  Tile cuts: RA=[{result['ra_min']}, {result['ra_max']}] "
        f"DEC=[{result['dec_min']}, {result['dec_max']}]"
    )
    print(f"  Source FITS: {result['source_fits']}")
    print(f"  Output FITS: {result['output_path']}")
    print(f"  Non-zero pixels: {result['non_zero_pixels']}")
    print(f"  Unique values: {result['unique_values']}")


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1.")
    if args.row_chunk_size < 1:
        raise SystemExit("--row-chunk-size must be at least 1.")

    tile_cuts = parse_tile_cuts(Path(args.tile_cuts))
    tiles = collect_tiles(args)
    skipped_tiles_log = Path(args.missing_tiles_log)
    problem_log = Path(args.problem_log)
    skipped_tiles_log.parent.mkdir(parents=True, exist_ok=True)
    problem_log.parent.mkdir(parents=True, exist_ok=True)
    skipped_tiles_log.write_text("", encoding="utf-8")
    problem_log.write_text("", encoding="utf-8")

    available_tiles: list[str] = []
    skipped_count = 0
    for tile in tiles:
        if tile not in tile_cuts:
            message = f"Tile {tile} not found in {args.tile_cuts}"
            print(f"WARNING: {message}")
            append_log(skipped_tiles_log, tile)
            skipped_count += 1
        else:
            available_tiles.append(tile)

    processed_count = 0
    failed_count = 0

    if args.jobs == 1:
        for tile in available_tiles:
            try:
                result = process_tile(tile, args, tile_cuts[tile])
                print_result(result)
                processed_count += 1
            except Exception as exc:
                failed_count += 1
                message = f"{tile} :: {type(exc).__name__}: {exc}"
                print(f"WARNING: {message}")
                append_log(problem_log, message)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_tile = {
                executor.submit(process_tile, tile, args, tile_cuts[tile]): tile
                for tile in available_tiles
            }
            for future in as_completed(future_to_tile):
                tile = future_to_tile[future]
                try:
                    result = future.result()
                    print_result(result)
                    processed_count += 1
                except Exception as exc:
                    failed_count += 1
                    message = f"{tile} :: {type(exc).__name__}: {exc}"
                    print(f"WARNING: {message}")
                    append_log(problem_log, message)

    if skipped_count or failed_count:
        print(
            f"Finished with warnings: processed {processed_count} tile(s), "
            f"skipped {skipped_count} tile(s) with missing cuts, "
            f"failed {failed_count} tile(s)."
        )
        print(f"Missing tile-cut log: {skipped_tiles_log}")
        print(f"Processing problem log: {problem_log}")
    else:
        print(f"Finished: processed {processed_count} tile(s), no skipped tiles.")
        print(f"Missing tile-cut log: {skipped_tiles_log} (empty)")
        print(f"Processing problem log: {problem_log} (empty)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
