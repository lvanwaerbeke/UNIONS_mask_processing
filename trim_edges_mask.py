#!/usr/bin/env python3
"""Build a trim-edge FITS mask from tile sky boundaries and image WCS.

For each tile, pixels outside the RA/Dec box listed in ``FullSky_tiles_cuts.txt``
receive the value 512 (2^9). Pixels inside the box receive 0.
The output FITS keeps the original header/WCS unchanged.
"""

from __future__ import annotations

import argparse
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


def process_tile(
    tile: str, args: argparse.Namespace, tile_cuts: dict[str, tuple[float, float, float, float]]
) -> None:
    if tile not in tile_cuts:
        raise KeyError(f"Tile {tile} not found in {args.tile_cuts}")

    ra_min, ra_max, dec_min, dec_max = tile_cuts[tile]
    source_fits = choose_source_fits(tile, args.source_templates)
    output_path = Path(args.output_template.format(tile=tile))
    if output_path.suffix != ".gz":
        output_path = Path(f"{output_path}.gz")

    with fits.open(source_fits, memmap=False) as hdul:
        header = hdul[0].header.copy()

    mask = build_trim_mask(header, ra_min, ra_max, dec_min, dec_max, args.row_chunk_size)
    write_mask_fits(output_path, header, mask)

    print(f"Tile {tile}")
    print(f"  Tile cuts: RA=[{ra_min}, {ra_max}] DEC=[{dec_min}, {dec_max}]")
    print(f"  Source FITS: {source_fits}")
    print(f"  Output FITS: {output_path}")
    print(f"  Non-zero pixels: {int(np.count_nonzero(mask))}")
    print(f"  Unique values: {sorted(int(v) for v in np.unique(mask))}")


def main() -> int:
    args = parse_args()
    tile_cuts = parse_tile_cuts(Path(args.tile_cuts))
    skipped_tiles: list[str] = []
    processed_count = 0

    for tile in collect_tiles(args):
        try:
            process_tile(tile, args, tile_cuts)
        except KeyError as exc:
            print(f"WARNING: {exc}")
            skipped_tiles.append(tile)
            continue
        processed_count += 1

    skipped_tiles_log = Path(args.missing_tiles_log)
    skipped_tiles_log.parent.mkdir(parents=True, exist_ok=True)
    skipped_tiles_log.write_text(
        "".join(f"{tile}\n" for tile in skipped_tiles),
        encoding="utf-8",
    )

    if skipped_tiles:
        print(
            f"Finished with warnings: processed {processed_count} tile(s), "
            f"skipped {len(skipped_tiles)} tile(s)."
        )
        print(f"Skipped tile log written to {skipped_tiles_log}")
    else:
        print(f"Finished: processed {processed_count} tile(s), no skipped tiles.")
        print(f"Skipped tile log written to {skipped_tiles_log} (empty)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
