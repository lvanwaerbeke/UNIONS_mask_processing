#!/usr/bin/env python3
"""Combine UNIONS mask layers into a final summed FITS mask."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile", help="Single tile ID, e.g. 122.316.")
    parser.add_argument("--tile-list-file", help="Text file containing one tile ID per line.")
    parser.add_argument(
        "--maskstars-dir",
        default=".",
        help="Directory containing UNIONS.{tile}_r_maskstars.mask.fits",
    )
    parser.add_argument(
        "--maskstars-template",
        help="Optional full path template for UNIONS.{tile}_r_maskstars.mask.fits(.gz).",
    )
    parser.add_argument(
        "--extpoly-dir",
        default=".",
        help="Directory containing UNIONS.{tile}.ext_polygon.mask.fits",
    )
    parser.add_argument(
        "--extpoly-template",
        help="Optional full path template for UNIONS.{tile}.ext_polygon.mask.fits(.gz).",
    )
    parser.add_argument(
        "--ugriz-dir",
        default=".",
        help="Directory containing UNIONS.{tile}_ugriz.mask.fits",
    )
    parser.add_argument(
        "--ugriz-template",
        help="Optional full path template for UNIONS.{tile}_ugriz.mask.fits(.gz).",
    )
    parser.add_argument(
        "--trim-dir",
        default=".",
        help="Directory containing UNIONS.{tile}_r_masktrim.fits",
    )
    parser.add_argument(
        "--trim-template",
        help="Optional full path template for UNIONS.{tile}_r_masktrim.fits(.gz).",
    )
    parser.add_argument(
        "--maximask-dir",
        default=".",
        help="Directory containing UNIONS.{tile}.r.maximask.fits",
    )
    parser.add_argument(
        "--maximask-template",
        help="Optional full path template for UNIONS.{tile}.r.maximask.fits(.gz).",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory where UNIONS.{tile}_final.mask.fits.gz will be written.",
    )
    parser.add_argument(
        "--problem-log",
        default="combine_masks_problems.log",
        help="Log file where missing, unreadable, or inconsistent tile inputs will be recorded.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes to use. Use 1 for serial combining.",
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


def input_paths_for_tile(tile: str, args: argparse.Namespace) -> list[tuple[str, Path]]:
    maskstars = (
        Path(args.maskstars_template.format(tile=tile))
        if args.maskstars_template
        else Path(args.maskstars_dir) / f"UNIONS.{tile}_r_maskstars.mask.fits"
    )
    extpoly = (
        Path(args.extpoly_template.format(tile=tile))
        if args.extpoly_template
        else Path(args.extpoly_dir) / f"UNIONS.{tile}.ext_polygon.mask.fits"
    )
    ugriz = (
        Path(args.ugriz_template.format(tile=tile))
        if args.ugriz_template
        else Path(args.ugriz_dir) / f"UNIONS.{tile}_ugriz.mask.fits"
    )
    trim = (
        Path(args.trim_template.format(tile=tile))
        if args.trim_template
        else Path(args.trim_dir) / f"UNIONS.{tile}_r_masktrim.fits"
    )
    maximask = (
        Path(args.maximask_template.format(tile=tile))
        if args.maximask_template
        else Path(args.maximask_dir) / f"UNIONS.{tile}.r.maximask.fits"
    )
    return [
        ("maskstars", maskstars),
        ("ext_polygon", extpoly),
        ("ugriz", ugriz),
        ("trim", trim),
        ("maximask", maximask),
    ]


def resolve_input_path(path: Path) -> Path | None:
    if path.exists():
        return path
    gz_path = Path(f"{path}.gz")
    if gz_path.exists():
        return gz_path
    return None


def append_problem(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")
        handle.flush()


def combine_one(tile: str, args: argparse.Namespace) -> dict[str, object]:
    requested_inputs = input_paths_for_tile(tile, args)
    resolved_inputs: list[tuple[str, Path]] = []
    missing: list[str] = []
    for label, path in requested_inputs:
        resolved = resolve_input_path(path)
        if resolved is None:
            missing.append(f"{label}: {path}(.gz)")
        else:
            resolved_inputs.append((label, resolved))
    if missing:
        raise FileNotFoundError(f"Missing input mask(s) for tile {tile}: {', '.join(missing)}")

    combined: np.ndarray | None = None
    output_header: fits.Header | None = None
    shape: tuple[int, int] | None = None

    for label, path in resolved_inputs:
        try:
            with fits.open(path, memmap=False) as hdul:
                header = hdul[0].header.copy()
                data = np.asarray(hdul[0].data)
        except Exception as exc:
            raise OSError(f"Unreadable FITS for tile {tile}, {label}: {path} :: {type(exc).__name__}: {exc}") from exc

        if shape is None:
            shape = data.shape
            combined = data.astype(np.uint16, copy=True)
            output_header = header
        else:
            if data.shape != shape:
                raise ValueError(
                    f"Shape mismatch for tile {tile}: {path} has {data.shape}, expected {shape}"
                )
            combined = combined + data.astype(np.uint16, copy=False)

    assert combined is not None
    assert output_header is not None

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"UNIONS.{tile}_final.mask.fits.gz"
    fits.PrimaryHDU(data=combined.astype(np.uint16, copy=False), header=output_header).writeto(
        output_path, overwrite=True
    )

    return {
        "tile": tile,
        "inputs": [(label, str(path)) for label, path in resolved_inputs],
        "output": str(output_path),
        "unique_values": sorted(int(v) for v in np.unique(combined)),
        "non_zero_pixels": int(np.count_nonzero(combined)),
    }


def print_success(result: dict[str, object]) -> None:
    print(f"Tile {result['tile']}")
    for label, path in result["inputs"]:
        print(f"  {label}: {path}")
    print(f"  Output: {result['output']}")
    print(f"  Unique values: {result['unique_values']}")
    print(f"  Non-zero pixels: {result['non_zero_pixels']}")


def main() -> int:
    args = parse_args()
    problem_log = Path(args.problem_log)
    problem_log.parent.mkdir(parents=True, exist_ok=True)
    problem_log.write_text("", encoding="utf-8")

    processed_count = 0
    skipped_count = 0
    tile_ids = collect_tiles(args)

    if args.jobs <= 1:
        for tile in tile_ids:
            try:
                result = combine_one(tile, args)
                print_success(result)
                processed_count += 1
            except Exception as exc:
                skipped_count += 1
                message = f"{tile} :: {type(exc).__name__}: {exc}"
                print(f"WARNING: {message}")
                append_problem(problem_log, message)
                continue
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_tile = {executor.submit(combine_one, tile, args): tile for tile in tile_ids}
            for future in as_completed(future_to_tile):
                tile = future_to_tile[future]
                try:
                    result = future.result()
                    print_success(result)
                    processed_count += 1
                except Exception as exc:
                    skipped_count += 1
                    message = f"{tile} :: {type(exc).__name__}: {exc}"
                    print(f"WARNING: {message}")
                    append_problem(problem_log, message)
                    continue

    if skipped_count:
        print(
            f"Finished with warnings: combined {processed_count} tile(s), "
            f"skipped {skipped_count} tile(s)."
        )
        print(f"Problem log: {problem_log}")
    else:
        print(f"Finished: combined {processed_count} tile(s), no skipped tiles.")
        print(f"Problem log: {problem_log} (empty)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
