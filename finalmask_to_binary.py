#!/usr/bin/env python3
"""Extract a binary 0/1 mask from a final 16-bit UNIONS mask FITS image."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits


BIT_INFO: dict[int, str] = {
    0: "faint star halo",
    1: "bright star halo",
    2: "star mask",
    3: "manual mask / external polygon",
    4: "u footprint",
    5: "g footprint",
    6: "r footprint",
    7: "i footprint",
    8: "z footprint",
    9: "tile RA/DEC trim",
    10: "maximask",
    11: "PanSTARRS z2 footprint",
    12: "free bit",
    13: "free bit",
    14: "free bit",
    15: "unused",
}


MASK_ALIASES: dict[str, list[int]] = {
    "faint_star_halo": [0],
    "bright_star_halo": [1],
    "star_mask": [2],
    "manual_mask": [3],
    "u": [4],
    "g": [5],
    "r": [6],
    "i": [7],
    "z": [8],
    "trim": [9],
    "maximask": [10],
    "z2": [11],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--list-masks",
        action="store_true",
        help="Print all accepted --masks keywords and exit.",
    )
    parser.add_argument("--input", help="Explicit input final mask FITS path for single-file mode.")
    parser.add_argument("--output", help="Explicit output binary FITS path for single-file mode.")
    parser.add_argument("--tile", help="Single tile ID, e.g. 122.316.")
    parser.add_argument("--tile-list-file", help="Text file containing one tile ID per line.")
    parser.add_argument(
        "--input-template",
        default="UNIONS.{tile}_final.mask.fits.gz",
        help="Input template used with --tile or --tile-list-file.",
    )
    parser.add_argument(
        "--output-template",
        default="UNIONS.{tile}_selected.mask.fits.gz",
        help="Output template used with --tile or --tile-list-file.",
    )
    parser.add_argument(
        "--select",
        nargs="*",
        default=[],
        help=(
            "Mask keywords or bit numbers to activate in the output binary mask. "
            "Unselected bits are ignored. Examples: --select trim manual_mask or "
            "--select 3 9 10"
        ),
    )
    parser.add_argument(
        "--bits",
        nargs="*",
        default=[],
        help="Legacy mode: bit numbers to select, e.g. --bits 3 9 10 or --bits 3,9,10",
    )
    parser.add_argument(
        "--masks",
        nargs="*",
        default=[],
        help=(
            "Legacy mode: named mask aliases to select, e.g. --masks manual_mask trim "
            "maximask or --masks ugriz"
        ),
    )
    parser.add_argument(
        "--require-all",
        action="store_true",
        help=(
            "Legacy mode: set output pixel to 1 only if all selected bits are present. "
            "Default is any selected bit."
        ),
    )
    parser.add_argument(
        "--dtype",
        choices=("uint8", "uint16"),
        default="uint8",
        help="Output FITS integer type for the binary image.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes to use with --tile-list-file.",
    )
    parser.add_argument(
        "--problem-log",
        default="finalmask_to_binary_problems.log",
        help="Log file for missing or unreadable inputs when processing multiple tiles.",
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


def split_csv_like(values: list[str]) -> list[str]:
    items: list[str] = []
    for value in values:
        for part in value.split(","):
            stripped = part.strip()
            if stripped:
                items.append(stripped)
    return items


def normalize_alias(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def parse_selected_bits(args: argparse.Namespace) -> list[int]:
    selected: set[int] = set()

    for token in split_csv_like(args.bits):
        try:
            bit = int(token)
        except ValueError as exc:
            raise SystemExit(f"Invalid bit value: {token}") from exc
        if bit < 0 or bit > 15:
            raise SystemExit(f"Bit number out of range [0, 15]: {bit}")
        selected.add(bit)

    for name in split_csv_like(args.masks):
        key = normalize_alias(name)
        if key not in MASK_ALIASES:
            known = ", ".join(sorted(MASK_ALIASES))
            raise SystemExit(f"Unknown mask alias: {name}. Known aliases: {known}")
        selected.update(MASK_ALIASES[key])

    if not selected:
        raise SystemExit("Select at least one bit with --bits and/or one alias with --masks.")

    return sorted(selected)


def expand_selection_key(key: str) -> list[int]:
    normalized = normalize_alias(key)
    if normalized in MASK_ALIASES:
        return MASK_ALIASES[normalized]
    try:
        bit = int(normalized)
    except ValueError as exc:
        known = ", ".join(sorted(MASK_ALIASES))
        raise SystemExit(
            f"Unknown bit or mask alias: {key}. Known aliases: {known}"
        ) from exc
    if bit < 0 or bit > 15:
        raise SystemExit(f"Bit number out of range [0, 15]: {bit}")
    return [bit]


def parse_select_bits(args: argparse.Namespace) -> list[int]:
    selected: set[int] = set()
    for token in split_csv_like(args.select):
        key = token.split("=", 1)[0].strip()
        if not key:
            continue
        selected.update(expand_selection_key(key))
    if not selected:
        raise SystemExit("Select at least one mask with --select, e.g. --select trim manual_mask.")
    return sorted(selected)


def describe_bits(bits: list[int]) -> str:
    return ", ".join(f"{bit}:{BIT_INFO.get(bit, 'unknown')}" for bit in bits)


def print_mask_aliases() -> None:
    print("Accepted --masks keywords:")
    for name in sorted(MASK_ALIASES, key=lambda key: (MASK_ALIASES[key][0], len(MASK_ALIASES[key]), key)):
        bits = ", ".join(str(bit) for bit in MASK_ALIASES[name])
        print(f"  {name:<18} -> {bits}")


def print_state_syntax() -> None:
    print("Selection syntax for --select:")
    print("  key    -> activate this mask bit in the output binary mask")
    print("  key=*  -> same as key; any '=...' suffix is ignored for compatibility")
    print("  unspecified bits are ignored")
    print("Examples:")
    print("  --select trim")
    print("  --select manual_mask trim maximask")
    print("  --select 3 9 10")


def append_problem(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")
        handle.flush()


def output_dtype(name: str) -> np.dtype[np.generic]:
    return np.uint16 if name == "uint16" else np.uint8


def default_output_from_input(input_path: Path) -> Path:
    name = input_path.name
    if name.endswith(".fits.gz"):
        return input_path.with_name(name[:-8] + ".binary.fits.gz")
    if name.endswith(".fits"):
        return input_path.with_name(name[:-5] + ".binary.fits")
    return input_path.with_name(name + ".binary.fits")


def resolve_single_input_output(args: argparse.Namespace) -> tuple[Path, Path]:
    if not args.input:
        raise SystemExit("Single-file mode requires --input.")
    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else default_output_from_input(input_path)
    return input_path, output_path


def collect_tiles(args: argparse.Namespace) -> list[str]:
    if args.tile_list_file:
        return read_tile_list_file(Path(args.tile_list_file))
    if args.tile:
        return [args.tile]
    raise SystemExit("Provide either --input, or --tile / --tile-list-file.")


def input_output_for_tile(tile: str, args: argparse.Namespace) -> tuple[Path, Path]:
    input_path = Path(args.input_template.format(tile=tile))
    output_path = Path(args.output_template.format(tile=tile))
    return input_path, output_path


def build_binary_mask(data: np.ndarray, bits: list[int], require_all: bool) -> np.ndarray:
    values = np.asarray(data, dtype=np.uint16)
    if require_all:
        selected = np.ones(values.shape, dtype=bool)
        for bit in bits:
            selected &= (values & np.uint16(1 << bit)) != 0
        return selected

    bitmask = np.uint16(sum(1 << bit for bit in bits))
    return (values & bitmask) != 0


def process_one_input(
    input_path: Path,
    output_path: Path,
    bits: list[int],
    dtype_name: str,
) -> dict[str, object]:
    if not input_path.exists():
        raise FileNotFoundError(f"Input final mask does not exist: {input_path}")

    try:
        with fits.open(input_path, memmap=False) as hdul:
            header = hdul[0].header.copy()
            data = np.asarray(hdul[0].data)
    except Exception as exc:
        raise OSError(
            f"Unreadable FITS input {input_path} :: {type(exc).__name__}: {exc}"
        ) from exc

    binary_bool = build_binary_mask(data, bits, require_all=False)
    binary = binary_bool.astype(output_dtype(dtype_name), copy=False)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=binary, header=header).writeto(output_path, overwrite=True)

    return {
        "input": str(input_path),
        "output": str(output_path),
        "ones": int(np.count_nonzero(binary)),
        "zeros": int(binary.size - np.count_nonzero(binary)),
        "shape": tuple(int(v) for v in binary.shape),
    }


def process_one_tile(
    tile: str,
    args: argparse.Namespace,
    bits: list[int],
) -> dict[str, object]:
    input_path, output_path = input_output_for_tile(tile, args)
    result = process_one_input(
        input_path,
        output_path,
        bits,
        args.dtype,
    )
    result["tile"] = tile
    return result


def print_result(result: dict[str, object]) -> None:
    label = f"Tile {result['tile']}" if "tile" in result else "File"
    print(label)
    print(f"  Input: {result['input']}")
    print(f"  Output: {result['output']}")
    print(f"  Shape: {result['shape']}")
    print(f"  Ones: {result['ones']}")
    print(f"  Zeros: {result['zeros']}")


def main() -> int:
    args = parse_args()
    if args.list_masks:
        print_mask_aliases()
        return 0
    bits: list[int]
    if args.select:
        bits = parse_select_bits(args)
        print(f"Selected bits: {describe_bits(bits)}")
        print("Selection mode: binary superposition of selected masks")
    else:
        bits = parse_selected_bits(args)
        print(f"Selected bits: {describe_bits(bits)}")
        print(
            "Selection mode: all selected bits"
            if args.require_all
            else "Selection mode: any selected bit"
        )

    if args.input:
        input_path, output_path = resolve_single_input_output(args)
        result = process_one_input(
            input_path,
            output_path,
            bits,
            args.dtype,
        )
        print_result(result)
        return 0

    problem_log = Path(args.problem_log)
    problem_log.parent.mkdir(parents=True, exist_ok=True)
    problem_log.write_text("", encoding="utf-8")

    tiles = collect_tiles(args)
    processed_count = 0
    skipped_count = 0

    if args.jobs <= 1:
        for tile in tiles:
            print(f"Checking tile {tile}")
            try:
                result = process_one_tile(tile, args, bits)
                print_result(result)
                processed_count += 1
            except Exception as exc:
                skipped_count += 1
                message = f"{tile} :: {type(exc).__name__}: {exc}"
                print(f"WARNING: {message}")
                append_problem(problem_log, message)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_tile = {
                executor.submit(process_one_tile, tile, args, bits): tile
                for tile in tiles
            }
            for future in as_completed(future_to_tile):
                tile = future_to_tile[future]
                try:
                    result = future.result()
                    print_result(result)
                    processed_count += 1
                except Exception as exc:
                    skipped_count += 1
                    message = f"{tile} :: {type(exc).__name__}: {exc}"
                    print(f"WARNING: {message}")
                    append_problem(problem_log, message)

    if skipped_count:
        print(
            f"Finished with warnings: wrote {processed_count} binary mask(s), "
            f"skipped {skipped_count} tile(s)."
        )
        print(f"Problem log: {problem_log}")
    else:
        print(f"Finished: wrote {processed_count} binary mask(s), no skipped tiles.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
