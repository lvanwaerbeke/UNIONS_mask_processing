#!/usr/bin/env python3
"""Compute RA/DEC cut boundaries from a UNIONS tile ID.

The convention implemented here is:

- dec_center = yyy / 2 - 90
- ra_center = (xxx / 2) / cos(dec_center)
- DECmin = dec_center - 0.25
- DECmax = dec_center + 0.25
- RAmin = ra_center - 0.25 / cos(dec_center)
- RAmax = ra_center + 0.25 / cos(dec_center)

RA bounds are wrapped into [0, 360).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile", help="Single tile ID, e.g. 122.316.")
    parser.add_argument("--tile-list-file", help="Text file containing one tile ID per line.")
    parser.add_argument(
        "--output",
        help="Optional output text file. If omitted, results are printed to stdout.",
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


def wrap_ra_deg(value: float) -> float:
    return value % 360.0


def compute_tile_cuts(tile: str) -> tuple[float, float, float, float]:
    try:
        xxx, yyy = tile.split(".")
        x = int(xxx)
        y = int(yyy)
    except ValueError as exc:
        raise ValueError(f"Invalid tile ID format: {tile}") from exc

    dec_center = y / 2.0 - 90.0
    cos_dec = math.cos(math.radians(dec_center))
    if abs(cos_dec) < 1.0e-12:
        raise ValueError(f"Tile center too close to pole for stable RA cut calculation: {tile}")

    ra_center = (x / 2.0) / cos_dec
    half_ra = 0.25 / cos_dec

    ra_min = wrap_ra_deg(ra_center - half_ra)
    ra_max = wrap_ra_deg(ra_center + half_ra)
    dec_min = dec_center - 0.25
    dec_max = dec_center + 0.25
    return ra_min, ra_max, dec_min, dec_max


def format_line(tile: str) -> str:
    ra_min, ra_max, dec_min, dec_max = compute_tile_cuts(tile)
    return f"{tile:>7s} {ra_min:14.8f} {ra_max:14.8f} {dec_min:14.8f} {dec_max:14.8f}"


def main() -> int:
    args = parse_args()
    lines = [format_line(tile) for tile in collect_tiles(args)]

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Wrote {len(lines)} tile cut(s) to {output_path}")
    else:
        for line in lines:
            print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
