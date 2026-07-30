#!/usr/bin/env python3
"""Clean DS9 maskstar region files and convert them to bit-mask FITS images."""

from __future__ import annotations

import argparse
import math
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from PIL import Image, ImageDraw


GLOBAL_RE = re.compile(r"^\s*global\b", re.IGNORECASE)
COORDSYS_RE = re.compile(r"^\s*(image|physical|fk5)\s*$", re.IGNORECASE)
POLYGON_RE = re.compile(r"^\s*(?:image\s*;\s*)?polygon\s*\(([^)]*)\)\s*(?:#\s*(.*))?$", re.IGNORECASE)
COLOR_RE = re.compile(r"(?:^|\s)color\s*=\s*([A-Za-z]+)", re.IGNORECASE)

COLOR_TO_MASK = {
    "cyan": np.uint16(1),
    "magenta": np.uint16(2),
    "green": np.uint16(4),
    "white": np.uint16(8),
}


@dataclass
class PolygonRegion:
    color: str
    vertices: np.ndarray


@dataclass
class CleanStats:
    total_polygons: int = 0
    kept_polygons: int = 0
    skipped_polygons: int = 0
    duplicate_vertices_removed: int = 0
    redundant_header_lines: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile", help="Single tile ID, e.g. 122.316.")
    parser.add_argument("--tile-list-file", help="Text file containing one tile ID per line.")
    parser.add_argument(
        "--input-template",
        default="UNIONS.{tile}_r_maskstars.reg",
        help="Input region filename template.",
    )
    parser.add_argument(
        "--cleaned-template",
        default="UNIONS.{tile}_r_maskstars.fixed.reg",
        help="Cleaned region filename template.",
    )
    parser.add_argument(
        "--output-template",
        default="UNIONS.{tile}_r_maskstars.mask.fits.gz",
        help="Output FITS mask filename template. The script always writes gzip-compressed FITS.",
    )
    parser.add_argument(
        "--source-templates",
        nargs="+",
        default=[
            "UNIONS.{tile}_ugriz.mask.fits",
            "UNIONS.{tile}_r.fits",
        ],
        help="Source FITS templates checked in order for WCS/header information.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes to use. Use 16 for a 16-CPU machine.",
    )
    parser.add_argument(
        "--problem-log",
        default="reg2fits_problems.log",
        help="Log file where missing or unreadable tiles will be recorded.",
    )
    return parser.parse_args()


def read_tile_list_file(path: Path) -> list[str]:
    tile_ids: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tile_ids.append(stripped.split()[0])
    return tile_ids


def extract_line_color(line_color: str | None, global_color: str) -> str:
    if line_color:
        match = COLOR_RE.search(line_color)
        if match:
            return match.group(1).lower()
    return global_color.lower()


def clean_vertices(vertices: np.ndarray) -> tuple[np.ndarray, int]:
    cleaned: list[tuple[float, float]] = []
    removed = 0
    for x, y in vertices:
        point = (float(x), float(y))
        if cleaned and point == cleaned[-1]:
            removed += 1
            continue
        cleaned.append(point)

    if len(cleaned) >= 2 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
        removed += 1

    if len(cleaned) < 3:
        return np.empty((0, 2), dtype=np.float64), removed

    unique_points = len({point for point in cleaned})
    if unique_points < 3:
        return np.empty((0, 2), dtype=np.float64), removed

    return np.asarray(cleaned, dtype=np.float64), removed


def parse_and_clean_region_file(path: Path) -> tuple[list[PolygonRegion], CleanStats]:
    polygons: list[PolygonRegion] = []
    stats = CleanStats()
    current_color = "green"

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#") or GLOBAL_RE.match(stripped) or COORDSYS_RE.match(stripped):
            stats.redundant_header_lines += 1
            if GLOBAL_RE.match(stripped):
                match = COLOR_RE.search(stripped)
                if match:
                    current_color = match.group(1).lower()
            continue

        match = POLYGON_RE.match(stripped)
        if match is None:
            stats.redundant_header_lines += 1
            continue

        values = [float(part.strip()) for part in match.group(1).split(",") if part.strip()]
        stats.total_polygons += 1
        if len(values) < 6 or len(values) % 2:
            stats.skipped_polygons += 1
            continue

        color = extract_line_color(match.group(2), current_color)
        if color not in COLOR_TO_MASK:
            stats.skipped_polygons += 1
            continue

        vertices = np.asarray(values, dtype=np.float64).reshape(-1, 2)
        cleaned_vertices, removed = clean_vertices(vertices)
        stats.duplicate_vertices_removed += removed
        if cleaned_vertices.size == 0:
            stats.skipped_polygons += 1
            continue

        polygons.append(PolygonRegion(color=color, vertices=cleaned_vertices))
        stats.kept_polygons += 1

    return polygons, stats


def format_float(value: float) -> str:
    return f"{value:.12f}".rstrip("0").rstrip(".")


def write_cleaned_region_file(path: Path, polygons: list[PolygonRegion]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii") as handle:
        handle.write("# Region file format: DS9 version 4.1\n")
        handle.write("image\n")
        for polygon in polygons:
            coords = ",".join(
                f"{format_float(x)},{format_float(y)}" for x, y in polygon.vertices
            )
            handle.write(f"polygon({coords}) # color={polygon.color}\n")


def choose_source_fits(tile: str, templates: list[str]) -> Path:
    for template in templates:
        candidate = Path(template.format(tile=tile))
        if candidate.exists():
            return candidate
    attempted = ", ".join(template.format(tile=tile) for template in templates)
    raise FileNotFoundError(f"No source FITS found for tile {tile}. Tried: {attempted}")


def polygon_bbox(vertices: np.ndarray, width: int, height: int) -> tuple[int, int, int, int] | None:
    xmin = max(1, int(math.floor(float(vertices[:, 0].min()))))
    xmax = min(width, int(math.ceil(float(vertices[:, 0].max()))))
    ymin = max(1, int(math.floor(float(vertices[:, 1].min()))))
    ymax = min(height, int(math.ceil(float(vertices[:, 1].max()))))
    if xmin > xmax or ymin > ymax:
        return None
    return xmin, xmax, ymin, ymax


def rasterize_polygon_mask(vertices: np.ndarray, bbox: tuple[int, int, int, int]) -> np.ndarray:
    xmin, xmax, ymin, ymax = bbox
    local_width = xmax - xmin + 1
    local_height = ymax - ymin + 1
    image = Image.new("L", (local_width, local_height), 0)
    draw = ImageDraw.Draw(image)
    shifted = [
        (float(x) - xmin + 1.0, float(y) - ymin + 1.0)
        for x, y in vertices
    ]
    draw.polygon(shifted, fill=1, outline=1)
    return np.asarray(image, dtype=np.uint8)


def build_mask(width: int, height: int, polygons: list[PolygonRegion]) -> np.ndarray:
    mask = np.zeros((height, width), dtype=np.uint16)
    for polygon in polygons:
        bbox = polygon_bbox(polygon.vertices, width, height)
        if bbox is None:
            continue
        xmin, xmax, ymin, ymax = bbox
        local_mask = rasterize_polygon_mask(polygon.vertices, bbox)
        bit_value = COLOR_TO_MASK[polygon.color]
        mask_slice = mask[ymin - 1 : ymax, xmin - 1 : xmax]
        mask_slice |= local_mask * bit_value
    return mask


def write_mask_fits(output_path: Path, source_fits: Path, mask: np.ndarray) -> None:
    with fits.open(source_fits, memmap=True) as hdul:
        header = hdul[0].header.copy()
    primary = fits.PrimaryHDU(data=mask.astype(np.uint16, copy=False), header=header)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    primary.writeto(output_path, overwrite=True)


def append_problem(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")
        handle.flush()


def process_tile(tile: str, args: argparse.Namespace) -> dict[str, object]:
    input_reg = Path(args.input_template.format(tile=tile))
    if not input_reg.exists():
        raise FileNotFoundError(f"Input region file not found: {input_reg}")

    cleaned_reg = Path(args.cleaned_template.format(tile=tile))
    output_fits = Path(args.output_template.format(tile=tile))
    if output_fits.suffix != ".gz":
        output_fits = Path(f"{output_fits}.gz")
    source_fits = choose_source_fits(tile, args.source_templates)

    polygons, stats = parse_and_clean_region_file(input_reg)
    write_cleaned_region_file(cleaned_reg, polygons)

    with fits.open(source_fits, memmap=True) as hdul:
        data_shape = hdul[0].data.shape
    height, width = int(data_shape[0]), int(data_shape[1])
    mask = build_mask(width, height, polygons)
    write_mask_fits(output_fits, source_fits, mask)

    return {
        "tile": tile,
        "input_reg": str(input_reg),
        "cleaned_reg": str(cleaned_reg),
        "source_fits": str(source_fits),
        "output_fits": str(output_fits),
        "kept_polygons": stats.kept_polygons,
        "total_polygons": stats.total_polygons,
        "skipped_polygons": stats.skipped_polygons,
        "duplicate_vertices_removed": stats.duplicate_vertices_removed,
        "unique_values": sorted(int(value) for value in np.unique(mask)),
        "non_zero_pixels": int(np.count_nonzero(mask)),
    }


def print_result(result: dict[str, object]) -> None:
    print(f"Tile {result['tile']}")
    print(f"  Input region file: {result['input_reg']}")
    print(f"  Cleaned region file: {result['cleaned_reg']}")
    print(f"  Source FITS: {result['source_fits']}")
    print(f"  Output FITS: {result['output_fits']}")
    print(
        "  Polygons: "
        f"{result['kept_polygons']}/{result['total_polygons']} kept, "
        f"{result['skipped_polygons']} skipped"
    )
    print(
        f"  Duplicate consecutive vertices removed: "
        f"{result['duplicate_vertices_removed']}"
    )
    print(f"  Unique mask values present: {result['unique_values']}")
    print(f"  Non-zero pixels: {result['non_zero_pixels']}")


def collect_tiles(args: argparse.Namespace) -> list[str]:
    if args.tile_list_file:
        return read_tile_list_file(Path(args.tile_list_file))
    if args.tile:
        return [args.tile]
    raise SystemExit("Provide either --tile or --tile-list-file.")


def main() -> int:
    args = parse_args()
    problem_log = Path(args.problem_log)
    problem_log.parent.mkdir(parents=True, exist_ok=True)
    problem_log.write_text("", encoding="utf-8")

    processed_count = 0
    skipped_count = 0
    tiles = collect_tiles(args)

    if args.jobs <= 1:
        for tile in tiles:
            try:
                result = process_tile(tile, args)
                print_result(result)
                processed_count += 1
            except Exception as exc:
                skipped_count += 1
                message = f"{tile} :: {type(exc).__name__}: {exc}"
                print(f"WARNING: {message}")
                append_problem(problem_log, message)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_tile = {executor.submit(process_tile, tile, args): tile for tile in tiles}
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
            f"Finished with warnings: processed {processed_count} tile(s), "
            f"skipped {skipped_count} tile(s)."
        )
        print(f"Problem log: {problem_log}")
    else:
        print(f"Finished: processed {processed_count} tile(s), no skipped tiles.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
