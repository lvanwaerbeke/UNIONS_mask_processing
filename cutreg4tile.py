#!/usr/bin/env python3
"""Create per-tile DS9 region files from a master sky-region catalogue.

This script mirrors the purpose of ``cutreg4tile.ipynb`` but is written as a
standalone program and tightens a few correctness details:

1. It uses the tile WCS from the FITS header to derive the tile sky footprint.
2. It keeps a polygon only if its projected pixel polygon intersects the tile
   image rectangle, so polygons grazing the image by even one pixel are kept.
3. It writes DS9 ``image`` coordinates in the correct 1-based convention.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


POLYGON_RE = re.compile(r"^\s*polygon\s*\(([^)]*)\)\s*(?:#.*)?$", re.IGNORECASE)
REGION_HEADERS = {
    "image",
    "physical",
    "fk5",
    "icrs",
    "j2000",
    "galactic",
    "ecliptic",
}
GRID_STEP_DEG = 1.0
EPS = 1.0e-12
def wrap_ra_deg(value: float) -> float:
    return value % 360.0


def wrap_delta_deg(delta: float) -> float:
    return ((delta + 180.0) % 360.0) - 180.0


def parse_float_card_value(card: str) -> float:
    return float(card.split("/")[0].strip().strip("'"))


def parse_str_card_value(card: str) -> str:
    return card.split("/")[0].strip().strip("'").strip()


def read_fits_header(path: Path) -> dict[str, str]:
    header: dict[str, str] = {}
    with path.open("rb") as handle:
        while True:
            block = handle.read(2880)
            if not block:
                break
            for offset in range(0, len(block), 80):
                card = block[offset : offset + 80].decode("ascii")
                key = card[:8].strip()
                if key == "END":
                    return header
                if card[8:10] == "= ":
                    header[key] = card[10:80]
    return header


@dataclass(frozen=True)
class Interval:
    low: float
    high: float


@dataclass
class TileInfo:
    tile_id: str
    ra_low_raw: float
    ra_high_raw: float
    dec_low: float
    dec_high: float

    @property
    def width_deg(self) -> float:
        return self.ra_high_raw - self.ra_low_raw

    @property
    def ra_intervals(self) -> list[Interval]:
        return continuous_interval_to_wrapped(self.ra_low_raw, self.ra_high_raw)

    @property
    def center_ra_deg(self) -> float:
        return wrap_ra_deg(self.ra_low_raw + 0.5 * self.width_deg)

    @property
    def center_dec_deg(self) -> float:
        return 0.5 * (self.dec_low + self.dec_high)


@dataclass
class TanWCS:
    naxis1: int
    naxis2: int
    crpix1: float
    crpix2: float
    crval1_deg: float
    crval2_deg: float
    cd11: float
    cd12: float
    cd21: float
    cd22: float

    @classmethod
    def from_fits(cls, path: Path) -> "TanWCS":
        header = read_fits_header(path)
        ctype1 = parse_str_card_value(header["CTYPE1"])
        ctype2 = parse_str_card_value(header["CTYPE2"])
        if ctype1 != "RA---TAN" or ctype2 != "DEC--TAN":
            raise ValueError(f"Unsupported WCS types: {ctype1} / {ctype2}")
        return cls(
            naxis1=int(parse_float_card_value(header["NAXIS1"])),
            naxis2=int(parse_float_card_value(header["NAXIS2"])),
            crpix1=parse_float_card_value(header["CRPIX1"]),
            crpix2=parse_float_card_value(header["CRPIX2"]),
            crval1_deg=parse_float_card_value(header["CRVAL1"]),
            crval2_deg=parse_float_card_value(header["CRVAL2"]),
            cd11=parse_float_card_value(header["CD1_1"]),
            cd12=parse_float_card_value(header["CD1_2"]),
            cd21=parse_float_card_value(header["CD2_1"]),
            cd22=parse_float_card_value(header["CD2_2"]),
        )

    def world_to_image(self, points: list[tuple[float, float]]) -> list[tuple[float, float]]:
        ra0 = math.radians(self.crval1_deg)
        dec0 = math.radians(self.crval2_deg)
        sin_dec0 = math.sin(dec0)
        cos_dec0 = math.cos(dec0)
        det = self.cd11 * self.cd22 - self.cd12 * self.cd21
        if abs(det) < EPS:
            raise ValueError("Singular CD matrix")

        inv11 = self.cd22 / det
        inv12 = -self.cd12 / det
        inv21 = -self.cd21 / det
        inv22 = self.cd11 / det

        projected: list[tuple[float, float]] = []
        for ra_deg, dec_deg in points:
            ra = math.radians(ra_deg)
            dec = math.radians(dec_deg)
            delta_ra = math.radians(wrap_delta_deg(ra_deg - self.crval1_deg))
            sin_dec = math.sin(dec)
            cos_dec = math.cos(dec)
            cosc = sin_dec0 * sin_dec + cos_dec0 * cos_dec * math.cos(delta_ra)
            if cosc <= 0.0:
                return []
            xi_deg = math.degrees(cos_dec * math.sin(delta_ra) / cosc)
            eta_deg = math.degrees(
                (cos_dec0 * sin_dec - sin_dec0 * cos_dec * math.cos(delta_ra)) / cosc
            )
            dx = inv11 * xi_deg + inv12 * eta_deg
            dy = inv21 * xi_deg + inv22 * eta_deg
            projected.append((self.crpix1 + dx, self.crpix2 + dy))
        return projected

    def image_to_world(self, points: list[tuple[float, float]]) -> list[tuple[float, float]]:
        ra0 = math.radians(self.crval1_deg)
        dec0 = math.radians(self.crval2_deg)
        sin_dec0 = math.sin(dec0)
        cos_dec0 = math.cos(dec0)

        world: list[tuple[float, float]] = []
        for x, y in points:
            dx = x - self.crpix1
            dy = y - self.crpix2
            xi_deg = self.cd11 * dx + self.cd12 * dy
            eta_deg = self.cd21 * dx + self.cd22 * dy
            xi = math.radians(xi_deg)
            eta = math.radians(eta_deg)

            denom = cos_dec0 - eta * sin_dec0
            ra = ra0 + math.atan2(xi, denom)
            dec = math.atan2(
                sin_dec0 + eta * cos_dec0,
                math.sqrt(denom * denom + xi * xi),
            )
            world.append((wrap_ra_deg(math.degrees(ra)), math.degrees(dec)))
        return world

    @property
    def image_rect(self) -> tuple[float, float, float, float]:
        return (0.5, self.naxis1 + 0.5, 0.5, self.naxis2 + 0.5)


@dataclass
class PolygonRecord:
    index: int
    sky_points: list[tuple[float, float]]
    ra_intervals: list[Interval]
    dec_low: float
    dec_high: float


def continuous_interval_to_wrapped(low_raw: float, high_raw: float) -> list[Interval]:
    width = high_raw - low_raw
    if width < 0:
        raise ValueError(f"Invalid interval: {low_raw} .. {high_raw}")
    low = wrap_ra_deg(low_raw)
    high = low + width
    if high <= 360.0 + EPS:
        return [Interval(low, min(high, 360.0))]
    return [Interval(low, 360.0), Interval(0.0, high - 360.0)]


def polygon_ra_intervals(points: list[tuple[float, float]]) -> list[Interval]:
    ras = [wrap_ra_deg(ra) for ra, _ in points]
    unwrapped = ras[:]
    if max(ras) - min(ras) > 180.0:
        unwrapped = [ra + 360.0 if ra < 180.0 else ra for ra in ras]
    return continuous_interval_to_wrapped(min(unwrapped), max(unwrapped))


def sky_bounds_from_points(points: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    ras = [wrap_ra_deg(ra) for ra, _ in points]
    decs = [dec for _, dec in points]
    unwrapped = ras[:]
    if max(ras) - min(ras) > 180.0:
        unwrapped = [ra + 360.0 if ra < 180.0 else ra for ra in ras]
    return min(unwrapped), max(unwrapped), min(decs), max(decs)


def parse_master_regions(path: Path) -> list[PolygonRecord]:
    polygons: list[PolygonRecord] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.lower() in REGION_HEADERS:
            continue
        match = POLYGON_RE.match(stripped)
        if match is None:
            continue
        values = [float(part.strip()) for part in match.group(1).split(",") if part.strip()]
        if len(values) < 6 or len(values) % 2:
            continue
        points = list(zip(values[0::2], values[1::2]))
        polygons.append(
            PolygonRecord(
                index=len(polygons),
                sky_points=points,
                ra_intervals=polygon_ra_intervals(points),
                dec_low=min(dec for _, dec in points),
                dec_high=max(dec for _, dec in points),
            )
        )
    return polygons


def sampled_border_points(wcs: TanWCS, samples_per_edge: int = 17) -> list[tuple[float, float]]:
    xmin, xmax, ymin, ymax = wcs.image_rect
    xs = np.linspace(xmin, xmax, samples_per_edge).tolist()
    ys = np.linspace(ymin, ymax, samples_per_edge).tolist()
    points: list[tuple[float, float]] = []
    points.extend((x, ymin) for x in xs)
    points.extend((x, ymax) for x in xs)
    points.extend((xmin, y) for y in ys[1:-1])
    points.extend((xmax, y) for y in ys[1:-1])
    return points


def tile_from_wcs(tile_id: str, wcs: TanWCS) -> TileInfo:
    border_world = wcs.image_to_world(sampled_border_points(wcs))
    ra_low_raw, ra_high_raw, dec_low, dec_high = sky_bounds_from_points(border_world)
    return TileInfo(
        tile_id=tile_id,
        ra_low_raw=ra_low_raw,
        ra_high_raw=ra_high_raw,
        dec_low=dec_low,
        dec_high=dec_high,
    )


def interval_cell_range(interval: Interval, step_deg: float) -> range:
    low_idx = int(math.floor(interval.low / step_deg))
    high_idx = int(math.floor(max(interval.high - EPS, interval.low) / step_deg))
    return range(low_idx, high_idx + 1)


def dec_cell_range(low: float, high: float, step_deg: float) -> range:
    low_idx = int(math.floor((low + 90.0) / step_deg))
    high_idx = int(math.floor((max(high - EPS, low) + 90.0) / step_deg))
    return range(low_idx, high_idx + 1)


def build_polygon_index(polygons: list[PolygonRecord], step_deg: float) -> dict[tuple[int, int], list[int]]:
    index: dict[tuple[int, int], list[int]] = defaultdict(list)
    for polygon in polygons:
        dec_cells = list(dec_cell_range(polygon.dec_low, polygon.dec_high, step_deg))
        for interval in polygon.ra_intervals:
            for ra_cell in interval_cell_range(interval, step_deg):
                for dec_cell in dec_cells:
                    index[(ra_cell % int(360.0 / step_deg), dec_cell)].append(polygon.index)
    return index


def candidate_polygon_ids(
    tile: TileInfo,
    spatial_index: dict[tuple[int, int], list[int]],
    step_deg: float,
) -> set[int]:
    candidates: set[int] = set()
    dec_cells = dec_cell_range(tile.dec_low, tile.dec_high, step_deg)
    n_ra_cells = int(360.0 / step_deg)
    for interval in tile.ra_intervals:
        for ra_cell in interval_cell_range(interval, step_deg):
            for dec_cell in dec_cells:
                candidates.update(spatial_index.get((ra_cell % n_ra_cells, dec_cell), ()))
    return candidates


def intervals_overlap(a: Interval, b: Interval) -> bool:
    return a.low <= b.high + EPS and b.low <= a.high + EPS


def sky_bbox_may_overlap(tile: TileInfo, polygon: PolygonRecord) -> bool:
    if polygon.dec_high < tile.dec_low - EPS or polygon.dec_low > tile.dec_high + EPS:
        return False
    return any(
        intervals_overlap(tile_interval, poly_interval)
        for tile_interval in tile.ra_intervals
        for poly_interval in polygon.ra_intervals
    )


def point_in_rect(point: tuple[float, float], rect: tuple[float, float, float, float]) -> bool:
    x, y = point
    xmin, xmax, ymin, ymax = rect
    return xmin - EPS <= x <= xmax + EPS and ymin - EPS <= y <= ymax + EPS


def orient(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> bool:
    return (
        min(a[0], c[0]) - EPS <= b[0] <= max(a[0], c[0]) + EPS
        and min(a[1], c[1]) - EPS <= b[1] <= max(a[1], c[1]) + EPS
    )


def segments_intersect(
    p1: tuple[float, float],
    p2: tuple[float, float],
    q1: tuple[float, float],
    q2: tuple[float, float],
) -> bool:
    o1 = orient(p1, p2, q1)
    o2 = orient(p1, p2, q2)
    o3 = orient(q1, q2, p1)
    o4 = orient(q1, q2, p2)

    def sgn(value: float) -> int:
        if value > EPS:
            return 1
        if value < -EPS:
            return -1
        return 0

    s1, s2, s3, s4 = map(sgn, (o1, o2, o3, o4))
    if s1 * s2 < 0 and s3 * s4 < 0:
        return True
    if s1 == 0 and on_segment(p1, q1, p2):
        return True
    if s2 == 0 and on_segment(p1, q2, p2):
        return True
    if s3 == 0 and on_segment(q1, p1, q2):
        return True
    if s4 == 0 and on_segment(q1, p2, q2):
        return True
    return False


def point_in_polygon(point: tuple[float, float], polygon: list[tuple[float, float]]) -> bool:
    x, y = point
    inside = False
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if abs(orient((x1, y1), (x2, y2), point)) <= EPS and on_segment((x1, y1), point, (x2, y2)):
            return True
        crosses = (y1 > y) != (y2 > y)
        if crosses:
            xinters = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if xinters >= x - EPS:
                inside = not inside
    return inside


def polygon_intersects_rect(
    polygon: list[tuple[float, float]], rect: tuple[float, float, float, float]
) -> bool:
    xmin, xmax, ymin, ymax = rect
    rect_points = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]

    if any(point_in_rect(point, rect) for point in polygon):
        return True
    if any(point_in_polygon(corner, polygon) for corner in rect_points):
        return True

    rect_edges = list(zip(rect_points, rect_points[1:] + rect_points[:1]))
    poly_edges = list(zip(polygon, polygon[1:] + polygon[:1]))
    for p1, p2 in poly_edges:
        for q1, q2 in rect_edges:
            if segments_intersect(p1, p2, q1, q2):
                return True
    return False


def format_region_polygon(pixel_points: list[tuple[float, float]]) -> str:
    coords = ",".join(f"{x:.6f},{y:.6f}" for x, y in pixel_points)
    return f"polygon({coords})"


def output_path_for_tile(template: str, tile_id: str) -> Path:
    return Path(template.format(tile=tile_id))


def fits_path_for_tile(template: str, tile_id: str) -> Path:
    return Path(template.format(tile=tile_id))


def build_tile_region_file(
    tile: TileInfo,
    wcs: TanWCS,
    polygons: list[PolygonRecord],
    spatial_index: dict[tuple[int, int], list[int]],
    output_path: Path,
) -> tuple[int, int]:
    candidates = candidate_polygon_ids(tile, spatial_index, GRID_STEP_DEG)
    kept_lines: list[str] = [
        "# Region file format: DS9 version 4.1",
        'global color=white',
        "image",
    ]
    kept_count = 0
    considered = 0

    for polygon_id in sorted(candidates):
        polygon = polygons[polygon_id]
        if not sky_bbox_may_overlap(tile, polygon):
            continue
        considered += 1
        pixel_points = wcs.world_to_image(polygon.sky_points)
        if not pixel_points:
            continue
        if polygon_intersects_rect(pixel_points, wcs.image_rect):
            kept_lines.append(format_region_polygon(pixel_points))
            kept_count += 1

    output_path.write_text("\n".join(kept_lines) + "\n", encoding="utf-8")
    return considered, kept_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--master-reg",
        default="large_galaxy.reg",
        help="Master DS9 region file in sky coordinates.",
    )
    parser.add_argument(
        "--fits-template",
        default="CFISimages/CFIS.{tile}.r.fits",
        help="Template for tile FITS filenames; the FITS WCS defines the tile sky footprint.",
    )
    parser.add_argument(
        "--output-template",
        default="UNIONS.{tile}.ext_polygon.reg",
        help="Template for output region filenames.",
    )
    parser.add_argument(
        "--tiles",
        nargs="*",
        help="Tile IDs to process, e.g. 145.287 001.247 234.311.",
    )
    parser.add_argument(
        "--tile-list-file",
        help="Text file containing one tile ID per line.",
    )
    parser.add_argument(
        "--missing-tiles-log",
        default="skipped_tiles.log",
        help="Path to a text file where skipped tile IDs will be written.",
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


def main() -> int:
    args = parse_args()
    master_path = Path(args.master_reg)
    polygons = parse_master_regions(master_path)
    spatial_index = build_polygon_index(polygons, GRID_STEP_DEG)

    if args.tile_list_file:
        tile_ids = read_tile_list_file(Path(args.tile_list_file))
    elif args.tiles:
        tile_ids = args.tiles
    else:
        raise SystemExit("Provide either --tiles or --tile-list-file.")

    skipped_tile_ids: list[str] = []
    processed_count = 0

    for tile_id in tile_ids:
        fits_path = fits_path_for_tile(args.fits_template, tile_id)
        if not fits_path.exists():
            print(f"WARNING: Missing FITS file for tile {tile_id}: {fits_path}")
            skipped_tile_ids.append(tile_id)
            continue
        wcs = TanWCS.from_fits(fits_path)
        tile = tile_from_wcs(tile_id, wcs)
        output_path = output_path_for_tile(args.output_template, tile_id)
        considered, kept = build_tile_region_file(
            tile=tile,
            wcs=wcs,
            polygons=polygons,
            spatial_index=spatial_index,
            output_path=output_path,
        )
        print(
            f"{tile_id}: considered {considered} candidate polygons, "
            f"kept {kept}, wrote {output_path} [fits WCS]"
        )
        processed_count += 1

    skipped_tiles_log = Path(args.missing_tiles_log)
    skipped_tiles_log.write_text(
        "".join(f"{tile_id}\n" for tile_id in skipped_tile_ids),
        encoding="utf-8",
    )

    if skipped_tile_ids:
        print(
            f"Finished with warnings: processed {processed_count} tile(s), "
            f"skipped {len(skipped_tile_ids)} tile(s)."
        )
        print(f"Skipped tile log written to {skipped_tiles_log}")
    else:
        print(f"Finished: processed {processed_count} tile(s), no skipped tiles.")
        print(f"Skipped tile log written to {skipped_tiles_log} (empty)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
