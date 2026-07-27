#!/usr/bin/env python3
"""Run maximask on UNIONS r-band tiles and write final UNIONS maximask FITS files.

This script replaces the old shell pipeline by:
1. copying the native r-band tile into a temporary work directory,
2. running ``/opt/conda/bin/maximask``,
3. copying the source header/WCS onto the generated mask,
4. dilating the mask with a 5x5 structuring element,
5. writing the final gzip-compressed UNIONS maximask with values {0, 1024}.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import binary_dilation


DEFAULT_INPUT_TEMPLATE = (
    "/arc/projects/unions/catalogues/unions/GAaP_photometry/UNIONS_DR6/"
    "UNIONS.{tile}/r/UNIONS.{tile}_r.fits"
)
DEFAULT_OUTPUT_TEMPLATE = (
    "/arc/projects/unions/catalogues/unions/GAaP_photometry/UNIONS_DR6_maximask/"
    "UNIONS.{tile}.r.maximask.fits.gz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile", help="Single tile ID, e.g. 122.316.")
    parser.add_argument("--tile-list-file", help="Text file containing one tile ID per line.")
    parser.add_argument(
        "--input-template",
        default=DEFAULT_INPUT_TEMPLATE,
        help="Template for the native UNIONS r-band FITS file.",
    )
    parser.add_argument(
        "--output-template",
        default=DEFAULT_OUTPUT_TEMPLATE,
        help="Template for the final UNIONS maximask FITS file (.fits.gz).",
    )
    parser.add_argument(
        "--maximask-bin",
        default="/opt/conda/bin/maximask",
        help="Path to the maximask executable.",
    )
    parser.add_argument(
        "--config-dir",
        default=str(Path(__file__).resolve().parent / "maximask_config"),
        help="Path to the maximask configuration directory.",
    )
    parser.add_argument(
        "--work-root",
        default=None,
        help="Parent directory for temporary work directories. Defaults to /scratch if it exists, else system tmp.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes to use. Use 16 for a 16-CPU machine.",
    )
    parser.add_argument(
        "--problem-log",
        default="make_unions_maximask_problems.log",
        help="Log file where missing or failed tiles will be recorded.",
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


def append_problem(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{message}\n")
        handle.flush()


def default_work_root() -> Path:
    scratch = Path("/scratch")
    if scratch.exists():
        return scratch
    return Path(tempfile.gettempdir())


def run_maximask_command(
    tile: str,
    maximask_bin: str,
    config_dir: str,
    input_path: Path,
    stream_prefix: str | None,
) -> Path:
    command = [
        maximask_bin,
        "--config_dir",
        config_dir,
        "--single_mask",
        "1",
        str(input_path),
    ]
    if stream_prefix is None:
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"maximask failed for {input_path} (exit code {exc.returncode})"
            ) from exc

        mask_path = input_path.with_name(input_path.name.replace(".fits", ".mask.fits"))
        if not mask_path.exists():
            raise FileNotFoundError(f"maximask did not produce {mask_path}")
        return mask_path

    output_lines: list[str] = []
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    ) as process:
        assert process.stdout is not None
        for line in process.stdout:
            text = line.rstrip()
            output_lines.append(text)
            print(f"{stream_prefix}{text}", flush=True)
        return_code = process.wait()

    if return_code != 0:
        suffix = f" (exit code {return_code})"
        if output_lines:
            suffix += f" :: {' | '.join(output_lines[-10:])}"
        raise RuntimeError(f"maximask failed for {input_path}{suffix}")

    mask_path = input_path.with_name(input_path.name.replace(".fits", ".mask.fits"))
    if not mask_path.exists():
        suffix = f" :: {' | '.join(output_lines[-10:])}" if output_lines else ""
        raise FileNotFoundError(f"maximask did not produce {mask_path}{suffix}")
    return mask_path


def build_output_mask(source_fits: Path, mask_fits: Path, output_fits: Path) -> dict[str, object]:
    with fits.open(source_fits, memmap=False) as hdul:
        source_header = hdul[0].header.copy()
        source_shape = hdul[0].data.shape

    with fits.open(mask_fits, memmap=False) as hdul:
        mask_data = np.asarray(hdul[0].data)

    if mask_data.shape != source_shape:
        raise ValueError(
            f"Shape mismatch: source {source_fits} has {source_shape}, mask {mask_fits} has {mask_data.shape}"
        )

    bool_mask = mask_data != 0
    dilated_mask = binary_dilation(
        bool_mask,
        structure=np.ones((5, 5), dtype=bool),
        iterations=1,
    )
    output_data = np.where(dilated_mask, 1024, 0).astype(np.uint16, copy=False)

    output_fits.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=output_data, header=source_header).writeto(output_fits, overwrite=True)

    return {
        "shape": output_data.shape,
        "unique_values": np.unique(output_data).tolist(),
        "non_zero_pixels": int(np.count_nonzero(output_data)),
    }


def process_tile(
    tile: str,
    input_template: str,
    output_template: str,
    maximask_bin: str,
    config_dir: str,
    work_root: str | None,
    stream_prefix: str | None,
) -> dict[str, object]:
    input_path = Path(input_template.format(tile=tile))
    if not input_path.exists():
        raise FileNotFoundError(f"Input r-band FITS not found: {input_path}")

    output_path = Path(output_template.format(tile=tile))
    work_parent = Path(work_root) if work_root else default_work_root()
    work_parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"maximask_{tile}_", dir=work_parent) as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        temp_input = temp_dir / input_path.name
        shutil.copy2(input_path, temp_input)

        temp_mask = run_maximask_command(
            tile=tile,
            maximask_bin=maximask_bin,
            config_dir=config_dir,
            input_path=temp_input,
            stream_prefix=stream_prefix,
        )
        stats = build_output_mask(
            source_fits=input_path,
            mask_fits=temp_mask,
            output_fits=output_path,
        )

    return {
        "tile": tile,
        "input": str(input_path),
        "output": str(output_path),
        "shape": stats["shape"],
        "unique_values": stats["unique_values"],
        "non_zero_pixels": stats["non_zero_pixels"],
    }


def print_success(result: dict[str, object]) -> None:
    print(f"Tile {result['tile']}")
    print(f"  Input: {result['input']}")
    print(f"  Output: {result['output']}")
    print(f"  Shape: {result['shape']}")
    print(f"  Unique values: {result['unique_values']}")
    print(f"  Non-zero pixels: {result['non_zero_pixels']}")


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1.")

    problem_log = Path(args.problem_log)
    problem_log.parent.mkdir(parents=True, exist_ok=True)
    problem_log.write_text("", encoding="utf-8")

    tiles = collect_tiles(args)
    processed_count = 0
    skipped_count = 0

    if args.jobs == 1 or len(tiles) == 1:
        for tile in tiles:
            try:
                result = process_tile(
                    tile=tile,
                    input_template=args.input_template,
                    output_template=args.output_template,
                    maximask_bin=args.maximask_bin,
                    config_dir=args.config_dir,
                    work_root=args.work_root,
                    stream_prefix=None,
                )
                print_success(result)
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
        else:
            print(f"Finished: processed {processed_count} tile(s), no skipped tiles.")
        print(f"Problem log: {problem_log}")
        return 0

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_tile = {
            executor.submit(
                process_tile,
                tile,
                args.input_template,
                args.output_template,
                args.maximask_bin,
                args.config_dir,
                args.work_root,
                f"[{tile}] ",
            ): tile
            for tile in tiles
        }
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

    if skipped_count:
        print(
            f"Finished with warnings: processed {processed_count} tile(s), "
            f"skipped {skipped_count} tile(s)."
        )
    else:
        print(f"Finished: processed {processed_count} tile(s), no skipped tiles.")
    print(f"Problem log: {problem_log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
