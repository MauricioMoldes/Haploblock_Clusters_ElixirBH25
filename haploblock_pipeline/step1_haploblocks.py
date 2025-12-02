#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib
import data_parser

logger = logging.getLogger(__name__)


def haploblocks_to_tsv(haploblock_boundaries: list[tuple[int, int]],
                       chrom: str,
                       out_dir: pathlib.Path):
    """Save haploblock boundaries to a TSV file."""
    output_file = out_dir / f"haploblock_boundaries_chr{chrom}.tsv"
    with output_file.open('w') as f:
        f.write("START\tEND\n")
        f.writelines(f"{start}\t{end}\n" for start, end in haploblock_boundaries)


def run_haploblocks(recombination_file: pathlib.Path,
                    chrom: str,
                    out_dir: pathlib.Path,
                    gpu: bool = False,
                    gpu_id: int = 0):
    """
    Generate haploblock boundaries.
    GPU is optional (currently placeholder for future acceleration).
    """

    logger.info(f"Parsing recombination file: {recombination_file}")
    logger.info(f"GPU mode: {gpu} (gpu_id={gpu_id})")

    # ------------------------------------------------------------------
    # GPU ACCELERATION Perspective
    # ------------------------------------------------------------------
    if gpu:
        try:
            import cupy as cp
            logger.info("CuPy detected. GPU acceleration is available.")
        except Exception:
            logger.warning("GPU was requested but CuPy is not available. Falling back to CPU.")
            gpu = False

    # ------------------------------------------------------------------
    # CURRENTLY: always use CPU parsing (fast enough)
    # ------------------------------------------------------------------
    haploblock_boundaries = data_parser.parse_recombination_rates(
        recombination_file,
        chrom
    )

    # Create output directory if it doesn't exist
    out_dir.mkdir(parents=True, exist_ok=True)

    haploblocks_to_tsv(haploblock_boundaries, chrom, out_dir)
    logger.info(f"Haploblock boundaries written to {out_dir}")


# -------------------------------------------------------------------------
# CLI wrapper for standalone execution
# -------------------------------------------------------------------------
if __name__ == "__main__":
    script_name = pathlib.Path(sys.argv[0]).name

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="Generate haploblock boundaries using a recombination map."
    )
    parser.add_argument('--recombination_file', type=pathlib.Path, required=True,
                        help='Path to recombination file (Halldorsson et al., 2019)')
    parser.add_argument('--chr', type=str, required=True, help='Chromosome number')
    parser.add_argument('--out', type=pathlib.Path, required=True, help='Output folder path')

    args = parser.parse_args()

    try:
        run_haploblocks(args.recombination_file, args.chr, args.out)
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)


# -------------------------------------------------------------------------
# Pipeline entry point
# -------------------------------------------------------------------------
def run(recombination_file, chr, out, threads=None, gpu=False, gpu_id=0):
    """
    Pipeline-compatible run function.

    Arguments:
        recombination_file: Path to recombination file
        chr: Chromosome number
        out: Output directory
        threads: Optional, not used in this step (kept for consistency)
        gpu: Enable GPU mode
        gpu_id: GPU device index
    """
    run_haploblocks(
        pathlib.Path(recombination_file),
        str(chr),
        pathlib.Path(out),
        gpu=gpu,
        gpu_id=gpu_id
    )
