#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib
import numpy as np
import data_parser

logger = logging.getLogger(__name__)


def haploblocks_to_tsv(haploblock_boundaries: list[tuple[int, int]], chrom: str, out_dir: pathlib.Path) -> None:
    """Save haploblock boundaries to TSV file."""
    output_file = out_dir / f"haploblock_boundaries_chr{chrom}.tsv"
    with output_file.open('w') as f:
        f.write("START\tEND\n")
        f.writelines(f"{start}\t{end}\n" for start, end in haploblock_boundaries)


def run_haploblocks(recombination_file: pathlib.Path, chrom: str, out_dir: pathlib.Path) -> None:
    """
    Modular function to generate haploblock boundaries.

    Can be imported and called from another script/pipeline.
    """
    logger.info(f"Parsing recombination file: {recombination_file}")
    haploblock_boundaries = data_parser.parse_recombination_rates(recombination_file, chrom)

    # Ensure output directory exists
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Saving haploblock boundaries")
    haploblocks_to_tsv(haploblock_boundaries, chrom, out_dir)


# -------------------------------------------------------------------------
# CUDA Perspective (Optional, Commented Out)
# -------------------------------------------------------------------------
# For extremely large datasets and GPU-enabled environments, this
# section shows how hashing could be performed on CUDA for massive speedups.
#
# from numba import cuda
#
# @cuda.jit
# def generate_hashes_cuda(output, width):
#     i = cuda.grid(1)
#     if i < output.size:
#         val = i
#         for j in range(width):
#             bit = (val >> (width - j - 1)) & 1
#             output[i, j] = bit
#
# def generate_haploblock_hashes_cuda(n_blocks):
#     import numpy as np
#     output = np.zeros((n_blocks, HAPLOBLOCK_HASH_LENGTH), dtype=np.uint8)
#     threads_per_block = 256
#     blocks_per_grid = (n_blocks + (threads_per_block - 1)) // threads_per_block
#     generate_hashes_cuda[blocks_per_grid, threads_per_block](output, HAPLOBLOCK_HASH_LENGTH)
#     cuda.synchronize()
#     return ["".join(str(bit) for bit in row) for row in output]
# -------------------------------------------------------------------------


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
        description="Generate haploblock boundaries and hashes from recombination data."
    )
    parser.add_argument('--recombination_file', type=pathlib.Path, required=True,
                        help='Path to recombination file (Halldorsson et al., 2019)')
    parser.add_argument('--chr', type=str, required=True, help='Chromosome name or number')
    parser.add_argument('--out', type=pathlib.Path, required=True, help='Output folder path')

    args = parser.parse_args()

    try:
        run_haploblocks(args.recombination_file, args.chr, args.out)
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)


# Alias for pipeline
def run(recombination_file, chr, out, threads=None):
    """
    Pipeline-compatible run function.

    Arguments:
        recombination_file: Path to recombination file
        chr: Chromosome number
        out: Output directory
        threads: Optional, not used in this step (kept for consistency)
    """
    run_haploblocks(pathlib.Path(recombination_file), str(chr), pathlib.Path(out))
