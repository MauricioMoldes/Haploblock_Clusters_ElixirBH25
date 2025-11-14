#!/usr/bin/env python3
import os
import sys
import logging
import pathlib
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import data_parser

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# Utility Functions
# ----------------------------------------------------------------------

def calculate_mmseq_params(variant_counts_file: pathlib.Path):
    haploblock2min_id = {}
    haploblock2cov_fraction = {}

    with open(variant_counts_file, "r") as f:
        header = f.readline().strip()
        if not header.startswith("START\t"):
            logger.error("Unexpected header in variant counts file: %s", header)
            raise ValueError("Variant counts file missing header")

        for line in f:
            split_line = line.strip().split("\t")
            if len(split_line) != 4:
                raise ValueError(f"Malformed line in {variant_counts_file}: {line}")

            start, end, mean, stdev = split_line
            hap_len = int(end) - int(start)

            haploblock2min_id[(start, end)] = 1 - (float(mean) / hap_len)
            haploblock2cov_fraction[(start, end)] = 1 - (682 / hap_len)

    return haploblock2min_id, haploblock2cov_fraction


def compute_clusters(input_fasta: str, out: str, min_seq_id: float, cov_fraction: float, cov_mode: int,
                     chrom: str, start: str, end: str):
    output_prefix = pathlib.Path(out) / "clusters" / f"chr{chrom}_{start}-{end}"
    tmp_dir = pathlib.Path(out) / "tmp"

    cmd = [
        "mmseqs", "easy-cluster",
        input_fasta,
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(cov_fraction),
        "--cov-mode", str(cov_mode),
        "--remove-tmp-files", "1"
    ]

    logger.debug("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


# ----------------------------------------------------------------------
# Main Workflow
# ----------------------------------------------------------------------

def run_clusters(boundaries_file: pathlib.Path, merged_consensus_dir: pathlib.Path,
                 variant_counts_file: pathlib.Path, chrom: str, out_dir: pathlib.Path,
                 cov_mode: int = 0, threads: int = None):

    for fpath in [boundaries_file, merged_consensus_dir, variant_counts_file]:
        if not fpath.exists():
            logger.error("Path not found: %s", fpath)
            raise FileNotFoundError(f"{fpath} does not exist")

    cluster_dir = out_dir / "clusters"
    if cluster_dir.exists():
        logger.error("Output directory %s already exists — remove it first", cluster_dir)
        raise FileExistsError("Output directory exists")

    cluster_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "tmp").mkdir(parents=True, exist_ok=True)

    logger.info("Parsing haploblock boundaries from %s", boundaries_file)
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    logger.info("Parsing variant counts file")
    haploblock2min_id, haploblock2cov_fraction = calculate_mmseq_params(variant_counts_file)

    if len(haploblock_boundaries) != len(haploblock2min_id):
        raise ValueError("Haploblock count mismatch between boundaries and variant counts")

    logger.info("Computing clusters with MMSeqs2 (parallel)...")

    # --- Parallel execution ---
    max_workers = threads or (os.cpu_count() - 1 or 1)
    futures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for (start, end) in haploblock_boundaries:
            chrom_str = str(chrom)
            start_str, end_str = str(start), str(end)
            input_fasta = merged_consensus_dir / f"chr{chrom_str}_region_{start_str}-{end_str}.fa"
            logger.info("Looking for FASTA: %s", input_fasta.resolve())
            if not input_fasta.exists():
                logger.warning("Skipping missing FASTA file: %s", input_fasta)
                continue
            futures.append(
                executor.submit(
                    compute_clusters,
                    str(input_fasta),
                    str(out_dir),
                    haploblock2min_id[(start, end)],
                    haploblock2cov_fraction[(start, end)],
                    cov_mode,
                    chrom, start, end
                )
            )

        # Wait for all futures
        for fut in as_completed(futures):
            fut.result()  # propagate exceptions if any

    logger.info("All clusters computed successfully.")


# ----------------------------------------------------------------------
# CUDA Perspective (Commented)
# ----------------------------------------------------------------------
# 💡 The current MMseqs2 command-line tool does not natively use CUDA.
#     For very large haploblocks, one could consider:
#     - Using GPU-accelerated alignment tools (e.g., cuBLAST, GPU-MMseqs)
#     - Offloading compute-intensive distance calculations to GPU via PyTorch or CuPy
#     - Preprocessing sequences in batches on GPU before clustering
#     These would require rewriting `compute_clusters` or replacing MMseqs2 backend.

# ----------------------------------------------------------------------
# Pipeline wrapper for main.py
# ----------------------------------------------------------------------

def run(boundaries_file, merged_consensus_dir, variant_counts, chr, out, cov_mode=0, threads=None):
    """
    Pipeline-compatible run function.
    Matches argument names used in main.py.
    """
    run_clusters(
        pathlib.Path(boundaries_file),
        pathlib.Path(merged_consensus_dir),
        pathlib.Path(variant_counts),
        str(chr),
        pathlib.Path(out),
        cov_mode=cov_mode,
        threads=threads
    )


# ----------------------------------------------------------------------
# CLI for standalone execution
# ----------------------------------------------------------------------

if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO
    )
    logger = logging.getLogger(script_name)

    import argparse
    parser = argparse.ArgumentParser(
        prog=script_name,
        description="Cluster haploblock consensus sequences using MMseqs2"
    )

    parser.add_argument("--boundaries_file", type=pathlib.Path, required=True)
    parser.add_argument("--merged_consensus_dir", type=pathlib.Path, required=True)
    parser.add_argument("--variant_counts", type=pathlib.Path, required=True)
    parser.add_argument("--chr", type=str, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--cov_mode", type=int, default=0)
    parser.add_argument("--threads", type=int, default=None)

    args = parser.parse_args()

    try:
        run(
            boundaries_file=args.boundaries_file,
            merged_consensus_dir=args.merged_consensus_dir,
            variant_counts=args.variant_counts,
            chr=args.chr,
            out=args.out,
            cov_mode=args.cov_mode,
            threads=args.threads
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {e}\n")
        sys.exit(1)

