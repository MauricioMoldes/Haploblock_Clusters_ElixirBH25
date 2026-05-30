#!/usr/bin/env python3
import os
import time
import logging
import pathlib
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import data_parser

logger = logging.getLogger(__name__)


# ------------------------------------------------------------
# MMSEQS PARAMS
# ------------------------------------------------------------
def calculate_mmseq_params(variant_counts_file: pathlib.Path):

    haploblock2min_id = {}
    haploblock2cov_fraction = {}

    with open(variant_counts_file, "r") as f:
        header = f.readline()

        if not header.startswith("START\t"):
            raise ValueError("Missing header in variant counts file")

        for line in f:
            start, end, mean, stdev = line.strip().split("\t")

            start = int(start)
            end = int(end)

            hap_len = end - start

            haploblock2min_id[(start, end)] = 1 - (float(mean) / hap_len)
            haploblock2cov_fraction[(start, end)] = 1 - (682 / hap_len)

    return haploblock2min_id, haploblock2cov_fraction


# ------------------------------------------------------------
# BLOCK CLASSIFICATION
# ------------------------------------------------------------
def classify_block(start, end, fasta_path: pathlib.Path):
    """
    Returns: (class, size_bytes)
    """

    try:
        size = fasta_path.stat().st_size
    except Exception:
        return "unknown", 0

    # HARD EXCLUSION: gigantic blocks
    if size > 5 * 1024**3:   # > 5 GB
        return "excluded", size

    if size > 2 * 1024**3:
        return "large", size

    if size > 200 * 1024**2:
        return "medium", size

    return "small", size


# ------------------------------------------------------------
# MMSEQS RUNNER
# ------------------------------------------------------------
def compute_clusters(
    input_fasta: str,
    out: str,
    min_seq_id: float,
    cov_fraction: float,
    cov_mode: int,
    chrom: str,
    start: int,
    end: int,
    mmseq_threads: int,
    tmp_base: str
):

    t0 = time.time()

    input_fasta = str(pathlib.Path(input_fasta).resolve())

    output_prefix = pathlib.Path(out) / "clusters" / f"chr{chrom}_{start}-{end}"
    output_prefix = str(output_prefix.resolve())

    # --------------------------------------------------------
    # ISOLATED MMSEQS TMP PER BLOCK
    # --------------------------------------------------------
    tmp_dir = pathlib.Path(tmp_base) / f"chr{chrom}_{start}_{end}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = str(tmp_dir)

    # skip existing
    if pathlib.Path(f"{output_prefix}_cluster.tsv").exists():
        return True

    cmd = [
        "mmseqs",
        "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(cov_fraction),
        "--cov-mode", str(cov_mode),
        "--remove-tmp-files", "1",
        "--threads", str(mmseq_threads)
    ]

    try:
        subprocess.run(cmd, check=True)

        logger.info("DONE chr%s:%s-%s", chrom, start, end)
        return True

    except Exception as e:
        logger.error("FAILED chr%s:%s-%s | %s", chrom, start, end, e)
        return False


# ------------------------------------------------------------
# ADAPTIVE SCHEDULER
# ------------------------------------------------------------
def worker_config(block_class):
    """
    Adaptive CPU allocation
    """

    if block_class == "small":
        return 24, 4

    if block_class == "medium":
        return 16, 8

    if block_class == "large":
        return 8, 16

    return None, None  # excluded or unknown


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def run_clusters(
    boundaries_file: pathlib.Path,
    merged_consensus_dir: pathlib.Path,
    variant_counts_file: pathlib.Path,
    chrom: str,
    out_dir: pathlib.Path,
    cov_mode: int,
    threads: int | None,
    max_retries: int = 2
):

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "clusters").mkdir(exist_ok=True)
    (out_dir / "tmp").mkdir(exist_ok=True)

    tmp_base = out_dir / "mmseqs_tmp"
    tmp_base.mkdir(exist_ok=True)

    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)

    haploblock2min_id, haploblock2cov_fraction = calculate_mmseq_params(variant_counts_file)

    logger.info("Blocks found: %d", len(haploblock_boundaries))

    remaining_blocks = list(haploblock_boundaries)

    retries = 0

    while remaining_blocks and retries <= max_retries:

        failed = []

        with ThreadPoolExecutor(max_workers=8) as executor:

            futures = {}

            for start, end in remaining_blocks:

                fasta = merged_consensus_dir / f"chr{chrom}_region_{start}-{end}.fa"

                if not fasta.exists():
                    continue

                block_class, size = classify_block(start, end, fasta)

                if block_class == "excluded":
                    logger.warning("SKIP HUGE BLOCK %s-%s (%d bytes)", start, end, size)
                    continue

                mmseq_threads, workers_hint = worker_config(block_class)

                if mmseq_threads is None:
                    continue

                # IO throttling: small blocks = more concurrency
                future = executor.submit(
                    compute_clusters,
                    str(fasta),
                    str(out_dir),
                    haploblock2min_id[(start, end)],
                    haploblock2cov_fraction[(start, end)],
                    cov_mode,
                    chrom,
                    start,
                    end,
                    mmseq_threads,
                    str(tmp_base)
                )

                futures[future] = (start, end, block_class)

            for fut in as_completed(futures):

                start, end, cls = futures[fut]

                try:
                    ok = fut.result()
                    if not ok:
                        failed.append((start, end))
                except Exception:
                    failed.append((start, end))

        remaining_blocks = failed
        retries += 1

        logger.warning("Retry round %d | failed=%d", retries, len(failed))

    logger.info("DONE all clusters")
