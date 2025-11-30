#!/usr/bin/env python3
import sys
import os
import yaml
import argparse
from pathlib import Path

import step1_haploblocks
import step2_phased_sequences
import step3_merge_fasta
import step4_clusters
import step5_variant_hashes

from utils.logging import setup_logger


def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Run the Haploblocks Pipeline")
    parser.add_argument("--config", type=Path, required=True,
                        help="Path to YAML configuration file")
    parser.add_argument("--step", type=str, help="Pipeline step to run ('all' or select one between 1 and 5)")
    parser.add_argument("--threads", type=int, help="Number of CPU threads")

    args = parser.parse_args()
    
    logger = setup_logger()

    logger.info(f"Loading configuration from {args.config}")
    cfg = load_config(args.config)

    step = args.step or cfg["pipeline"]["step"]
    threads = args.threads if args.threads else cfg["pipeline"]["threads"]

    if threads == "auto":
        threads = max(1, (os.cpu_count() or 2) - 1)

    logger.info(f"Starting Haploblocks pipeline (Step: {step}, Threads: {threads})")

    try:
        # ---- STEP 1 ---------------------------------------------------------
        if step in ["1", "all"]:
            step1_haploblocks.run(
                recombination_file=cfg["data"]["recombination_file"],
                chr=cfg["chromosome"]["number"],
                out=cfg["outputs"]["out_dir"],
                threads=threads
            )

        # ---- STEP 2 ---------------------------------------------------------
        if step in ["2", "all"]:
            samples_file = cfg["data"].get("samples_file")
            samples_file = Path(samples_file) if samples_file else None

            step2_phased_sequences.run(
                boundaries_file=cfg["outputs"]["boundaries_file"],
                vcf=cfg["data"]["vcf"],
                ref=cfg["data"]["ref"],
                chr_map=cfg["data"]["chr_map"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["step2_out"]),
                samples_file=samples_file,
                threads=threads
            )

        # ---- STEP 3 ---------------------------------------------------------
        if step in ["3", "all"]:
            step3_merge_fasta.run(
                input_dir=Path(cfg["outputs"]["step3_input"]),
                output_dir=Path(cfg["outputs"]["step3_output"]),
                threads=threads
            )

        # ---- STEP 4 ---------------------------------------------------------
        if step in ["4", "all"]:
            step4_clusters.run(
                boundaries_file=cfg["outputs"]["boundaries_file"],
                merged_consensus_dir=cfg["outputs"]["merged_consensus_dir"],
                variant_counts=cfg["outputs"]["variant_counts"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["out_dir"]),
                threads=threads
            )

        # ---- STEP 5 ---------------------------------------------------------
        if step in ["5", "all"]:
            variants_file = cfg["data"].get("variants")
            variants_file = Path(variants_file) if variants_file else None
            vcf=Path(cfg["data"]["vcf"]) if variants_file else None

            samples_file = cfg["data"].get("samples_file")
            samples_file = Path(samples_file) if samples_file else None

            step5_variant_hashes.run(
                boundaries_file=cfg["outputs"]["boundaries_file"],
                clusters_file=cfg["outputs"]["clusters"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["out_dir"]),
                variants_file=variants_file,
                vcf=vcf,
                samples_file=samples_file if samples_file else None,
                threads=threads,
            )

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)

    logger.info("Pipeline finished successfully!")


if __name__ == "__main__":
    main()

