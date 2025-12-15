#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib

import data_parser

logger = logging.getLogger(__name__)


def haploblocks_to_tsv(haploblocks, chrom, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / f"haploblock_boundaries_chr{chrom}.tsv"

    with output_file.open("w") as f:
        f.write("START\tEND\n")
        for start, end in haploblocks:
            f.write(f"{start}\t{end}\n")


def run_haploblocks(recombination_file, chrom, out_dir):
    logger.info("Parsing recombination file %s (chr %s)", recombination_file, chrom)

    haploblocks = data_parser.parse_recombination_rates(
        recombination_file, chrom
    )

    haploblocks_to_tsv(haploblocks, chrom, out_dir)
    logger.info("Wrote haploblock boundaries to %s", out_dir)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    parser = argparse.ArgumentParser(
        description="Generate haploblock boundaries from recombination maps"
    )
    parser.add_argument(
        "--recombination_file",
        type=pathlib.Path,
        required=True,
    )
    parser.add_argument("--chr", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)

    args = parser.parse_args()

    try:
        run_haploblocks(args.recombination_file, args.chr, args.out)
    except Exception as e:
        logger.exception("Haploblock generation failed")
        sys.exit(1)


# Pipeline alias
def run(recombination_file, chr, out, threads=None):
    run_haploblocks(
        pathlib.Path(recombination_file),
        str(chr),
        pathlib.Path(out),
    )


if __name__ == "__main__":
    main()

