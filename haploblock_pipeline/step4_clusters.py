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

# ----------------------------------------------------------------------
# MMseqs2 params
# ----------------------------------------------------------------------
def calculate_mmseq_params(variant_counts_file: pathlib.Path):

    haploblock2min_id = {}
    haploblock2cov_fraction = {}

    with open(variant_counts_file, "r") as f:
        header = f.readline()

        if not header.startswith("START\t"):
            raise ValueError(
                f"Variant counts file missing header: {header.strip()}"
            )

        for line in f:
            start, end, mean, stdev = line.strip().split("\t")

            start = int(start)
            end = int(end)

            hap_len = end - start

            haploblock2min_id[(start, end)] = (
                1 - (float(mean) / hap_len)
            )

            haploblock2cov_fraction[(start, end)] = (
                1 - (682 / hap_len)
            )

    return haploblock2min_id, haploblock2cov_fraction


# ----------------------------------------------------------------------
# Run clustering per FASTA
# ----------------------------------------------------------------------
def compute_clusters(
    input_fasta: str,
    out: str,
    min_seq_id: float,
    cov_fraction: float,
    cov_mode: int,
    chrom: str,
    start: int,
    end: int,
    mmseq_threads: int
):

    t0 = time.time()

    input_fasta = str(pathlib.Path(input_fasta).resolve())

    output_prefix = (
        pathlib.Path(out)
        / "clusters"
        / f"chr{chrom}_{start}-{end}"
    )

    output_prefix = str(output_prefix.resolve())

    # ------------------------------------------------------------------
    # CRITICAL FIX:
    # Dedicated tmp directory per haploblock
    # Avoids MMseqs tmp contention / filesystem collapse
    # ------------------------------------------------------------------
    tmp_dir = (
        pathlib.Path(out)
        / "tmp"
        / f"chr{chrom}_{start}_{end}"
    )

    tmp_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = str(tmp_dir.resolve())

    # ------------------------------------------------------------------
    # Skip existing results
    # ------------------------------------------------------------------
    if pathlib.Path(f"{output_prefix}_cluster.tsv").exists():

        logger.info(
            "Skipping existing cluster chr%s:%s-%s",
            chrom,
            start,
            end
        )

        return True

    cmd = [
        "mmseqs",
        "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(cov_fraction),
        "--cov-mode",
        str(cov_mode),
        "--remove-tmp-files",
        "1",
        "--threads",
        str(mmseq_threads)
    ]

    logger.debug("Running: %s", " ".join(cmd))

    try:

        subprocess.run(cmd, check=True)

        runtime = time.time() - t0

        logger.info(
            "Finished clustering chr%s:%s-%s | runtime=%.1fs",
            chrom,
            start,
            end,
            runtime
        )

        return True

    except subprocess.CalledProcessError as e:

        runtime = time.time() - t0

        logger.error(
            (
                "MMseqs failed chr%s:%s-%s "
                "| exit_code=%s "
                "| runtime=%.1fs "
                "| FASTA=%s"
            ),
            chrom,
            start,
            end,
            e.returncode,
            runtime,
            input_fasta
        )

        return False

    except Exception as e:

        runtime = time.time() - t0

        logger.error(
            (
                "Unexpected error chr%s:%s-%s "
                "| runtime=%.1fs "
                "| FASTA=%s "
                "| error=%s"
            ),
            chrom,
            start,
            end,
            runtime,
            input_fasta,
            str(e)
        )

        return False


# ----------------------------------------------------------------------
# Main workflow
# ----------------------------------------------------------------------
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

    haploblock_boundaries = (
        data_parser.parse_haploblock_boundaries(
            boundaries_file
        )
    )

    (
        haploblock2min_id,
        haploblock2cov_fraction
    ) = calculate_mmseq_params(variant_counts_file)

    logger.info(
        "Found %d haploblocks",
        len(haploblock_boundaries)
    )

    # ------------------------------------------------------------------
    # SMART CPU SCHEDULING
    # ------------------------------------------------------------------
    total_cpus = os.cpu_count() or 1

    cpu_budget = threads or total_cpus

    # ------------------------------------------------------------------
    # PERFORMANCE FIX:
    #
    # OLD:
    #   50 parallel MMseqs jobs × 4 threads
    #
    # NEW:
    #   fewer stronger jobs
    #
    # This dramatically reduces:
    #   - filesystem contention
    #   - tmp collisions
    #   - MMseqs DB overhead
    #   - IO saturation
    # ------------------------------------------------------------------
    mmseq_threads = 16

    max_workers = max(
        1,
        cpu_budget // mmseq_threads
    )

    # safety cap
    max_workers = min(max_workers, 12)

    logger.info(
        (
            "CPU budget=%d "
            "| mmseq_threads=%d "
            "| max_workers=%d"
        ),
        cpu_budget,
        mmseq_threads,
        max_workers
    )

    retries = 0

    remaining_blocks = list(haploblock_boundaries)

    while remaining_blocks and retries <= max_retries:

        logger.info(
            (
                "Clustering attempt %d "
                "for %d haploblocks"
            ),
            retries + 1,
            len(remaining_blocks)
        )

        futures = {}

        failed_blocks = []

        with ThreadPoolExecutor(
            max_workers=max_workers
        ) as executor:

            for start, end in remaining_blocks:

                input_fasta = (
                    merged_consensus_dir
                    / f"chr{chrom}_region_{start}-{end}.fa"
                )

                output_file = (
                    out_dir
                    / "clusters"
                    / f"chr{chrom}_{start}-{end}_cluster.tsv"
                )

                if output_file.exists():

                    logger.debug(
                        "Already clustered: chr%s_%s-%s",
                        chrom,
                        start,
                        end
                    )

                    continue

                future = executor.submit(
                    compute_clusters,
                    str(input_fasta),
                    str(out_dir),
                    haploblock2min_id[(start, end)],
                    haploblock2cov_fraction[(start, end)],
                    cov_mode,
                    chrom,
                    start,
                    end,
                    mmseq_threads
                )

                futures[future] = (start, end)

            for fut in as_completed(futures):

                start, end = futures[fut]

                try:

                    success = fut.result()

                    if not success:
                        failed_blocks.append((start, end))

                except Exception as e:

                    logger.error(
                        (
                            "Cluster job crashed "
                            "for chr%s:%s-%s "
                            "| error=%s"
                        ),
                        chrom,
                        start,
                        end,
                        str(e)
                    )

                    failed_blocks.append((start, end))

        if failed_blocks:

            retries += 1

            logger.warning(
                (
                    "Retrying %d failed haploblocks "
                    "(attempt %d)"
                ),
                len(failed_blocks),
                retries + 1
            )

            # ----------------------------------------------------------
            # SAFETY CHECK:
            # prevents endless IO-heavy retry storms
            # ----------------------------------------------------------
            if len(failed_blocks) == len(remaining_blocks):

                logger.error(
                    (
                        "All remaining blocks failed in this round. "
                        "Possible filesystem saturation or MMseqs issue."
                    )
                )

            remaining_blocks = failed_blocks

        else:
            break

    if remaining_blocks:

        logger.error(
            "Some haploblocks failed after %d retries:",
            max_retries
        )

        for start, end in remaining_blocks:

            logger.error(
                "FAILED chr%s:%s-%s",
                chrom,
                start,
                end
            )

    else:

        logger.info(
            (
                "All haploblocks clustered successfully "
                "after %d attempts."
            ),
            retries + 1
        )


# ----------------------------------------------------------------------
# Pipeline wrapper
# ----------------------------------------------------------------------
def run(
    boundaries_file,
    merged_consensus_dir,
    variant_counts,
    chr,
    out,
    cov_mode=2,
    threads=None,
    max_retries=2
):

    run_clusters(
        pathlib.Path(boundaries_file),
        pathlib.Path(merged_consensus_dir),
        pathlib.Path(variant_counts),
        str(chr),
        pathlib.Path(out),
        cov_mode=cov_mode,
        threads=threads,
        max_retries=max_retries
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=logging.INFO
    )

    parser = argparse.ArgumentParser(
        description=(
            "Cluster haploblock consensus sequences "
            "using MMseqs2"
        )
    )

    parser.add_argument(
        "--boundaries_file",
        type=pathlib.Path,
        required=True,
        help="TSV file with header (START END)"
    )

    parser.add_argument(
        "--merged_consensus_dir",
        type=pathlib.Path,
        required=True,
        help="Folder with merged phased FASTA files"
    )

    parser.add_argument(
        "--variant_counts",
        type=pathlib.Path,
        required=True,
        help="TSV file with START END MEAN STDEV"
    )

    parser.add_argument(
        "--chr",
        type=str,
        required=True,
        help="Chromosome"
    )

    parser.add_argument(
        "--out",
        type=pathlib.Path,
        required=True,
        help="Output folder"
    )

    parser.add_argument(
        "--cov_mode",
        type=int,
        default=2,
        help="MMSeqs2 coverage mode"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Total CPU budget"
    )

    parser.add_argument(
        "--max_retries",
        type=int,
        default=2,
        help=(
            "Maximum number of retries "
            "for failed haploblocks"
        )
    )

    args = parser.parse_args()

    run(
        boundaries_file=args.boundaries_file,
        merged_consensus_dir=args.merged_consensus_dir,
        variant_counts=args.variant_counts,
        chr=args.chr,
        out=args.out,
        cov_mode=args.cov_mode,
        threads=args.threads,
        max_retries=args.max_retries
    )
