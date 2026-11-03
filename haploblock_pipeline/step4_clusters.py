#!/usr/bin/env python3
import os
import time
import logging
import pathlib
import argparse
import subprocess
import signal
from concurrent.futures import ThreadPoolExecutor, as_completed

import data_parser

logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
# Clustering mode -> MMseqs2 workflow
#   cluster  : mmseqs easy-cluster  (cascaded clustering; production default)
#   linclust : mmseqs easy-linclust (linear-time, k-mer based heuristic)
# Both write <prefix>_cluster.tsv / _rep_seq.fasta / _all_seqs.fasta, so
# step 5 reads either unchanged. Cluster membership can DIFFER between
# the two modes - never point linclust at a production out_dir.
# ----------------------------------------------------------------------
CLUSTER_MODES = {
    "cluster": "easy-cluster",
    "linclust": "easy-linclust",
}
DEFAULT_CLUSTER_MODE = "cluster"


def _resolve_cluster_mode(cluster_mode):
    """Normalize and validate the mode; None means the default."""
    if cluster_mode is None:
        return DEFAULT_CLUSTER_MODE
    mode = str(cluster_mode).strip().lower()
    if mode not in CLUSTER_MODES:
        raise ValueError(
            f"Unknown cluster_mode {cluster_mode!r}; "
            f"expected one of {sorted(CLUSTER_MODES)}"
        )
    return mode


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
# Run clustering per FASTA (STABLE v2 FIXED)
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
    mmseq_threads: int,
    cluster_mode: str = DEFAULT_CLUSTER_MODE
):

    t0 = time.time()

    cluster_mode = _resolve_cluster_mode(cluster_mode)
    mmseqs_workflow = CLUSTER_MODES[cluster_mode]

    input_fasta = str(pathlib.Path(input_fasta).resolve())

    output_prefix = (
        pathlib.Path(out)
        / "clusters"
        / f"chr{chrom}_{start}-{end}"
    )
    output_prefix = str(output_prefix.resolve())

    # ------------------------------------------------------------------
    # isolated MMseqs workspace
    # ------------------------------------------------------------------
    tmp_dir = (
        pathlib.Path(out)
        / "mmseqs_tmp"
        / f"chr{chrom}_{start}_{end}"
    )
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = str(tmp_dir.resolve())

    # ------------------------------------------------------------------
    # skip if done
    # (the skip check does not know which mode produced the existing
    # file - use a separate out_dir per mode)
    # ------------------------------------------------------------------
    if pathlib.Path(f"{output_prefix}_cluster.tsv").exists():
        logger.info(
            "Skipping chr%s:%s-%s (existing output; requested mode=%s)",
            chrom, start, end, cluster_mode
        )
        return True

    # ------------------------------------------------------------------
    # force isolated env
    # ------------------------------------------------------------------
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir
    env["MMSEQS_TMPDIR"] = tmp_dir
    env["APPTAINER_TMPDIR"] = tmp_dir
    env["SINGULARITY_TMPDIR"] = tmp_dir

    cmd = [
        "mmseqs",
        mmseqs_workflow,
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

    logger.debug("Running (%s): %s", cluster_mode, " ".join(cmd))

    # ================================================================
    # 🔥 STABLE EXECUTION (PROCESS GROUP CONTROL)
    # ================================================================
    try:

        p = subprocess.Popen(
            cmd,
            env=env,
            preexec_fn=os.setsid
        )

        timeout_sec = 3600  # 1h per haploblock safety cap

        while True:

            if p.poll() is not None:
                break

            if time.time() - t0 > timeout_sec:
                logger.error("TIMEOUT chr%s:%s-%s", chrom, start, end)
                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                return False

            time.sleep(5)

        returncode = p.poll()

        if returncode != 0:
            logger.error(
                "MMseqs failed chr%s:%s-%s | mode=%s | exit=%s",
                chrom, start, end, cluster_mode, returncode
            )
            return False

        runtime = time.time() - t0

        logger.info(
            "Finished chr%s:%s-%s | mode=%s | runtime=%.1fs",
            chrom, start, end, cluster_mode, runtime
        )

        return True

    except Exception as e:

        runtime = time.time() - t0

        logger.error(
            "Unexpected error chr%s:%s-%s | mode=%s | runtime=%.1fs | error=%s",
            chrom, start, end, cluster_mode, runtime, str(e)
        )

        return False


# ----------------------------------------------------------------------
# Main workflow (unchanged logic)
# ----------------------------------------------------------------------
def run_clusters(
    boundaries_file: pathlib.Path,
    merged_consensus_dir: pathlib.Path,
    variant_counts_file: pathlib.Path,
    chrom: str,
    out_dir: pathlib.Path,
    cov_mode: int,
    threads: int | None,
    max_retries: int = 2,
    cluster_mode: str = DEFAULT_CLUSTER_MODE
):

    cluster_mode = _resolve_cluster_mode(cluster_mode)

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "clusters").mkdir(exist_ok=True)

    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)

    haploblock2min_id, haploblock2cov_fraction = calculate_mmseq_params(variant_counts_file)

    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    cpu_budget = threads or (os.cpu_count() or 1)
    mmseq_threads = 16

    max_workers = max(1, cpu_budget // mmseq_threads)
    max_workers = min(max_workers, 12)

    logger.info(
        "CPU=%d mmseq_threads=%d workers=%d cluster_mode=%s",
        cpu_budget, mmseq_threads, max_workers, cluster_mode
    )

    remaining_blocks = sorted(haploblock_boundaries, key=lambda x: x[1] - x[0])

    retries = 0

    while remaining_blocks and retries <= max_retries:

        futures = {}
        failed = []

        with ThreadPoolExecutor(max_workers=max_workers) as ex:

            for start, end in remaining_blocks:

                inp = merged_consensus_dir / f"chr{chrom}_region_{start}-{end}.fa"

                fut = ex.submit(
                    compute_clusters,
                    str(inp),
                    str(out_dir),
                    haploblock2min_id[(start, end)],
                    haploblock2cov_fraction[(start, end)],
                    cov_mode,
                    chrom,
                    start,
                    end,
                    mmseq_threads,
                    cluster_mode
                )

                futures[fut] = (start, end)

            for f in as_completed(futures):

                start, end = futures[f]

                try:
                    if not f.result():
                        failed.append((start, end))
                except Exception:
                    failed.append((start, end))

        if failed:
            retries += 1
            remaining_blocks = failed
        else:
            break

    if remaining_blocks:
        logger.error("FAILED BLOCKS:")
        for s, e in remaining_blocks:
            logger.error("%s:%s-%s", chrom, s, e)
    else:
        logger.info("ALL DONE")


# ----------------------------------------------------------------------
# Wrapper for pipeline
# ----------------------------------------------------------------------
def run(boundaries_file, merged_consensus_dir, variant_counts, chr, out,
        threads=None, cov_mode=0, gpu=None, gpu_id=None, max_retries=2,
        cluster_mode=DEFAULT_CLUSTER_MODE):
    """
    gpu/gpu_id accepted-and-ignored: MMseqs2 has no GPU-accelerated
    clustering path (only its search/profile-alignment functions do),
    so this step stays CPU-only regardless of pipeline.gpu in the config.

    cluster_mode: "cluster" (mmseqs easy-cluster, default, production) or
    "linclust" (mmseqs easy-linclust). Use a separate `out` directory for
    linclust - cluster membership may differ from easy-cluster, and the
    skip-if-done check would otherwise reuse existing output.
    """
    run_clusters(
        boundaries_file=pathlib.Path(boundaries_file),
        merged_consensus_dir=pathlib.Path(merged_consensus_dir),
        variant_counts_file=pathlib.Path(variant_counts),
        chrom=str(chr),
        out_dir=pathlib.Path(out),
        cov_mode=cov_mode,
        threads=threads,
        max_retries=max_retries,
        cluster_mode=cluster_mode,
    )


# ----------------------------------------------------------------------
# CLI (standalone use / comparison runs without main.py)
# ----------------------------------------------------------------------
if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    parser = argparse.ArgumentParser(
        prog="step4_clusters",
        description="Cluster merged haploblock consensus FASTAs with MMseqs2 (Step 4)"
    )
    parser.add_argument('--boundaries_file', type=pathlib.Path, required=True)
    parser.add_argument('--merged_consensus_dir', type=pathlib.Path, required=True)
    parser.add_argument('--variant_counts', type=pathlib.Path, required=True)
    parser.add_argument('--chr', type=str, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True,
                        help="general out_dir (a clusters/ subdirectory is created under it)")
    parser.add_argument('--threads', type=int, default=None)
    parser.add_argument('--cov_mode', type=int, default=0)
    parser.add_argument('--cluster_mode', choices=sorted(CLUSTER_MODES),
                        default=DEFAULT_CLUSTER_MODE)
    parser.add_argument('--max_retries', type=int, default=2)

    args = parser.parse_args()

    run(
        boundaries_file=args.boundaries_file,
        merged_consensus_dir=args.merged_consensus_dir,
        variant_counts=args.variant_counts,
        chr=args.chr,
        out=args.out,
        threads=args.threads,
        cov_mode=args.cov_mode,
        max_retries=args.max_retries,
        cluster_mode=args.cluster_mode,
    )
