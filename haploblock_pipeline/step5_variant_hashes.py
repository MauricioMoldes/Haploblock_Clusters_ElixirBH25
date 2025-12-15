#!/usr/bin/env python3
import os
import logging
import pathlib
from typing import Dict, List, Optional

import numpy as np
import cupy as cp  # for GPU
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count

import data_parser

logger = logging.getLogger(__name__)

CLUSTER_HASH_LENGTH = 20
HAPLOBLOCK_HASH_LENGTH = 20
PARALLEL_THRESHOLD = 1000  # parallelize > 1000 haploblocks

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def _make_hash_cpu(i: int, width: int) -> str:
    return np.binary_repr(i, width=width)

def _make_hash_gpu(indices: cp.ndarray, width: int) -> cp.ndarray:
    """Convert integer indices to binary arrays on GPU (0/1)."""
    bin_array = ((indices[:, None] & (1 << cp.arange(width)[::-1])) > 0).astype(cp.uint8)
    return bin_array

# ---------------------------------------------------------------------
# Hash generators
# ---------------------------------------------------------------------
def generate_haploblock_hashes(haploblock_boundaries: list[tuple[int,int]], use_gpu=False):
    n = len(haploblock_boundaries)
    logger.info(f"Generating hashes for {n:,} haploblocks (GPU={use_gpu})")

    if use_gpu:
        indices = cp.arange(n, dtype=cp.uint32)
        bin_array = _make_hash_gpu(indices, HAPLOBLOCK_HASH_LENGTH)
        haploblock2hash = {hap: "".join(map(str, cp.asnumpy(bin_array[i]))) for i, hap in enumerate(haploblock_boundaries)}
    else:
        if n > PARALLEL_THRESHOLD:
            with Pool(cpu_count()) as pool:
                hashes = pool.starmap(_make_hash_cpu, [(i, HAPLOBLOCK_HASH_LENGTH) for i in range(n)])
        else:
            hashes = [_make_hash_cpu(i, HAPLOBLOCK_HASH_LENGTH) for i in range(n)]
        haploblock2hash = dict(zip(haploblock_boundaries, hashes))
    return haploblock2hash

def generate_cluster_hashes(clusters: list[int], use_gpu=False):
    n = len(clusters)
    if use_gpu:
        indices = cp.arange(n, dtype=cp.uint32)
        bin_array = _make_hash_gpu(indices, CLUSTER_HASH_LENGTH)
        cluster2hash = {cluster: "".join(map(str, cp.asnumpy(bin_array[i]))) for i, cluster in enumerate(clusters)}
    else:
        cluster2hash = {cluster: np.binary_repr(i, width=CLUSTER_HASH_LENGTH) for i, cluster in enumerate(clusters)}
    return cluster2hash

# ---------------------------------------------------------------------
# Individual hash generation (CPU / GPU)
# ---------------------------------------------------------------------
def generate_individual_hash(individual,
                             individual2cluster,
                             cluster2hash,
                             haploblock2hash,
                             chr_hash,
                             variant2hash=None):
    strand = individual[-1]
    strand_hash = "0001" if strand == "0" else "0010"

    # Extract start-end robustly
    individual_split = individual.split("_")
    region_str = individual_split[3].replace(".fa", "").replace(".fasta", "").replace(".vcf", "")
    start, end = map(int, region_str.split("-"))
    haploblock_hash = haploblock2hash[(start, end)]   
    cluster_hash = cluster2hash[individual2cluster[individual]]
    hash_val = strand_hash + chr_hash + haploblock_hash + cluster_hash
    if variant2hash:
        hash_val += variant2hash[individual]

    return individual, hash_val

def generate_individual_hashes_parallel(individual2cluster, cluster2hash, haploblock2hash,
                                        chr_hash, variant2hash=None, max_workers=None):
    max_workers = max_workers or (os.cpu_count() - 1 or 1)
    individual2hash = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(generate_individual_hash, ind, individual2cluster, cluster2hash,
                                   haploblock2hash, chr_hash, variant2hash)
                   for ind in individual2cluster]
        for fut in as_completed(futures):
            ind, h = fut.result()
            individual2hash[ind] = h
    return individual2hash

# ---------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------
def run_hashes(boundaries_file: pathlib.Path,
               clusters_dir: pathlib.Path,
               chrom: str,
               out: pathlib.Path,
               variants_file: Optional[pathlib.Path] = None,
               vcf: Optional[pathlib.Path] = None,
               samples_file: Optional[pathlib.Path] = None,
               threads: Optional[int] = None,
               gpu: bool = False,
               gpu_id: Optional[int] = None):

    if gpu:
        try:
            cp.cuda.Device(gpu_id or 0).use()
            logger.info(f"GPU enabled (ID={gpu_id or 0})")
        except Exception as e:
            logger.warning(f"GPU not available: {e}. Falling back to CPU.")
            gpu = False

    chr_hash = np.binary_repr(int(chrom))

    # Haploblock hashes
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries, use_gpu=gpu)

    # Save haploblock hashes
    out.mkdir(parents=True, exist_ok=True)
    haploblock_hashes_file = out / "haploblock_hashes.tsv"
    with haploblock_hashes_file.open("w") as f:
        f.write("START\tEND\tHASH\n")
        for (start, end), h in haploblock2hash.items():
            f.write(f"{start}\t{end}\t{h}\n")

    # Variant hashes
    variant2hash = None
    if variants_file:
        samples = data_parser.parse_samples(samples_file) if samples_file else data_parser.parse_samples_from_vcf(vcf)
        variants = data_parser.parse_variants_of_interest(variants_file)
        variant2hash = data_parser.generate_variant_hashes(variants, vcf, chrom, haploblock_boundaries, samples)

    # Cluster and individual hashes
    for (start, end) in haploblock_boundaries:
        cluster_file = clusters_dir / f"chr{chrom}_{start}-{end}_cluster.tsv"
        try:
            individual2cluster, clusters = data_parser.parse_clusters(cluster_file)
            cluster2hash = generate_cluster_hashes(clusters, use_gpu=gpu)

            # --- Save cluster hashes ---
            cluster_hash_file = out / f"cluster_hashes_{start}-{end}.tsv"
            with cluster_hash_file.open("w") as f:
                f.write("CLUSTER\tHASH\n")
                for cl, h in cluster2hash.items():
                    f.write(f"{cl}\t{h}\n")

            max_workers = threads or (os.cpu_count() - 1 or 1)
            individual2hash = generate_individual_hashes_parallel(
                individual2cluster, cluster2hash, haploblock2hash, chr_hash,
                variant2hash=variant2hash, max_workers=max_workers
            )

            # Save individual hashes
            out_file = out / f"individual_hashes_{start}-{end}.tsv"
            with out_file.open("w") as f:
                f.write("INDIVIDUAL\tHASH\n")
                for ind, h in individual2hash.items():
                    f.write(f"{ind}\t{h}\n")

        except Exception as e:
            logger.error(f"Failed processing haploblock {start}-{end}: {e}")

    logger.info("All hashes generated successfully")

# ---------------------------------------------------------------------
# Wrapper for pipeline
# ---------------------------------------------------------------------
def run(boundaries_file, clusters_dir, chr, out,
        variants_file=None, vcf=None, samples_file=None,
        threads=None, gpu=False, gpu_id=None):
    run_hashes(
        boundaries_file=pathlib.Path(boundaries_file),
        clusters_dir=pathlib.Path(clusters_dir),
        chrom=str(chr),
        out=pathlib.Path(out),
        variants_file=pathlib.Path(variants_file) if variants_file else None,
        vcf=pathlib.Path(vcf) if vcf else None,
        samples_file=pathlib.Path(samples_file) if samples_file else None,
        threads=threads,
        gpu=gpu,
        gpu_id=gpu_id
    )

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )
    logger = logging.getLogger(__file__)

    parser = argparse.ArgumentParser(
        prog="step5_variant_hashes",
        description="Generate haploblock, cluster, and individual hashes (Step 5)"
    )
    parser.add_argument('--boundaries_file', type=pathlib.Path, required=True)
    parser.add_argument('--clusters_dir', type=pathlib.Path, required=True)
    parser.add_argument('--chr', type=str, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument('--variants', type=pathlib.Path, default=None)
    parser.add_argument('--vcf', type=pathlib.Path, default=None)
    parser.add_argument('--samples', type=pathlib.Path, default=None)
    parser.add_argument('--threads', type=int, default=None)
    parser.add_argument('--gpu', action='store_true', help="Use GPU if available")
    parser.add_argument('--gpu_id', type=int, default=0, help="GPU ID (if multiple GPUs)")

    args = parser.parse_args()

    run(
        boundaries_file=args.boundaries_file,
        clusters_dir=args.clusters_dir,
        chr=args.chr,
        out=args.out,
        variants_file=args.variants,
        vcf=args.vcf,
        samples_file=args.samples,
        threads=args.threads,
        gpu=args.gpu,
        gpu_id=args.gpu_id
    )

