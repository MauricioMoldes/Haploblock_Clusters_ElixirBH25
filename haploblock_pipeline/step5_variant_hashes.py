#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import pathlib
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count

import data_parser

logger = logging.getLogger(__name__)

CLUSTER_HASH_LENGTH = 20

HAPLOBLOCK_HASH_LENGTH = 20
#PARALLEL_THRESHOLD = 1_000_000  # Enable multiprocessing if > 1M haploblocks
PARALLEL_THRESHOLD = 1000  # DEV only option Trigger parallelization when >1000 haploblock


def _make_hash(i: int) -> str:
    """Helper function for multiprocessing."""
    return np.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH)


def generate_haploblock_hashes(haploblock_boundaries: list[tuple[int, int]]) -> dict[tuple[int, int], str]:
    """
    Generate unique binary hashes for each haploblock boundary.

    Uses multiprocessing for large datasets to speed up hash generation.

    Args:
        haploblock_boundaries: List of (start, end) tuples.

    Returns:
        Dict mapping (start, end) -> binary hash string.
    """
    n_blocks = len(haploblock_boundaries)
    logger.info(f"Generating hashes for {n_blocks:,} haploblocks")

    if n_blocks == 0:
        logger.warning("No haploblocks provided.")
        return {}

    # Decide execution strategy
    if n_blocks > PARALLEL_THRESHOLD:
        logger.info(f"Large dataset detected (> {PARALLEL_THRESHOLD:,} blocks). Using parallel generation.")
        with Pool(cpu_count()) as pool:
            hashes = pool.map(_make_hash, range(n_blocks))
    else:
        # Vectorized version for smaller datasets
        hashes = [np.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH) for i in range(n_blocks)]

    return dict(zip(haploblock_boundaries, hashes))


def generate_cluster_hashes(clusters):
    """
    Generate a unique binary hash for each cluster.
    Parallelization is not needed here — small overhead.
    """
    cluster2hash = {cluster: np.binary_repr(idx, width=CLUSTER_HASH_LENGTH)
                    for idx, cluster in enumerate(clusters)}
    return cluster2hash


def generate_individual_hash(individual, individual2cluster, cluster2hash,
                             variant2hash, haploblock2hash, chr_hash):
    """Generate hash string for a single individual."""
    strand = individual[-1]
    strand_hash = "0001" if strand == "0" else "0010"

    # Extract start-end robustly
    individual_split = individual.split("_")
    region_str = individual_split[3].replace(".fa", "").replace(".fasta", "").replace(".vcf", "")
    start, end = region_str.split("-")
    haploblock_hash = haploblock2hash[(start, end)]
    cluster_hash = cluster2hash[individual2cluster[individual]]
    variant_hash = variant2hash[individual]

    return individual, strand_hash + chr_hash + haploblock_hash + cluster_hash + variant_hash


def generate_individual_hashes_parallel(individual2cluster, cluster2hash, variant2hash,
                                        haploblock2hash, chr_hash, max_workers=None):
    """
    Parallelized version: generate hashes for all individuals using ThreadPoolExecutor.
    """
    haploblock2hash = {(str(s), str(e)): h for (s, e), h in haploblock2hash.items()}
    individual2hash = {}

    max_workers = max_workers or (os.cpu_count() - 1 or 1)
    futures = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for individual in individual2cluster:
            futures.append(
                executor.submit(
                    generate_individual_hash,
                    individual, individual2cluster, cluster2hash,
                    variant2hash, haploblock2hash, chr_hash
                )
            )
        for fut in as_completed(futures):
            ind, h = fut.result()
            individual2hash[ind] = h

    return individual2hash


def haploblock_hashes_to_tsv(haploblock2hash: dict[tuple[int, int], str], chrom: str, out_dir: pathlib.Path) -> None:
    """Save haploblock hashes to TSV file."""
    output_file = out_dir / f"haploblock_hashes_chr{chrom}.tsv"
    with output_file.open('w') as f:
        f.write("START\tEND\tHASH\n")
        f.writelines(f"{start}\t{end}\t{hash_}\n" for (start, end), hash_ in haploblock2hash.items())


def hashes_to_TSV(individual2hash, out):
    output_path = os.path.join(out, "individual_hashes.tsv")
    with open(output_path, 'w') as f:
        f.write("INDIVIDUAL\tHASH\n")
        for individual, hash_val in individual2hash.items():
            f.write(f"{individual}\t{hash_val}\n")
    logger.info(f"Individual hashes written to {output_path}")


def run_variant_hashes(boundaries_file, chrom: str, clusters_file, variant_hashes, chr, out, threads: int = None):
    logger.info("Generating haploblock hashes")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries)
    haploblock_hashes_to_tsv(haploblock2hash, chrom, out)

    logger.info("Parsing clusters")
    individual2cluster, clusters = data_parser.parse_clusters(clusters_file)
    logger.info(f"Found {len(clusters)} clusters with {len(individual2cluster)} individuals")

    logger.info("Parsing variant hashes")
    variant2hash = data_parser.parse_variant_hashes(variant_hashes)

    logger.info("Generating individual hashes (parallelized)")
    chr_hash = np.binary_repr(int(chr))
    cluster2hash = generate_cluster_hashes(clusters)
    max_workers = threads or (os.cpu_count() - 1 or 1)
    individual2hash = generate_individual_hashes_parallel(
        individual2cluster, cluster2hash, variant2hash, haploblock2hash, chr_hash, max_workers=max_workers
    )

    logger.info("Saving individual hashes")
    hashes_to_TSV(individual2hash, out)


# ----------------------------------------------------------------------
# CUDA perspective
# ----------------------------------------------------------------------
#    Bitwise hash computation is perfect for GPU acceleration:
#    - Represent hashes as uint8/uint32 arrays on GPU
#    - Use CuPy or PyTorch for batch concatenation of strand, chr, haploblock, cluster, variant hashes
#    - Can scale to millions of individuals efficiently
#    - Final conversion to string can be done on CPU before TSV export


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name, description="Generate individual variant hashes per haploblock"
    )
    parser.add_argument('--boundaries_file', type=pathlib.Path, required=True,
                        help='boundaries_file')
    parser.add_argument('--clusters', required=True, type=pathlib.Path,
                        help='Path to clusters file generated with MMSeqs2')
    parser.add_argument('--variant_hashes', required=True, type=pathlib.Path,
                        help='Path to file with variant hashes')
    parser.add_argument('--chr', required=True, type=str, help='Chromosome number')
    parser.add_argument('--out', required=True, type=pathlib.Path,
                        help='Output folder path')
    parser.add_argument(
    "--threads",
    type=int,
    default=None,
    help="Number of worker threads for generating individual hashes (default: auto-detect CPU cores)"
    )

    args = parser.parse_args()

    try:
        run_variant_hashes(
            args.boundaries_file,
            clusters_file=args.clusters,
            variant_hashes=args.variant_hashes,
            chr=args.chr,
            out=args.out,
            threads=args.threads
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)

def run(boundaries_file, clusters, variant_hashes, chr, out, threads=None):
    run_variant_hashes(
        pathlib.Path(boundaries_file),
        pathlib.Path(clusters),
        pathlib.Path(variant_hashes),
        str(chr),
        pathlib.Path(out),
        threads
    )

