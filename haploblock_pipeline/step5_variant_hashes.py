#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import pathlib
import subprocess
from typing import Dict, List, Optional

import numpy

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
    hash = numpy.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH)
    return(hash)


def generate_haploblock_hashes(haploblock_boundaries: list[tuple[int, int]]) -> dict[tuple[int, int], str]:
    """
    Generate unique binary hashes for haploblocks.

    Uses multiprocessing for large datasets to speed up hash generation.

    arguments:
    - haploblock_boundaries: List of (start, end) tuples

    returns:
    - haploblock2hash: dict, key=(start, end), value=binary hash
    """
    n_haploblocks = len(haploblock_boundaries)
    if n_haploblocks == 0:
        logger.warning("No haploblocks provided.")
        return {}
    logger.info(f"Generating hashes for {n_haploblocks:,} haploblocks")

    # Decide execution strategy
    if n_haploblocks > PARALLEL_THRESHOLD:
        logger.info(f"Large dataset detected (> {PARALLEL_THRESHOLD:,} blocks). Using parallel generation.")
        with Pool(cpu_count()) as pool:
            hashes = pool.map(_make_hash, range(n_haploblocks))
    else:
        # Vectorized version for smaller datasets
        hashes = [numpy.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH) for i in range(n_haploblocks)]
    
    haploblock2hash = dict(zip(haploblock_boundaries, hashes))

    return(haploblock2hash)


def generate_cluster_hashes(clusters: pathlib.Path):
    """
    Generate unique binary hashes for clusters.

    Parallelization is not needed here — small overhead.

    arguments:
    - clusters: pathlib.Path, path to clusters file generated with MMSeqs2

    returns:
    - cluster2hash: dict, key=clusterID, values=hash
    """
    cluster2hash = {cluster: numpy.binary_repr(idx, width=CLUSTER_HASH_LENGTH)
                    for idx, cluster in enumerate(clusters)}

    return(cluster2hash)


def generate_variant_hashes(variants: List[str],
                            vcf: pathlib.Path,
                            chrom: str,
                            haploblock_boundaries: List[tuple],
                            samples: Optional[List[str]]) -> Dict[str, str]:
    """
    Generate unique binary hashes for variants of interest (if provided by user).
    
    arguments:
    - variants: list of string variant positions
    - vcf: 
    - chrom
    - haploblock_boundaries:
    - samples:

    returns:
    - variant2hash: dict, key=individual, values=hash

    """
    if not samples:
        logger.warning("No samples provided for variant hash generation. Returning empty dict.")
        return {}

    first_pos = variants[0]

    # find region containing first variant
    start = end = None
    for (s, e) in haploblock_boundaries:
        if int(s) <= int(first_pos) <= int(e):
            start, end = s, e
            break

    if start is None:
        raise ValueError(f"Variant {first_pos} not found in any haploblock")

    idx_map = {str(v): i for i, v in enumerate(variants)}
    variant2hash = {
        f"{sample}_chr{chrom}_region_{start}-{end}_hap{h}": ["0"] * len(variants)
        for sample in samples for h in (0, 1)
    }

    region = f"{chrom}:{first_pos}-{end}"
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n",
         "-s", ",".join(samples), "--force-samples", "-r", region, str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )

    for line in result.stdout.splitlines():
        _, pos, *gts = line.split("\t")
        if pos not in idx_map:
            continue

        i = idx_map[pos]
        for sample_idx, gt in enumerate(gts):
            if "|" not in gt:
                continue
            a0, a1 = gt.split("|")
            sample = samples[sample_idx]
            if a0 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap0"][i] = "1"
            if a1 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap1"][i] = "1"

    variant2hash = {k: "".join(v) for k, v in variant2hash.items()}

    return(variant2hash)


def generate_individual_hash(individual,
                             individual2cluster,
                             cluster2hash,
                             haploblock2hash,
                             chr_hash,
                             variant2hash=None):
    """Generate hash string for a single individual."""
    strand = individual[-1]
    strand_hash = "0001" if strand == "0" else "0010"

    # Extract start-end robustly
    individual_split = individual.split("_")
    region_str = individual_split[3].replace(".fa", "").replace(".fasta", "").replace(".vcf", "")
    start, end = region_str.split("-")
    haploblock_hash = haploblock2hash[(start, end)]
    cluster_hash = cluster2hash[individual2cluster[individual]]
    hash = strand_hash + chr_hash + haploblock_hash + cluster_hash
    if variant2hash:
        variant_hash = variant2hash[individual]
        hash += variant_hash

    return(individual, hash)


def generate_individual_hashes_parallel(individual2cluster, cluster2hash, haploblock2hash,
                                        chr_hash, variant2hash=None, max_workers=None):
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
                    haploblock2hash, chr_hash, variant2hash
                )
            )
        for fut in as_completed(futures):
            ind, h = fut.result()
            individual2hash[ind] = h

    return(individual2hash)


def haploblock_hashes_to_tsv(haploblock2hash: dict[tuple[int, int], str], out_dir: pathlib.Path) -> None:
    """Save haploblock hashes to TSV file."""
    out = out_dir / f"haploblock_hashes.tsv"
    with out.open('w') as f:
        f.write("START\tEND\tHASH\n")
        f.writelines(f"{start}\t{end}\t{hash_}\n" for (start, end), hash_ in haploblock2hash.items())
    logger.info(f"Haploblock hashes written to {out}")


def cluster_hashes_to_tsv(cluster2hash: dict[tuple[int, int], str],
                          start: int,
                          end: int,
                          out_dir: pathlib.Path) -> None:
    """Save cluster hashes to TSV file."""
    out = out_dir / f"cluster_hashes_{str(start)}-{str(end)}.tsv"
    with out.open("w") as f:
        f.write("CLUSTER\tHASH\n")
        for ind, h in cluster2hash.items():
            f.write(f"{ind}\t{h}\n")
    logger.info(f"Cluster hashes written to {out}")


def variant_hashes_to_TSV(variant2hash: Dict[str, str], out_dir: pathlib.Path) -> None:
    out = out_dir / "variant_hashes.tsv"
    with out.open("w") as f:
        f.write("VARIANT\tHASH\n")
        for ind, h in variant2hash.items():
            f.write(f"{ind}\t{h}\n")
    logger.info(f"Variant hashes written to {out}")


def individual_hashes_to_TSV(individual2hash: Dict[str, str],
                             start: int,
                             end: int,
                             out_dir: pathlib.Path):
    out = out_dir / f"individual_hashes_{str(start)}-{str(end)}.tsv"
    with open(out, 'w') as f:
        f.write("INDIVIDUAL\tHASH\n")
        for individual, hash_val in individual2hash.items():
            f.write(f"{individual}\t{hash_val}\n")
    logger.info(f"Individual hashes written to {out}")


def run_hashes(boundaries_file: pathlib.Path,
               clusters_dir: pathlib.Path,
               chrom: str,
               out: pathlib.Path,
               variants_file: Optional[pathlib.Path] = None,
               vcf: Optional[pathlib.Path] = None,
               samples_file: Optional[pathlib.Path] = None,
               threads: Optional[int] = None):

    chr_hash = numpy.binary_repr(int(chrom))

    logger.info("Generating haploblock hashes")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries)
    haploblock_hashes_to_tsv(haploblock2hash, out)

    variant2hash = None
    if variants_file:
        logger.info("Generating variant hashes")
        samples = (data_parser.parse_samples(samples_file)
               if samples_file else data_parser.parse_samples_from_vcf(vcf))
        variants = data_parser.parse_variants_of_interest(variants_file)
        variant2hash = generate_variant_hashes(variants, vcf, chrom, haploblock_boundaries, samples)
        variant_hashes_to_TSV(variant2hash, out)

    logger.info("Generating cluster hashes and individual hashes")
    for (start, end) in haploblock_boundaries:
        cluster_file = out / "clusters" / f"chr{chrom}_{start}-{end}_cluster.tsv"
        (individual2cluster, clusters) = data_parser.parse_clusters(cluster_file, len(clusters))
        cluster2hash = generate_cluster_hashes(clusters)
        cluster_hashes_to_tsv(cluster2hash, start, end, out)

        max_workers = threads or (os.cpu_count() - 1 or 1)
        individual2hash = generate_individual_hashes_parallel(
            individual2cluster, cluster2hash, haploblock2hash, chr_hash, variant2hash, max_workers=max_workers
        )

        individual_hashes_to_TSV(individual2hash, start, end, out)


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
        prog=script_name, description="Generate individual hashes"
    )
    parser.add_argument('--boundaries_file', type=pathlib.Path, required=True,
                        help="TSV file with header (START\tEND) and 2 columns: start end")
    parser.add_argument('--clusters_dir', required=True, type=pathlib.Path,
                        help='Path to directory with cluster files generated with MMSeqs2')
    parser.add_argument('--chr', required=True, type=str, help='Chromosome number')
    parser.add_argument('--out', required=True, type=pathlib.Path,
                        help='Output folder path')
    parser.add_argument("--variants", type=pathlib.Path, default=None,
                        help="File with one variant of interest per line")
    parser.add_argument("--vcf", type=pathlib.Path, default=None,
                        help="Phased VCF file")
    parser.add_argument("--samples", type=pathlib.Path, default=None,
                        help="TSV file with samples from 1000Genomes (optional)")
    parser.add_argument("--threads", type=int, default=None,
                        help="Number of worker threads for generating individual hashes (default: auto-detect CPU cores)")

    args = parser.parse_args()

    try:
        run_hashes(
            boundaries_file=args.boundaries_file,
            clusters_dir=args.clusters_dir,
            chrom=args.chr,
            out=args.out,
            variants_file=args.variants,
            vcf=args.vcf,
            samples_file=args.samples,
            threads=args.threads
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)

def run(boundaries_file, clusters_dir, chr, out, variants_file, vcf, samples_file, threads):
    run_hashes(
        pathlib.Path(boundaries_file),
        pathlib.Path(clusters_dir),
        str(chr),
        pathlib.Path(out),
        pathlib.Path(variants_file) if variants_file else None,
        pathlib.Path(vcf) if vcf else None,
        pathlib.Path(samples_file) if samples_file else None,
        threads
    )

