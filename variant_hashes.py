#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import pathlib
import numpy

import data_parser

# Logger setup — inherits root config if used as a module
logger = logging.getLogger(__name__)

CLUSTER_HASH_LENGTH = 20


def generate_cluster_hashes(clusters):
    """
    Generate a unique binary hash for each cluster.

    Args:
        clusters (list): Cluster IDs.
    Returns:
        dict: clusterID -> binary hash (string of 0/1s)
    """
    cluster2hash = {}
    for idx, cluster in enumerate(clusters):
        cluster2hash[cluster] = numpy.binary_repr(idx, width=CLUSTER_HASH_LENGTH)
    return cluster2hash


def generate_individual_hashes(individual2cluster, cluster2hash, variant2hash, haploblock2hash, chr_hash):
    """
    Generate a 64-character binary hash per individual.

    Hash structure:
      - Strand hash: 4 bits
      - Chromosome hash: 10 bits
      - Haploblock hash: variable length
      - Cluster hash: 20 bits
      - Variant hash: variable length
    """
    individual2hash = {}

    for individual, cluster in individual2cluster.items():
        # Strand hash
        strand = individual[-1]
        strand_hash = "0001" if strand == "0" else "0010"

        # Haploblock hash
        parts = individual.split("_")
        start, end = parts[3].split("-")
        haploblock_hash = haploblock2hash[(start, end)]

        # Cluster and variant hashes
        cluster_hash = cluster2hash[cluster]
        variant_hash = variant2hash[individual]

        # Full individual hash
        full_hash = strand_hash + chr_hash + haploblock_hash + cluster_hash + variant_hash
        individual2hash[individual] = full_hash

    return individual2hash


def hashes_to_TSV(individual2hash, out):
    """
    Save individual → hash mapping to a TSV file.
    """
    output_path = os.path.join(out, "individual_hashes.tsv")
    with open(output_path, 'w') as f:
        f.write("INDIVIDUAL\tHASH\n")
        for individual, hash_val in individual2hash.items():
            f.write(f"{individual}\t{hash_val}\n")
    logger.info(f"Individual hashes written to {output_path}")


def main(clusters_file, variant_hashes, haploblock_hashes, chr, out):
    logger.info("Parsing clusters")
    individual2cluster, clusters = data_parser.parse_clusters(clusters_file)
    logger.info(f"Found {len(clusters)} clusters with {len(individual2cluster)} individuals")

    logger.info("Parsing variant hashes")
    variant2hash = data_parser.parse_variant_hashes(variant_hashes)

    logger.info("Parsing haploblock hashes")
    haploblock2hash = data_parser.parse_haploblock_hashes(haploblock_hashes)

    logger.info("Generating individual hashes")
    chr_hash = numpy.binary_repr(int(chr))
    cluster2hash = generate_cluster_hashes(clusters)
    individual2hash = generate_individual_hashes(
        individual2cluster, cluster2hash, variant2hash, haploblock2hash, chr_hash
    )

    logger.info("Saving individual hashes")
    hashes_to_TSV(individual2hash, out)

    # Perspective:
    # ─────────────
    # # Save also as JSON for programmatic access
    # import json
    # with open(os.path.join(out, "individual_hashes.json"), "w") as jf:
    #     json.dump(individual2hash, jf, indent=2)


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG,
    )
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name, description="Generate individual variant hashes per haploblock"
    )
    parser.add_argument('--clusters', required=True, type=pathlib.Path,
                        help='Path to clusters file generated with MMSeqs2')
    parser.add_argument('--variant_hashes', required=True, type=pathlib.Path,
                        help='Path to file with variant hashes')
    parser.add_argument('--haploblock_hashes', required=True, type=pathlib.Path,
                        help='Path to file with haploblock hashes')
    parser.add_argument('--chr', required=True, type=str, help='Chromosome number')
    parser.add_argument('--out', required=True, type=pathlib.Path,
                        help='Output folder path')

    args = parser.parse_args()

    try:
        main(
            clusters_file=args.clusters,
            variant_hashes=args.variant_hashes,
            haploblock_hashes=args.haploblock_hashes,
            chr=args.chr,
            out=args.out,
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {script_name}: {repr(e)}\n")
        sys.exit(1)

