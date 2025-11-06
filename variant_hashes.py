import os
import sys
import logging
import argparse
import pathlib

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def generate_individual_hashes(individual2cluster, variant2hash, haploblock2hash, chr_hash):
    """
    Generate individual hashes, ie 64-character strings of 0/1s, each contains:
        strand hash: 4 chars
        chromosome hash: 10 chars
        haploblock hash: 20 chars
        cluster hash: 20 chars
        individual hash: 10 chars

    arguments:
    - individual2cluster: dict, key=individual, value=unique clusterID
    - variant2hash: dict, key=individual, key=hash
    - haploblock2hash: dict, key=(start, end), key=hash
    - chr_hash: 10-digit
    """
    individual2hash = {}
    # generate cluster hashes

    for individual in individual2cluster:

    pass


def main(clusters_file, variant_hashes, haploblock_hashes):

    logger.info("Parsing clusters")
    # this should be fixed, beacuse we will have one clusters file per haploblock
    (individual2cluster, num_clusters) = data_parser.parse_clusters(clusters_file)
    logger.info("Found %i clusters with %i individuals in total", num_clusters, len(individual2cluster))

    logger.info("Parsing variant hashes")
    variant2hash = data_parser.parse_variant_hashes(variant_hashes)
    
    logger.info("Parsing haploblock hashes")
    haploblock2hash = data_parser.parse_haploblock_hashes(haploblock_hashes)

    # for now hardcode chromosome hash
    chr_hash = "0000100000"  # chr6
    variant2hash = generate_individual_hashes(individual2cluster, variant2hash, haploblock2hash, chr_hash)


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="TODO"
    )
    
    parser.add_argument('--clusters',
                        help='Path to clusters file generated with MMSeqs2',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--variant_hashes',
                        help='Path to file with variant hashes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--haploblock_hashes',
                        help='Path to file with haploblock hashes',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(clusters_file=args.clusters,
             variant_hashes=args.variant_hashes,
             haploblock_hashes=args.haploblock_hashes)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)