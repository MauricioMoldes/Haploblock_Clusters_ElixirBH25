import os
import sys
import logging
import argparse
import pathlib

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def parse_individual_hashes(individual_hashes):
    pass

def parse_haploblock_hashes(haploblock_hashes):
    pass

def generate_variant_hases(individual2cluster):
    # generate cluster hashes
    pass



def main(clusters_file, individual_hashes, haploblock_hashes):

    logger.info("Parsing clusters")
    # this should be fixed, beacuse we will have one clusters file per haploblock
    (individual2cluster, num_clusters) = data_parser.parse_clusters(clusters_file)
    logger.info("Found %i clusters with %i individuals in total", num_clusters, len(individual2cluster))




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
    parser.add_argument('--individual_hashes',
                        help='Path to file with individual hashes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--haploblock_hashes',
                        help='Path to file with haploblock hashes',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(clusters_file=args.clusters,
             individual_hashes=args.individual_hashes,
             haploblock_hashes=args.haploblock_hashes)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)