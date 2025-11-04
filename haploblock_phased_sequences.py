import os
import sys
import logging
import argparse
import pathlib
import subprocess

import numpy

import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)


def count_variants(vcf):
    """
    Count the number of variants in the vcf file
    """
    # extract genotype (GT) strings
    GTs = subprocess.run(["bcftools", "query",
                          "-f", "[ %GT]",
                          vcf],
                          capture_output=True,
                          text=True,
                          check=True)

    count_0 = 0  # left haplotype
    count_1 = 0  # right haplotype

    for line in GTs.stdout.splitlines():
        line = line.strip()
        if not line or line == "." or line == "./." or line == ".|.":
            continue  # skip missing genotypes

        # Handle phased and unphased cases
        if "|" in line:
            parts = line.split("|")
        elif "/" in line:
            parts = line.split("/")
        else:
            continue  # skip malformed entries

        if len(parts) != 2:
            continue  # unexpected format, skip

        (left, right) = parts

        if left == "1":
            count_0 += 1
        if right == "1":
            count_1 += 1

    return(count_0, count_1)


def generate_consensus_fasta(fasta, vcf, out):
    """
    Apply variants from VCF to reference sequence

    Generates the following files in out/ :
    - ref_chr6.fa.gz.fai
    - ref_chr6.fa.gz.gzi
    - chr{chr}_region_{start}-{end}.fa.gz
    """
    output_hap1 = os.path.join(out, pathlib.Path(vcf.stem).stem + "_hap1.fa")  # removes .vcf.gz
    output_hap2 = os.path.join(out, pathlib.Path(vcf.stem).stem + "_hap2.fa")  # removes .vcf.gz
    
    # create a consensus sequence (fasta) from reference and variants extracted from VCF
    # haploid sequence 1
    subprocess.run(["bcftools", "consensus",
                    "-H", "1",
                    "-f", fasta,
                    vcf],
                    stdout=open(output_hap1, "w"),
                    check=True)
    
    # haploid sequence 2
    subprocess.run(["bcftools", "consensus",
                    "-H", "2",
                    "-f", fasta,
                    vcf],
                    stdout=open(output_hap2, "w"),
                    check=True)

    return(output_hap1, output_hap2)


def main(boundaries_file, samples_file, vcf, ref, chr_map, chr, out, variant_counts_file):
    # sanity check
    if not os.path.exists(boundaries_file):
        logger.error(f"File {boundaries_file} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(samples_file):
        logger.error(f"File {samples_file} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(vcf):
        logger.error(f"File {vcf} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(ref):
        logger.error(f"File {ref} does not exist.")
        raise Exception("File does not exist")
    if not os.path.exists(chr_map):
        logger.error(f"File {chr_map} does not exist.")
        raise Exception("File does not exist")

    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %i haploblocks", len(haploblock_boundaries))

    logger.info("Parsing samples")
    samples = data_parser.parse_samples(samples_file)
    # samples = data_parser.parse_samples_from_vcf(vcf)
    logger.info("Found %i samples", len(samples))

    # dict for variant counts, key=(start, end), value=list(mean, stdev)
    haploblock2count = {}

    for (start, end) in haploblock_boundaries:
        logger.info(f"Generating phased VCF for haploblock {start}-{end}")
        region_vcf = data_parser.extract_region_from_vcf(vcf, chr, chr_map, start, end, out)

        logger.info(f"Generating phased fasta for haploblock {start}-{end}")
        region_fasta = data_parser.extract_region_from_fasta(ref, chr, start, end, out)

        # list for the number variants (separate for each haplotype)
        haploblock_counts = []

        logger.info(f"Generating consensus fasta files for haploblock {start}-{end}")
        for sample in samples:
            # logger.info(f"Generating phased VCF for haploblock {start}-{end} for sample %s", sample)
            sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out)

            # calculate of the number variants
            (count_0, count_1) = count_variants(sample_vcf)
            haploblock_counts.append(count_0)
            haploblock_counts.append(count_1)

            (sample_hap1, sample_hap2) = generate_consensus_fasta(region_fasta, sample_vcf, out)

        # calculate mean and stdev of the number variants
        mean = sum(haploblock_counts) / len(haploblock_counts)
        stdev = numpy.std(haploblock_counts)

        haploblock2count[(start, end)] = [mean, stdev]
    
    data_parser.save_haploblock_variant_counts(haploblock2count, variant_counts_file)


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
    
    parser.add_argument('--boundaries_file',
                        help='Path to boundaries file generated from Halldorsson et al., 2019',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--samples_file',
                        help='Path to samples file from 1000Genomes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--vcf',
                        help='Path to phased VCF file (bgzipped) from 1000Genomes',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--ref',
                        help='Path to reference sequence (bgzipped)',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--chr_map',
                        help='Path to chr_map: one mapping per line, ie "6 chr6"',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--chr',
                        help='chromosome',
                        type=str,
                        required=True)
    parser.add_argument('--out',
                        help='Path to output folder',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--variant_counts_file',
                        help='Path to a file to variant counts (mean, stdev)',
                        type=pathlib.Path,
                        required=True)

    args = parser.parse_args()

    try:
        main(boundaries_file=args.boundaries_file,
             samples_file=args.samples_file,
             vcf=args.vcf,
             ref=args.ref,
             chr_map=args.chr_map,
             chr=args.chr,
             out=args.out,
             variant_counts_file=args.variant_counts_file)

    except Exception as e:
        # details on the issue should be in the exception name, print it to stderr and die
        sys.stderr.write("ERROR in " + script_name + " : " + repr(e) + "\n")
        sys.exit(1)