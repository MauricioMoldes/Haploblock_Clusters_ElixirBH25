import os
import logging
import pathlib

import subprocess

import numpy as np
import scipy.ndimage
import cupy as cp

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

def parse_recombination_rates(recombination_file, chromosome, use_gpu: bool = False, gpu_id: int = 0):
    """
    Parses recombination map from Halldorsson et al., 2019,
    use Gaussian smoothing to find high recombination rates.

    arguments:
    - recombination_file with header lines (lines starting with '#')
      and 5 columns per data line: Chr Begin End cMperMb cM

    - use_gpu: try to perform the smoothing and local-maximum detection on GPU (CuPy).
    - gpu_id: integer GPU id (if CuPy is available).

    returns:
    - haploblock_boundaries: list of tuples with haploblock boundaries (start, end)
    """
    if not str(chromosome).startswith("chr"):
        chromosome = "chr" + str(chromosome)

    positions = []
    rates = []

    # Read file - line by line is memory friendly and matches original behavior
    try:
        with open(recombination_file, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                split_line = line.rstrip().split('\t')
                if len(split_line) != 5:
                    logger.error("Recombination file %s has bad line (not 5 tab-separated fields): %s",
                                 recombination_file, line)
                    raise Exception("Bad line in the recombination file")
                (chr_field, start, end, rate, cM) = split_line
                if chr_field == chromosome:
                    try:
                        positions.append(int(start))
                        rates.append(float(rate))
                    except ValueError:
                        logger.error("Bad numeric values in line: %s", line)
                        raise

    except Exception as e:
        logger.error("Opening or parsing provided recombination file %s: %s", recombination_file, e)
        raise Exception("Cannot open/parse provided recombination file")

    if len(positions) == 0:
        logger.warning("No positions found for chromosome %s in %s", chromosome, recombination_file)
        return []

    # Convert to numpy arrays
    positions = np.array(positions, dtype=np.int64)
    rates = np.array(rates, dtype=np.float64)
    n = rates.size

    # If data is too short, avoid smoothing artifacts
    if n < 5:
        logger.info("Very short recombination map (n=%d) - returning single haploblock", n)
        return [(1, int(positions[-1]))]

    # Smoothing parameters (match prior behavior: sigma=5 in index units)
    sigma = 5.0

    # Try GPU smoothing + local maxima detection when requested
    if use_gpu:
        try:            
            # optional cupyx for ndimage; but we'll implement a safe GPU smoothing using convolution
            try:
                cp.cuda.Device(gpu_id).use()
            except Exception:
                logger.warning("Could not set GPU device to id %s; using default device.", gpu_id)

            # Build Gaussian kernel on GPU
            kernel_radius = int(max(1, np.ceil(3 * sigma)))
            kernel_x = cp.arange(-kernel_radius, kernel_radius + 1)
            kernel = cp.exp(-0.5 * (kernel_x / float(sigma)) ** 2)
            kernel = kernel / kernel.sum()

            # Move rates to GPU and convolve (mode='same' equivalent)
            rates_gpu = cp.asarray(rates)
            smoothed_gpu = cp.convolve(rates_gpu, kernel, mode='same')

            # find local maxima on GPU
            # create shifted arrays and compare
            left = cp.empty_like(smoothed_gpu)
            left[1:] = smoothed_gpu[:-1]
            left[0] = -cp.inf
            right = cp.empty_like(smoothed_gpu)
            right[:-1] = smoothed_gpu[1:]
            right[-1] = -cp.inf

            maxima_mask = (smoothed_gpu > left) & (smoothed_gpu > right)

            # bring maxima indices back to host
            maxima_idx = cp.asnumpy(cp.nonzero(maxima_mask)[0])

            # positions of maxima
            high_rates_positions = positions[maxima_idx].tolist()

            # if no peaks found, fallback to CPU smoothing (more robust in edge-cases)
            if len(high_rates_positions) == 0:
                logger.warning("No peaks detected on GPU smoothing - falling back to CPU smoothing.")
                raise RuntimeError("No peaks on GPU smoothing")

            # proceed to haploblocks creation below with high_rates_positions (already host list)
        except Exception as e:
            # if anything fails with GPU path, fallback gracefully to CPU
            logger.warning("GPU path failed (%s). Falling back to CPU.", e)
            use_gpu = False  # ensure CPU path below

    # CPU fallback path (or initial path if not using GPU)
    if not use_gpu:
        # Use scipy's gaussian_filter1d which is robust and matches original code
        try:
            smoothed = scipy.ndimage.gaussian_filter1d(rates, sigma=sigma)
        except Exception as e:
            # As last resort, use numpy convolution with computed kernel
            logger.warning("scipy.ndimage.gaussian_filter1d failed: %s. Using numpy convolution fallback.", e)
            kernel_radius = int(max(1, np.ceil(3 * sigma)))
            kernel_x = np.arange(-kernel_radius, kernel_radius + 1)
            kernel = np.exp(-0.5 * (kernel_x / float(sigma)) ** 2)
            kernel = kernel / kernel.sum()
            smoothed = np.convolve(rates, kernel, mode='same')

        # detect strict local maxima (same logic as before)
        # avoid endpoints as before
        maxima_idx = []
        # iterate in python to be explicit and safe for small arrays (arithmetic done index by index)
        for i in range(1, n - 1):
            if smoothed[i] > smoothed[i - 1] and smoothed[i] > smoothed[i + 1]:
                maxima_idx.append(i)

        high_rates_positions = positions[np.array(maxima_idx, dtype=int)].tolist()

    # Safety: if we still don't have any peak, create coarse breakpoints:
    if len(high_rates_positions) == 0:
        logger.warning("No recombination peaks detected for chromosome %s; returning single haploblock.", chromosome)
        return [(1, int(positions[-1]))]

    # Build haploblock boundaries using peaks as boundaries (mirrors original behavior)
    haploblock_boundaries = []

    # ensure sorted
    high_rates_positions = sorted(high_rates_positions)

    # add first haploblock from 1 to first peak
    haploblock_boundaries.append((1, int(high_rates_positions[0])))

    # intermediate haploblocks
    for i in range(1, len(high_rates_positions)):
        start = int(high_rates_positions[i - 1])
        end = int(high_rates_positions[i])
        haploblock_boundaries.append((start, end))

    # last haploblock from last peak to last position
    haploblock_boundaries.append((int(high_rates_positions[-1]), int(positions[-1])))

    logger.info("Found %i haploblocks", len(haploblock_boundaries))

    return haploblock_boundaries


def parse_haploblock_boundaries(boundaries_file):
    """
    Parses haploblock boundaries with header and 2 columns: start end
    Use Gaussian smoothing to find high recombination rates

    arguments:
    - boundaries_file
    
    returns:
    - haploblock_boundaries: list of tuples with haploblock boundaries (start, end)
    """
    haploblock_boundaries = []

    try:
        f = open(boundaries_file, 'r')
    except Exception as e:
        logger.error("Opening provided boundaries file %s: %s", boundaries_file, e)
        raise Exception("Cannot open provided boundaries file")
    
    # skip header
    line = f.readline()
    if not line.startswith("START\t"):
        logging.error("boundaries file %s is headerless? expecting headers but got %s",
                      boundaries_file, line)
        raise Exception("boundaries file problem")
    
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 2:
            logger.error("Boundaries file %s has bad line (not 2 tab-separated fields): %s",
                         boundaries_file, line)
            raise Exception("Bad line in the boundaries file")
        
        (start, end) = split_line

        haploblock_boundaries.append((start, end))

    return(haploblock_boundaries)


def parse_samples(samples_file):
    """
    Parse samples file for a population from 1000Genomes

    arguments:
    - samples file with 9 columns:
    Sample name, Sex, Biosample ID, Population code, Population name, Superpopulation code,
    Superpopulation name, Population elastic ID, Data collections

    returns:
    - list of sample names
    """
    samples = []

    try:
        f = open(samples_file, 'r')
    except Exception as e:
        logger.error("Opening provided samples file %s: %s", samples_file, e)
        raise Exception("Cannot open provided samples file")
    
    # skip header
    line = f.readline()
    if not line.startswith("Sample name\t"):
        logging.error("samples file %s is headerless? expecting headers but got %s",
                      samples_file, line)
        raise Exception("samples file problem")
    
    for line in f:
        split_line = line.rstrip().split('\t')
        
        (sample, *_) = split_line

        if not (sample.startswith("HG") or sample.startswith("NA")):
            logger.error("Samples file %s has bad line (not 9 tab-separated fields): %s",
                         samples_file, line)
            raise Exception("Bad line in the samples file")

        samples.append(sample)
    
    logger.info("Found %d samples", len(samples))
    return(samples)


def parse_samples_from_vcf(vcf):
    """
    Parse samples file for a population VCF from 1000Genomes

    arguments:
    - vcf file

    returns:
    - list of sample names
    """
    samples = []

    try:
        f = open(vcf, 'r')
    except Exception as e:
        logger.error("Opening provided vcf file %s: %s", vcf, e)
        raise Exception("Cannot open provided vcf file")
    
    samples = subprocess.run(["bcftools", "query",
                              "-l",
                              vcf],
                              check=True,
                              capture_output=True,
                              text=True).stdout.splitlines()

    logger.info("Found %d samples", len(samples))
    return(samples)


def parse_variants_of_interest(variants_file):
    """
    Parses variants of interest file with one variant per line in format chr:pos

    arguments:
    - variants_file: file with one variant per line (--variants),
        all must be in the same haploblock and formatatted as: "chr(number only):position"
    
    returns:
    - variants: list of string variant positions
    """
    variants = []

    try:
        f = open(variants_file, 'r')
    except Exception as e:
        logger.error("Opening provided variants file %s: %s", variants_file, e)
        raise Exception("Cannot open provided variants file")
    
    for line in f:
        line_split = line.rstrip().split(":")

        if len(line_split) != 2:
            logger.error("Variants file %s has bad line: %s",
                         variants_file, line)
            raise Exception("Bad line in the variants file")
        
        chr, position = line_split

        variants.append(position)

    return(variants)


def extract_region_from_vcf(vcf, chr, chr_map, start, end, out):
    """
    Extract variants from a specific region from a VCF file,
    if VCF has 6 instead of chr6, which will be required by bcftools consensus
    create file chr_map: "6 chr6" one mapping per line and provide it with --chr_map

    Generates the following files in out/tmp/:
    - {chr}_region_{start}-{end}.vcf.gz
    - {chr}_region_{start}-{end}.vcf.gz.csi
    - chr{chr}_region_{start}-{end}.vcf
    - chr{chr}_region_{start}-{end}.vcf.gz
    - chr{chr}_region_{start}-{end}.vcf.gz.csi

    returns:
    - output_vcf: pathlib.Path to bgzipped vcf for region
    """
    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    # extract region start-end from VCF and index
    temporary_vcf = os.path.join(out, "tmp", f"{chr}_region_{start}-{end}.vcf.gz")

    subprocess.run(["bcftools", "view",
                    "-r", f"{chr}:{start}-{end}",
                    "--min-af", "0.05",
                    vcf,
                    "-o", temporary_vcf],
                    check=True)

    subprocess.run(["bcftools", "index",
                    temporary_vcf],
                    check=True)
    
    # map chr6 to 6 in VCF, bgzip and index
    output_vcf = os.path.join(out, "tmp", f"chr{chr}_region_{start}-{end}.vcf")
    output_index = os.path.join(out, "tmp", f"chr{chr}_region_{start}-{end}.vcf.gz.csi")

    subprocess.run(["bcftools", "annotate",
                    "--rename-chrs", chr_map,
                    temporary_vcf],
                    stdout=open(output_vcf, "w"),
                    check=True)
    
    subprocess.run(["bgzip",
                    output_vcf],
                    check=True)

    output_vcf_bgzip = output_vcf + ".gz"

    subprocess.run(["bcftools", "index",
                    "-c",
                    "-o", output_index,
                    output_vcf_bgzip],
                    check=True)
    
    return(pathlib.Path(output_vcf_bgzip))


def extract_sample_from_vcf(vcf, sample, out):
    """
    Extract a specific sample from a VCF file
    
    returns:
    - output_vcf: pathlib.Path to bgzipped VCF
    """

    output_vcf = os.path.join(out, "tmp", sample + "_" + vcf.stem + ".gz")

    # extract sample from VCF and index
    subprocess.run(["bcftools", "view",
                    "--force-samples",  # only warn about unknown subset samples
                    "-s", sample,
                    "-o", output_vcf,
                    vcf],
                    check=True)

    subprocess.run(["bcftools", "index",
                    output_vcf],
                    check=True)
    
    return(pathlib.Path(output_vcf))


def extract_region_from_fasta(fasta, chr, start, end, out):
    """
    Extract a specific region from a fasta file
    
    returns:
    - output_fasta: pathlib.Path to fasta with region start-end
    """
    # index reference
    subprocess.run(["samtools", "faidx",
                    fasta],
                    check=True)
    
    output_fasta = os.path.join(out, "tmp", f"chr{chr}_region_{start}-{end}.fa")
    # extract region start-end from reference fasta
    subprocess.run(["samtools", "faidx",
                    fasta,
                    f"chr{chr}:{start}-{end}"],
                    stdout=open(output_fasta, "w"),
                    check=True)

    return(output_fasta)


def parse_clusters(clusters_file):
    """
    Parses clusters TSV file from MMSeqs2 (no header) with 2 columns: representative, individual
    assign unique ids for each cluster.
    We want to create unique cluster ID (starting at 0) based on cluster representatives,
    and match individual to cluster IDs.
    individual is a string formatted as:

    arguments:
    - clusters_file

    returns:
    - individual2cluster: dict, key=individual, value=unique clusterID
    - clusters: list of clusterIDs
    """
    try:
        f = open(clusters_file, 'r')
    except Exception as e:
        logger.error("Opening provided clusters file %s: %s", clusters_file, e)
        raise Exception("Cannot open provided clusters file")
    
    representative2cluster = {}
    individual2representative = {}
    individual2cluster = {}
    clusters = []
    num_clusters = 0
    for line in f:
        split_line = line.rstrip().split('\t')

        if len(split_line) != 2:
            logger.error("Clusters file %s has bad line (not 2 tab-separated fields): %s",
                         clusters_file, line)
            raise Exception("Bad line in the clusters file")
        
        (representative, individual) = split_line
        
        if not representative in representative2cluster:
            representative2cluster[representative] = num_clusters
            clusters.append(num_clusters)
            num_clusters += 1

        individual2representative[individual] = representative
    
    for individual in individual2representative:
        representative = individual2representative[individual]
        cluster = representative2cluster[representative]
        individual2cluster[individual] = cluster

    return(individual2cluster, clusters)