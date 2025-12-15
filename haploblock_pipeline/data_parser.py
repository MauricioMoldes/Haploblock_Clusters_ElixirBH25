#!/usr/bin/env python3
import logging
import pathlib
import subprocess
import os

import numpy as np
from scipy.ndimage import gaussian_filter1d

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Recombination map parsing (MAIN CPU HOTSPOT – optimized)
# ---------------------------------------------------------------------
def parse_recombination_rates(recombination_file, chromosome):
    """
    Fast CPU-optimized parsing of Halldorsson2019 recombination map.
    Returns list of (start, end) tuples.
    """
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"

    data = np.fromiter(
        ((int(line[1]), float(line[3]))
         for line in map(str.split, open(recombination_file))
         if not line[0].startswith("#") and line[0] != "Chr" and line[0] == chromosome),
        dtype=[('start', 'i8'), ('rate', 'f8')]
    )

    positions = data['start']
    rates = data['rate']

    if len(rates) < 3:
        raise ValueError(f"Not enough data points for {chromosome}")

    smoothed = gaussian_filter1d(rates, sigma=5)
    peaks = np.where(
        (smoothed[1:-1] > smoothed[:-2]) &
        (smoothed[1:-1] > smoothed[2:])
    )[0] + 1

    if len(peaks) == 0:
        raise ValueError(f"No recombination peaks detected for {chromosome}")

    peak_positions = positions[peaks]

    haploblocks = [(1, peak_positions[0])]
    haploblocks.extend(zip(peak_positions[:-1], peak_positions[1:]))
    haploblocks.append((peak_positions[-1], positions[-1]))

    logger.info("Found %d haploblocks for %s", len(haploblocks), chromosome)
    return haploblocks


# ---------------------------------------------------------------------
# Simple TSV parsers
# ---------------------------------------------------------------------
def parse_haploblock_boundaries(boundaries_file):
    """
    Parses haploblock boundaries TSV (header: START\tEND)
    Returns integer tuples.
    """
    with open(boundaries_file) as f:
        header = next(f)
        if not header.startswith("START\t"):
            raise ValueError("Boundaries file missing header")
        return [tuple(map(int, line.rstrip().split("\t"))) for line in f]


def parse_samples(samples_file):
    with open(samples_file) as f:
        header = next(f)
        if not header.startswith("Sample name\t"):
            raise ValueError("Samples file missing header")

        samples = []
        for line in f:
            sample = line.split("\t", 1)[0]
            if not (sample.startswith("HG") or sample.startswith("NA")):
                raise ValueError(f"Invalid sample line: {line}")
            samples.append(sample)

    logger.info("Found %d samples", len(samples))
    return samples


def parse_samples_from_vcf(vcf):
    samples = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()

    logger.info("Found %d samples", len(samples))
    return samples


def parse_variants_of_interest(variants_file):
    variants = []
    with open(variants_file) as f:
        for line in f:
            chr_pos = line.rstrip().split(":")
            if len(chr_pos) != 2:
                raise ValueError(f"Bad variant line: {line}")
            variants.append(chr_pos[1])
    return variants


# ---------------------------------------------------------------------
# VCF / FASTA extraction (unchanged, I/O bound)
# ---------------------------------------------------------------------
def extract_region_from_vcf(vcf, chr, chr_map, start, end, out):
    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    temp_vcf = tmp_dir / f"{chr}_region_{start}-{end}.vcf.gz"

    subprocess.run(
        ["bcftools", "view",
         "-r", f"{chr}:{start}-{end}",
         "--min-af", "0.05",
         vcf,
         "-o", temp_vcf],
        check=True,
    )
    subprocess.run(["bcftools", "index", temp_vcf], check=True)

    output_vcf = tmp_dir / f"chr{chr}_region_{start}-{end}.vcf"
    subprocess.run(
        ["bcftools", "annotate", "--rename-chrs", chr_map, temp_vcf],
        stdout=open(output_vcf, "w"),
        check=True,
    )
    subprocess.run(["bgzip", output_vcf], check=True)
    subprocess.run(
        ["bcftools", "index", "-c", output_vcf.with_suffix(".vcf.gz")],
        check=True,
    )

    return output_vcf.with_suffix(".vcf.gz")


def extract_sample_from_vcf(vcf, sample, out):
    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    output_vcf = tmp_dir / f"{sample}_{vcf.stem}.gz"

    subprocess.run(
        ["bcftools", "view",
         "--force-samples",
         "-s", sample,
         "-o", output_vcf,
         str(vcf)],
        check=True,
    )
    subprocess.run(["bcftools", "index", output_vcf], check=True)
    return output_vcf


def extract_region_from_fasta(fasta, chr, start, end, out):
    subprocess.run(["samtools", "faidx", fasta], check=True)

    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    output_fasta = tmp_dir / f"chr{chr}_region_{start}-{end}.fa"
    subprocess.run(
        ["samtools", "faidx", fasta, f"chr{chr}:{start}-{end}"],
        stdout=open(output_fasta, "w"),
        check=True,
    )
    return output_fasta


# ---------------------------------------------------------------------
# Cluster parsing (single-pass optimized)
# ---------------------------------------------------------------------
def parse_clusters(clusters_file):
    representative2cluster = {}
    individual2cluster = {}
    clusters = []
    next_cluster = 0

    with open(clusters_file) as f:
        for line in f:
            rep, indiv = line.rstrip().split("\t")
            if rep not in representative2cluster:
                representative2cluster[rep] = next_cluster
                clusters.append(next_cluster)
                next_cluster += 1
            individual2cluster[indiv] = representative2cluster[rep]

    return individual2cluster, clusters


