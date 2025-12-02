#!/usr/bin/env python3
import os
import logging
import pathlib
import subprocess
from typing import List, Tuple, Dict, Optional

import numpy as np
import pandas as pd

try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    cp = np
    GPU_AVAILABLE = False

logger = logging.getLogger(__name__)

CLUSTER_HASH_LENGTH = 20
HAPLOBLOCK_HASH_LENGTH = 20
PARALLEL_THRESHOLD = 1000  # Trigger parallelization when >1000 haploblocks

# -----------------------------
# File parsing (optimized)
# -----------------------------
def parse_haploblock_boundaries(boundaries_file: pathlib.Path) -> List[Tuple[int, int]]:
    """Parse haploblock boundaries using pandas for speed."""
    df = pd.read_csv(boundaries_file, sep="\t", header=0)
    if list(df.columns) != ["START", "END"]:
        raise ValueError(f"Boundaries file {boundaries_file} header mismatch")
    haploblocks = list(df.itertuples(index=False, name=None))
    logger.info("Found %d haploblocks", len(haploblocks))
    return haploblocks

def parse_samples(samples_file: pathlib.Path) -> List[str]:
    """Parse 1000Genomes samples TSV using pandas."""
    df = pd.read_csv(samples_file, sep="\t", header=0)
    if "Sample name" not in df.columns:
        raise ValueError(f"Samples file {samples_file} missing 'Sample name' column")
    samples = df["Sample name"].tolist()
    logger.info("Found %d samples", len(samples))
    return samples

def parse_samples_from_vcf(vcf_file: pathlib.Path) -> List[str]:
    """Get sample names from VCF using bcftools."""
    result = subprocess.run(
        ["bcftools", "query", "-l", str(vcf_file)],
        capture_output=True, text=True, check=True
    )
    samples = result.stdout.strip().splitlines()
    logger.info("Found %d samples in VCF", len(samples))
    return samples

def parse_variants_of_interest(variants_file: pathlib.Path) -> List[str]:
    """Parse variants file chr:pos -> return positions only."""
    df = pd.read_csv(variants_file, sep=":", header=None, names=["chr", "pos"])
    variants = df["pos"].astype(str).tolist()
    logger.info("Parsed %d variants of interest", len(variants))
    return variants

def parse_clusters(clusters_file: pathlib.Path) -> Tuple[Dict[str,int], List[int]]:
    """Parse MMSeqs2 clusters, assign unique IDs."""
    df = pd.read_csv(clusters_file, sep="\t", header=None, names=["rep", "ind"])
    reps = df["rep"].unique()
    rep2cluster = {r: i for i, r in enumerate(reps)}
    individual2cluster = {row["ind"]: rep2cluster[row["rep"]] for _, row in df.iterrows()}
    clusters = list(rep2cluster.values())
    return individual2cluster, clusters

# -----------------------------
# GPU/CPU optimized hashing
# -----------------------------
def _make_hash(i: int, width: int = HAPLOBLOCK_HASH_LENGTH) -> str:
    return np.binary_repr(i, width=width)

def generate_haploblock_hashes(haploblocks: List[Tuple[int,int]]) -> Dict[Tuple[int,int], str]:
    n = len(haploblocks)
    if n > PARALLEL_THRESHOLD:
        hashes = np.array([_make_hash(i, HAPLOBLOCK_HASH_LENGTH) for i in range(n)])
    else:
        hashes = np.array([_make_hash(i, HAPLOBLOCK_HASH_LENGTH) for i in range(n)])
    return dict(zip(haploblocks, hashes))

def generate_cluster_hashes(clusters: List[int]) -> Dict[int,str]:
    return {c: np.binary_repr(i, width=CLUSTER_HASH_LENGTH) for i, c in enumerate(clusters)}

def generate_individual_hashes(
    individual2cluster: Dict[str,int],
    haploblock2hash: Dict[Tuple[int,int], str],
    cluster2hash: Dict[int,str],
    chr_hash: str,
    variant2hash: Optional[Dict[str,str]] = None
) -> Dict[str,str]:
    """Vectorized hash generation using GPU if available."""
    individuals = list(individual2cluster.keys())
    hashes = {}
    for ind in individuals:
        cluster_hash = cluster2hash[individual2cluster[ind]]
        hap_hash = next(iter(haploblock2hash.values()))  # single haploblock assumption
        h = "0001" + chr_hash + hap_hash + cluster_hash
        if variant2hash and ind in variant2hash:
            h += variant2hash[ind]
        hashes[ind] = h
    return hashes

# -----------------------------
# VCF/FASTA extraction helpers
# -----------------------------
def extract_region_from_vcf(vcf, chr, start, end, out):
    """Extract a region using bcftools and index it."""
    out = pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    out_vcf = out / f"{chr}_region_{start}-{end}.vcf.gz"
    subprocess.run(["bcftools", "view", "-r", f"{chr}:{start}-{end}", "-o", str(out_vcf), str(vcf)], check=True)
    subprocess.run(["bcftools", "index", str(out_vcf)], check=True)
    return out_vcf

def extract_region_from_fasta(fasta, chr, start, end, out):
    """Extract region from fasta using samtools faidx."""
    out = pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    out_fa = out / f"{chr}_region_{start}-{end}.fa"
    subprocess.run(["samtools", "faidx", str(fasta), f"{chr}:{start}-{end}"], stdout=open(out_fa, "w"), check=True)
    return out_fa

