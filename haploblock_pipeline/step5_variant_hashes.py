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
PARALLEL_THRESHOLD = 1000  # DEV only: trigger parallelization when >1000 haploblocks

def _make_hash(i: int) -> str:
    return numpy.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH)

def generate_haploblock_hashes(haploblock_boundaries: list[tuple[int,int]]) -> dict:
    n = len(haploblock_boundaries)
    if n == 0:
        logger.warning("No haploblocks provided.")
        return {}
    logger.info(f"Generating hashes for {n:,} haploblocks")

    if n > PARALLEL_THRESHOLD:
        logger.info(f"Large dataset (> {PARALLEL_THRESHOLD:,}) detected. Using multiprocessing.")
        with Pool(cpu_count()) as pool:
            hashes = pool.map(_make_hash, range(n))
    else:
        hashes = [numpy.binary_repr(i, width=HAPLOBLOCK_HASH_LENGTH) for i in range(n)]
    return dict(zip(haploblock_boundaries, hashes))

def generate_cluster_hashes(clusters: List[str]) -> dict:
    return {c: numpy.binary_repr(i, width=CLUSTER_HASH_LENGTH) for i,c in enumerate(clusters)}

def generate_variant_hashes(variants: List[str], vcf: pathlib.Path, chrom: str,
                            haploblock_boundaries: List[tuple], samples: Optional[List[str]]) -> Dict[str,str]:
    if not samples:
        logger.warning("No samples provided. Returning empty variant hash dict.")
        return {}

    first_pos = variants[0]
    start = end = None
    for s,e in haploblock_boundaries:
        if int(s) <= int(first_pos) <= int(e):
            start, end = s, e
            break
    if start is None:
        raise ValueError(f"Variant {first_pos} not in any haploblock")

    idx_map = {str(v): i for i,v in enumerate(variants)}
    variant2hash = {f"{s}_chr{chrom}_region_{start}-{end}_hap{h}":["0"]*len(variants)
                    for s in samples for h in (0,1)}

    region = f"{chrom}:{first_pos}-{end}"
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n",
         "-s", ",".join(samples), "--force-samples", "-r", region, str(vcf)],
        capture_output=True, text=True, check=True
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
            if a0=="1": variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap0"][i]="1"
            if a1=="1": variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap1"][i]="1"

    return {k:"".join(v) for k,v in variant2hash.items()}

def generate_individual_hash(individual, individual2cluster, cluster2hash, haploblock2hash,
                             chr_hash, variant2hash=None):
    strand = individual[-1]
    strand_hash = "0001" if strand=="0" else "0010"
    region_str = individual.split("_")[3].replace(".fa","").replace(".fasta","").replace(".vcf","")
    start, end = region_str.split("-")
    hap_hash = haploblock2hash[(start,end)]
    clus_hash = cluster2hash[individual2cluster[individual]]
    h = strand_hash + chr_hash + hap_hash + clus_hash
    if variant2hash:
        h += variant2hash[individual]
    return individual,h

def generate_individual_hashes_parallel(individual2cluster, cluster2hash, haploblock2hash,
                                        chr_hash, variant2hash=None, max_workers=None,
                                        gpu: bool=False, gpu_id: Optional[int]=None):
    haploblock2hash = {(str(s), str(e)): h for (s,e),h in haploblock2hash.items()}
    individual2hash = {}

    if gpu:
        import cupy as cp
        if gpu_id is not None:
            cp.cuda.Device(gpu_id).use()  # set specific GPU device
        chr_hash_arr = cp.array([int(b) for b in chr_hash], dtype=cp.uint8)
        hap_hash_arr = {k: cp.array([int(b) for b in v], dtype=cp.uint8) for k,v in haploblock2hash.items()}
        clus_hash_arr = {k: cp.array([int(b) for b in v], dtype=cp.uint8) for k,v in cluster2hash.items()}
        variant_hash_arr = {k: cp.array([int(b) for b in v], dtype=cp.uint8) for k,v in (variant2hash or {}).items()}

        for ind in individual2cluster:
            strand = ind[-1]
            strand_hash = cp.array([0,0,0,1] if strand=="0" else [0,0,1,0], dtype=cp.uint8)
            region_str = ind.split("_")[3].replace(".fa","").replace(".fasta","").replace(".vcf","")
            start, end = region_str.split("-")
            h = cp.concatenate([strand_hash, chr_hash_arr, hap_hash_arr[(start,end)], clus_hash_arr[individual2cluster[ind]]])
            if variant_hash_arr: h = cp.concatenate([h, variant_hash_arr[ind]])
            individual2hash[ind] = "".join(h.get().astype(str))
        return individual2hash

    # CPU fallback
    max_workers = max_workers or (os.cpu_count()-1 or 1)
    futures=[]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for ind in individual2cluster:
            futures.append(executor.submit(generate_individual_hash, ind, individual2cluster, cluster2hash,
                                           haploblock2hash, chr_hash, variant2hash))
        for fut in as_completed(futures):
            ind,h = fut.result()
            individual2hash[ind]=h
    return individual2hash

def run_hashes(boundaries_file: pathlib.Path, clusters_dir: pathlib.Path, chrom: str, out: pathlib.Path,
               variants_file: Optional[pathlib.Path]=None, vcf: Optional[pathlib.Path]=None,
               samples_file: Optional[pathlib.Path]=None, threads: Optional[int]=None,
               gpu: bool=False,gpu_id: Optional[int]=None):

    chr_hash = numpy.binary_repr(int(chrom))
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries)
    out.mkdir(parents=True, exist_ok=True)
    with (out/"haploblock_hashes.tsv").open("w") as f:
        f.write("START\tEND\tHASH\n")
        for (s,e),h in haploblock2hash.items():
            f.write(f"{s}\t{e}\t{h}\n")

    variant2hash=None
    if variants_file:
        samples = data_parser.parse_samples(samples_file) if samples_file else data_parser.parse_samples_from_vcf(vcf)
        variants = data_parser.parse_variants_of_interest(variants_file)
        variant2hash = generate_variant_hashes(variants, vcf, chrom, haploblock_boundaries, samples)
        with (out/"variant_hashes.tsv").open("w") as f:
            f.write("VARIANT\tHASH\n")
            for ind,h in variant2hash.items(): f.write(f"{ind}\t{h}\n")

    for (start,end) in haploblock_boundaries:
        cluster_file = clusters_dir / f"chr{chrom}_{start}-{end}_cluster.tsv"
        individual2cluster, clusters = data_parser.parse_clusters(cluster_file)
        cluster2hash = generate_cluster_hashes(clusters)
        with (out/f"cluster_hashes_{start}-{end}.tsv").open("w") as f:
            f.write("CLUSTER\tHASH\n")
            for ind,h in cluster2hash.items(): f.write(f"{ind}\t{h}\n")

        individual2hash = generate_individual_hashes_parallel(individual2cluster, cluster2hash,
                                                              haploblock2hash, chr_hash, variant2hash,
                                                              max_workers=threads, gpu=gpu,gpu_id=gpu_id)
        with (out/f"individual_hashes_{start}-{end}.tsv").open("w") as f:
            f.write("INDIVIDUAL\tHASH\n")
            for ind,h in individual2hash.items(): f.write(f"{ind}\t{h}\n")

def run(boundaries_file, clusters_dir, chr, out, variants_file=None, vcf=None, samples_file=None, threads=None, gpu=False,gpu_id: Optional[int]=None):
    run_hashes(
        pathlib.Path(boundaries_file),
        pathlib.Path(clusters_dir),
        str(chr),
        pathlib.Path(out),
        pathlib.Path(variants_file) if variants_file else None,
        pathlib.Path(vcf) if vcf else None,
        pathlib.Path(samples_file) if samples_file else None,
        threads,
        gpu=gpu,
        gpu_id=gpu_id
    )

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Generate individual hashes")
    parser.add_argument('--boundaries_file', type=pathlib.Path, required=True)
    parser.add_argument('--clusters_dir', type=pathlib.Path, required=True)
    parser.add_argument('--chr', type=str, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument("--variants", type=pathlib.Path, default=None)
    parser.add_argument("--vcf", type=pathlib.Path, default=None)
    parser.add_argument("--samples", type=pathlib.Path, default=None)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--gpu", action="store_true", help="Use GPU acceleration if available")
    args=parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    try:
        run(
            boundaries_file=args.boundaries_file,
            clusters_dir=args.clusters_dir,
            chr=args.chr,
            out=args.out,
            variants_file=args.variants,
            vcf=args.vcf,
            samples_file=args.samples,
            threads=args.threads,
            gpu=args.gpu
        )
    except Exception as e:
        sys.stderr.write(f"ERROR: {repr(e)}\n")
        sys.exit(1)