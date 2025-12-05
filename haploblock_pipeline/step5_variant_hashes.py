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

# CPU AND GPU sorted 
"""
def generate_cluster_hashes(clusters: List[str]) -> dict:
    # Ensure deterministic ordering
    clusters_sorted = sorted(clusters)

    return {
        c: numpy.binary_repr(i, width=CLUSTER_HASH_LENGTH)
        for i, c in enumerate(clusters_sorted)
    }
"""
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
    """
    Build individual hashes in parallel; GPU-optimized batch implementation when gpu=True.

    Returns: dict {individual_name: bitstring}
    """
    # Normalize haploblock keys to string pairs, to match how keys are looked up
    haploblock2hash = { (str(s), str(e)) if not isinstance(k, tuple) else (str(k[0]), str(k[1])): v
                        for k, v in (haploblock2hash.items() if hasattr(haploblock2hash, "items") else haploblock2hash) }

    # Quick return on empty
    individuals = list(individual2cluster.keys())
    N = len(individuals)
    if N == 0:
        return {}

    # Decide CPU path early if not using GPU
    if not gpu:
        # CPU threaded fallback (keeps previous generate_individual_hash semantics)
        max_workers = max_workers or (os.cpu_count()-1 or 1)
        individual2hash = {}
        futures = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for ind in individuals:
                futures.append(
                    executor.submit(generate_individual_hash, ind, individual2cluster, cluster2hash,
                                    haploblock2hash, chr_hash, variant2hash)
                )
            for fut in as_completed(futures):
                ind, h = fut.result()
                individual2hash[ind] = h
        return individual2hash

    # ---- GPU path (batched) ----
    try:
        import cupy as cp
        # set device if provided
        if gpu_id is not None:
            try:
                cp.cuda.Device(gpu_id).use()
            except Exception as e:
                logger.warning("Could not set CuPy device %s: %s (using default)", gpu_id, e)

        # Convert chr_hash to bit array
        chr_bits = [int(b) for b in str(chr_hash)]
        chr_len = len(chr_bits)
        chr_arr = cp.array(chr_bits, dtype=cp.uint8)  # shape (chr_len,)

        # Prepare haploblock hashes map -> cp arrays
        # Ensure all hap hashes have same length
        sample_hap = next(iter(haploblock2hash.values()))
        hap_len = len(sample_hap)
        hap_hash_map = {}
        for k, v in haploblock2hash.items():
            if len(v) != hap_len:
                raise ValueError("Inconsistent haploblock hash lengths")
            hap_hash_map[k] = cp.array([int(b) for b in v], dtype=cp.uint8)

        # Prepare cluster hash map
        sample_cl = next(iter(cluster2hash.values())) if len(cluster2hash) > 0 else "0"*CLUSTER_HASH_LENGTH
        clus_len = len(sample_cl)
        clus_hash_map = {k: cp.array([int(b) for b in v], dtype=cp.uint8) for k, v in cluster2hash.items()}

        # Variant hashes (may be None). Determine variant length if present
        var_len = 0
        var_hash_map = {}
        if variant2hash:
            # variant2hash keys may not map 1:1 to individuals; treat missing per-individual as zeros
            sample_var = next(iter(variant2hash.values()))
            var_len = len(sample_var)
            for k, v in variant2hash.items():
                if len(v) != var_len:
                    raise ValueError("Inconsistent variant hash lengths")
                var_hash_map[k] = cp.array([int(b) for b in v], dtype=cp.uint8)

        # total length and allocation
        strand_len = 4
        total_len = strand_len + chr_len + hap_len + clus_len + var_len

        # allocate full array on GPU
        full_gpu = cp.zeros((N, total_len), dtype=cp.uint8)

        # Build strand array (N x 4)
        # rule from original: strand = individual[-1]; strand_hash = "0001" if strand=="0" else "0010"
        strand_arr = cp.zeros((N, strand_len), dtype=cp.uint8)
        # fill by indices
        strand0_idx = [i for i, ind in enumerate(individuals) if ind[-1] == "0"]
        strand1_idx = [i for i, ind in enumerate(individuals) if ind[-1] != "0"]
        if len(strand0_idx) > 0:
            strand_arr[cp.array(strand0_idx, dtype=cp.int64), :] = cp.array([0,0,0,1], dtype=cp.uint8)
        if len(strand1_idx) > 0:
            strand_arr[cp.array(strand1_idx, dtype=cp.int64), :] = cp.array([0,0,1,0], dtype=cp.uint8)

        # place strand
        off = 0
        full_gpu[:, off:off+strand_len] = strand_arr
        off += strand_len

        # place chr hash repeated for all individuals
        # tile chr_arr to (N, chr_len)
        full_gpu[:, off:off+chr_len] = cp.tile(chr_arr[cp.newaxis, :], (N, 1))
        off += chr_len

        # Build hap hash per individual: need hap key for each individual
        # parse region like original: region_str = individual.split("_")[3] ...
        hap_matrix = cp.zeros((N, hap_len), dtype=cp.uint8)
        for i, ind in enumerate(individuals):
            try:
                region_str = ind.split("_")[3].replace(".fa","").replace(".fasta","").replace(".vcf","")
                s, e = region_str.split("-")
                key = (str(s), str(e))
                if key not in hap_hash_map:
                    # attempt reverse-int keys maybe present as ints
                    key_alt = (str(int(s)), str(int(e)))
                    if key_alt in hap_hash_map:
                        key = key_alt
                    else:
                        raise KeyError(key)
                hap_matrix[i, :] = hap_hash_map[key]
            except Exception as ex:
                raise RuntimeError(f"Cannot determine haploblock hash for individual '{ind}': {ex}")

        full_gpu[:, off:off+hap_len] = hap_matrix
        off += hap_len

        # Build cluster hash per individual
        clus_matrix = cp.zeros((N, clus_len), dtype=cp.uint8)
        for i, ind in enumerate(individuals):
            clus_name = individual2cluster[ind]
            if clus_name not in clus_hash_map:
                raise KeyError(f"Cluster '{clus_name}' missing from cluster2hash")
            clus_matrix[i, :] = clus_hash_map[clus_name]

        full_gpu[:, off:off+clus_len] = clus_matrix
        off += clus_len

        # Variant hash: fill from variant2hash map if provided, otherwise leave zeros
        if var_len > 0:
            var_matrix = cp.zeros((N, var_len), dtype=cp.uint8)
            # variant2hash keys were previously formatted as f"{sample}_chr{chrom}_region_{start}-{end}_hap{h}"
            # we expect variant2hash keys to match individual names, otherwise try match by prefix
            for i, ind in enumerate(individuals):
                # prefer exact key match
                if ind in var_hash_map:
                    var_matrix[i, :] = var_hash_map[ind]
                else:
                    # fallback: try to find any variant key that startswith the sample id
                    # sample id is first token before first '_'
                    sample = ind.split("_")[0]
                    matched = None
                    for k in var_hash_map:
                        if k.startswith(sample + "_"):
                            matched = k
                            break
                    if matched:
                        var_matrix[i, :] = var_hash_map[matched]
                    else:
                        # leave zeros if no variant vector found for individual
                        pass
            full_gpu[:, off:off+var_len] = var_matrix
            off += var_len

        # Sanity: off should equal total_len
        if off != total_len:
            raise RuntimeError("Internal error: computed offsets mismatch")

        # Transfer to host once
        full_numpy = cp.asnumpy(full_gpu)  # shape (N, total_len), dtype uint8

        # Convert each row to bitstring (list comprehension is reasonably fast)
        individual2hash = {}
        for i, ind in enumerate(individuals):
            row = full_numpy[i]
            # faster join: map to characters and join
            s = "".join(map(str, row.tolist()))
            individual2hash[ind] = s

        return individual2hash

    except Exception as e:
        logger.warning("GPU batch hashing failed (%s). Falling back to CPU threaded implementation.", e)
        # Fall back to CPU threaded implementation
        max_workers = max_workers or (os.cpu_count()-1 or 1)
        individual2hash = {}
        futures = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for ind in individuals:
                futures.append(
                    executor.submit(generate_individual_hash, ind, individual2cluster, cluster2hash,
                                    haploblock2hash, chr_hash, variant2hash)
                )
            for fut in as_completed(futures):
                ind, h = fut.result()
                individual2hash[ind] = h
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