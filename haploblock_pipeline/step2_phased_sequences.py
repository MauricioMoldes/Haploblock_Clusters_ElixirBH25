#!/usr/bin/env python3
import sys
import os
import logging
import argparse
import pathlib
import subprocess
from multiprocessing import Pool, cpu_count
from typing import Tuple, Dict, Optional

import numpy

import data_parser

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# Variant counting
# -------------------------------------------------------------------------
def count_variants(vcf: pathlib.Path) -> Tuple[int, int]:
    """
    Count the number of variants in the VCF file for each haplotype separately.
    Returns:
    - count_hap0: int
    - count_hap1: int
    """
    result = subprocess.run(
        ["bcftools", "query", "-f", "[ %GT]", str(vcf)],
        capture_output=True,
        text=True,
        check=True)

    count_0 = count_1 = 0
    for token in result.stdout.split():
        if token in {".", "./.", ".|."}:
            continue

        if "|" in token:
            left, right = token.split("|")
        elif "/" in token:
            left, right = token.split("/")
        else:
            continue

        count_0 += left == "1"
        count_1 += right == "1"

    return(count_0, count_1)


# -------------------------------------------------------------------------
# VCF processing
# -------------------------------------------------------------------------
def normalize_vcf(vcf: pathlib.Path, ref: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    """
    Normalize indels

    returns:
    - normalized_vcf: pathlib.Path
    """
    vcf_name = vcf.stem
    if vcf_name.endswith(".vcf"):
        vcf_name = vcf_name[:-4]

    normalized_vcf = out_dir / "tmp" / f"{vcf_name}.norm.vcf.gz"
    bcftools_log_file = out_dir / "tmp" / "log_bcftools.txt"  # file for BCFtools log messages

    subprocess.run(["bcftools", "norm", "-f", str(ref), str(vcf), "-Ob", "-o", normalized_vcf],
                    stdout=subprocess.PIPE,
                    stderr=open(bcftools_log_file, "a+"),
                    text=True,
                    check=True)
    
    # index normalized VCF
    subprocess.run(["bcftools", "index", normalized_vcf],
                    check=True)
    
    return(normalized_vcf)


def filter_vcf(vcf: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    """
    Filter adjacent indels within 5bp and variants with quality < 40

    returns:
    - filtered_vcf: pathlib.Path
    """
    vcf_name = vcf.stem
    if vcf_name.endswith(".vcf"):
        vcf_name = vcf_name[:-4]

    filtered_vcf = out_dir / "tmp" / f"{vcf_name}.flt.vcf.gz"
    bcftools_log_file = out_dir / "tmp" / "log_bcftools.txt"  # file for BCFtools log messages
    
    subprocess.run(["bcftools", "filter", "--IndelGap", "5", "-e", "QUAL<40", str(vcf), "-Ob", "-o", filtered_vcf],
                    stdout=subprocess.PIPE,
                    stderr=open(bcftools_log_file, "a+"),
                    text=True,
                    check=True)
    
    # index filtered VCF
    subprocess.run(["bcftools", "index", filtered_vcf],
                    check=True)

    return(filtered_vcf)


# -------------------------------------------------------------------------
# Consensus FASTA generation
# -------------------------------------------------------------------------
def generate_consensus_fasta(ref: pathlib.Path,
                             vcf: pathlib.Path,
                             out_dir: pathlib.Path) -> Tuple[pathlib.Path, pathlib.Path]:
    """
    Apply variants from VCF to the reference, create phased consensus FASTA files.
    Use bcftools consensus, as described here:
    https://samtools.github.io/bcftools/howtos/consensus-sequence.html
    """
    # Create temporary directory for phased sequences if it doesn't exist
    tmp_dir = out_dir / "tmp" / "consensus_fasta"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    vcf_name = vcf.stem
    if vcf_name.endswith(".vcf"):
        vcf_name = vcf_name[:-4]
    if vcf_name.endswith(".norm.flt"):
        vcf_name = vcf_name[:-9]

    output_hap0 = tmp_dir / f"{vcf_name}_hap0.fa"
    output_hap1 = tmp_dir / f"{vcf_name}_hap1.fa"

    bcftools_log_file = os.path.join(out_dir, "tmp", "log_bcftools.txt")  # file for BCFtools log messages

    for (hap, out_file) in [(1, output_hap0), (2, output_hap1)]:
        subprocess.run(
                ["bcftools", "consensus", "-H", str(hap), "-f", str(ref), str(vcf)],
                stdout=out_file.open("w"),
                stderr=open(bcftools_log_file, "a+"),
                check=True)

    return(output_hap0, output_hap1)


# -------------------------------------------------------------------------
# Write outputs
# -------------------------------------------------------------------------
def variant_counts_to_TSV(haploblock2count: Dict[tuple, Tuple[float, float]], out_dir: pathlib.Path):
    out_file = out_dir / "variant_counts.tsv"
    with open(out_file, "w") as f:
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (s, e), (mean, stdev) in haploblock2count.items():
            f.write(f"{s}\t{e}\t{mean:.3g}\t{stdev:.3g}\n")


# -------------------------------------------------------------------------
# Worker function for Pool
# -------------------------------------------------------------------------

def _process_sample_for_region(args, gpu=False, gpu_id=0):
    """Worker: extract sample, normalize, filter, count variants, generate FASTAs."""
    if args is None:
        return None

    sample, region_vcf, region_fasta, ref, out_dir = args
    logger.info("Processing sample %s on %s", sample, "GPU" if gpu else "CPU")

    try:
        sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out_dir)
        normalized_vcf = normalize_vcf(sample_vcf, ref, out_dir)
        filtered_normalized_vcf = filter_vcf(normalized_vcf, out_dir)

        # ----------------------------
        # Variant counting (GPU or CPU)
        # ----------------------------
        if gpu:
            import cupy as cp
            try:
                # Load VCF GT field
                result = subprocess.run(
                    ["bcftools", "query", "-f", "[ %GT]", str(filtered_normalized_vcf)],
                    capture_output=True, text=True, check=True
                )
                gt_tokens = [t for t in result.stdout.split() if t not in {".", "./.", ".|."}]

                # Convert to GPU array
                gt_array = cp.array(gt_tokens, dtype=object)

                # Separate haplotypes (assume phased "|")
                hap0 = cp.array([int(t.split("|")[0]) for t in gt_array])
                hap1 = cp.array([int(t.split("|")[1]) for t in gt_array])

                count_0 = int(cp.sum(hap0 == 1))
                count_1 = int(cp.sum(hap1 == 1))

            except Exception as e:
                logger.warning("GPU variant counting failed: %s. Falling back to CPU.", e)
                count_0, count_1 = count_variants(filtered_normalized_vcf)
        else:
            count_0, count_1 = count_variants(filtered_normalized_vcf)

        # ----------------------------
        # Consensus FASTA (CPU for bcftools)
        # ----------------------------
        generate_consensus_fasta(region_fasta, filtered_normalized_vcf, out_dir)

        return sample, count_0, count_1

    except Exception as e:
        logger.exception("Error processing sample %s: %s", sample, e)
        return sample, 0, 0


# -------------------------------------------------------------------------
# Main function (called by pipeline)
# -------------------------------------------------------------------------
def run_phased_sequences(boundaries_file: pathlib.Path,
                         vcf: pathlib.Path,
                         ref: pathlib.Path,
                         chr_map: pathlib.Path,
                         chrom: str,
                         out_dir: pathlib.Path,
                         samples_file: Optional[pathlib.Path] = None,
                         workers: Optional[int] = None,
                         gpu: bool = False,
                         gpu_id: int = 0):
    """Generate haploblock phased sequences (CPU + GPU)."""
    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    samples = (data_parser.parse_samples(samples_file)
               if samples_file else data_parser.parse_samples_from_vcf(vcf))
    if not samples:
        logger.warning("No samples found; proceeding with empty sample list.")

    cpu_cores = cpu_count()
    if not workers or workers <= 0:
        workers = cpu_cores
    logger.info("Using %d worker(s) (available cores: %d)", workers, cpu_cores)
    logger.info("GPU mode: %s (gpu_id=%d)", gpu, gpu_id)

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "tmp"
    (tmp_dir / "consensus_fasta").mkdir(parents=True, exist_ok=True)

    haploblock2count = {}
    for start, end in haploblock_boundaries:
        logger.info("Processing haploblock %s-%s", start, end)

        region_vcf = data_parser.extract_region_from_vcf(vcf, chrom, chr_map, start, end, out_dir)
        region_fasta = data_parser.extract_region_from_fasta(ref, chrom, start, end, out_dir)

        work = [(s, region_vcf, region_fasta, ref, out_dir) for s in samples]
        if workers == 1:
            results = [_process_sample_for_region(w, gpu=gpu, gpu_id=gpu_id) for w in work]
        else:
            from functools import partial
            with Pool(workers) as pool:
                results = pool.map(partial(_process_sample_for_region, gpu=gpu, gpu_id=gpu_id), work)

        results = [r for r in results if r is not None]
        counts = [c for (_, c0, c1) in results for c in (c0, c1)]
        haploblock2count[(start, end)] = (
            float(numpy.mean(counts)) if counts else 0.0,
            float(numpy.std(counts)) if counts else 0.0
        )

    variant_counts_to_TSV(haploblock2count, out_dir)
    logger.info(f"Variant counts written to {out_dir}")



# -------------------------------------------------------------------------
# Pipeline entrypoint
# -------------------------------------------------------------------------
def run(boundaries_file, vcf, ref, chr_map, chr, out, samples_file, threads=None, gpu=False, gpu_id=0):
    run_phased_sequences(
        pathlib.Path(boundaries_file),
        pathlib.Path(vcf),
        pathlib.Path(ref),
        pathlib.Path(chr_map),
        str(chr),
        pathlib.Path(out),
        pathlib.Path(samples_file) if samples_file else None,
        threads,
        gpu=gpu,
        gpu_id=gpu_id
    )


# -------------------------------------------------------------------------
# CLI for standalone execution
# -------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG,
    )
    logger = logging.getLogger(pathlib.Path(sys.argv[0]).name)

    parser = argparse.ArgumentParser(description="Generate haploblock phased sequences")
    parser.add_argument("--boundaries_file", type=pathlib.Path, required=True,
                        help="TSV file with header (START\tEND) and 2 columns: start end")
    parser.add_argument("--vcf", type=pathlib.Path, required=True,
                        help="Phased VCF file")
    parser.add_argument("--ref", type=pathlib.Path, required=True,
                        help="FASTA file with reference sequence")
    parser.add_argument("--chr_map", type=pathlib.Path, required=True,
                        help="File with one chromosome number-to-name mapping per line (e.g., '6 chr6')")
    parser.add_argument("--chr", type=str, required=True,
                        help="Chromosome number")
    parser.add_argument("--out", type=pathlib.Path, required=True,
                        help="Output folder path")
    parser.add_argument("--samples", type=pathlib.Path, default=None,
                        help="TSV file with samples from 1000Genomes (optional)")
    parser.add_argument("--workers", type=int, default=None,
                        help="TODO")

    args = parser.parse_args()

    try:
        run_phased_sequences(
            args.boundaries_file,
            args.vcf,
            args.ref,
            args.chr_map,
            args.chr,
            args.out,
            args.samples,
            args.workers,
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {pathlib.Path(sys.argv[0]).name}: {repr(e)}\n")
        sys.exit(1)

