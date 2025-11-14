#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib
import subprocess
import numpy as np
from multiprocessing import Pool, cpu_count
from typing import Tuple, List, Dict, Optional

from haploblock_pipeline import data_parser

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# Variant counting
# -------------------------------------------------------------------------
def count_variants(vcf: pathlib.Path) -> Tuple[int, int]:
    """
    Count the number of variants in a bgzipped VCF file for each haplotype.
    Returns: (count_hap0, count_hap1)
    """
    result = subprocess.run(
        ["bcftools", "query", "-f", "[ %GT]", str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )

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

    return count_0, count_1


# -------------------------------------------------------------------------
# Consensus FASTA generation
# -------------------------------------------------------------------------
def generate_consensus_fasta(ref_fasta: pathlib.Path,
                             vcf: pathlib.Path,
                             out_dir: pathlib.Path) -> Tuple[pathlib.Path, pathlib.Path]:
    """
    Apply variants from VCF to the reference to create phased consensus FASTA files.
    Uses bcftools consensus (external command).
    """
    tmp_dir = out_dir / "tmp" / "consensus_fasta"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    vcf_name = vcf.stem
    if vcf_name.endswith(".vcf"):
        vcf_name = vcf_name[:-4]

    output_hap0 = tmp_dir / f"{vcf_name}_hap0.fa"
    output_hap1 = tmp_dir / f"{vcf_name}_hap1.fa"

    for hap, outfile in [(1, output_hap0), (2, output_hap1)]:
        with outfile.open("w") as f_out:
            subprocess.run(
                ["bcftools", "consensus", "-H", str(hap), "-f", str(ref_fasta), str(vcf)],
                stdout=f_out,
                check=True,
            )

    return output_hap0, output_hap1


# -------------------------------------------------------------------------
# Variant hash generation
# -------------------------------------------------------------------------
def generate_variant_hashes(variants: List[str],
                            vcf: pathlib.Path,
                            chrom: str,
                            haploblock_boundaries: List[tuple],
                            samples: Optional[List[str]]) -> Dict[str, str]:
    """Generate binary variant presence hashes for all samples and haplotypes."""
    if not samples:
        logger.warning("No samples provided for variant hash generation. Returning empty dict.")
        return {}

    first_pos = variants[0]

    # find region containing first variant
    start = end = None
    for (s, e) in haploblock_boundaries:
        if int(s) <= int(first_pos) <= int(e):
            start, end = s, e
            break

    if start is None:
        raise ValueError(f"Variant {first_pos} not found in any haploblock")

    idx_map = {str(v): i for i, v in enumerate(variants)}
    variant2hash = {
        f"{sample}_chr{chrom}_region_{start}-{end}_hap{h}": ["0"] * len(variants)
        for sample in samples for h in (0, 1)
    }

    region = f"{chrom}:{first_pos}-{end}"
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n",
         "-s", ",".join(samples), "--force-samples", "-r", region, str(vcf)],
        capture_output=True,
        text=True,
        check=True,
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
            if a0 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap0"][i] = "1"
            if a1 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap1"][i] = "1"

    return {k: "".join(v) for k, v in variant2hash.items()}


# -------------------------------------------------------------------------
# Write outputs
# -------------------------------------------------------------------------
def write_tsv_variant_hashes(variant2hash: Dict[str, str], out_dir: pathlib.Path) -> None:
    out = out_dir / "variant_hashes.tsv"
    with out.open("w") as f:
        f.write("INDIVIDUAL\tHASH\n")
        for ind, h in variant2hash.items():
            f.write(f"{ind}\t{h}\n")


def write_tsv_variant_counts(haploblock2count: Dict[tuple, Tuple[float, float]], out_dir: pathlib.Path) -> None:
    out = out_dir / "variant_counts.tsv"
    with out.open("w") as f:
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (s, e), (mean, stdev) in haploblock2count.items():
            f.write(f"{s}\t{e}\t{mean:.3g}\t{stdev:.3g}\n")


# -------------------------------------------------------------------------
# Worker function for Pool
# -------------------------------------------------------------------------
def _process_sample_for_region(args):
    """Worker: extract sample, count variants, generate FASTAs."""
    if args is None:
        return None  # skip empty worker

    sample, region_vcf, region_fasta, out_dir = args
    try:
        sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out_dir)
        c0, c1 = count_variants(sample_vcf)
        generate_consensus_fasta(region_fasta, sample_vcf, out_dir)
        return sample, c0, c1
    except Exception as e:
        logger.exception("Error processing sample %s: %s", sample, e)
        return sample, 0, 0


# -------------------------------------------------------------------------
# Main function (called by pipeline)
# -------------------------------------------------------------------------
def run_phased_sequences(boundaries_file: pathlib.Path,
                         samples_file: Optional[pathlib.Path],
                         vcf: pathlib.Path,
                         ref: pathlib.Path,
                         chr_map: pathlib.Path,
                         chrom: str,
                         variants_file: pathlib.Path,
                         out_dir: pathlib.Path,
                         workers: Optional[int] = None) -> None:

    # Create output paths
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "tmp"
    (tmp_dir / "consensus_fasta").mkdir(parents=True, exist_ok=True)

    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    logger.info("Parsing samples")
    samples = (data_parser.parse_samples(samples_file)
               if samples_file else data_parser.parse_samples_from_vcf(vcf))
    if not samples:
        logger.warning("No samples found; proceeding with empty sample list.")
    logger.info("Found %d samples", len(samples))

    logger.info("Parsing variants of interest")
    variants = data_parser.parse_variants_of_interest(variants_file)

    # Sequential variant-hash generation
    variant2hash = generate_variant_hashes(variants, vcf, chrom, haploblock_boundaries, samples)
    write_tsv_variant_hashes(variant2hash, out_dir)

    # Worker count
    cpu_cores = cpu_count()
    if not workers or workers <= 0:
        workers = cpu_cores
    logger.info("Using %d worker(s) (available cores: %d)", workers, cpu_cores)

    # Per-haploblock processing
    haploblock2count = {}

    for (start, end) in haploblock_boundaries:
        logger.info("Processing haploblock %s-%s", start, end)

        region_vcf = data_parser.extract_region_from_vcf(vcf, chrom, chr_map, start, end, out_dir)
        region_fasta = data_parser.extract_region_from_fasta(ref, chrom, start, end, out_dir)

        work = [(s, region_vcf, region_fasta, out_dir) for s in samples]

        if workers == 1:
            results = [_process_sample_for_region(w) for w in work]
        else:
            with Pool(workers) as pool:
                results = pool.map(_process_sample_for_region, work)

        # filter out None results if no samples
        results = [r for r in results if r is not None]

        counts = [c for (_, c0, c1) in results for c in (c0, c1)]
        haploblock2count[(start, end)] = (
            float(np.mean(counts)) if counts else 0.0,
            float(np.std(counts)) if counts else 0.0,
        )

    write_tsv_variant_counts(haploblock2count, out_dir)


# -------------------------------------------------------------------------
# Pipeline entrypoint
# -------------------------------------------------------------------------
def run(boundaries_file, samples_file, vcf, ref, chr_map, chr, variants, out, threads=None):
    run_phased_sequences(
        pathlib.Path(boundaries_file),
        pathlib.Path(samples_file) if samples_file else None,
        pathlib.Path(vcf),
        pathlib.Path(ref),
        pathlib.Path(chr_map),
        str(chr),
        pathlib.Path(variants),
        pathlib.Path(out),
        threads,
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

    parser = argparse.ArgumentParser(description="Generate haploblock phased sequences and hashes.")
    parser.add_argument("--boundaries_file", type=pathlib.Path, required=True)
    parser.add_argument("--vcf", type=pathlib.Path, required=True)
    parser.add_argument("--ref", type=pathlib.Path, required=True)
    parser.add_argument("--chr_map", type=pathlib.Path, required=True)
    parser.add_argument("--chr", type=str, required=True)
    parser.add_argument("--variants", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--samples_file", type=pathlib.Path, default=None)  # <-- optional now
    parser.add_argument("--workers", type=int, default=None)

    args = parser.parse_args()

    try:
        run_phased_sequences(
            args.boundaries_file,
            args.samples_file,
            args.vcf,
            args.ref,
            args.chr_map,
            args.chr,
            args.variants,
            args.out,
            args.workers,
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {pathlib.Path(sys.argv[0]).name}: {repr(e)}\n")
        sys.exit(1)

