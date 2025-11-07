import sys
import logging
import argparse
import pathlib
import subprocess
import numpy as np
import data_parser
# from multiprocessing import Pool, cpu_count  # optional parallel version

logger = logging.getLogger(__name__)


def count_variants(vcf: pathlib.Path) -> tuple[int, int]:
    """
    Count the number of variants in a bgzipped VCF file.

    Returns:
        (count_hap0, count_hap1)
    """
    result = subprocess.run(
        ["bcftools", "query", "-f", "[ %GT]", vcf],
        capture_output=True,
        text=True,
        check=True
    )

    count_0 = count_1 = 0
    for line in result.stdout.split():
        if line in {".", "./.", ".|."}:
            continue

        if "|" in line:
            left, right = line.split("|")
        elif "/" in line:
            left, right = line.split("/")
        else:
            continue

        count_0 += left == "1"
        count_1 += right == "1"

    return count_0, count_1


def generate_consensus_fasta(ref_fasta: pathlib.Path, vcf: pathlib.Path, out_dir: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
    """
    Apply variants from VCF to the reference to create phased consensus FASTA files.
    """
    tmp_dir = out_dir / "tmp" / "consensus_fasta"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    output_hap0 = tmp_dir / f"{vcf.stem}_hap0.fa"
    output_hap1 = tmp_dir / f"{vcf.stem}_hap1.fa"

    for hap, outfile in [(1, output_hap0), (2, output_hap1)]:
        with outfile.open("w") as f_out:
            subprocess.run(
                ["bcftools", "consensus", "-H", str(hap), "-f", str(ref_fasta), str(vcf)],
                stdout=f_out,
                check=True
            )

    return output_hap0, output_hap1


def generate_variant_hashes(variants: list[str], vcf: pathlib.Path, chrom: str, haploblock_boundaries, samples: list[str]) -> dict[str, str]:
    """
    Generate binary variant presence hashes for all samples and haplotypes.
    """
    first_variant_pos = variants[0]
    start = end = None
    for (bound_start, bound_end) in haploblock_boundaries:
        if int(bound_start) <= int(first_variant_pos) <= int(bound_end):
            start, end = bound_start, bound_end
            break

    if start is None:
        raise ValueError(f"Variant {first_variant_pos} not found in any haploblock")

    variant_indices = {str(pos): i for i, pos in enumerate(variants)}
    variant2hash = {
        f"{sample}_chr{chrom}_region_{start}-{end}_hap{h}": ["0"] * len(variants)
        for sample in samples for h in (0, 1)
    }

    query_region = f"{chrom}:{first_variant_pos}-{end}"
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n",
         "-s", ",".join(samples), "--force-samples", "-r", query_region, str(vcf)],
        capture_output=True, text=True, check=True
    )

    for line in result.stdout.splitlines():
        chrom, pos, *genotypes = line.split("\t")
        if pos not in variant_indices:
            continue
        idx = variant_indices[pos]

        for sample_idx, genotype in enumerate(genotypes):
            if "|" not in genotype:
                continue
            hap0, hap1 = genotype.split("|")
            sample = samples[sample_idx]
            if hap0 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap0"][idx] = "1"
            if hap1 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap1"][idx] = "1"

    # convert hash lists to strings
    return {k: "".join(v) for k, v in variant2hash.items()}

    # -------------------------------------------------------------------------
    # 💡 Perspective: Parallelize bcftools queries across haploblocks
    # Using Pool to distribute haploblocks or samples.
    #
    # def process_region(region):
    #     return generate_variant_hashes_for_region(...)
    #
    # with Pool(cpu_count()) as pool:
    #     results = pool.map(process_region, haploblock_boundaries)
    # -------------------------------------------------------------------------


def write_tsv_variant_hashes(variant2hash: dict[str, str], out_dir: pathlib.Path) -> None:
    """Write variant hashes to TSV."""
    out_path = out_dir / "variant_hashes.tsv"
    with out_path.open("w") as f:
        f.write("INDIVIDUAL\tHASH\n")
        f.writelines(f"{ind}\t{hash_}\n" for ind, hash_ in variant2hash.items())


def write_tsv_variant_counts(haploblock2count: dict[tuple[int, int], list[float]], out_dir: pathlib.Path) -> None:
    """Write mean and stdev variant counts to TSV."""
    out_path = out_dir / "variant_counts.tsv"
    with out_path.open("w") as f:
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (start, end), (mean, stdev) in haploblock2count.items():
            f.write(f"{start}\t{end}\t{mean:.3g}\t{stdev:.3g}\n")


def main(boundaries_file: pathlib.Path, samples_file: pathlib.Path | None,
         vcf: pathlib.Path, ref: pathlib.Path, chr_map: pathlib.Path,
         chrom: str, variants_file: pathlib.Path, out_dir: pathlib.Path) -> None:

    # Ensure output dirs exist
    tmp_dir = out_dir / "tmp"
    if tmp_dir.exists():
        logger.error(f"Output directory {out_dir} already exists; remove it first.")
        raise FileExistsError(f"{out_dir} exists")
    (tmp_dir / "consensus_fasta").mkdir(parents=True, exist_ok=True)

    logger.info("Parsing haploblock boundaries")
    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    logger.info("Parsing samples")
    samples = (data_parser.parse_samples(samples_file)
               if samples_file else data_parser.parse_samples_from_vcf(vcf))
    logger.info("Found %d samples", len(samples))

    logger.info("Parsing variants of interest")
    variants = data_parser.parse_variants_of_interest(variants_file)

    variant2hash = generate_variant_hashes(variants, vcf, chrom, haploblock_boundaries, samples)
    write_tsv_variant_hashes(variant2hash, out_dir)

    haploblock2count = {}
    for (start, end) in haploblock_boundaries:
        logger.info("Processing haploblock %s-%s", start, end)
        region_vcf = data_parser.extract_region_from_vcf(vcf, chrom, chr_map, start, end, out_dir)
        region_fasta = data_parser.extract_region_from_fasta(ref, chrom, start, end, out_dir)

        counts = []
        for sample in samples:
            sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out_dir)
            count_0, count_1 = count_variants(sample_vcf)
            counts.extend([count_0, count_1])
            generate_consensus_fasta(region_fasta, sample_vcf, out_dir)

        haploblock2count[(start, end)] = [np.mean(counts), np.std(counts)]

    write_tsv_variant_counts(haploblock2count, out_dir)


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
    parser.add_argument("--samples_file", type=pathlib.Path, required=False)

    args = parser.parse_args()
    try:
        main(args.boundaries_file, args.samples_file, args.vcf, args.ref,
             args.chr_map, args.chr, args.variants, args.out)
    except Exception as e:
        sys.stderr.write(f"ERROR in {pathlib.Path(sys.argv[0]).name}: {repr(e)}\n")
        sys.exit(1)

