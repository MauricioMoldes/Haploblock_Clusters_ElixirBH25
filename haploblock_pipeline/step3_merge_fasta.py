#!/usr/bin/env python3
"""
Fast, deterministic merge of per-sample per-haplotype FASTAs into one FASTA per region.

Replaces the old bash pipeline. Behavior:
 - Input: directory with files like:
     HG00096_HG00096_chr6_region_31480875-31598421.vcf.norm_flt_hap0.fa
 - Output: one FASTA per region: <region>.fa (e.g. chr6_region_31480875-31598421.fa)
 - Order inside merged FASTA: lexicographic by sample name then hap (hap0 before hap1)
 - Streams files; does not load entire files into memory.
 - Parallel across regions using ProcessPoolExecutor.
"""

from __future__ import annotations
import os
import sys
import re
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple

# ---------- logging ----------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger("step3_merge_fasta")

# ---------- filename parsing ----------
# Example filename:
# HG00096_HG00096_chr6_region_31480875-31598421.vcf.norm_flt_hap0.fa
FNAME_RE = re.compile(
    r"""^(?P<sample>[^_]+)_[^_]*_           # sample prefix (first field); second field may duplicate sample
        (?P<region>chr[0-9XYM]+_region_[0-9]+-[0-9]+)  # region token
        .*_hap(?P<hap>[01])                 # hap number (0 or 1)
        \.(?:fa|fasta)(?:\.gz)?$            # extension
    """,
    re.VERBOSE
)

# ---------- helpers ----------
def discover_region_files(input_dir: Path) -> Dict[str, List[Path]]:
    """
    Scan input_dir for fasta files and group them by region token.
    Returns dict: region -> list of Path
    """
    files_by_region: Dict[str, List[Path]] = {}
    for p in sorted(input_dir.iterdir()):
        if not p.is_file():
            continue
        if not (p.suffix.lower() in {".fa", ".fasta"} or p.name.endswith(".fa.gz") or p.name.endswith(".fasta.gz")):
            continue
        m = FNAME_RE.match(p.name)
        if not m:
            # fallback: attempt to find region token inside filename
            region_search = re.search(r"(chr[0-9XYM]+_region_[0-9]+-[0-9]+)", p.name)
            region = region_search.group(1) if region_search else "unknown_region"
        else:
            region = m.group("region")
        files_by_region.setdefault(region, []).append(p)
    return files_by_region

def canonical_key_for_fname(p: Path) -> Tuple[str,int]:
    """
    Return a sorting key for deterministic ordering:
    (sample_name, hap)
    If we can't parse the sample, use full filename as sample.
    """
    m = FNAME_RE.match(p.name)
    if m:
        sample = m.group("sample")
        hap = int(m.group("hap"))
        return (sample, hap)
    # fallback: try to extract hap number
    hap_search = re.search(r"_hap([01])", p.name)
    hap = int(hap_search.group(1)) if hap_search else 0
    # fallback sample: filename without suffix
    sample = p.name.rsplit(".", 1)[0]
    return (sample, hap)

def merge_region(region: str, files: List[Path], output_dir: Path, buffer_size: int = 1 << 20) -> Path:
    """
    Merge the list of FASTA files into a single output file named <region>.fa in output_dir.
    Writes atomically (tmp file then rename). Returns path to output file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{region}.fa"
    tmp_path = output_dir / f".{region}.fa.tmp"

    # sort files deterministically by sample then hap
    files_sorted = sorted(files, key=canonical_key_for_fname)

    with tmp_path.open("w", buffering=buffer_size) as out_f:
        for p in files_sorted:
            m = FNAME_RE.match(p.name)
            if m:
                sample = m.group("sample")
                hap = m.group("hap")
                header = f"{sample}_{region}_hap{hap}"
            else:
                # fallback header: filename without extension
                header = p.name.rsplit(".", 1)[0]

            out_f.write(f">{header}\n")

            # stream sequence lines skipping any existing headers
            # support gz compressed input if needed
            if p.suffixes and p.suffixes[-1] == ".gz":
                import gzip
                opener = lambda path: gzip.open(path, "rt")
            else:
                opener = lambda path: open(path, "r")

            with opener(p) as in_f:
                for line in in_f:
                    if not line:
                        continue
                    if line.startswith(">"):
                        continue
                    out_f.write(line.rstrip("\n") + "\n")

    # atomic replace
    tmp_path.replace(out_path)
    try:
        os.chmod(out_path, 0o644)
    except Exception:
        pass
    return out_path

# ---------- main ----------
def main(input_dir: str, output_dir: str, threads: int = None, clean: bool = False):
    input_dir_p = Path(input_dir)
    output_dir_p = Path(output_dir)

    if not input_dir_p.exists():
        logger.error("Input directory does not exist: %s", input_dir)
        sys.exit(1)

    if clean and output_dir_p.exists():
        logger.info("Cleaning output directory %s", output_dir)
        for f in output_dir_p.iterdir():
            try:
                f.unlink()
            except Exception:
                pass

    output_dir_p.mkdir(parents=True, exist_ok=True)

    # determine worker count
    if threads is None or threads <= 0:
        threads = max(1, (os.cpu_count() or 2) - 1)

    logger.info("Scanning input directory: %s", input_dir_p)
    files_by_region = discover_region_files(input_dir_p)
    n_regions = len(files_by_region)
    if n_regions == 0:
        logger.warning("No FASTA files found in %s", input_dir_p)
        return

    logger.info("Found %d regions to merge (threads=%d)", n_regions, threads)

    # Create jobs: list of (region, filelist)
    jobs = [(region, files) for region, files in sorted(files_by_region.items())]

    # Parallel merge
    with ProcessPoolExecutor(max_workers=threads) as ex:
        future_to_region = {
            ex.submit(merge_region, region, files, output_dir_p): region for region, files in jobs
        }
        for fut in as_completed(future_to_region):
            region = future_to_region[fut]
            try:
                out_path = fut.result()
                logger.info("Merged region %s -> %s", region, out_path)
            except Exception as e:
                logger.exception("Failed to merge region %s: %s", region, e)

    logger.info("All merges completed. Output dir: %s", output_dir_p)


def run(input_dir, output_dir, threads=None, gpu=False):
    """
    Entry point used by main.py for Step 3.
    Matches the call signature exactly.
    """
    # We do not use GPU in step3, but we accept it to avoid errors
    logger.info(f"[step3] GPU flag received (ignored): {gpu}")

    # Ensure paths are strings for our main()
    input_dir = str(input_dir)
    output_dir = str(output_dir)

    logger.info(f"[step3] input_dir={input_dir}")
    logger.info(f"[step3] output_dir={output_dir}")
    logger.info(f"[step3] threads={threads}")

    # Call the actual merging logic
    main(input_dir, output_dir, threads=threads, clean=False)


# ---------- CLI ----------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Merge per-sample haplotype FASTAs into one FASTA per region")
    ap.add_argument("input_dir", help="Directory with per-sample FASTAs (step2 output)")
    ap.add_argument("output_dir", help="Directory to write merged FASTAs")
    ap.add_argument("threads", nargs="?", type=int, default=None, help="Number of parallel workers (optional)")
    ap.add_argument("--clean", action="store_true", help="Remove output files before merging")
    args = ap.parse_args()
    main(args.input_dir, args.output_dir, threads=args.threads, clean=args.clean)

