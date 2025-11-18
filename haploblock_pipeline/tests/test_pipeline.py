import subprocess
import hashlib
from pathlib import Path


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def test_pipeline_md5():
    project_root = Path(__file__).resolve().parents[1]

    # Run the pipeline
    cmd = [
        "python3",
        str(project_root / "main.py"),
        "--config",
        str(project_root / "config" / "default.yaml"),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    print(result.stdout)
    print(result.stderr)
    assert result.returncode == 0, "Pipeline did not finish successfully"

    # ---- Define expected outputs ----
    expected = {
        "out_dir/haploblock_boundaries_chr6.tsv":
            "070267da8da74793d7758456de11ced2",

        "out_dir/haploblock_hashes_chr6.tsv":
            "7f85da31ec608a86879318de94b1e314",

        "out_dir/TNFa/haploblock_phased_seq_merged/chr6_region_31480875-31598421.fa":
            "9f54e03acfb35f54e2a9fc3b6b197ad6",        
    }

    # ---- Validate each file MD5 ----
    for rel_path, expected_md5 in expected.items():
        file_path = project_root / rel_path
        assert file_path.exists(), f"Output file missing: {file_path}"

        computed = md5sum(file_path)
        assert computed == expected_md5, (
            f"MD5 mismatch for {rel_path}:\n"
            f"  expected: {expected_md5}\n"
            f"  got:      {computed}"
        )

