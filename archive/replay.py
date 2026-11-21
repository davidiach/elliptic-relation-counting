"""Replay the historical checks in a temporary copy of their original layout."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile


ARCHIVE = Path(__file__).resolve().parent
CHECKPOINTS = [
    f"hillclimb_{i:03}/" + ("audit_package.py" if i == 4 else "run_all.py")
    for i in range(1, 10)
]
OTHER_CHECKS = (
    "c38_arithmetic_join/verify_join.py", "short_hash_count/verify.py",
    "slot_arity/verify.py", "mixed_coefficient/verify.py", "mixed_coefficient/check_frontier.py",
    "joint_parity/verify.py", "joint_parity/check_fresh_windows.py",
    "joint_parity/check_seed_arithmetic.py", "joint_parity/check_target_incidence.py",
    "joint_parity/check_moment_decoding.py", "joint_parity/probe_planted_moments.py",
    "scalar_norm/verify.py", "implicit_resultant/verify.py", "second_digit/verify.py",
    "section_product/verify.py", "single_target/verify.py", "single_target/covers_check.py",
    "single_target/verify_block_unit.py", "single_target/verify_exception_product.py",
    "single_target/verify_prime_projection.py",
)


def prepare(workspace, scratch):
    research = scratch / "research"
    shutil.copytree(ARCHIVE / "research", research, ignore=shutil.ignore_patterns("runs", "__pycache__"))
    for note in (ARCHIVE.parent / "derivations").glob("*/*.tex"):
        if (research / note.parent.name).is_dir():
            shutil.copyfile(note, research / note.parent.name / note.name)
    first = research / "hillclimb_001"
    manifest = json.loads((first / "source_manifest.json").read_text(encoding="utf-8"))
    (first / "vendor").mkdir()
    for relative, source in manifest["files"].items():
        with ZipFile(workspace / source["archive"]) as z:
            data = z.read(source["member"])
        if hashlib.sha256(data).hexdigest() != source["sha256"]:
            raise ValueError(f"Historical module differs from the recorded source: {relative}")
        (first / relative).write_bytes(data)
    # The package audit verifies its original inputs, including the private rerun records.
    for relative in ("_review/compressed_relation_search_package",
                     "_review/reruns/compressed_relation_search_20261031",
                     "_review/sharper-complete-ecdl-v6-reproducibility/sharper-complete-ecdl-v6/reproducibility/paper"):
        shutil.copytree(workspace / relative, scratch / relative,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copyfile(workspace / "compressed_relation_search_package.zip",
                    scratch / "compressed_relation_search_package.zip")
    return research


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path, help="directory containing the supplied ZIPs and original _review inputs")
    args = parser.parse_args()
    with TemporaryDirectory() as tmp:
        research = prepare(args.workspace.resolve(), Path(tmp))
        env = {**os.environ, "PYTHONUTF8": "1", "PYTHONHASHSEED": "0"}
        env.pop("PYTHONOPTIMIZE", None)
        env.pop("PYTHONPATH", None)
        for script in (*CHECKPOINTS, *OTHER_CHECKS):
            path = research / script
            subprocess.run([sys.executable, "-B", str(path)], cwd=path.parent, env=env, check=True)
            print(f"Passed {script}", flush=True)


if __name__ == "__main__":
    main()
