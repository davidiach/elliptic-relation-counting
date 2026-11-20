"""Run each paper check and compare its output with the recorded evidence."""
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
CHECKS = (
    "check_accounting",
    "check_pooled_zero",
    "probe_pooled_zero",
    "verify_span",
    "check_strided_isolation",
    "check_block_support",
    "check_seed_cuts",
    "check_grid_values",
    "check_address_parity",
    "probe_signed_moment",
)


def main():
    for name in CHECKS:
        run = subprocess.run([sys.executable, "-B", str(ROOT / "checks" / (name+".py"))],
                             cwd=ROOT, capture_output=True, text=True, check=False)
        sys.stdout.write(run.stdout)
        if run.returncode:
            sys.stderr.write(run.stderr)
            raise SystemExit(run.returncode)
        expected = (ROOT / "evidence" / (name+".txt")).read_text(encoding="utf-8-sig")
        if expected.splitlines() != run.stdout.splitlines():
            raise SystemExit(f"{name}: computed output differed from the recorded evidence.")
        print(f"{name}: matches recorded evidence.", flush=True)


if __name__ == "__main__":
    main()
