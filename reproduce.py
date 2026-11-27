"""Run the reproducibility checks and compare their output with the recorded evidence."""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile


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
    "check_nonlinear_branches",
    "check_nonlinear_readout",
    "check_signed_spectrum",
    "check_hidden_weights",
    "check_character_fourth",
    "probe_nonlinear_moments",
    "probe_pinned_moments",
    "probe_dense_residues",
    "check_projective_coefficient",
    "check_direct_algebra",
    "check_scalar_join",
    "check_direct_fiber",
    "check_direct_projection",
    "check_section_families",
    "check_section_deformation",
    "check_multiscale",
    "check_orbit_parity",
    "check_chord_difference",
    "check_difference_algebra",
    "check_difference_fiber",
    "check_difference_projection",
    "check_chord_sections",
    "check_implicit_chord_algebra",
    "check_implicit_chord_fiber",
    "check_quarter_cosets",
    "check_implicit_chord_projection",
    "check_candidate_targets",
    "check_candidate_offsets",
    "check_candidate_splices",
    "check_candidate_arithmetic",
    "check_candidate_fingerprint",
    "check_batch",
    "check_geometry",
    "check_sign_identity",
    "check_target_jets",
    "probe_target_jets",
    "check_quartic_decoder",
    "check_opposite_targets",
    "check_dense_weighted",
    "check_seed_rank",
    "check_weak_recovery",
    "probe_weak_gate",
)

C38_CHECKS = ("check_strided_compiler", "check_reciprocal", "check_sparse", "check_weighted")
# These reference implementations come from the supplied archives listed in the README.
C38_SOURCES = (
    ("ECDLP_v6_complete.zip", "b213c4b3ffb4ef9301fda5f2ac277a5bbfb70ec95cd35d76f9679b4f37922d80",
     ("reusable-profiled-ecdl-v6/reproducibility/new_code/arithmetic_alternating_chord.py",)),
    ("simplex-alphabet-ecdl-v6.zip", "b4ed0d2ab9bf5dae521e3e63634ccbfee8a85af234b709053c10abfe60901f6c",
     ("simplex-alphabet-ecdl/reproducibility/code/simplex_contraction.py",
      "simplex-alphabet-ecdl/reproducibility/code/native_jet_contraction.py",
      "simplex-alphabet-ecdl/reproducibility/code/weighted_matrix_join.py")),
)


def compare_checks(names, reference_dir=None):
    """Run every check and return the names whose output or exit status differed."""
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    env.pop("PYTHONOPTIMIZE", None)
    if reference_dir is not None:
        env["PYTHONPATH"] = str(reference_dir)
    failures = []
    for name in names:
        run = subprocess.run([sys.executable, "-B", str(ROOT / "checks" / (name+".py"))],
                             cwd=ROOT, capture_output=True, encoding="utf-8", env=env)
        sys.stdout.write(run.stdout)
        if run.returncode:
            sys.stderr.write(run.stderr)
            print(f"{name}: FAILED with exit status {run.returncode}.", flush=True)
            failures.append(name)
            continue
        transcript = ROOT / "evidence" / (name + ".txt")
        try:
            expected = transcript.read_text(encoding="utf-8-sig")
        except (OSError, UnicodeError) as error:
            print(f"{name}: FAILED to read {transcript.name}: {error}", flush=True)
            failures.append(name)
            continue
        if expected.splitlines() != run.stdout.splitlines():
            print(f"{name}: computed output DIFFERED from the recorded evidence.", flush=True)
            failures.append(name)
            continue
        print(f"{name}: matches recorded evidence.", flush=True)
    return failures


def finish(failures, total):
    if failures:
        raise SystemExit(f"{len(failures)} of {total} checks failed or differed: "
                         + ", ".join(failures))
    print(f"All {total} checks matched the recorded evidence.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-archives", type=Path, metavar="DIR",
                        help="also run the four C38 integration checks with the original ZIPs in DIR")
    args = parser.parse_args()
    if args.source_archives is None:
        failures = compare_checks(CHECKS)
        print("Four C38 integration checks are available with --source-archives DIR.")
        finish(failures, len(CHECKS))
        return
    with TemporaryDirectory() as tmp:
        for filename, digest, members in C38_SOURCES:
            path = args.source_archives / filename
            if not path.is_file():
                parser.error(f"Missing source archive: {path}")
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                parser.error(f"Source archive has an unexpected SHA-256 digest: {path}")
            with ZipFile(path) as archive:
                for member in members:
                    (Path(tmp) / Path(member).name).write_bytes(archive.read(member))
        failures = compare_checks(C38_CHECKS, tmp)
    failures += compare_checks(CHECKS)
    finish(failures, len(C38_CHECKS)+len(CHECKS))


if __name__ == "__main__":
    main()
