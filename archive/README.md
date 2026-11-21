# Historical reports and checks

| Date | Report | Scope |
|---|---|---|
| 31 October 2026 | [Validated reuse and coordinate-dependent tensor contraction (PDF)](2026-10-31-research-synthesis.pdf) | Version 6.0 synthesis of three supplied research packages, with their assumptions and finite reproducibility results. |
| 31 October 2026 | [Compressed relation-search audit](research/hillclimb_004/README.txt) · [Recorded audit result](research/hillclimb_004/runs/20261031T180405.300821Z/independent_audit.json) | Review of the pilot-reuse argument, unit selectors and spectral join, with independent finite checks. |
| 17 November 2026 | [Progress toward the square-root ECDLP goal (PDF)](2026-11-17-research-progress.pdf) | A dated account of the conditional bounds and open problems. |

The [research directory](research/) preserves the earlier checking programs, checkpoint reports, and the latest successful run records from each of the nine October checkpoints. Its filenames and source hashes follow the original workspace layout. The [v2 paper](../paper.pdf) and [derivations index](../derivations/README.md) give the current results and their hypotheses.

## Replay the historical checks

The [replay script](replay.py) runs the nine checkpoint entry points and twenty additional historical checks in a temporary directory. It restores the published TeX notes alongside the original code and verifies the hashes of the five reference modules against the [source manifest](research/hillclimb_001/source_manifest.json). Earlier scripts write new run records inside that temporary copy.

Install the dependencies in the [main README](../README.md#reproduce-the-experiments). From the repository root, pass the path to the original research workspace.

```sh
python archive/replay.py /path/to/original/workspace
```

The replay uses the three ZIPs named in the source manifest, `compressed_relation_search_package.zip`, and the original `_review` inputs for the package audit and spectrum check. These inputs include the compressed-package review manifest, its independently rerun records under `reruns/compressed_relation_search_20261031`, and the sharper-package paper sources. For the current paper and supplement, use [the release runner](../README.md#reproduce-the-experiments).
