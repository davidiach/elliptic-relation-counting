# Research ledger, 27 November 2026

## Revisions

The initial checkout used a separate clone and fetched all branch heads and pull request heads. The second fetch during proof verification found two new submissions. The manuscript on the initial main revision is **A polynomial certificate for deterministic elliptic-curve relation systems**. The earlier v4 paper remains a historical comparison.

| Revision | Exact SHA |
|---|---|
| Initial main / PR #4 merge | 7e77dd0fcd222dd39780d4e1d5eb1c464c3b3a87 |
| Initial reviewed PR #7 research commit from the request | 801bee85e89785ef0ef6afde2a46f2a59a7be3bb |
| PR #7 head, including matrix-construction repair | 914dded075b823e6fc4573db87eae0926d918c41 |
| PR #7 merge | ce89d5adff8e712b3bed45a5f1f785970604851a |
| PR #8 head | 75f684bed85e831c1bec709378ca0f5eb1710ee8 |
| Main used for submission / PR #8 merge | 04a8b3d8bb904af8a30aa996f1af4027b52f6dd4 |
| PR #9 open research head inspected | c88236d2416f315012614924ac58a17b4e9124b6 |
| PR #3 head / merge | 89211f7652076a4cde9d89bcb9f9e797750268fb / d155b9468f386f55a41a20a5cc7d9d505a04f769 |
| PR #4 head | 3c47f0969ffc23f1af7274f99b8e58dbc993b453 |
| PR #1 historical branch head | 87077af75cd188ebb2306247bfcac4a0e2facccb |
| PR #2 historical branch head / v4 merge | a4c74d4efa568ddacffa8d4f0d57bb00f89c3e89 / fc3a8019f1763bd15250b34cc9c15cfabda8131c |
| PR #5 branch head | 7511f5cc72918066c57a7e3b395979199f784ef3 |
| PR #6 branch head | 81ec1e0c6c851e3e54b29345e78e28d79a4f6cee |
| Related Long Mathematics manuscript | b11be05092e0f95862d9c7a8aef57e27e9a231c2 |

The author-hosted LZR PDF hash, inaccessible ePrint artifact hash, formal scope, and primary literature records appear in the [audit](INDEPENDENT_AUDIT.md) and [literature review](LITERATURE_REVIEW.md). Research notes inspected included the nine top-level requested files, the fixed-shift Typst manuscript, every note, verifier and transcript in short_relations_adaptive_rank, the PR #7 note and both verifiers, and the additional PR #8 and #9 research.

## Attempts and decisions

| Proposition tested | Proof or falsification | Outcome |
|---|---|---|
| Every target has an identifiable single triple relative to a fixed reference. | A \(C_7\) fixture fails at every nonzero target with its identity reference fixed. On actual sufficiently large Stage I sets, remove one reference per block and construct the corank-one shift subsystem. | False for the chosen baseline and genuine subsystems. The latter family is rationally deficient. |
| A nonaffine score must vary on a prescribed triple fiber, with an identical third tuple available. | Nonconstant finite differences split the admissible support into two levels; Cauchy–Davenport counts mixed centers. Integer interval extremizers check sharpness. | Proved for \(3s>2r+1\) in arbitrary abelian codomains. |
| The prescribed fiber has the exact integral value-zero lattice. | Use its own quotient as the score group and change one coordinate at a time. Check for torsion by mapping the quotient back onto the prime group. | Proved, including common-third generators. |
| Quotient spanning automatically returns an element of the given row space. | Affine spanning permits combinations but can miss zero among individual points. Listing fibers takes \(\Omega(s(2s-r))\) records. | Small witness exists, construction generally open. |
| Class-restricted representations can be found by ordinary Stage I membership calls. | A private lookup may violate the degree limit. Prove uniform restricted density by Fourier mixing, then replace the lookup by a signed finite-order count with literal pins. | Proved at fixed positive joint retention. Construct both endpoints and certify in the supplied rows. |
| First-order trimming should work generally at constant product mass. | Exact \(C_5\) test gives minimum masses \(-267,-6,36\) at odd orders \(1,3,5\). | False for this signed first-order rule. Higher odd orders and a conditional tail estimate are justified. |
| Unweighted flatness protects a single Gram. | Count first and second incidences for \(C_r\setminus\{0\}\) on the prescribed zero three-sum. | The single Gram loses half the raw rank; a vertical stack of common-third cells restores it in the dense proper one-block setting. |
| Inserting a quotient discrepancy might remove the exponential coefficient cut. | Compute each value block's rank exactly. Mixed centers force a nonzero block. | Cut remains exponential for nonaffine scores in the specified encoding. |
| The full integral shift presentation or \(d+1\) bound can be attributed here. | Refetched PR #9 and compared Long Mathematics. | Those general claims have independent prior formulations. Restricted generators and a constructive restricted-label selector are the claims submitted for review. |

## Reproducibility

Run the following commands in this directory. The new verifiers use Python 3.9+ and SymPy. They use synthetic cyclic-group labels, exact integer HNF, rational arithmetic and modular row reduction. They do not implement a cryptographic solver or tensor-contraction oracle.

~~~sh
python3 -B verify.py
python3 -B verify_restrictions.py
tectonic manuscript.tex
~~~

Expected outputs are in [verification.txt](verification.txt) and [restriction_verification.txt](restriction_verification.txt). The common-third HNF was checked independently of the HNF of all same-target contrasts. The restricted-label program recomputes prefix masses using a direct unweighted point-pin counter; a separate exhaustive loop checks the small order-five fixture.

**Inherited verifiers rerun.** The fixed-shift Gram and deterministic-point verifiers completed. All ten programs from short_relations_adaptive_rank completed; eight matched their individual recorded transcripts, and the other two matched their sections of the combined transcript. Both PR #7 verifiers completed. After the second fetch the expanded PR #8 pinning verifier required Python ≥3.10 for strict zip, and passed on Python 3.11. Its integer minor verifier and the PR #9 lattice checks also passed. The test outcomes include 18 integer stacks and the PR #9 rational-basis target obstruction. No asymptotic theorem is inferred from the finite checks. No production curve inputs or third-party cryptographic infrastructure were used.

**Final strengthening.** Re-examining the deletion-resilient seed showed that the two endpoint searches may both omit entire blocks. A single collision comparison certifies the required untrimmed law on an arbitrary subset. On the simultaneous resilient seed sorting class fractions suffices. The padded-shift coordinate-subspace example on these very masks separates this sufficient condition from the all-block joint product condition; it has unpriced basis construction. This refinement still leaves unconcentrated residual labels on all large subsets unresolved.

**Final revision refresh and positive response to PR #10.** Main and PR #9 were unchanged on the refresh. The new open pointwise-target-obstruction head was febe120cc3e69e248a1dd3472f87fa51c9ab8dd2. I checked its basis/offset mechanism, inspected its proof and verifier, and reran its recorded finite checks. Its offsets supply a short integral solution of the inhomogeneous equations while the homogeneous residual direction has diffuse coordinates. This led to the short-lift theorem and a certified response on the same synthetic instance. Run **python3 -B verify_height.py** in this directory; the supplied sparse JSON is the independently selected square basis and public indexed blocks from its \(C_{101}\) dictionary. Each row is checked for exact doubled-tuple realizability. No hidden offset or obstruction direction is included in that fixture. [Output](height_verification.txt) · [Source PR](https://github.com/davidiach/elliptic-relation-counting/pull/10)
