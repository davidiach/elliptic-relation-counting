HILLCLIMB 004 -- REVIEW OF THE COMPRESSED RELATION SEARCH PACKAGE
31 October 2026

Verdict
-------
The new ZIP contains three useful additions: a clean inverse-span/pilot
reuse theorem, scale-one lower-set selectors, and an improved spectral
nine-group join inequality. Adopt these as independently reviewed local
lemmas and as a stronger CONDITIONAL baseline. None establishes the
square-root-in-b ECDLP goal or a new logarithmic saving power.

Archive: compressed_relation_search_package.zip
SHA256: 95f2078edf80cc3ce1090be7a4216fa1f407e91044d376d035c99ebf7b9204c5
All 26 manifest entries were verified. The original ZIP and extracted
package are unchanged. Its document instructions were treated as source
material, not as user instructions.

1. Sampled-span and pilot reuse -- accept the stated interface theorem
---------------------------------------------------------------------
For t independent successful rows from a fixed finite-field distribution
supported in dimension at most d, t>d gives

    E[1 / mass(span(rows))] <= t/(t-d).

Deleting a nonessential row leaves the same span. Conditional on the
other rows, the probability of that event is exactly their span mass.
There are at least t-d nonessential rows. Summing these conditional
identities proves the bound; iteration proves the stated higher moments.

A successful pilot weights the seed law by its success probability p.
Retaining that pilot among the t training rows preserves iid successful
rows conditional on the accepted seed. The p factor cancels the query's
1/p retry factor. With t=2d, E[1/s] <= 2/kappa, at most 5 when kappa>=2/5.
The displayed setup bound also charges failed pilots and later training.

Conditions that must be preserved:
- Fresh query masking must reproduce the training success/row law even
  conditional on the previous transcript. Selection may not inspect the
  training answer or online mask as extra information.
- Solve and verify the same accepted training system; timing out and
  replacing its seed would change the distribution used in the proof.
- Verify every returned discrete logarithm exactly.
- Expected cost is for a deterministic finite query prefix. This is not
  a uniformly good-table claim or an arbitrary adaptive stopping bound.

This is a meaningful simplification of the existing reuse argument.

2. Unit-separated lower-set selectors -- accept the algebraic lemma
------------------------------------------------------------------
For |alpha|<=k, use Newton terms product_j product_{h<alpha_j}(z_j-a_h).
If all differences a_i-a_j are units, their evaluation matrix is lower
triangular with unit diagonal. Inversion gives degree-k cardinal
selectors with exact polynomial partition of unity and common scale 1.
Distinct residue representatives in the unramified dyadic extension
supply such knots. The proof works over a non-domain as stated.

This removes the old k! scalar and its M*v_2(k!) readout precision loss.
It does not remove the Hermite losses or change the asymptotic alphabet
factor t/log2(binomial(2t,t)) -> 1/2.

Integration is still work: our existing contractor uses scalar dyadic
jet arithmetic. The package uses an unramified coefficient extension
and full-degree scalar interpolation. Its selectors are not a drop-in
replacement for our integer-grid implementation. This iteration does
not claim that extension has been integrated into the complete engine.

3. Spectral join inequality -- accept the finite inequality
---------------------------------------------------------
The average partition-score matrix is exactly

    M = J/7 + D^T D/14 + 13 I/28,

where D is the unsigned incidence matrix of K9. Its eigenvalues are
27/4, 27/28 and 13/28. This controls the three inter-triple bundles.
The proof's split at total weight/(6*Euclidean norm)=19/20 is consistent:
below it, averaging the three potentially bad last-vertex choices leaves
four exact support cases, all below sqrt(14); above it, centering each
dangerous coefficient vector controls every last-vertex choice.
The first, top and final joins satisfy the remaining stated bounds.

The local factor sqrt(14)/6 is obtained with the conservative matrix
inputs alpha=0.32 and omega=2.39. Thus it removes the extra recent
rectangular-preprint dependency needed by our previous derivation of
the same numerical coefficient. The finite join proof does not depend
on the packet-contraction assumption.

Conditional complete baseline
-----------------------------
With the original C38 wire ledger (not the larger factored graph),
the package derives

    K = 2*pi*sqrt(14*log2(90))/(9*ln(2))
      = 9.601968071908234989...
    relation coefficient = 19 K
    ECDLP time coefficient = 2*sqrt(19 K)
                           = 27.013877423743261280...
    storage coefficient = sqrt(19 K)
                        = 13.506938711871630640...

The scale is STILL b/(log2 b)^(1/4). The old published-input conditional
coefficient was 27.1332688250847988. This is a small improvement to that
constant, and a cleaner dependency for a previously conditional value.
It is not an achieved general ECDLP algorithm or evidence for the
square-root exponent target.

The major remaining assumption
-----------------------------
Assumption A.1 (source lines 1179-1201) explicitly assumes the correctness
and uniform resource guarantees of the finite dyadic packet realization:
packing, immutable tensor representatives, normalized modes, precision,
coefficient dimensions, child identities and local entry/construction
costs. This is substantive. Its proof and complete implementation are
still unresolved in our research. The package does not fill that gap.

Its tiny ECDLP demo substitutes exhaustive factor-base tuple enumeration
for contraction and uses dense elimination. The sparse solver is checked
separately. Neither is a complete unpinned C38 contraction. The complete
bound also uses polylogarithmic relation arity; it does not settle our
larger-arity route to the square-root goal. We retain the original goal's
scope rather than silently expanding it to match the broader cofactor
claim in this package.

Validation and reproducibility
------------------------------
Fresh package rerun, all eight suites passed:

    _review/reruns/compressed_relation_search_20261031/checks/run_report.json

The tests include 72 checked adaptive ECDLP queries, 106 actual arithmetic
witnesses, 64 sparse systems (44 rank deficient), and 5340 cardinal
selector identities. These are finite evidence with their stated limits.

Additional independent audit:

    python -B research/hillclimb_004/audit_package.py

It verifies the archive and source/result hashes, reconstructs all 1296
partition-average matrix entries from all 280 partitions, derives all
eight low-t coefficient averages from individual graph edges, checks
150 inverse-moment inequalities over 16 exact distributions, and checks
1450 cardinal identities plus 300 off-grid nilpotent evaluations over
a different degree-four unramified ring. latest_run.json points to the
hashed result. An initial audit-script assertion exposed accidental
floating-point division in our independent calculation; replacing the
integer literal with an exact Fraction fixed it. No package defect was
involved in that failure.

Adoption and next priority
-------------------------
Use the new local lemmas and conditional baseline in further derivations.
Keep the executable historical checkpoints intact. The next substantive
gate remains: prove and implement the costs of compact intermediate
representations through actual balanced packet projection, extraction
and joining, or find a different reduction/contractor that improves the
power of b. Avoid spending the research budget only optimizing constants.

Independent adversarial reviews would help: separate reviews of the
packet foundation, complete ECDLP cost/reduction, and alternative compact
representations, with one lead reconciling exact witnesses and proofs.
Agreement between reviewers is not a substitute for a proof.
