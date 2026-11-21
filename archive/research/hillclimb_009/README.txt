HILLCLIMB 009 -- ACTUAL PACKET SPACES AND PARITY RECOVERY
31 October 2026

Outcome
-------
The smaller-polynomial-space hypothesis fails on the tested complete
generic packets. Their legal indicator-tensor deformations span every
polynomial allowed by the degree bound:

  packet size L     degree-space dimension     certified rank
       2                       7                      7
       4                      29                     29
       8                     121                    121

At L=8 this remains true with a fixed actual C38 source-predicate
restriction: its 256 source assignments span the full 121-dimensional
space. The predicate is regenerated from the public arithmetic instance
and checked on every surviving Boolean assignment. It is not a newly
invented equality tensor. The local packing is complete and legal, but
does not claim to satisfy the asymptotic balanced schedule.

Products of k independent copies of that size-eight packet span every
polynomial of degree at most 120*k. This is a proved fixed-width family,
not a theorem for all growing packet sizes or all complete C38 graphs.

The parity-only loophole was checked separately
---------------------------------------------
Exact linear readout on a rational span differs from recovering parity
on integer combinations of the actual generators. A new finite solver
optimizes the latter problem. It retains both reconstruction weights
and an input combination with odd target value but zero observations at
the next smaller precision profile.

For jet precisions 1+floor(lambda/2^v2(tau)), the certified scales are:

  L   family       full Hermite   optimal exact full space   family parity
  2   generic          10                   9                     9
  2   C38 copy         10                   9                     4
  4   generic          45                  27                    27
  4   C38 copy         45                  27                     8
  8   C38 copy        128                 102                   102

The small copy cases have proper rational spans, so their gains cannot
be credited to coefficient integrality alone. At L=8 the full span and
the parity optimum both match the general degree-limited model.

The L=8 lower certificate has odd value at zero, but zero jets modulo
2^102 at odd centers, 2^51 at valuation-one centers, and 2^26 at
valuation-two centers. Those observations therefore cannot determine
the parity on the full integer-generated family. This is an observation
model lower bound. It does not exclude reading the source tensor itself
or using stronger restrictions on the source weights.

Validation
----------
The complete audit passed in about 17 seconds and retains:

* Six exact rational-rank certificates via nonzero prime-field minors.
  The complete generic families at L=2,4 also pass exact elimination.
* 435 native-kernel coefficient comparisons.
* 180 native selector-jet comparisons using the actual C38 predicate.
* 75 coefficients checked by complete original-assignment deformation
  sums for generic indicator networks.
* Five optimal parity-profile certificates, each with upper weights
  and a lower-precision collision, plus checked exact-readout controls.
* Two controls with the same rational span and different coefficient
  lattices, confirming that the modular-only solver treats them differently.
* 100 random integer-input trials with independently changed observation
  representatives, and 11 corrupted-certificate negative controls.

The generic size-eight rank search found a full minor in about 54
seconds. Reproduction rebuilds its 121 saved legal generators and verifies
the determinant, rather than repeating the search. This is an exact
certificate, not a sampled numerical rank claim.

Scope
-----
The variable source weights can be represented by a joint source tensor
adjacent to the pinned C38 predicate. They are not asserted to equal the
source predicates of a complete unrestricted C38 relation instance.
The source families and matrices are explicit: construction uses 4^L
generic entries or 2^L copy entries. Their time, storage and coefficient
sizes must be charged before using this approach in a contractor.

Neither a complete unpinned C38 contraction nor an improved ECDLP
asymptotic exponent has been obtained. The conditional baseline remains
27.0138774*b/(log2 b)^(1/4), and the square-root-in-b goal remains active.

Next research decision
----------------------
Do not infer a global compression method from a small local polynomial
space. The next gate returns to the actual arithmetic network: establish
the exact fibers of its source projection, including every auxiliary
wire and state. That can support complete structured marginalization of
a small unpinned C38 graph and expose which arithmetic structure survives
elimination. It must be distinguished from the earlier one-way witness
construction and from a fast generic packet contractor.

Reproduce
---------
    python -B research/hillclimb_009/run_all.py

latest_run.json selects the successful run. The run directory contains
source hashes, exact certificates, stdout, stderr and execution status.
Original archives, PDF files, vendor modules, and prior executable
checkpoints were preserved.
