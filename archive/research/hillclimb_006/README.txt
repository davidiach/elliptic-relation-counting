HILLCLIMB 006 -- NONUNIFORM GRIDS AND THEIR TRUE PRECISION COST
31 October 2026

Outcome
-------
Tested the geometric coordinate set {0,1,2,4,...,2^(n-2)} as a way to
replace linear mode weights by quadratic ones. The quadratic diagonal
gain is real. The naive interpolation grid is invalid: a local constant
coefficient has a negative 2-adic valuation of cubic magnitude in n.
Changing the jet-variable scale cannot repair a constant term.

Rather than stopping there, constructed a polynomial-size repair by
adding one nearby partner of each nonzero interpolation node. Derived
the minimal common separation and, for n>=4, the minimal common local
jet scale for that paired construction. The repaired coefficients are
integral, but the readout cost erases the apparent gain for the tested
odd-node, unrestricted-enumeration route.

Concrete findings
-----------------
- The distinct difference grid has (n-2)^2+(n-2)+3 nodes. Its derivative
  valuation at zero is ((n-2)^3+5(n-2))/3; at one it is
  (n^2-n-4)/2. For n=8, the uncorrected local constant has valuation -56.
- Pairing nodes restores integrality with O(n^2) nodes, O(n^3)-bit
  coordinates, and a common jet scale exponent of order n^3. It does
  not require exponentially many interpolation nodes.
- The even-block mode increments are k(k+1), checked by independent
  exact Cauchy elimination. Normalized triangular coefficients are
  integral; their off-diagonal valuations are not all zero.
- A realizable network with all original variables pinned to one has
  original count one but highly divisible local jets. It forces both
  sufficient jet length and sufficient precision at a retained node.
- For n>=7, reconstruction using only odd-node jet residues followed by
  unrestricted feasible-mode enumeration must encounter at least
  2^(M/4) candidates at one required node, where M is the edge count
  and M is divisible by 4n. Thus this particular repaired route does
  not produce a sublinear exponent in M.
- Uniform scaling of every coordinate and interpolation node leaves
  the deformation exactly unchanged after scaling T. Uniformly larger
  gaps alone cannot improve the algorithm.

This is not a lower bound for all nonuniform grids or all ECDLP solvers.
The candidate-universe result does not say every candidate contributes.
Support-aware methods, coefficient-dependent precision, even-node
schedules, and different repairs are not ruled out. The all-one witness
is itself easy to contract; it is a necessary readout test, not a hard
instance for arbitrary algorithms. No new complete runtime is claimed.

Evidence
--------
Independent standard-library code checks:
  16 coordinate sizes, n=3 through 18;
  5024 derivative valuations plus direct derivative-product controls;
  16 exact Cauchy matrix eliminations and reconstruction identities;
  22 directly expanded packet jets, including the n=7 threshold case;
  54 exact uniform-dilation identities;
  15 initial balanced-cluster grid fixtures.

The asymptotic statements are proved in grid_audit.tex. The finite checks
are supporting evidence. One initial diagnostic incorrectly expected
all triangular entries to be units; elimination caught this and the
correct integral valuation formula was derived and checked instead.

Next branch
-----------
The simple geometric set is poorly balanced across parity classes. A
more relevant candidate is a union of full low-bit clusters:

  X_(u,n) = {2^u*x+r : x in X_n, 0 <= r < 2^u}.

It fills every u-bit orientation equally. The pilot shows that checking
only central nodes would miss failures: at n=8, u=1,2,3, some retained
nodes have constant valuation -16 even though the central node is
integral. The next gate is to choose admissible nodes explicitly and
derive their full Hermite multiplicities, precision and mode spectrum.
The original precision schedule cannot simply be reused.

Reproduce
---------
    python -B research/hillclimb_006/run_all.py

Each run retains source hashes, exact result JSON, stdout and stderr in
a fresh runs/ directory. latest_run.json identifies the latest checkpoint.
Original archives and earlier executable checkpoints remain unchanged.

The square-root-in-b ECDLP goal remains active and unachieved. The latest
conditional baseline is still 27.0138774*b/(log2 b)^(1/4), with its stated
packet-framework assumptions.
