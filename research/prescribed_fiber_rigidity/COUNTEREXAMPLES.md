# Exact examples and their scope

### Sharp integrality threshold

In \(C_r\), take the interval \([0,s-1]\), labeling each point by its integer representative. Its ordered double totals run through \([0,2s-2]\), giving exactly \(2s-r-1\) mixed residues when \(3s>2r\). At \(s=\lfloor(2r+1)/3\rfloor,\ r\ge5\), every triple summing to \(r-1\) has integer sum \(r-1\). This is a nonaffine integer score constant on the complete prescribed fiber. Translation gives an interval avoiding identity. Over \(\mathbb F_r\) this particular score is affine. The example establishes the optimal **abelian-codomain** density threshold and the need for its strict inequality.

### A single decomposition can fail at every nonzero target

The verifier takes \(C_7\) with blocks \((1,2,3)\) and \((2,4,6)\), deleting the first coordinate of each. Let \(W\) be the span of the incidence change from \((3,4)\) at zero to \((2,6)\) at one. The reference triple at zero is the lexicographically first reduced incidence in the sorted fiber, namely \((0,0,0,3)\). All six nonzero targets have every single translated incidence outside \(W\). At every target the affine hull in the three-dimensional quotient contains zero. The program verifies a combination explicitly by field arithmetic. This fixture concerns a chosen subsystem and fixed identity reference. It supplies no lower bound on the size of certificates under all choices of reference or subsystem.

### A flat unweighted fiber with a large averaging radical

Take one block \(F=C_r\setminus\{0\}\) for any prime \(r\ge5\), and all ordered triples of sum zero. Each fixed slot and value has mass \(r-2\). Let \(R\) reflect \(x\) to \(-x\) and \(J\) be all ones. Counting positions gives \(n=(r-1)(r-2)\), first moment \(v=3(r-2)\mathbf1\), and second moment \(H=3(r-2)I+6(J-R)\). A mixed-slot pair of values \(i,j\) completes uniquely precisely when \(j\ne-i\). Hence the Gram of **all ordered pairs of triples** is
\[
 2(nH-vv^\top)\equiv-24(I+R)-48J\pmod r .
\]
It kills every odd score, including the nonaffine \(u(x)=x^3\), and every constant score. Its action on mean-zero even scores is invertible. Its rank is \((r-3)/2\), compared with raw rank \(r-3\). Pinning the common third slot on both sides restores raw rank for this example and for every dense proper one-block set in the pinning theorem. A single aggregate and a stack of cell Grams have different kernels. The sum of within-cell Grams also omits cross-cell terms from the all-pairs Gram.

The single-slot law on \(F\) has chi-square \(1/(r-1)\) relative to uniform, which tends to zero. The example has one large block. It gives no counterexample to the parity-affine Stage I construction with growing arity, to PR #7's two-pin fixed-shift problem, or to the existence of cheap selectors for that construction. [Reproducible arithmetic](verification.txt)

### A corank-one family on the actual Stage I blocks

For sufficiently large inputs, remove just the reference point of each block. In the manuscript we show that the surviving short sumset still has size greater than \(r/2\), using \(M=\prod_i |F_i|=O(r)\) and \(\min_i|F_i|\ge r^{1/m}2^{-O(m)}\). Take the shift rows on fully nonreference doubled tuples and add shifts exchanging the single reference occurrence between two block positions along a tree. Their modular row space is exactly the hyperplane of vectors whose entries sum to zero. Every row is genuine, and a basis exists. Its efficient construction with all right sides by this argument is open.

Fix the all-reference triple, whose group value is explicitly known and need not be the identity. Every triple at a different value has \(1,\dots,3m\) nonreference occurrences, so every individual target functional relative to that baseline is ambiguous for this single missing direction. Allow both endpoints to move: the joint nonreference fraction is \(1-o(1)\), so the given hyperplane with known right sides does certify every target by the first-order selector. This is evidence that efficient choice of representation can matter even on the prescribed blocks. The constructed genuine subsystem has rational corank one and is not the full-rational basis from integer extraction. PR #9 gives a separate finite obstruction for such a rational basis. [Proof](manuscript.pdf) · [Concurrent rational-basis obstruction](https://github.com/davidiach/elliptic-relation-counting/pull/9)
