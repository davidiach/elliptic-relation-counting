# Direct Newton certificates and sharpness at active weights four and five

27 November 2026. Follow-up to [PR #19](https://github.com/davidiach/elliptic-relation-counting/pull/19), `a2dc9e57d2d27659b152d0d1e4744c495e6506ee`. All results are about abstract integer row families and exact moments; none establishes a faster algorithm for, or a lower bound on, the project's actual elliptic relation systems. Classical Boolean Fourier expansion, Newton divided differences and orthogonal arrays are not claimed novel. Independent priority and proof review remain open.

## 1. Direct local coefficients for the previous certificate

Write d for the number of nonzero coordinates of a field direction u over F_r, r an odd prime. Discard inactive coordinates. Let X_W={x in Z^d: ||x||_1<=W}, where 1<=W<r. A target c is any residue. We construct a rational polynomial P whose restriction to X_W is nonnegative and vanishes exactly where x dot u=c mod r. All coordinates of u are nonzero after projection.

For each magnitude vector alpha in the nonnegative integer simplex |alpha|<=W, let S=supp(alpha), and let v_alpha count the violating sign patterns on S. A nonzero alpha always has at least half its patterns violating: pair patterns by flipping one nonzero coordinate. The score changes by 2 alpha_i u_i, a nonzero residue. At alpha=0, there is one pattern, legal if c=0 and violating otherwise.

Use B(alpha)=sum_i alpha_i^2 if c=0 and B(alpha)=1 otherwise. On each magnitude pattern define F_alpha to be zero on legal signs and B(alpha)*2^|S|/v_alpha on violating signs. When v_alpha=0 set F_alpha=0 (this occurs only at the zero pattern for c=0). This function has sign-average B(alpha). In the homogeneous case it is even under simultaneous sign reversal, so its nonconstant Boolean Fourier coefficients are indexed by nonempty even subsets R of S. In the affine case all nonempty subsets can occur.

Let f_R(alpha) denote the ordinary normalized Fourier coefficient of F_alpha at the sign character R. Explicitly it equals B(alpha)/v_alpha times the sum of that character over violating signs. Fix a nonempty R, eta=1_R, and set

    g_R(gamma)=f_R(gamma+eta)/product_{i in R}(gamma_i+1),
    |gamma| <= W-|R|.

Define the local nodes y_{i,j}=(j+eta_i)^2. For beta in the simplex |beta|<=W-|R|, define the divided difference

    C_{R,beta} = sum_{0<=gamma<=beta} g_R(gamma) /
                 product_i product_{0<=j<=beta_i, j!=gamma_i}
                      (y_{i,gamma_i}-y_{i,j}).

Then the polynomial is

    P(x) = B(x) + sum_R product_{i in R} x_i *
           sum_{|beta|<=W-|R|} C_{R,beta} *
              product_i product_{j=0}^{beta_i-1}(x_i^2-y_{i,j}),

where B(x)=sum x_i^2 for c=0 and B(x)=1 for c!=0, and for c=0 retain only even R. Empty sums/products have their usual meanings.

**Proof of interpolation.** For each R, the ordinary tensor-product Newton interpolation formula holds on every Cartesian prefix grid of index gamma; its coefficients are the displayed product divided differences. At a simplex point gamma, every basis element with some beta_i>gamma_i vanishes because its factor includes node y_{i,gamma_i}. Thus only beta<=gamma contribute, all contained in the simplex. This proves the same interpolation formula on the lower set. Multiplication by product_{i in R} x_i gives zero where R is not contained in the support. Otherwise it restores the product of magnitudes and the sign character. Summing the Fourier expansion yields exactly F_alpha on every x in X_W, proving nonnegativity and the exact zero criterion. QED.

The degree of a nonconstant term is |R|+2|beta|<=2W-|R|. Thus the homogeneous bound max(2,2W-2) from PR #19 is retained, and the new affine bound is 2W-1 for c!=0. In particular a nonnegative exact rational moment certificate works for nonzero right sides too, without count or multiplicity restrictions. If every active coordinate is zero, use the constant zero or one certificate according to the target.

## 2. Representation size and exact construction work

There is one potential Newton summand per pair (R,beta). The map (R,beta)->(alpha=beta+1_R,R) identifies these with magnitude patterns |alpha|<=W and nonempty subsets of their supports (even subsets only in the homogeneous case). Hence the number of Newton terms is at most 2^W binom(d+W,W); each involves at most W distinct variables. Each term expands into at most product_i(beta_i+1)<=2^W monomials. Including B, the expanded term count is at most 4^W binom(d+W,W)+d+1. For fixed W, this gives O_W(d^W) rather than the crude degree-only O_W(d^(2W-2)) bound for W>=2.

Each coefficient is local: gamma<=beta has at most 2^W choices, and g_R uses at most 2^W sign evaluations, all on at most W coordinates. The denominators of each divided-difference summand divide the common positive integer

    L_W=(2^W)! * ((W^2)!)^W.

Indeed the sign count v_alpha divides (2^W)! in the integer-divisibility sense; the remaining denominator factors are the |R| magnitude factors and |beta| squared-node differences, each of absolute value at most W^2, and there are at most W such factors. Thus each rational coefficient has numerator and denominator bit length 2^{O(W)}: its absolute value is at most 2^W max(1,W^2), since all nonzero integer factors have magnitude at least one. Exact integer arithmetic, with a common denominator if desired, gives a total construction bound 2^{O(W)} binom(d+W,W) poly(log r,log(d+1)) bit operations in a sparse Newton representation, allowing ordinary polynomial-time integer operations and enumerating sparse multi-indices. This count includes reading u. No global linear system in d is needed. The expanded monomial bound follows by distributing the products; it is not necessary to expand them for the algebraic result.

This is a bound for explicit abstract certificates; it does not prove that obtaining or evaluating the requisite exact moments for the actual implicit relation systems fits their target resource budget. In particular d may be large and W may grow. No new end-to-end ECDLP complexity claim is made.

## 3. Sharp homogeneous bounds at W=4 and W=5

The following exact compressed moment trades show that the homogeneous degree 2W-2 is optimal for W=4,5 in the unrestricted integer-multiset model, with r>W, full rational rank and modular rank separation. Both use r=7 and u=(1,2,3). The signed permutation orbit O(a) of a in Z^3 is formed by multiplying each residue label i by k in F_7^*, writing ki=sign*j with j in {1,2,3}, and sending coefficient a_i to sign*a_i at j. Include each distinct vector in the orbit once. This multiplies the modular score by k and preserves l1 norm. An orbit specification (a,n) gives multiplicity n to every element of O(a).

For W=4 the two compressed populations are:

    H: (-2,-1,-1):1, (-2,1,0):2
    J: (-2,0,-1):2, (-2,1,1):1.

Each has 18 rows counting multiplicity, 12 distinct rows, and is centrally symmetric. Each H row is legal; every J row violates u. Their exact moments agree through degree five and differ at degree six. For example the even monomial equality can be checked by the orbit definition (only degrees 0,2,4 need checking); the odd moments vanish. Their ranks over Q are both three, over F_7 both two, and the union has rank three. Thus H+H versus H+J gives full rational rank and respective modular ranks two and three, with equal fifth moments.

For W=5, a finite orbit table with 949,116 rows counting multiplicity in each population is included as H5,J5 in `verify.py`. There are 32 and 144 distinct rows respectively. Every H row is legal, every J row violates; the first seven degrees agree exactly, and degree eight differs. Both have full rational rank, and the respective ranks over F_7 are two and three. The orbit weights were found with numerical linear optimization and integer-lattice search, but are fixed integers and the submitted proof certificate is checked entirely with exact arithmetic (enumerating all monomials). The displayed weight-four certificate and the weight-five table are finite constructive counterexamples, not a theorem about all larger weights. No optimizer is required to reproduce the checks, and no minimality of the row counts is claimed.

Both examples have counts above r. They do NOT give no-wrap sharpness at W=4 or W=5; the no-wrap and signed-unit cases there remain open. They use entries of magnitude larger than one.

## 4. Distinct-row lifting without changing active weight

For any two integer multisets agreeing in exact moments through t>=2, choose h>=1 with 2^h at least the largest multiplicity. Append h groups of t+1 inactive binary tag coordinates. Assign to each copy of a base row a different tag class, specified by one parity bit per group; replace the copy by all tags in that class. Each class has size 2^{ht}, and every assignment of at most t selected coordinates has identical counts in every class: at least one coordinate in each parity group remains free. Consequently every mixed monomial of total degree <=t has the original moment multiplied by a common tag-marginal count. Rows are distinct, can be ordered lexicographically, and their active weight and modular legality are unchanged on extending u by zeros. Counts grow and do not in general satisfy no-wrap.

Over any field of characteristic not two each parity class affinely spans the entire tag space: the even parity set includes 0 and e_i+e_j for i!=j, which span the group space for group length t+1>=3. Odd parity is its image under a coordinate complement, an invertible affine map. Products of classes preserve this affine-span fact. After appending an inactive constant coordinate 1, tag differences at a single base row span every tag coordinate; rank is therefore the tag dimension plus the rank of the augmented base rows. Central symmetry makes the constant coordinate add exactly one rank to the original linear rank over Q and over F_7. The rank gap is preserved. This supplies distinct-row versions of section 3 with the same active weight and full rational rank; no-wrap is not preserved, and no full large tagged instance is enumerated.

For an affine example at W=2, r=3, u=1, c=1, take H=(-2,1,1), J=(-1,-1,2). The zeroth, first and second moments are (3,0,6) in both cases. Only H satisfies the affine equation, and third moments differ. Thus the affine degree 2W-1 bound is sharp at W=2 for unrestricted multisets; at W=1 degree one is needed since counts alone cannot distinguish x=0 from x=1. Larger affine sharpness is open.

## 5. Verification and limitations

Run `python -B verify.py` with Python 3.10+ (standard library only, assertions required). The recorded transcript checks 54 homogeneous/affine certificates on complete small integer balls, 4,464 exact rational point comparisons, the two compressed homogeneous trades, the affine W=2 example, and 1,096 tag marginal comparisons with exact affine ranks. Rational rank of the homogeneous trade supports is certified by full rank modulo the separate prime 1009. These are finite checks, not a substitute for the Newton proof or a wider novelty review. The signed-unit and no-wrap sharpness cases for W>=4, and the representation/evaluation problem for the project's actual structured row family, are unresolved.
