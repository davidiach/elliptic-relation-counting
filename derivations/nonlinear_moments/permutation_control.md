# A control in which every address is legal and each grid has one zero

Consider $D=2^d\le p$ legal binary addresses and a public set $S$ of $D$ distinct field labels. Draw a uniform bijection $\pi:\{0,1\}^d\to S$ and independently choose uniform $u\in S$. Give address $a$ the residual $\pi(a)-u$. Exactly one address $v_*=\pi^{-1}(u)$ has residual zero. The decoder knows $S$ and may even know $u$. Its only access to the assignment $\pi$ is through the specified sketches and, optionally, $A$ bits of advice computed from $\pi$ before choosing $u$. This is a comparison model for the role of the residual values' arrangement. Actual curve seeds give a structured assignment whose distribution is different.

For each uniform binary challenge $r$ permit $t$ preselected Walsh masks, chosen using $r$, independent coins and the advice, before receiving residual values or using $u$ to select a mask. At every such mask $h$ return all signed field moments through order $k$,

$$
\sum_a(-1)^{h\cdot a}(\pi(a)-u)^j,\qquad 1\le j\le k .
$$

The output $f(r)\in[-1,1]$ may use any nonlinear computation on the values, the challenge, advice and shift. It may query adaptively within a preselected menu if all masks in the menu are charged. Define

$$
B=A+2tk\log_2p,\qquad e_D=D-\log_2\binom{D}{D/2}.
$$

Then its expected absolute Walsh coefficient at the zero address satisfies

$$
\mathbb E|\widehat f(v_*)|
\le \left(\frac1D+\sqrt{\frac{2\ln2}{D}(B+e_D)}\right)^{1/2}.
$$

For $D=2^{b-o(b)}$ and $tk,A=2^{o(b)}$ this is at most $2^{-b/4+o(b)}$. The bound allows arbitrary polarity on each realization and free postprocessing. The advice may encode any function of the assignment within its stated bit length.

To prove it, square the Walsh coefficient and introduce independent challenges $r,r'$. When $r=r'$ bound the integrand by one. Otherwise $s=r+r'\ne0$, with addition in the binary cube. Condition on these challenges and all independent coins. For each of the at most $2t$ masks, a single message can give the **unshifted** moments through order $k$. Binomial translation answers all shifted observations for *any* $u$ from those values, since the constant moment is public. Together with the advice, this message $Z$ carries at most $B$ bits of entropy. Its choice is independent of $u$.

Index the labels by $x\in S$. The signs $Y_x=(-1)^{s\cdot\pi^{-1}(x)}$ form a uniform balanced sign string, of entropy $D-e_D$ bits. Hence $H(Y\mid Z)\ge D-e_D-B$. Subadditivity gives $\sum_x H(Y_x\mid Z)\ge H(Y\mid Z)$. For a biased sign of mean $a$, the binary entropy satisfies $1-h_2((1+a)/2)\ge a^2/(2\ln 2)$ by differentiation. Summing over labels yields $\mathbb E_{x,Z}|\mathbb E(Y_x\mid Z)|\le\sqrt{2\ln2(B+e_D)/D}$. The product of the two readouts is bounded by one and determined by $Z,u$ at the fixed challenges and coins. Its expected correlation with $Y_u$ is at most that expression. Average over challenges, then take a square root. Also $e_D\le\log_2(D+1)$ by the largest binomial coefficient bound. The entropy argument is the classical random-access-code counting principle. [Nayak, 1999](https://arxiv.org/abs/quant-ph/9904093)

There is an exact finite contrast. The field values $M_0(r)-M_{p-1}(r)$ give $(-1)^{r\cdot v_*}$ on every such input. The shifted functions $x\mapsto1_{x=u}$ span $D$ dimensions on $S$, so the translation argument would have to retain $D$ independent weights to simulate these zero tests at one mask for all shifts. More generally the same proof charges the dimension of the span of all translated per-value functions on $S$ in place of $k$ (with any known constant component free). A family of weighted sums of low powers has small translation dimension even when its final readout is unrestricted.

The control shows why isolation, a large signed spectrum, and unrestricted nonlinear arithmetic on the returned low moments alone do not guarantee success for arbitrary residual assignments. Transferring the bound to the elliptic law would require a new argument about structured curve rows. [Actual mixed-curve spectrum calculation](signed_spectrum.md)
