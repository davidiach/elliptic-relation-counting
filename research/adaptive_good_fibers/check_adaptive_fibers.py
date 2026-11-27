"""Exact finite controls for affine-mask good-fiber and translation lemmas.
Only small synthetic cyclic groups and toy prime-order elliptic curves.
Standard library; assertions required. No asymptotic contractor is implemented.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from functools import lru_cache
from math import prod, comb
from random import Random
import subprocess
import sys
from unittest import TestCase

if not __debug__:
    raise RuntimeError("This verifier requires Python assertions.")


def conv(a,b):
    r=len(a); z=[0]*r
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:z[(i+j)%r]+=x*y
    return z


def counts(r,blocks):
    z=[1]+[0]*(r-1)
    for B in blocks:
        a=[0]*r
        for x in B:a[x%r]+=1
        z=conv(z,a)
    assert sum(z)==prod(map(len,blocks))
    return z


def curve(p,a,b):
    def add(P,R):
        if P is None:return R
        if R is None:return P
        x,y=P;z,w=R
        if x==z and (y+w)%p==0:return None
        t=((3*x*x+a)*pow(2*y,-1,p) if P==R else (w-y)*pow(z-x,-1,p))%p
        X=(t*t-x-z)%p
        return X,(t*(x-X)-y)%p
    pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-b)%p==0]
    G=[None]; P=pts[0]
    for _ in pts:G.append(add(G[-1],P))
    assert len(set(G))==len(G) and set(G)=={None,*pts} and add(G[-1],P) is None
    return G,add


def retained(m,k):
    return [I for a in range(m-k,m+1) for I in combinations(range(m),a)]


def pins(I,h):
    return [J for a in range(min(h,len(I))+1) for J in combinations(I,a)]


def energy(r,B,I,J):
    """r^2 E_{I,J}(g), by independent pinned completion counts."""
    R=[i for i in I if i not in J]
    d=counts(r,[B[i] for i in R]);W=sum(d)
    n=counts(r,[B[i] for i in J])
    return [sum(n[x]*(r*d[(g-x)%r]-W)**2 for x in range(r)) for g in range(r)]


def pair_conv(r,A,B):
    out={}
    for (x,y),a in A.items():
        for (z,w),b in B.items():
            key=((x+z)%r,(y+w)%r)
            out[key]=out.get(key,0)+a*b
    return out


def shared_count_energy(r,tables,I,J):
    z={(0,0):1};den=1
    for i in I:
        H,D=tables[i];den*=D
        A={xy:v for xy,v in H.items() if xy[0]==xy[1]} if i in J else H
        z=pair_conv(r,z,A)
    rows=[0]*r;cols=[0]*r
    for (x,y),v in z.items():rows[x]+=v;cols[y]+=v
    W=sum(z.values())
    return [Q(r*r*z.get((g,g),0)-r*(rows[g]+cols[g])+W,den) for g in range(r)]


def fixed_tables(B):
    return [({(x,y):1 for x in F for y in F},1) for F in B]


def phi_arrays(r,B,q,h,k):
    result=[]
    for I in retained(len(B),k):
        mu=q**len(I)/r;z=[Q(0)]*r
        for J in pins(I,h):
            w=q**(2*len(J))/(r*r*mu*mu)
            z=[a+w*b for a,b in zip(z,energy(r,B,I,J))]
        result.append((I,z))
    return result


def small_patterns(m,limit,coeffs):
    """At most two distinct labels per block. Ordered variables overcount."""
    out=[]
    for s in range(1,limit+1):
        for slots in product(range(m),repeat=s):
            if tuple(sorted(slots))!=slots or any(slots.count(i)>2 for i in set(slots)):continue
            for cs in product(coeffs,repeat=s):out.append(tuple(zip(slots,cs)))
    return out


def pattern_expect(r,tables,pattern):
    den=1;D=[1]+[0]*(r-1)
    for i in sorted(set(i for i,c in pattern)):
        cs=[c for j,c in pattern if j==i]
        H,d=tables[i];den*=d;a=[0]*r
        if len(cs)==1:
            for (x,y),v in H.items():
                if x==y:a[cs[0]*x%r]+=v
        else:
            for (x,y),v in H.items():
                if x!=y:a[(cs[0]*x+cs[1]*y)%r]+=v
        D=conv(D,a)
    return Q(D[0],den)


def direct_pattern(r,B,pat):
    n=0
    for X in product(*(B[i] for i,c in pat)):
        if any(pat[i][0]==pat[j][0] and X[i]==X[j] for i in range(len(pat)) for j in range(i)):continue
        if sum(c*x for (_,c),x in zip(pat,X))%r==0:n+=1
    return n


def mask_prefix_case(p,a,b,T,m,h,k,do_prefix):
    G,add=curve(p,a,b);r=len(G);q=Q((r-1)*T,p);theta=Q(T,p)
    opts=tuple(product(range(p),repeat=2));domains=[[opts,opts] for _ in range(m)]
    @lru_cache(None)
    def parity(par,S):
        sets=[{i for i in range(1,r) if G[i][1]%2==par and (u*G[i][0]+v)%p<T} for u,v in S]
        one={x:sum(x in z for z in sets) for x in range(1,r)}
        pair={(x,y):sum(x in z and y in z for z in sets) for x,y in product(range(1,r),repeat=2)}
        return one,pair
    @lru_cache(None)
    def table(A,B):
        a,aa=parity(0,A);b,bb=parity(1,B);H={}
        for x,y in product(range(1,r),repeat=2):
            xp=G[x][1]%2;yp=G[y][1]%2
            w=aa[x,y]*len(B) if xp==yp==0 else bb[x,y]*len(A) if xp==yp==1 else a[x]*b[y] if xp==0 else b[x]*a[y]
            if w:H[x,y]=w
        return H,len(A)*len(B)
    Isets=retained(m,k); IJ=[(I,J) for I in Isets for J in pins(I,h)]
    # Two-variable finite girth penalties. These are weaker than asymptotic guards.
    pats=small_patterns(m,2,(-1,1))
    anchors=(0,1);weight=2
    root_checks=0;pattern_checks=0
    root=[table(*s) for s in domains]
    H,den=root[0]
    for x,y in product(range(1,r),repeat=2):assert Q(H.get((x,y),0),den)==(theta if x==y else theta**2)
    for I,J in IJ:
        mu=q**len(I)/r
        e=shared_count_energy(r,root,I,J)
        assert all(0<=v<=2*mu*r*r for v in e);root_checks+=r
    for pat in pats:
        v=pattern_expect(r,root,pat)
        assert v<=Q(r*theta)**len(pat)/r;pattern_checks+=1
    def expect():
        tables=[table(*s) for s in domains];v=Q(0)
        for H,d in tables:
            n=Q(sum(a for (x,y),a in H.items() if x==y),d);n2=Q(sum(H.values()),d)
            v+=m**4*(n2-2*q*n+q*q)/q**2
        for I,J in IJ:
            mu=q**len(I)/r;e=shared_count_energy(r,tables,I,J)
            assert min(e)>=0
            v+=weight*q**(2*len(J))/(r*r*mu*mu)*(sum(e)/r+sum(e[g] for g in anchors))
        return v+sum(pattern_expect(r,tables,pat) for pat in pats)
    if not do_prefix:
        print(f'Actual affine law E/F{p}, order {r}, m={m}, h={h}, deletion {k}: {root_checks} higher-pin mean inequalities and {pattern_checks} relation-pattern bounds.')
        return root_checks,pattern_checks,0,None
    old=initial=expect();steps=0;bits=(p-1).bit_length();width=2*bits
    for i,par in product(range(m),range(2)):
        for bit in reversed(range(width)):
            S=domains[i][par]
            children=[tuple(x for x in S if (((x[0]<<bits)|x[1])>>bit)&1==v) for v in (0,1)]
            C=[x for x in children if x]
            if len(C)<2:continue
            vals=[]
            for child in C:domains[i][par]=child;vals.append(expect())
            assert old*len(S)==sum(len(child)*v for child,v in zip(C,vals))
            j=min(range(len(vals)),key=lambda j:vals[j]);domains[i][par]=C[j];old=vals[j];steps+=1
    B=[]
    for pair in domains:
        seed=[s[0] for s in pair]
        B.append([i for i in range(1,r) if (seed[G[i][1]%2][0]*G[i][0]+seed[G[i][1]%2][1])%p<T])
    direct=m**4*sum((Q(len(F))-q)**2/q**2 for F in B)
    for I,J in IJ:
        mu=q**len(I)/r;e=energy(r,B,I,J)
        assert [Q(x) for x in e]==shared_count_energy(r,fixed_tables(B),I,J)
        direct+=weight*q**(2*len(J))/(r*r*mu*mu)*(Q(sum(e),r)+sum(e[g] for g in anchors))
    for pat in pats:
        n=direct_pattern(r,B,pat);assert n==pattern_expect(r,fixed_tables(B),pat);direct+=n
    assert direct==old<=initial
    d=[0]*r
    for X in product(*B):
        z=None
        for x in X:z=add(z,G[x])
        d[G.index(z)]+=1
    assert d==counts(r,B)
    print(f'Actual affine prefix E/F{p}, order {r}: {root_checks} mean inequalities, {pattern_checks} pattern bounds, {steps} exact bit-prefix identities. Sizes {[len(F) for F in B]}, selected toy penalty {float(old):.5g}; asymptotic guard {old<Q(1,256)}.')
    return root_checks,pattern_checks,steps,B


def check_tuple_girth(r, blocks, t0):
    for x, y in combinations(product(*blocks), 2):
        if sum(x) % r == sum(y) % r:
            distance = sum(a != b for a, b in zip(x, y))
            assert 2*distance > t0, 'short equal-sum collision'


def two_slot_collision_check():
    r, t0 = 1009, 4
    good = [(1, 5), (25, 125)]
    bad = [(1, 5), (25, 29)]
    x, y = (1, 29), (5, 25)
    assert x in product(*bad) and y in product(*bad)
    assert sum(x) % r == sum(y) % r and all(a != b for a, b in zip(x, y))
    patterns = small_patterns(2, t0, (-2, -1, 1, 2))
    for blocks, expected_lengths in [(good, set()), (bad, {4})]:
        tables = fixed_tables(blocks)
        found = set()
        for pattern in patterns:
            n = direct_pattern(r, blocks, pattern)
            assert n == pattern_expect(r, tables, pattern)
            if n:
                found.add(len(pattern))
        assert found == expected_lengths, 'incorrect signed collision lengths'
    check_tuple_girth(r, good, t0)
    with TestCase().assertRaisesRegex(AssertionError, '^short equal-sum collision$'):
        check_tuple_girth(r, bad, t0)
    print('Two-slot collision at sum 30 rejected. Only length-four signed patterns detected it. The collision-free control passed.')


def shared_checks():
    rng=Random(1703);checks=patterns=0;girth_pass=0
    for r,m,h,k in [(5,5,3,1),(7,5,3,1),(11,6,4,1)]:
        for _ in range(8):
            B=[rng.sample(range(1,r),rng.randrange(0,min(r,5))) for _ in range(m)]
            T=fixed_tables(B)
            for I in retained(m,k):
                for J in pins(I,h):
                    assert energy(r,B,I,J)==shared_count_energy(r,T,I,J);checks+=r
            for pat in rng.sample(small_patterns(m,3,(-2,-1,1,2)),30):
                assert direct_pattern(r,B,pat)==pattern_expect(r,T,pat);patterns+=1
    # Sparse fixed blocks on a toy curve; direct verification of the consequence
    # that equal-total tuples cannot differ in at most floor(t0/2) slots.
    G,add=curve(13,0,2);r=len(G)
    for _ in range(100):
        B=[rng.sample(range(1,r),2) for i in range(3)];t0=2
        if any(direct_pattern(r,B,pat) for pat in small_patterns(3,t0,(-1,1))):continue
        check_tuple_girth(r, B, t0)
        girth_pass+=1
    print(f'Direct shared-pin checks: {checks} integer energy entries, including strength four and empty blocks; {patterns} signed relation-pattern comparisons; {girth_pass} finite girth controls.')


def select_translations(r,B,fs,s,encoding=None):
    # Each translate is an actual doubled tuple of the MAIN indexed blocks.
    slots=[tuple(F) for F in B for _ in range(2)];m=len(B)
    first=[(F[0],) for F in slots]
    domains=[first]+[list(slots) for _ in range(s-1)]
    assert all(slots)
    @lru_cache(None)
    def moment(j,D):
        d=counts(r,D);W=sum(d)
        return tuple(tuple(sum(Q(d[t],W)*f[(g-t)%r] for t in range(r)) for g in range(r)) for f in fs)
    def expected():
        Z=[moment(i,tuple(D)) for i,D in enumerate(domains)]
        return sum(prod(Z[j][i][g] for j in range(s)) for i in range(len(fs)) for g in range(r))
    base=counts(r,B);W=sum(base);chi=Q(r*sum(x*x for x in base)-W*W,W*W)
    means=[sum(f)/r for f in fs]
    assert all(Q(x,W*W)<=Q(1+chi,r) for x in counts(r,slots))
    old=initial=expected();bound=sum(r*v*((1+chi)*v)**(s-1) for v in means)
    assert old<=bound<1
    if encoding is None:encoding=lambda v:v
    bits=max(map(encoding,range(1,r))).bit_length();steps=0
    for j in range(1,s):
        for i in range(len(slots)):
            for bit in reversed(range(bits)):
                S=domains[j][i];C=[tuple(x for x in S if (encoding(x)>>bit)&1==v) for v in (0,1)];C=[x for x in C if x]
                if len(C)<2:continue
                vals=[]
                for c in C:domains[j][i]=c;vals.append(expected())
                assert old*len(S)==sum(v*len(c) for v,c in zip(vals,C))
                ii=min(range(len(C)),key=lambda x:vals[x]);domains[j][i]=C[ii];old=vals[ii];steps+=1
    T=[tuple(F[0] for F in D) for D in domains];a=[sum(t)%r for t in T]
    assert old==sum(prod(f[(g-z)%r] for z in a) for f in fs for g in range(r))<1
    assert all(any(f[(g-z)%r]<1 for z in a) for f in fs for g in range(r))
    # Verify incidence representation on the original columns (scalar coordinates
    # are used only to check toy ground truth, not to perform any cryptographic work).
    for t,z in zip(T,a):
        incidence=[sum(t[2*i+j]==x for j in (0,1)) for i,F in enumerate(B) for x in F]
        assert sum(incidence)==2*m
        assert sum(c*x for c,x in zip(incidence,[x for F in B for x in F]))%r==z
    return T,a,initial,old,steps


def checked_shift(r, blocks, witness, preset=()):
    assert len(witness) == 2*len(blocks), 'wrong translation length'
    assert all(type(x) is int and x in blocks[j//2] for j, x in enumerate(witness)), 'wrong translation block'
    assert all(0 <= j < len(witness) and witness[j] == x for j, x in preset), 'wrong translation pin'
    return sum(witness) % r


def check_fixed_family(r, blocks, witnesses, shifts, size):
    assert len(witnesses) == len(shifts) == size, 'wrong translation family size'
    for witness, shift in zip(witnesses, shifts):
        assert checked_shift(r, blocks, witness) == shift, 'wrong translation sum'


def translation_witness_checks():
    r = 11
    original = [(1, 2), (7, 8), (4, 6)]
    retained_blocks = [original[i] for i in (0, 2)]
    witness = (1, 2, 4, 6)
    shift = sum(witness) % r
    pins = [(0, 1)]
    assert checked_shift(r, retained_blocks, witness, pins) == shift
    check_fixed_family(r, retained_blocks, [witness], [shift], 1)
    for invalid, reason in [((shift,), 'wrong translation length'),
                            ((1, 2, 7, 8, 4, 6), 'wrong translation length'),
                            ((1, 2, 7, 6), 'wrong translation block'),
                            ((4, 2, 1, 6), 'wrong translation block')]:
        with TestCase().assertRaisesRegex(AssertionError, reason):
            check_fixed_family(r, retained_blocks, [invalid], [shift], 1)
        with TestCase().assertRaisesRegex(AssertionError, reason):
            checked_shift(r, retained_blocks, invalid, pins)
    with TestCase().assertRaisesRegex(AssertionError, 'wrong translation pin'):
        checked_shift(r, retained_blocks, (2, 1, 4, 6), pins)
    with TestCase().assertRaisesRegex(AssertionError, 'wrong translation sum'):
        check_fixed_family(r, retained_blocks, [witness], [shift+1], 1)
    with TestCase().assertRaisesRegex(AssertionError, 'wrong translation family size'):
        check_fixed_family(r, retained_blocks, [], [shift], 1)
    print('Translation validation rejected scalar sums, deleted or misplaced slots, a changed pin, an incorrect shift and an incomplete family.')


def translation_checks():
    rng=Random(1804);steps=coverage=actual=0
    # Synthetic nonnegative penalties and retained subsets, with actual main-tuple translates.
    for r in (13,19,31,47):
        B=[rng.sample(range(1,r),max(3,r//3)) for _ in range(3)]
        d=counts(r,B);W=sum(d);chi=Q(r*sum(x*x for x in d)-W*W,W*W)
        fs=[]
        for _ in range(3):
            bad=set(rng.sample(range(r),2));fs.append([Q(3,2) if x in bad else Q(0) for x in range(r)])
        s=2
        while sum(r*(sum(f)/r)*((1+chi)*sum(f)/r)**(s-1) for f in fs)>=1:s+=1
        T,a,ini,end,n=select_translations(r,B,fs,s);steps+=n;coverage+=len(fs)*r
        check_fixed_family(r, B, T, a, s)
        assert all(any(f[(g-z)%r] < 1 for z in a) for f in fs for g in range(r))
        print(f'  Synthetic translation cover C{r}: {s} doubled-tuple shifts {a}; exact prefix identities {n}; all targets covered.')
    # Real affine-mask-generated tiny curve blocks, with relaxed finite thresholds.
    # The initial search for a strong h=2, k=1 fixture at T=8 did not find a
    # sufficiently low mean/peak ratio. Keep a positive h=1 deletion fixture
    # and a separate h=2 nondeletion fixture, using explicit recorded masks.
    G,add=curve(13,0,2);r=len(G);p=13
    fixtures=[(5,1,1,5,[[(4,2),(9,4)],[(1,1),(0,10)],[(10,11),(10,6)],[(0,10),(6,1)],[(0,7),(3,10)]]),
              (4,2,0,3,[[(0,8),(8,8)],[(11,9),(11,4)],[(9,7),(10,10)],[(6,7),(6,0)]])]
    for m,h,k,T0,masks in fixtures:
        q=Q((r-1)*T0,p)
        B=[[i for i in range(1,r) if (mask[G[i][1]%2][0]*G[i][0]+mask[G[i][1]%2][1])%p<T0] for mask in masks]
        arrays=phi_arrays(r,B,q,h,k);peak=max(max(v) for I,v in arrays);cut=peak*Q(99,100)
        subset_steps = translations = proper_targets = 0
        for I, values in arrays:
            blocks = [B[i] for i in I]
            f = [x/cut for x in values]
            main_counts = counts(r, blocks)
            mean = sum(f)/r
            density_bound = Q(r*sum(x*x for x in main_counts), sum(main_counts)**2)
            assert density_bound*mean < Q(4, 5)
            s = 2
            while r*mean*(density_bound*mean)**(s-1) >= 1:
                s += 1
            T, shifts, _, _, n = select_translations(
                r, blocks, [f], s, lambda x: (G[x][0]<<4) | G[x][1])
            check_fixed_family(r, [B[i] for i in I], T, shifts, s)
            for witness, shift in zip(T, shifts):
                point = None
                for x in witness:
                    point = add(point, G[x])
                assert point == G[shift]
            for g in range(r):
                j = next(j for j, z in enumerate(shifts) if f[(g-z)%r] < 1)
                target = (g-shifts[j]) % r
                ct = 0
                for X in product(*blocks):
                    if sum(X)%r != target:
                        continue
                    point = None
                    for x in X+T[j]:
                        point = add(point, G[x])
                    assert point == G[g]
                    ct += 1
                assert ct == main_counts[target] and ct > 0, 'empty or incomplete selected fiber'
                actual += ct
                coverage += 1
            if len(I) < m:
                proper_targets += r
            subset_steps += n
            translations += s
        if k == 1:
            assert proper_targets == 95, 'missing proper-subset cover checks'
        steps += subset_steps
        print(f'  Actual affine masks E/F13, order {r}, h={h}, k={k}, families={len(arrays)}, block sizes {[len(F) for F in B]}. Checked {translations} translations, {subset_steps} exact prefix identities and {proper_targets} proper-subset targets (relaxed finite threshold).')
    print(f'Translation selection total: {steps} exact prefix identities; {coverage} retained-subset/target cover checks; {actual} translated tuples verified by curve group law.')


def adaptive_select(r,B,requests,preset=(),encoding=None):
    """One doubled tuple, several targets; requests are (g, f_array)."""
    m=len(B);slots=[tuple(F) for F in B for _ in range(2)]
    domains=list(slots)
    for j,x in preset:
        assert x in slots[j];domains[j]=(x,)
    def expect():
        d=counts(r,domains);W=sum(d)
        return sum(Q(d[t],W)*f[(g-t)%r] for g,f in requests for t in range(r))
    old=initial=expect();assert initial<1
    if encoding is None:encoding=lambda x:x
    bits=max(map(encoding,range(1,r))).bit_length();steps=0
    for j in range(2*m):
        for bit in reversed(range(bits)):
            S=domains[j];C=[tuple(x for x in S if (encoding(x)>>bit)&1==v) for v in (0,1)];C=[c for c in C if c]
            if len(C)<2:continue
            vals=[]
            for c in C:domains[j]=c;vals.append(expect())
            assert old*len(S)==sum(len(c)*v for c,v in zip(C,vals))
            i=min(range(len(vals)),key=lambda x:vals[x]);domains[j]=C[i];old=vals[i];steps+=1
    T=tuple(D[0] for D in domains);a=sum(T)%r
    assert all(T[j]==x for j,x in preset)
    assert old==sum(f[(g-a)%r] for g,f in requests)<1
    return T,initial,old,steps


def adaptive_checks():
    rng=Random(3267);queries=steps=cases=0;maxpins=0
    for r in (19,31,47):
        B=[rng.sample(range(1,r),max(5,r//3)) for _ in range(4)]
        for _ in range(8):
            I=list(range(len(B)))
            if rng.randrange(2):I.pop(rng.randrange(len(I)))
            F=[B[i] for i in I];m=len(F)
            preset=((0,F[0][0]),(3,F[1][0])) if m>1 else ()
            slots=[tuple(A) for A in F for _ in range(2)]
            for j,x in preset:slots[j]=(x,)
            d=counts(r,slots);W=sum(d)
            # Arbitrary sparse bad sets, small perturbations of the penalties.
            req=[]
            for j in range(3):
                g=rng.randrange(r);bad=set(rng.sample(range(r),2));req.append((g,[Q(3,2) if x in bad else Q(1,100) for x in range(r)]))
            if sum(Q(d[t],W)*f[(g-t)%r] for g,f in req for t in range(r))>=1:continue
            T,ini,end,n=adaptive_select(r,F,req,preset);steps+=n;queries+=len(req);cases+=1;maxpins=max(maxpins,len(preset))
            shift = checked_shift(r, F, T, preset)
            assert all(f[(g-shift)%r]<1 for g,f in req)
    # Exact count-expansion check for the adaptive mean at a specified target,
    # on actual small affine masks, using the shared-pin energy identity.
    G,add=curve(13,0,2);r=len(G);p=13;T0=3;m=4;h=2;k=0
    masks=[[(0,8),(8,8)],[(11,9),(11,4)],[(9,7),(10,10)],[(6,7),(6,0)]]
    B=[[i for i in range(1,r) if (mask[G[i][1]%2][0]*G[i][0]+mask[G[i][1]%2][1])%p<T0] for mask in masks]
    q=Q((r-1)*T0,p);I=tuple(range(m));phi=phi_arrays(r,B,q,h,k)[0][1];mu=q**m/r
    slots=[tuple(F) for F in B for _ in range(2)];d=counts(r,slots);W=sum(d);cut=max(phi)*Q(99,100);f=[x/cut for x in phi];checks=0;branches=0
    for g in range(r):
        direct=sum(Q(d[t],W)*phi[(g-t)%r] for t in range(r));expanded=Q(0)
        for J in pins(I,h):
            E=shared_count_energy(r,fixed_tables(B),I,J)
            expanded+=q**(2*len(J))/(r*r*mu**2)*sum(Q(d[t],W)*E[(g-t)%r] for t in range(r));checks+=1
        assert direct==expanded
        if direct/cut>=1:continue
        T,ini,end,n=adaptive_select(r,B,[(g,f)],encoding=lambda x:(G[x][0]<<4)|G[x][1]);steps+=n;queries+=1;cases+=1
        a = checked_shift(r, B, T)
        assert f[(g-a)%r] == end < 1
        N=counts(r,B)[(g-a)%r];ct=0
        for X in product(*B):
            if sum(X)%r!=(g-a)%r:continue
            S=None
            for x in X+T:S=add(S,G[x])
            assert S==G[g];ct+=1
        assert ct==N;branches+=1
    assert branches>0
    print(f'Adaptive known-incidence selector: {cases} successful finite selections, {queries} target requests, {steps} exact prefix identities; up to {maxpins} pre-pinned translate slots; {checks} expanded shared-pin mean terms checked, {branches} real-curve branches (relaxed threshold).')


def rectangle_count(q,allowed,weight=None):
    z=[1]+[0]*(q-1)
    for j,A in enumerate(allowed):
        a=[0]*q
        for x in A:a[x]=(weight(x) if j==0 and weight is not None else 1)
        z=conv(z,a)
    return z[0]


def one_word(q,allowed):
    # Find a word of prescribed zero sum, only for the synthetic parity code.
    suf=[[1]+[0]*(q-1)]
    for A in reversed(allowed):
        a=[int(x in A) for x in range(q)];suf.append(conv(a,suf[-1]))
    if not suf[-1][0]:return None
    S=0;X=[]
    for i,A in enumerate(allowed):
        x=next(x for x in A if suf[len(allowed)-i-1][(-S-x)%q]);X.append(x);S+=x
    assert S%q==0
    return X


def check_rectangle_mass(num, den, independent):
    assert den > 0 and 0 <= num <= den, 'invalid rectangle counts'
    assert 2*num >= independent*den, 'rectangle probability below the lower bound'


def check_score_pair(q, labels, left, right, X, Y):
    scores = []
    for allowed, word in [(left, X), (right, Y)]:
        assert word is not None and len(word) == len(allowed) == len(labels), 'wrong witness length'
        assert all(type(x) is int and 0 <= x < q and x in A
                   for x, A in zip(word, allowed)), 'witness outside its allowed set'
        assert sum(word) % q == 0, 'witness outside the parity code'
        scores.append(sum(u[x] for u, x in zip(labels, word)) % q)
    assert scores[0] != scores[1], 'equal modular scores'


def cylinder_regressions():
    # Include equality at the lower bound and a positive count just below it.
    check_rectangle_mass(3, 12, Q(1, 2))
    with TestCase().assertRaisesRegex(AssertionError, 'rectangle probability below the lower bound'):
        check_rectangle_mass(2, 12, Q(1, 2))

    q = 3
    labels = [[0, 1, 2], [0, 0, 0]]
    left, right = [{0}, {0}], [{1}, {2}]
    X, Y = (0, 0), (1, 2)
    check_score_pair(q, labels, left, right, X, Y)
    for invalid, reason in [((0,), 'wrong witness length'),
                            ((0, 0, 0), 'wrong witness length'),
                            ((-3, 0), 'witness outside its allowed set'),
                            ((2, 1), 'witness outside its allowed set')]:
        with TestCase().assertRaisesRegex(AssertionError, reason):
            check_score_pair(q, labels, left, right, invalid, Y)
    with TestCase().assertRaisesRegex(AssertionError, 'witness outside the parity code'):
        check_score_pair(q, labels, left, [{1}, {0}], X, (1, 0))
    # These admissible parity words have integer scores 0 and 3, both zero in F3.
    with TestCase().assertRaisesRegex(AssertionError, 'equal modular scores'):
        check_score_pair(q, [range(q), range(q)], left, right, X, Y)
    print('The rectangle boundary passed and a smaller positive mass failed. Malformed, off-code and equal-residue score witnesses were rejected.')


def trimming_checks():
    rng=Random(1707);rects=modal=diffuse=0
    for h,q in [(16,3),(32,5),(64,7)]:
        m=h+1;U=set(range(q))
        tau = Q(1, 2**(h+14)*m*m)
        delta = 8*tau
        def weight(x):
            return delta.denominator + int(x==0)
        # Measure both possible one-coordinate distortions against the theorem's tolerance.
        full = rectangle_count(q, [U]*m, weight)
        for label in (0, 1):
            pinned = rectangle_count(q, [{label}] + [U]*(m-1), weight)
            assert abs(Q(q*pinned, full)-1) <= 8*tau, 'perturbation exceeds the theorem tolerance'
        for rep in range(40):
            p=rng.randrange(h//4+1);P=set(rng.sample(range(m),p))
            allowed=[{rng.randrange(q)} if i in P else set(U) for i in range(m)]
            E=Q(0)
            for i in rng.sample([i for i in range(m) if i not in P],m-p):
                lost=rng.randrange(q//2+1)
                if E+Q(lost,q)>Q(h,16):continue
                allowed[i]-=set(rng.sample(range(q),lost));E+=Q(lost,q)
            pins_only=[A if i in P else U for i,A in enumerate(allowed)]
            # Positive perturbation of the h-wise uniform parity law. Every
            # cylinder's relative probability changes by <=delta.
            den=rectangle_count(q,pins_only,weight);num=rectangle_count(q,allowed,weight)
            # Check exact independent Bonferroni bracket and the claimed lower bound.
            es=[Q(q-len(allowed[i]),q) for i in range(m) if i not in P]
            s=[Q(1)]+[Q(0)]*len(es)
            for e in es:
                for j in reversed(range(1,len(s))):s[j]+=e*s[j-1]
            t=h-p if (h-p)%2 else h-p-1
            assert t>=1
            p0=prod(1-e for e in es)
            check_rectangle_mass(num, den, p0)
            low=sum((-1)**j*s[j] for j in range(min(t,len(s)-1)+1))
            assert low<=p0 and low>=p0-(s[t+1] if t+1<len(s) else 0)
            assert low-3*delta*sum(s[:t+1])>0
            # The code is exactly h-wise uniform; perturbation weights obey the
            # claimed relative distortion, so a rectangle cannot be lost.
            rects+=1
        # Witnesses against low-modal score directions, including diffuse blocks.
        for rep in range(40):
            labels=[[0]*q for _ in range(m)];E=Q(0)
            for i in rng.sample(range(m),m):
                if rng.random()<.18:
                    k=rng.randrange(1,q);loss=Q(k,q)
                    if E+loss>Q(h,16):continue
                    for x in rng.sample(range(q),k):labels[i][x]=x+1
                    E+=loss
            # Modal classes and the chosen differing labels are classes in F_q.
            labels = [[v % q for v in row] for row in labels]
            if all(len(set(row)) == 1 for row in labels):
                labels[0][0] = 1
                E = max(E, Q(1, q))
            # Actual modal class chosen independently per block, including ties.
            classes=[];nonmodal=[]
            for u in labels:
                v=max(set(u),key=lambda z:(u.count(z),-z));C={x for x in U if u[x]==v}
                classes.append(C);nonmodal.append(Q(q-len(C),q))
            assert sum(nonmodal)<=E<=Q(h,16)
            bad=[i for i,e in enumerate(nonmodal) if e>Q(1,2)];diffuse+=bool(bad)
            j=next(i for i,C in enumerate(classes) if len(C)<q)
            x=next(iter(classes[j]));y=next(v for v in U if labels[j][v]!=labels[j][x])
            A=[{min(C)} if i in bad else set(C) for i,C in enumerate(classes)]
            left, right = A.copy(), A.copy()
            left[j], right[j] = {x}, {y}
            X, Y = one_word(q, left), one_word(q, right)
            check_score_pair(q, labels, left, right, X, Y)
            modal += 1
    # h-wise uniform abstract arrays do not force a single score relation:
    # two independent parity constraints on a two-coordinate alphabet.
    q=3;m=4;alphabet=list(product(range(q),repeat=2));X=[]
    for prefix in product(alphabet,repeat=m-1):X.append(prefix+(((-sum(x for x,y in prefix))%q,(-sum(y for x,y in prefix))%q),))
    for I in combinations(range(m),m-1):assert len({tuple(z[i] for i in I) for z in X})==len(alphabet)**(m-1)
    assert all(sum(z[0] for z in x)%q==sum(z[1] for z in x)%q==0 for x in X)
    print(f'All {rects} measured rectangle probabilities attained the lower bound or exceeded it. All {modal} admissible witness pairs had distinct scores in F_q ({diffuse} with diffuse blocks). The two-constraint OA control passed.')


def score_moments(q,allowed,labels,weight=None):
    N=[1]+[0]*(q-1);S=[0]*q;T=[0]*q
    for j,A in enumerate(allowed):
        NN=[0]*q;SS=[0]*q;TT=[0]*q
        for y in range(q):
            for x in A:
                z=(x+y)%q;v=labels[j][x];w=weight(x) if j==0 and weight else 1
                NN[z]+=w*N[y];SS[z]+=w*(S[y]+v*N[y]);TT[z]+=w*(T[y]+2*v*S[y]+v*v*N[y])
        N,S,T=NN,SS,TT
    return N[0],S[0],T[0]


def trimmed_rank_checks():
    rng=Random(5428);checks=varchecks=0
    for h,q in [(96,5),(128,7)]:
        m=h+1;U=set(range(q));delta=Q(1,2**(h+11)*m*m);err=Q(8,2**(h//2))
        assert h>=16+8*m.bit_length()
        for _ in range(12):
            P=set(rng.sample(range(m),rng.randrange(h//32+1)));A=[{rng.randrange(q)} if i in P else set(U) for i in range(m)];E=Q(0)
            for i in rng.sample([j for j in range(m) if j not in P],m-len(P)):
                k=rng.randrange(q//2+1)
                if E+Q(k,q)>Q(h,16)-len(P):continue
                A[i]-=set(rng.sample(range(q),k));E+=Q(k,q)
            # Exactly h-wise uniform additive parity law; positive tiny weight
            # perturbation preserves all cylinder ratios within 1 +/- delta.
            weight=lambda x:delta.denominator+int(x==0)
            den=rectangle_count(q,A,weight);assert den>0
            for _ in range(40):
                J=rng.sample([i for i in range(m) if i not in P],rng.choice((1,2)))
                B=[set(x) for x in A]
                for j in J:B[j]={rng.choice(sorted(A[j]))}
                ratio=Q(rectangle_count(q,B,weight),den)*prod(len(A[j]) for j in J)
                assert abs(ratio-1)<=err;checks+=1
            U0=[[rng.randrange(-5,6) for _ in range(q)] for _ in range(m)]
            N,S,T=score_moments(q,A,U0,weight);assert N==den
            V=Q(T,N)-Q(S,N)**2
            VI=sum(Q(sum(U0[i][x]**2 for x in A[i]),len(A[i]))-Q(sum(U0[i][x] for x in A[i]),len(A[i]))**2 for i in range(m))
            assert VI/2<=V<=2*VI;varchecks+=1
    # All labels distinct, and the original integer and circular width sums
    # exceed the score modulus. A low-loss trim makes the width sum < modulus.
    h=96;m=h+1;q=5;r=1009
    labels=[list(range(q)) for _ in range(m)];A=[set(range(q)) for _ in range(m)]
    for i in range(10):labels[i][-1]=r//2;A[i].remove(q-1)
    def width(u):
        z=sorted(set(x%r for x in u));return r-max([b-a for a,b in zip(z,z[1:])]+[r+z[0]-z[-1]])
    original=sum(width(u) for u in labels);trim=sum(max(u[x] for x in B)-min(u[x] for x in B) for u,B in zip(labels,A))
    E=sum(Q(q-len(B),q) for B in A);modal=m*Q(q-1,q)
    assert original>=r and trim<r and E<=Q(h,16) and modal>Q(h,16)
    N,S,T=score_moments(q,A,labels);assert N>0 and N*T>S*S
    # Independent min/max exact DP under the sum-zero code constraint.
    mn=[0]+[None]*(q-1);mx=mn[:]
    for i,B in enumerate(A):
        n=[None]*q;x=[None]*q
        for y in range(q):
            if mn[y] is None:continue
            for a in B:
                z=(y+a)%q;v=labels[i][a];lo=mn[y]+v;hi=mx[y]+v
                n[z]=lo if n[z] is None else min(n[z],lo);x[z]=hi if x[z] is None else max(x[z],hi)
        mn,mx=n,x
    assert 0<mx[0]-mn[0]<=trim<r
    print(f'Trimmed-cylinder rank controls: {checks} exact one/two-pin marginal ratios, {varchecks} exact variance bounds. Distinct-label example: original circular width {original}, trimmed width {trim} < {r}, loss {E}, modal distance {modal}; actual score extrema {mn[0]}, {mx[0]}.')


def numeric_guards():
    # Conservative bit bounds using b, not enormous-prime enumeration.
    # They test algebraic inequalities only; no finite fixture is an asymptotic seed.
    for b,m,h,k in [(2**16,128,16,2),(2**20,128,16,2),(2**21,256,24,3),(2**22,256,32,4),(2**24,512,48,6)]:
        a=(h+k+1)*b//(128*m);den=m-k-2*h;t0=3*den//4
        s=(4*b+a-1)//a;lo=Q(b-1+2*a,den);hi=Q(b+2*a,den)+3
        # m^5/q + 2^(2m+1)(1+2)/sqrt K + 21^m q'^t/r.
        size=5*m.bit_length()-lo;mom=2*m+3-a;girth=5*m+t0*(hi+1)-(b-1)
        extra=1+2*(h+14+2*m.bit_length());cover=b+m+s*(extra-a)
        good=max(size,mom,girth)<-12 and cover<0 and a>2*extra
        print(f'  b={b}, m={m}, h={h}, k={k}, shift count {s}: conservative seed/cover guard {good}.')


def optimized_mode_check():
    # Import the file in a fresh optimized interpreter to exercise its startup guard.
    result = subprocess.run([sys.executable, '-I', '-B', '-O', '-c',
                             'import runpy, sys; runpy.run_path(sys.argv[1])', __file__],
                            capture_output=True, text=True)
    assert result.returncode != 0 and not result.stdout
    assert result.stderr.splitlines()[-1] == 'RuntimeError: This verifier requires Python assertions.'
    print('The verifier rejected optimized Python before producing output.')


if __name__=='__main__':
    optimized_mode_check()
    two_slot_collision_check()
    translation_witness_checks()
    cylinder_regressions()
    shared_checks()
    mask_prefix_case(5,2,1,3,5,2,1,True)
    mask_prefix_case(7,0,3,4,6,3,1,False)
    translation_checks()
    adaptive_checks()
    trimming_checks()
    trimmed_rank_checks()
    numeric_guards()
