"""Exact checks of affine-mask energy and common translation covering.
Synthetic cyclic groups and two tiny elliptic curves; no DLP solver or contractor.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from functools import lru_cache
from math import prod
from random import Random


def conv(A, B):
    r = len(A)
    return [sum(A[j] * B[(i-j) % r] for j in range(r)) for i in range(r)]


def multi_conv(r, arrays):
    z = [1] + [0]*(r-1)
    for a in arrays:
        z = conv(z, a)
    return z


def counts(r, blocks):
    z = [1] + [0]*(r-1)
    for B in blocks:
        nxt = [0]*r
        for a in B:
            for i, x in enumerate(z):
                nxt[(i+a) % r] += x
        z = nxt
    assert sum(z) == prod(map(len, blocks))
    return z


def auto(r, B):
    h = [0]*r
    for x in B:
        for y in B:
            h[(y-x) % r] += 1
    return h


def curve(p, a, b):
    def add(P, Q0):
        if P is None: return Q0
        if Q0 is None: return P
        x, y = P; z, w = Q0
        if x == z and (y+w) % p == 0: return None
        t = ((3*x*x+a)*pow(2*y,-1,p) if P == Q0 else (w-y)*pow(z-x,-1,p)) % p
        X = (t*t-x-z) % p
        return X, (t*(x-X)-y) % p
    pts = [(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-b) % p == 0]
    G = [None]
    P = pts[0]
    for _ in pts: G.append(add(G[-1],P))
    assert set(G) == {None,*pts} and add(G[-1],P) is None
    return G, add


def bernoulli_checks():
    checked=0; selections=0
    for r, m, num, den in [(3,3,1,2),(5,2,1,3),(5,3,1,2),(7,2,2,3)]:
        theta=Q(num,den); q=(r-1)*theta; a=q-(r-2)*theta**2
        opts=[]
        for mask in range(1<<(r-1)):
            B=[i+1 for i in range(r-1) if mask>>i&1]
            w=num**len(B)*(den-num)**(r-1-len(B))
            opts.append((B,w))
        tot=0;denom=den**((r-1)*m)
        for selections0 in product(opts,repeat=m):
            B=[x[0] for x in selections0];w=prod(x[1] for x in selections0)
            d=counts(r,B);W=sum(d)
            v=r*sum(x*x for x in d)-W*W
            assert v>=0
            tot+=w*v;selections+=1
        assert Q(tot,denom)==(r-1)*a**m
        checked+=1
    print(f'Exact autocorrelation energy identity: {checked} independent Bernoulli cases, {selections} complete mask arrays.')


def affine_case(p,a,b,T,m,do_prefix):
    G,add=curve(p,a,b);r=len(G)
    allseeds=tuple(product(range(p),repeat=2));q=Q((r-1)*T,p)
    @lru_cache(None)
    def parity_data(par, seeds):
        ss=[{i for i in range(1,r) if G[i][1]%2==par and (A*G[i][0]+C)%p<T} for A,C in seeds]
        one=[sum(x in S for S in ss) for x in range(r)]
        pair=[[sum(x in S and y in S for S in ss) for y in range(r)] for x in range(r)]
        return one,pair
    @lru_cache(None)
    def one_block(sa,sb):
        A,C=parity_data(0,sa);B,D=parity_data(1,sb);den=len(sa)*len(sb)
        h=[0]*r; n1=0;pairchecks=0
        for x in range(1,r):
            px=G[x][1]%2
            n1+=A[x]*len(sb) if px==0 else B[x]*len(sa)
            for y in range(1,r):
                py=G[y][1]%2
                w=C[x][y]*len(sb) if px==py==0 else D[x][y]*len(sa) if px==py==1 else A[x]*B[y] if px==0 else B[x]*A[y]
                h[(y-x)%r]+=w
        assert h[0]==n1
        return [Q(x,den) for x in h],Q(n1,den),Q(sum(h),den)
    h,n1,n2=one_block(allseeds,allseeds)
    theta=Q(T,p);c=(r-2)*theta**2; aa=q-c
    assert h==[q]+[c]*(r-1)
    assert n1==q and n2-q*q==q*(1-theta)
    for k in range(1,9):
        Z=multi_conv(r,[h]*k)
        assert r*Z[0]-sum(Z)==(r-1)*aa**k
    if not do_prefix:return 8,0
    subsets=[I for k in range(m-1,m+1) for I in combinations(range(m),k)];J=len(subsets)
    domains=[[allseeds,allseeds] for _ in range(m)]
    def expectation():
        data=[one_block(*pair) for pair in domains]
        v=J*m**4*sum((n2-2*q*n1+q*q)/(q*q) for h,n1,n2 in data)
        for I in subsets:
            Z=multi_conv(r,[data[i][0] for i in I])
            e=r*Z[0]-sum(Z)
            assert e>=0
            v+=e/(r*q**len(I))
        return v
    initial=expectation();old=initial;steps=0
    assert initial==J*m**5*(1-theta)/q+sum(Q(r-1,r)*(aa/q)**len(I) for I in subsets)
    width=(p-1).bit_length()
    forced=0
    for i,par in product(range(m),range(2)):
        for coefficient in range(2):  # slope, then intercept, most significant bit first
            prefix=0
            for bit in reversed(range(width)):
                parent=domains[i][par]
                children=[tuple(seed for seed in parent if seed[coefficient] >> bit == 2*prefix+value)
                          for value in (0,1)]
                sizes=tuple(map(len,children))
                assert sum(sizes)==len(parent)
                if (i,par,coefficient,bit)==(0,0,0,width-1):
                    first_split=sizes
                    zero_slopes=1 << bit
                    assert sizes==(p*zero_slopes, p*(p-zero_slopes))
                if 0 in sizes:
                    chosen=int(sizes[0]==0)
                    forced+=1
                else:
                    domains[i][par]=children[0];left=expectation()
                    domains[i][par]=children[1];right=expectation()
                    assert old*len(parent)==left*sizes[0]+right*sizes[1]
                    chosen=int(right<left)
                    old=min(left,right)
                    steps+=1
                domains[i][par]=children[chosen]
                prefix=2*prefix+chosen
        assert len(domains[i][par])==1
    assert steps+forced==4*m*width
    blocks=[]
    for pair in domains:
        seeds=[S[0] for S in pair]
        blocks.append([i for i in range(1,r) if (seeds[G[i][1]%2][0]*G[i][0]+seeds[G[i][1]%2][1])%p<T])
    direct=J*m**4*sum((Q(len(B))-q)**2/(q*q) for B in blocks)
    for I in subsets:
        d=counts(r,[blocks[i] for i in I]);W=sum(d)
        direct+=Q(r*sum(x*x for x in d)-W*W,r)/q**len(I)
    assert direct==old<=initial
    # Verify the retained full-tuple values using the actual group law.
    D=[0]*r
    for X in product(*blocks):
        S=None
        for i in X:S=add(S,G[i])
        D[G.index(S)]+=1
    assert D==counts(r,blocks)
    print(f'Affine-seed coefficient bits E/F{p}, order {r}, m={m}: first split {first_split[0]}/{first_split[1]}, {steps} exact parent/child identities, {forced} forced bits; block sizes {list(map(len,blocks))}; direct potential equals final {float(old):.6g} (initial {float(initial):.6g}).')
    return 8,steps


def chi2(r,d):
    W=sum(d);assert W>0
    return Q(r*sum(x*x for x in d)-W*W,W*W)


def expected_potential(r,ds,domains,expand=False):
    h=len(domains);v=Q(0);den=prod(map(len,domains))
    for d in ds:
        W=sum(d)
        if not expand:
            z=0
            for g in range(r):
                z+=prod(sum((r*d[(g-a)%r]-W)**2 for a in A) for A in domains)
            v+=Q(z,den*W**(2*h))
        else:
            # 3^h terms, each an ordinary count of alpha_j independent
            # tuple occurrences with the same translated target per shift.
            coeff=[W*W,-2*r*W,r*r]
            z=0
            for alpha in product(range(3),repeat=h):
                count=sum(prod(sum(d[(g-a)%r]**k for a in A) for A,k in zip(domains,alpha)) for g in range(r))
                z+=prod(coeff[k] for k in alpha)*count
            v+=Q(z,den*W**(2*h))
    return v


def shift_select(r,ds,h,encoding=None,check_expand=False):
    # First shift is P (=1); all other shifts uniform nonidentity points.
    domains=[[1]]+[list(range(1,r)) for _ in range(h-1)]
    if encoding is None:encoding=lambda x:x
    bits=max(encoding(x) for x in range(1,r)).bit_length()
    initial=expected_potential(r,ds,domains);old=initial
    es=[chi2(r,d) for d in ds]
    assert initial <= sum(r*e**h*Q(r,r-1)**(h-1) for e in es)
    expansion_checks=0
    if check_expand:
        assert initial==expected_potential(r,ds,domains,True);expansion_checks+=1
    steps=0
    for j in range(1,h):
        for bit in reversed(range(bits)):
            S=domains[j]
            children=[[x for x in S if (encoding(x)>>bit)&1==v] for v in (0,1)]
            nonempty=[A for A in children if A]
            if len(nonempty)==1:continue
            vals=[]
            for A in nonempty:
                domains[j]=A
                q=expected_potential(r,ds,domains);vals.append(q)
                if check_expand:
                    assert q==expected_potential(r,ds,domains,True);expansion_checks+=1
            assert old*len(S)==sum(q*len(A) for q,A in zip(vals,nonempty))
            jmin=min(range(len(vals)),key=lambda i:vals[i]);domains[j]=nonempty[jmin];old=vals[jmin];steps+=1
        assert len(domains[j])==1
    shifts=sorted(set(A[0] for A in domains))
    assert old<=initial
    holes0=[sum(x==0 for x in d) for d in ds]
    holes1=[sum(all(d[(g-a)%r]==0 for a in shifts) for g in range(r)) for d in ds]
    if old<1:assert not any(holes1)
    return shifts,initial,old,holes0,holes1,steps,expansion_checks


def nontrivial_pair(r,rng):
    # Forbid zero by taking disjoint negative supports.  The other
    # coefficients are sufficiently flat for the two-shift guard.
    for _ in range(200):
        A=rng.sample(range(1,r),(r-1)//2)
        B=sorted(set(range(1,r))-{-x%r for x in A})
        d=counts(r,[A,B]); assert d[0]==0
        e=chi2(r,d)
        if Q(r*r,r-1)*e*e<1: return A,B,d
    raise AssertionError('No toy with the two-shift guard')


def translated_support(r, support, shifts):
    return {(g+t)%r for g in support for t in shifts}


def check_returned_cover(r, blocks, histogram, shifts, before, after):
    direct=[0]*r
    for values in product(*blocks):
        direct[sum(values)%r]+=1
    assert direct==histogram, f'C{r}: stale tuple counts for the returned blocks'
    assert 1 in shifts and all(0<t<r for t in shifts)
    support={g for g,count in enumerate(direct) if count}
    covered=translated_support(r,support,shifts)
    assert before==r-len(support)
    assert after==r-len(covered)==0, f'C{r}: returned shifts leave holes'
    return support


def check_common_cover(r, blocks, deleted, h, check_expand=False):
    m=len(blocks)
    subsets=[I for k in range(m-deleted,m+1) for I in combinations(range(m),k)]
    retained=[[blocks[i] for i in I] for I in subsets]
    ds=[counts(r,group) for group in retained]
    shifts,ini,end,before,after,steps,checks=shift_select(r,ds,h,check_expand=check_expand)
    assert end<=ini<1 and sum(before)>0
    supports={}
    double=[t+s for t,s in product(shifts,repeat=2)]
    for I,group,d,b,a in zip(subsets,retained,ds,before,after,strict=True):
        support=check_returned_cover(r,group,d,shifts,b,a)
        supports[I]=support
        # Check the case where both indexed copies of the translating block remain.
        assert len(translated_support(r,support,double))==r
    print(f'Common translation C{r}, m={m}, deletions={deleted}: {len(subsets)} retained subsets, total holes {sum(before)} -> {sum(after)}, shifts {shifts}, potential {float(ini):.6g} -> {float(end):.6g}, {steps} exact prefix identities.')
    return supports,steps,checks


def translation_checks():
    rng=Random(1864);st=0;ex=0
    for r in (19,31,101,257,1009):
        A,B,d=nontrivial_pair(r,rng)
        G=add=None
        if r==19:
            G,add=curve(13,0,2)
            assert len(G)==r
            encoding=lambda x:(G[x][0]<<4)|G[x][1]
        else:encoding=None
        T,ini,end,h0,h1,n,c=shift_select(r,[d],2,encoding,check_expand=(r<=31))
        check_returned_cover(r,[A,B],d,T,h0[0],h1[0])
        assert ini<1 and end<1 and h0[0]>=1
        assert len(T)*len(A)*len(B)<r*r
        if G is not None:
            values={add(add(G[x],G[y]),G[t]) for x,y,t in product(A,B,T)}
            assert values==set(G)
        print(f'Two-shift coverage C{r}: holes {h0[0]} -> 0, shifts {T}, tuple product {len(T)*len(A)*len(B)} < {r*r}, potential {float(ini):.6g} -> {float(end):.6g}, {n} exact prefix identities.')
        if r==257:
            # Same block sizes as the fixture, paired with its stale histogram.
            replacement=[137*i % r for i in range(1,129)]
            try:
                check_returned_cover(r,[replacement,replacement],d,T,h0[0],h1[0])
            except AssertionError as error:
                assert str(error)=='C257: stale tuple counts for the returned blocks'
            else:
                raise AssertionError('Accepted stale counts for replacement blocks')
            print('C257 stale-count control: rejected the original histogram for two replacement blocks of size 128.')
        st+=n;ex+=c
    # Common translating block for all retained subsets.  The toy full
    # product is not required to fit the asymptotic r^2 guard here.
    for r,m,deleted,h in [(19,3,1,4),(31,4,2,4),(101,4,2,4)]:
        A,B,_=nontrivial_pair(r,rng);blocks=[A,B]+[rng.sample(range(1,r),(r-1)//2) for _ in range(m-2)]
        _,n,c=check_common_cover(r,blocks,deleted,h,check_expand=(r==19))
        st+=n;ex+=c
    r=23
    group=set(range(r))
    blocks=[sorted(group-{0,1,2}), sorted(group-{0,1,3,4}), [1,2]]
    supports,n,c=check_common_cover(r,blocks,deleted=1,h=2,check_expand=True)
    left,right=supports[0,2],supports[1,2]
    assert group-left=={2,3} and group-right=={2,5}
    assert group-translated_support(r,left,[1,2])=={4}
    assert translated_support(r,right,[1,2])==group
    assert translated_support(r,left,[1,4])==group
    assert group-translated_support(r,right,[1,4])=={6}
    print('C23 simultaneous repair: initial holes {2,3} and {2,5}; shifts [1,2] miss 4 in the first subset, [1,4] miss 6 in the second.')
    st+=n;ex+=c
    print(f'Translation selection totals: {st} exact parent/child identities; {ex} direct-vs-expanded count comparisons.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Run with assertions enabled')
    bernoulli_checks()
    k,s=affine_case(5,2,1,4,4,True)
    for T in (1,3,7):
        k0,s0=affine_case(13,0,2,T,4,False);k+=k0;s+=s0
    print(f'Actual affine masks: {k} exact autocorrelation/energy identities across two curves, {s} seed prefix checks.')
    translation_checks()
