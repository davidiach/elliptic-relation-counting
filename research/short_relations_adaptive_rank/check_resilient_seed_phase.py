"""Exact finite checks; no implementation of the assumed tensor contractor."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, factorial, prod
from random import Random


def rank(a, p):
    a = [list(row) for row in a]
    if not a:
        return 0
    n = len(a[0]); k = 0
    for j in range(n):
        piv = next((i for i in range(k, len(a)) if a[i][j] % p), None)
        if piv is None:
            continue
        a[k], a[piv] = a[piv], a[k]
        z = pow(a[k][j] % p, -1, p)
        a[k] = [(x*z) % p for x in a[k]]
        for i in range(k+1, len(a)):
            z = a[i][j] % p
            a[i] = [(x-z*y) % p for x,y in zip(a[i],a[k])]
        k += 1
        if k == len(a):
            break
    return k


def partition(xs):
    seen = {}
    return tuple(seen.setdefault(x, len(seen)) for x in xs)


def partition_rank(patterns, h, p):
    columns = []
    for pat in patterns:
        for j in range(1+max(pat)):
            columns.append([int(pat[i] == j)-int(pat[-1] == j) for i in range(h-1)])
    return rank(list(zip(*columns)), p)


def polynomial_masks():
    tests = 0; evals = 0
    for p,h in [(5,2),(5,3),(7,3),(5,4)]:
        seeds = list(product(range(p), repeat=h))
        for xs in combinations(range(p), h):
            # For each fixed nonzero leading coefficient, the lower h coefficients
            # are a bijection onto h evaluations. Check all possible lower words.
            for lead in range(1,p):
                seen = Counter(tuple((lead*x**h+sum(c*x**j for j,c in enumerate(cs)))%p
                                     for x in xs) for cs in seeds)
                assert len(seen) == p**h and set(seen.values()) == {1}
                tests += 1; evals += len(seeds)
    print(f'Polynomial hash: {tests} fixed-leading-coefficient evaluation maps, {evals} seeds checked, all bijections.')


def collision_partitions():
    tested=0; equal_cases=0
    for p,m,h in [(5,2,2),(5,3,2),(5,2,3),(5,3,3),(5,2,4),(7,2,3)]:
        counts=Counter()
        tuples=list(product(range(1,p),repeat=m))
        by=defaultdict(list)
        for t in tuples:
            by[sum(t)%p].append(t)
        for ts in by.values():
            for rows in product(ts,repeat=h):
                pats=tuple(partition(x) for x in zip(*rows))
                counts[pats]+=1
        for pats,c in counts.items():
            sizes=[1+max(x) for x in pats]
            R=partition_rank(pats,h,p)
            assert R>=max(sizes)-1 and c<=p**(sum(sizes)-R)
            tested+=1;equal_cases+=c
        # Compare exact expectation by inclusion-count expansion to rank bound.
        theta=Q(2,5); qprime=p*theta
        exact=sum(c*theta**sum(1+max(t) for t in pats) for pats,c in counts.items())
        upper=sum(qprime**sum(1+max(t) for t in pats)/p**partition_rank(pats,h,p) for pats in counts)
        assert exact<=upper
    print(f'Collision partitions: {tested} patterns, {equal_cases} equal-sum tuple arrays, all rank/count bounds passed.')


def conv(a,b):
    r=len(a);out=[0]*r
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:out[(i+j)%r]+=x*y
    return out


def counts(blocks,r):
    d=[1]+[0]*(r-1)
    for B in blocks:
        a=[0]*r
        for x in B:a[x]+=1
        d=conv(d,a)
    return d


def curve_cycle(p,a,c):
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        s=((3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
        X=(s*s-x-u)%p;return X,(s*(x-X)-y)%p
    P=next((x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-c)%p==0)
    G=[None]
    while True:
        z=add(G[-1],P)
        if z is None:break
        G.append(z)
    return G


def seed_selection():
    rng=Random(4483);fixtures=0; checks=0
    # Restricted, explicit seed libraries only. Conditional minimization is exact
    # for these libraries; the asymptotic theorem uses all polynomial coefficients.
    G=curve_cycle(13,0,2); assert len(G)==19
    for r,m,h,k,libs in [(19,4,3,1,4),(31,4,4,1,4),(43,5,3,1,3)]:
        options=[]
        for i in range(m):
            B=[]
            for _ in range(libs):
                if r==19:
                    p=13;T=4
                    cs=[[rng.randrange(p) for _ in range(h)]+[rng.randrange(1,p)] for _ in range(2)]
                    s=tuple(j for j,(x,y) in enumerate(G[1:],1) if sum(v*x**k for k,v in enumerate(cs[y%2]))%p<T)
                else:s=tuple(sorted(rng.sample(range(1,r),rng.randrange(4,9))))
                B.append(s)
            options.append(B)
        subsets=[I for size in range(m-k,m+1) for I in combinations(range(m),size)]
        mu={I:prod(Q(sum(len(B) for B in options[i]),libs) for i in I) for I in subsets}
        def potential(seed):
            S=Q(0)
            for I in subsets:
                d=counts([options[i][seed[i]] for i in I],r);W=sum(d);C=sum(x**h for x in d)
                ex=r**(h-1)*C-W**h
                assert ex>=0 and max(d)**h<=C
                S+=(Q(W,1)/mu[I]-1)**2+Q(ex,1)/mu[I]**h
            return S
        vals={s:potential(s) for s in product(range(libs),repeat=m)}
        total=len(vals);mean=sum(vals.values(),Q(0))/total
        pref=()
        for i in range(m):
            def avg(j):
                vv=[v for s,v in vals.items() if s[:i+1]==pref+(j,)]
                return sum(vv,Q(0))/len(vv)
            pref+=min(range(libs),key=avg),
        assert vals[pref]<=mean
        for I in subsets:
            d=counts([options[i][pref[i]] for i in I],r);W=sum(d)
            if W:
                chi=Q(r*sum(x*x for x in d),W*W)-1
                assert chi>=0
                # Support deficit <= chi by Cauchy-Schwarz, and by zero contributions.
                assert Q(d.count(0),r)<=chi
            checks+=1
        fixtures+=1
        print(f'Joint seed fixture r={r}: {total} seeds, {len(subsets)} retained subsets, mean {float(mean):.6g}, selected {float(vals[pref]):.6g}.')
    print(f'Joint potential: {fixtures} finite seed libraries, {checks} retained-subset checks.')


def phase_tests():
    # Rational approximation to 2*pi within 5e-10 (proof uses certified pi interval).
    omega=Q(6283185307,10**9); delta=Q(1,4096); k=4
    def poly(x):
        return sum((-1)**(j+1)*(omega*x)**(2*j)/factorial(2*j) for j in range(1,k+1))-delta
    assert all(poly(Q(i))<=0 for i in range(-20,21))
    rng=Random(2239); cases=0; pairs=0; extracted=0
    for r in [19,31,43,61]:
        for _ in range(8):
            blocks=[rng.sample(range(1,r),4) for _ in range(2)]
            labels=[[rng.randrange(2) for _ in B] for B in blocks]
            loc=[0,0,1,1]
            tuples=[]
            for t in product(range(4),repeat=4):
                g=sum(blocks[i][v] for i,v in zip(loc,t))%r
                y=sum(labels[i][v] for i,v in zip(loc,t))
                tuples.append((t,g,y))
            s=rng.randrange(r); c=r//2
            rows=[(X+Y, y-x-c) for X,g,x in tuples for Y,h,y in tuples if (h-g)%r==s]
            if not rows:continue
            weights=[poly(Q(d,r)) for _,d in rows]
            # Direct exact sum vs moment expansion (coefficients applied after counts).
            moments=[sum(d**(2*j) for _,d in rows) for j in range(1,k+1)]
            value=sum((-1)**(j+1)*omega**(2*j)*moments[j-1]/(factorial(2*j)*r**(2*j)) for j in range(1,k+1))-len(rows)*delta
            assert value==sum(weights,Q(0)) and value>0
            assert all(w<=0 for (_,d),w in zip(rows,weights) if d%r==0)
            rr=list(zip(rows,weights))
            for i in range(8):
                rr=next([v for v in rr if v[0][0][i]==a] for a in range(4)
                        if sum((v[1] for v in rr if v[0][0][i]==a),Q(0))>0)
            assert len(rr)==1 and rr[0][0][1]%r
            cases+=1;pairs+=len(rows);extracted+=1
    print(f'Rational phase polynomial: {cases} low-variance toy fixtures, {pairs} relations, {extracted} positive-mass witnesses extracted.')


def concentrated_obstruction():
    checks=0
    for r,d in [(101,4),(257,8),(521,16),(1009,32),(65537,64)]:
        nodes=[Q(k,2)-Q(k%2,2*r)-d//4 for k in range(d+1)]
        w=[1/prod(abs(x-y) for j,y in enumerate(nodes) if j!=i) for i,x in enumerate(nodes)]
        even=sum(w[::2]);odd=sum(w[1::2]);assert even==odd
        pe={i:w[i]/even for i in range(0,d+1,2)};po={i:w[i]/odd for i in range(1,d+1,2)}
        for i,x in enumerate(nodes):
            assert (x.denominator==1)==(i%2==0)
            assert (r*x).denominator==1
            b=Q(comb(d,i),2**(d-1))
            a=(pe if i%2==0 else po)[i]
            assert b/2<=a<=2*b
        for j in range(d):
            assert sum(p*nodes[i]**j for i,p in pe.items())==sum(p*nodes[i]**j for i,p in po.items())
            checks+=1
        print(f'Concentrated obstruction: r={r}, d={d}, exact moment agreement through degree {d-1}, parity-binomial density ratio in [1/2,2].')
    print(f'Concentrated obstruction: {checks} exact moment identities.')

def punctured_ranks(B,r):
    n=sum(map(len,B)); ix=[];vals=[]; off=0
    for block in B:
        ix.append(range(off,off+len(block)));off+=len(block);vals+=block
    def get(I):
        d=[0]*r;v=[[0]*n for _ in range(r)]
        for X in product(*(ix[i] for i in I)):
            g=sum(vals[j] for j in X)%r;d[g]+=1
            for j in X:v[g][j]+=1
        return d,v
    m=len(B);d,v=get(range(m));pun=[get([j for j in range(m) if j!=i]) for i in range(m)]
    if any(min(a)==0 or max(a)>=r for a in [d]+[a for a,_ in pun]):return None
    rows=[];sums=[]
    for i,J in enumerate(ix):
        a,b=pun[i]
        for j in J:
            R=[]
            for g in range(r):
                x=(g-vals[j])%r; rr=[(d[g]*(b[x][k]+a[x]*int(k==j))-a[x]*v[g][k])%r for k in range(n)]
                assert sum(x*y for x,y in zip(rr,vals))%r==0
                assert all(sum(rr[k] for k in J1)%r==0 for J1 in ix)
                R.append(rr)
            rows+=R;sums.append([sum(x)%r for x in zip(*R)])
    rk=rank(rows,r);sr=rank(sums,r);assert rk==n-m-1
    return rk,sr,max(d),min(d),[min(a) for a,b in pun]


def punctured_mean_checks():
    rng=Random(774); tested=0
    for r,m,k in [(19,3,6),(31,3,8),(43,3,9),(61,3,10),(101,4,8)]:
        for _ in range(120):
            blocks=[sorted(rng.sample(range(1,r),k)) for j in range(m)]
            a=punctured_ranks(blocks,r)
            if a is not None:
                tested+=1
    B=[[2,14,39,61,72,82,85,89],[3,4,17,24,33,42,44,63],[1,3,7,18,32,41,42,44],[11,37,45,54,58,72,81,87]]
    a=punctured_ranks(B,101)
    assert a==(27,26,53,28,[1,1,1,1])
    print(f'Punctured means: {tested} full-support no-wrap fixtures; complete rank passed. C101 regression: full rank {a[0]}, summed rank {a[1]}, full counts {a[3]}..{a[2]}, minimum puncture counts {a[4]}.')

if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    polynomial_masks();collision_partitions();seed_selection();phase_tests();concentrated_obstruction();punctured_mean_checks()
