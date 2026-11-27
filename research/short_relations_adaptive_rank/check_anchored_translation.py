"""Finite exact tests for anchored translation certificates, not a contractor.
Python 3.10+. Small cyclic groups, and an actual order-19 curve.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random
from math import prod


def rank(rows,r,n):
    B={}
    for a in rows:
        u=[x%r for x in a]
        for j in sorted(B):
            if u[j]:
                c=u[j];u=[(x-c*y)%r for x,y in zip(u,B[j])]
        k=next((j for j,x in enumerate(u) if x),None)
        if k is not None:
            c=pow(u[k],-1,r);B[k]=[c*x%r for x in u]
        if len(B)==n:break
    return len(B)


def dot(x,y):return sum(a*b for a,b in zip(x,y))


def curve_cycle(p,a,c,r):
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        s=(3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p)
        X=(s*s-x-u)%p;return X,(s*(x-X)-y)%p
    pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-c)%p==0]
    P=pts[0];G=[None]
    while True:
        Q=add(G[-1],P)
        if Q is None:break
        G.append(Q)
    assert len(G)==r
    return G,add


def curve13():return curve_cycle(13,0,2,19)


def fixture(r,blocks,target=4,elliptic=False):
    values=[a for B in blocks for a in B]+[1,target];N=len(values)
    ix=[];i=0
    for B in blocks:
        ix.append(range(i,i+len(B)));i+=len(B)
    if elliptic:
        p=13 if r==19 else 29
        G,add=curve13() if r==19 else curve_cycle(29,4,9,37);lookup={x:i for i,x in enumerate(G)}
        def value(X):
            P=None
            for j in X:P=add(P,G[values[j]])
            return lookup[P]
        def disp(y,x):
            P=G[x];return lookup[add(G[y],None if P is None else (P[0],-P[1]%p))]
    else:
        value=lambda X:sum(values[j] for j in X)%r
        disp=lambda y,x:(y-x)%r
    T=[(X,value(X)) for X in product(*ix)];W=len(T);d=Counter(g for X,g in T);C=sum(c*c for c in d.values())
    H=[]
    # All positive-slot pins, for the pointwise upper-bound and support selector.
    for k,B in enumerate(ix):
        for j in B:
            n=[0]*r;v=[[0]*N for _ in range(r)]
            for Y,y in T:
                if Y[k]!=j:continue
                for X,x in T:
                    t=disp(y,x);n[t]+=1
                    for h in X:v[t][h]-=1
                    for h in Y:v[t][h]+=1
            assert sum(n)==W*W//len(B)
            assert len(B)*max(n)**2<=C*C
            assert all((dot(row,values)-t*c)%r==0 for t,(row,c) in enumerate(zip(v,n)))
            H.append((n,v,k,j))
    return values,H,W,C


def row_at(n,v,r,t,s,point=None):
    q=(t+s)%r;c=n[t]*n[q];row=[n[t]*y-n[q]*x for x,y in zip(v[t],v[q])]
    if point is not None:row[point]-=c;rhs=0
    else:rhs=c*s
    return row,rhs


def tests():
    rng=Random(982);cases=0;directions=0;minimum=Fraction(1);sum_ranks=[];density_min=Fraction(1)
    fixtures=[(19,[(1,2,3),(1,4,7,10,13,16)],True),
              (37,[list(range(1,7)),[1+6*j for j in range(6)]],False),
              (101,[list(range(1,11)),[1+10*j for j in range(10)]],False),
              (127,[list(range(1,8)),[1+7*j for j in range(6)],[1+42*j for j in range(3)]],False),
              (257,[list(range(1,17)),[1+16*j for j in range(16)]],False)]
    for r,B,ec in fixtures:
        vals,H,W,C=fixture(r,B,r-2,ec);N=len(vals)
        eps=Fraction(r*C,W*W)-1
        for n,v,_,_ in H:
            mu=Fraction(sum(n),r);chi=sum((Fraction(c,1)/mu-1)**2 for c in n)/r
            assert chi<=eps
        n,v,_,_=H[0];good=[0<c<r for c in n];assert all(c<r for c in n)
        delta=Fraction(good.count(False),r);assert delta<=Fraction(1,256)
        den=sum(n);chi=Fraction(r*sum(c*c for c in n),den*den)-1
        density_min=min(density_min,1-delta)
        # All shifts for finite checks: zero bias. This is NOT the asymptotic list.
        shifts=list(range(r));families=[(s,None) for s in shifts]+[(s,j) for j,s in enumerate(vals)]
        def equations():
            for s,j in families:
                for t in range(r):
                    a,c=row_at(n,v,r,t,s,j);assert (dot(a,vals)-c)%r==0;yield a
        assert rank(equations(),r,N)==N
        # Is the simplest all-displacement sum enough? Record, not assert.
        sums=[]
        for s,j in families:
            z=[0]*N
            for t in range(r):
                a,c=row_at(n,v,r,t,s,j);z=[x+y for x,y in zip(z,a)]
            sums.append(z)
        sr=rank(sums,r,N);sum_ranks.append((r,sr,N))
        short=sorted({sgn*pow(2,k,r)%r for k in range(r.bit_length()) for sgn in (-1,1)})
        assert rank((row_at(n,v,r,0,s,j)[0] for s,j in [(s,None) for s in short]+[(s,j) for j,s in enumerate(vals)]),r,N)==N
        # Include wrong dedicated columns, near-correct and broadly wrong tables.
        candidates=[]
        for j in range(N):
            y=vals.copy();y[j]=(y[j]+1)%r;candidates.append(y)
        candidates += [[rng.randrange(r) for _ in range(N)] for _ in range(25)]
        candidates += [[(k*x+c)%r for x in vals] for k,c in [(0,0),(2,0),(1,1),(r-1,0)]]
        for y in candidates:
            if y==vals:continue
            f=[dot(a,y)*pow(c,-1,r)%r if ok else 0 for a,c,ok in zip(v,n,good)]
            errors=[(z-t)%r for t,z in enumerate(f)]
            def bad(s,expected):
                return sum(good[t] and good[(t+s)%r] and (errors[(t+s)%r]-errors[t]-expected)%r!=0 for t in range(r))
            best=max([bad(s,0) for s in shifts]+[bad(s,(y[j]-s)%r) for j,s in enumerate(vals)])
            assert 4*best>=r
            # Cross-check a direct row against the error test.
            for _ in range(3):
                s,j=rng.choice(families);t=rng.randrange(r);a,c=row_at(n,v,r,t,s,j)
                expected=0 if j is None else (y[j]-s)%r
                assert bool((dot(a,y)-c)%r)==bool(good[t] and good[(t+s)%r] and (errors[(t+s)%r]-errors[t]-expected)%r)
            minimum=min(minimum,Fraction(best,r));directions+=1
        cases+=1
        print(f'r={r}: full rank {N}/{N}, summed-family rank {sr}/{N}, anchor chi-square={float(chi):.6g}, all {len(H)} pin bounds passed.')
    print(f'Anchored translations: {cases} fixtures (1 elliptic), {directions} wrong candidates, minimum best violation fraction {float(minimum):.6f}, minimum anchor support {float(density_min):.6f}.')
    print('All-displacement sums are an experimental fast path, not a theorem:',sum_ranks)


def support_selector(counts,bundles,r):
    # Each constraint is (count-array index, translation). Use exact rational averages.
    means=[Fraction(sum(a),r) for a in counts]
    def potential(t,b):
        return sum((Fraction(counts[i][(t+s)%r],1)/means[i]-1)**2 for i,s in b)
    def choose(rem):
        S=list(range(r)); vals={t:sum(potential(t,b) for b in rem) for t in S}
        # Prefix fixing for canonical integer encodings; direct enumeration only in toy test.
        for bit in range((r-1).bit_length()-1,-1,-1):
            B=[[t for t in S if (t>>bit)&1==k] for k in (0,1)]
            B=[x for x in B if x]
            S=min(B,key=lambda x:sum(vals[t] for t in x)/len(x))
        assert len(S)==1 and vals[S[0]]<=sum(vals.values())/r
        return S[0]
    rem=list(bundles);chosen=[]
    while rem:
        t=choose(rem);chosen.append(t)
        new=[b for b in rem if not all(0<counts[i][(t+s)%r]<r for i,s in b)]
        assert len(new)<len(rem)
        rem=new
    assert all(any(all(0<counts[i][(t+s)%r]<r for i,s in b) for t in chosen) for b in bundles)
    return chosen


def selector_tests():
    r=257;counts=[]
    # Deliberate zeros: abstract dense count arrays, not claimed elliptic masks.
    for j in range(4):
        a=[17]*r;a[(13*j+7)%r]=0;counts.append(a)
    bundles=[[(j,s),(0,(s+9)%r)] for j in range(4) for s in range(r)]
    L=support_selector(counts,bundles,r);assert len(L)<=3
    print(f'Support selector: {len(bundles)} translated two-count bundles with genuine zeros covered by {len(L)} starts: {L}.')
    r=19;vals,H,W,C=fixture(r,[(1,2,3),(1,4,7,10,13,16)],17,True)
    counts=[h[0] for h in H];bundles=[[(i,s),(0,(s+1)%r)] for i in range(len(H)) for s in range(r)]
    L=support_selector(counts,bundles,r);assert L
    print(f'Elliptic support selector: {len(bundles)} bundles covered by {len(L)} starts.')


def summed_obstruction():
    r=37;B=[[1,3,5,10,18,19,27,32,34,36],[18,19,21,31]]
    vals,H,W,C=fixture(r,B,r-2,True);N=len(vals);n,v,*_=H[0]
    F=[(s,None) for s in range(r)]+[(s,j) for j,s in enumerate(vals)]
    sums=[]
    for s,j in F:
        z=[0]*N
        for t in range(r):
            a,c=row_at(n,v,r,t,s,j);z=[x+y for x,y in zip(z,a)]
        sums.append(z)
    assert 0<min(n)<=max(n)<r
    assert rank(sums,r,N)==N-1
    assert rank((row_at(n,v,r,t,s,j)[0] for s,j in F for t in range(r)),r,N)==N
    # Choose one valid doubled decomposition per group point, in indexed-log order.
    G,add=curve_cycle(29,4,9,37);lookup={x:i for i,x in enumerate(G)}
    ix=[range(len(B[0])),range(len(B[0]),len(B[0])+len(B[1]))]
    T=[]
    for X in product(*ix):T.append((add(G[vals[X[0]]],G[vals[X[1]]]),X))
    by={}
    for g,X in T:
        for h,Y in T:
            s=lookup[add(g,h)]
            if s not in by:
                row=[0]*N
                for j in X+Y:row[j]+=1
                by[s]=row
    assert len(by)==r
    E=[[int(i==N-2) for i in range(N)]]
    for j,s in enumerate(vals):
        row=[-x for x in by[s]];row[j]+=1;E.append(row)
    assert rank(sums+E,r,N)==N-1
    short=sorted({sgn*pow(2,k,r)%r for k in range(r.bit_length()) for sgn in (-1,1)})
    rr=rank((row_at(n,v,r,0,s,j)[0] for s,j in [(s,None) for s in short]+[(s,j) for j,s in enumerate(vals)]),r,N)
    print(f'Obstruction on E/F29, order 37: counts {min(n)}..{max(n)}, summed rank {N-1}/{N}, unchanged with chosen explicit rows, complete stack {N}/{N}, one-start short-list rank {rr}/{N}.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    tests();selector_tests();summed_obstruction()
