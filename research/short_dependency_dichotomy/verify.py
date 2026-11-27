#!/usr/bin/env python3
"""Exact finite checks of short-dependency certificates and planted sparse spans.
These are abstract vector-space experiments, not cryptographic instances.
"""
from itertools import product, combinations
from math import comb
from random import Random
from collections import Counter

if not __debug__:
    raise RuntimeError('Assertions must be enabled.')


def add(x,y,q): return tuple((a+b)%q for a,b in zip(x,y))
def scale(x,c,q): return tuple(a*c%q for a in x)
def sub(x,y,q): return add(x,scale(y,-1,q),q)
def zero(d): return (0,)*d


def rank(rows,q):
    piv={}
    for a in rows:
        a=[x%q for x in a]
        for j,b in sorted(piv.items()):
            if a[j]:
                c=a[j];a=[(v-c*z)%q for v,z in zip(a,b)]
        i=next((i for i,v in enumerate(a) if v),None)
        if i is not None:
            c=pow(a[i],-1,q);piv[i]=[v*c%q for v in a]
    return len(piv)


def weight(x): return sum(a!=0 for a in x)
def ball(K,w,q): return sum(comb(K,j)*(q-1)**j for j in range(min(K,w)+1))


def sparse_patterns(K,w,q):
    for j in range(min(K,w)+1):
        for I in combinations(range(K),j):
            for a in product(range(1,q),repeat=j): yield tuple(zip(I,a))


def dense(c,K):
    x=[0]*K
    for i,a in c:x[i]=a
    return tuple(x)


def eval_sparse(H,c,q):
    d=len(H[0]);out=[0]*d
    for i,a in c:
        for j,x in enumerate(H[i]):out[j]=(out[j]+a*x)%q
    return tuple(out)


def lincomb(B,a,q):
    return tuple(sum(t*x for t,x in zip(a,v))%q for v in zip(*B))


def witnesses(H,w,q):
    W={}
    for c in sparse_patterns(len(H),w,q): W.setdefault(eval_sparse(H,c,q),dense(c,len(H)))
    return W


def short_rank(H,w,q):
    return rank((dense(c,len(H)) for c in sparse_patterns(len(H),3*w,q) if eval_sparse(H,c,q)==zero(len(H[0]))),q)


def checked_witness(H, W, target, w, q):
    word = tuple(v % q for v in W[target])
    assert len(word) == len(H), 'witness has the wrong length'
    terms = [(i, v) for i, v in enumerate(word) if v]
    assert len(terms) <= w, 'witness exceeds the weight limit'
    assert eval_sparse(H, terms, q) == target, 'witness has the wrong syndrome'
    return word


def defect_certificate(H,B,W,w,q):
    ell=len(B);K=len(H);T=list(product(range(q),repeat=ell));Z=zero(ell)
    units=[tuple(int(i==j) for i in range(ell)) for j in range(ell)]
    f={a:checked_witness(H,W,lincomb(B,a,q),w,q) for a in T};assert f[Z]==zero(K)
    defects=[]
    for a in T:
        supp=[i for i,v in enumerate(a) if v]
        if not supp:continue
        i=supp[-1];part=scale(units[i],a[i],q)
        if len(supp)==1:
            c=sub(f[a],scale(f[units[i]],a[i],q),q)
        else:
            rest=sub(a,part,q)
            c=sub(f[a],add(f[rest],f[part],q),q)
        assert weight(c)<=3*w
        assert eval_sparse(H,tuple((j,v) for j,v in enumerate(c) if v),q)==zero(len(B[0]))
        defects.append(c)
    e=rank(defects,q);f_rank=rank(f.values(),q)
    assert f_rank==ell+e
    assert q**ell<=ball(f_rank,w,q)
    assert e<=short_rank(H,w,q) if K<=8 else True
    return e,len(defects)


def projective_defect_certificate(H,B,W,w,q):
    ell=len(B);K=len(H);P=canonical_lines(q,ell)
    units=[tuple(int(i==j) for i in range(ell)) for j in range(ell)]
    f={a:checked_witness(H,W,lincomb(B,a,q),w,q) for a in P}
    def lift(a):
        if not any(a):return zero(K)
        c=next(v for v in a if v)
        return scale(f[scale(a,pow(c,-1,q),q)],c,q)
    defects=[]
    for a in P:
        i=next(i for i,v in enumerate(a) if v);b=sub(a,units[i],q)
        if not any(b):continue
        c=sub(f[a],add(f[units[i]],lift(b),q),q)
        assert weight(c)<=3*w
        assert eval_sparse(H,tuple((j,v) for j,v in enumerate(c) if v),q)==zero(len(B[0]))
        defects.append(c)
    e=rank(defects,q);fr=rank(f.values(),q)
    assert fr==ell+e and q**ell<=ball(fr,w,q)
    assert len(defects)==len(P)-ell
    return len(defects)


def canonical_lines(q,d):
    return [a for a in product(range(q),repeat=d) if any(a) and next(x for x in a if x)==1]


def check_invalid_witnesses():
    basis = [(1, 0), (0, 1)]
    plane = {x: x for x in product(range(2), repeat=2)}
    # The repeated column cannot generate the second basis vector.
    wrong_syndrome = ([(1, 0), (1, 0)], basis, plane, 2, 2)
    # The extra nonzero coefficient is on a zero column, so the syndrome is correct.
    overweight = (basis + [(0, 0)], basis[:1], {(0, 0): (0, 0, 0), (1, 0): (1, 0, 1)}, 1, 2)
    for certificate in (defect_certificate, projective_defect_certificate):
        for args, reason in ((wrong_syndrome, 'witness has the wrong syndrome'),
                             (overweight, 'witness exceeds the weight limit')):
            try:
                certificate(*args)
            except AssertionError as error:
                assert str(error) == reason, error
            else:
                raise AssertionError(f'{certificate.__name__} accepted an invalid witness')
    print('Full and projective certificates both rejected a wrong syndrome and a weight-two witness at radius one.')


def subspaces(q,d,ell):
    lines=canonical_lines(q,d);found={}
    for B in combinations(lines,ell):
        if rank(B,q)!=ell:continue
        T=frozenset(lincomb(B,a,q) for a in product(range(q),repeat=ell));found.setdefault(T,B)
    return list(found.items())


def check_defects():
    print('SUBSPACE PROBES AND SHORT-DEPENDENCY CERTIFICATES')
    probes=covered=defects=projective=halfchecks=0;maxe=0
    for q,d,allcases in [(2,3,True),(3,3,False)]:
        lines=canonical_lines(q,d);rng=Random(271126)
        Hs=[list(H) for k in range(1,len(lines)+1) for H in combinations(lines,k)] if allcases else [rng.sample(lines,rng.randrange(2,9)) for _ in range(75)]
        Ts=[item for ell in (2,3) for item in subspaces(q,d,ell)]
        pc=cc=dc=0
        for H in Hs:
            for w in (1,2):
                W=witnesses(H,w,q);eall=short_rank(H,w,q);e6=short_rank(H,2*w,q)
                for T,B in Ts:
                    pc+=1;ell=len(B);ok=T<=W.keys()
                    if q**ell>ball(ell+eall,w,q):assert not ok
                    if q**ell>ball(ell+e6,2*w,q):
                        assert 2*len(T&W.keys())<=q**ell;halfchecks+=1
                    if ok:
                        e,ds=defect_certificate(H,B,W,w,q);cc+=1;dc+=ds;maxe=max(e,maxe);assert e<=eall
                        projective+=projective_defect_certificate(H,B,W,w,q)
        probes+=pc;covered+=cc;defects+=dc
        print(f'q={q}, d={d}, dictionaries={len(Hs)}, subspaces={len(Ts)}, tests={pc}, covered={cc}, defects checked={dc}')
    print(f'{probes} probe tests, {covered} covering certificates, {defects} exact short defects; max extracted rank {maxe}')
    print(f'{projective} projective defects and {halfchecks} half-density bound checks')


def check_sharpness():
    print('\nDIRECT-SUM PROJECTIVE COVER CONTROLS')
    for q,w,h in [(2,1,4),(3,1,3),(2,2,3),(3,2,2),(2,3,2)]:
        ell=w*h;lines=canonical_lines(q,h);H=[]
        for i in range(w):H += [zero(i*h)+a+zero((w-i-1)*h) for a in lines]
        W=witnesses(H,w,q);assert len(W)==q**ell
        B=[tuple(int(i==j) for i in range(ell)) for j in range(ell)]
        e,_=defect_certificate(H,B,W,w,q);assert e==len(H)-ell
        lb=0
        while ball(ell+lb,w,q)<q**ell:lb+=1
        if w==1:assert e==lb
        print(f'q={q}, w={w}, ell={ell}, K={len(H)}, exact e={e}, information-set lower bound={lb}')


def peel(H,q,s):
    K=len(H);S=[{j for j,v in enumerate(a) if v} for a in H]
    # Precompute all support intersections. A zero entry is irrelevant.
    overlap=[{j:len(a&b) for j,b in enumerate(S) if i!=j and a&b} for i,a in enumerate(S)]
    I=set(range(K));rounds=[]
    while True:
        safe=[]
        for i in sorted(I):
            v=sorted((a for j,a in overlap[i].items() if j in I),reverse=True)
            if sum(v[:s-1])<len(S[i]):safe.append(i)
        if not safe:break
        I.difference_update(safe);rounds.append(tuple(safe))
    J=sorted(I);E=len(J)-rank((H[j] for j in J),q)
    return J,E,rounds


def affine_lines(t,a):
    X=list(product(range(t),repeat=a));index={v:i for i,v in enumerate(X)}
    H=[]
    for v in canonical_lines(t,a):
        pending=set(X)
        while pending:
            x=min(pending);B={add(x,scale(v,c,t),t) for c in range(t)}
            assert B<=pending;pending.difference_update(B)
            supp={index[y] for y in B};H.append(tuple(int(j in supp) for j in range(len(X))))
    assert len(H)==t**(a-1)*(t**a-1)//(t-1)
    return H


def check_peel():
    print('\nCERTIFIED SHORT-DEPENDENCY PEELING')
    rng=Random(2711);checks=rescues=0
    for q,d in [(2,4),(3,3)]:
        X=list(product(range(q),repeat=d))
        for _ in range(50):
            H=rng.sample(X,rng.randrange(3,9));w=rng.choice((1,2))
            J,E,rounds=peel(H,q,3*w)
            for c in sparse_patterns(len(H),3*w,q):
                if eval_sparse(H,c,q)==zero(d):assert all(i in J for i,v in c)
            assert short_rank(H,w,q)<=E
            checks+=1;rescues+=len(rounds)>1
    print(f'{checks} exact finite peeling certificates; {rescues} multi-round rescues')
    for q,w,ln in [(2,1,7),(3,2,7)]:
        d=ln+3;core=[(1,0)+zero(d-2),(0,1)+zero(d-2),(1,1)+zero(d-2)]
        H=core+[tuple(int(j in (i,i+1)) for j in range(d)) for i in range(2,d-1)]
        J,E,R=peel(H,q,3*w);assert J==[0,1,2] and E==1 and len(R)>1
        ell=next(i for i in range(1,d+1) if q**i>ball(i+E,w,q))
        T=[a+zero(d-ell) for a in canonical_lines(q,ell)]
        W=witnesses(H,w,q);assert any(a not in W for a in T)
        print(f'path/core q={q}, w={w}, K={len(H)}, peeling rounds={list(map(len,R))}, E={E}, probes={len(T)}, absent={sum(a not in W for a in T)}')
    for t,a,q,w in [(3,3,2,1),(7,2,2,2)]:
        H=affine_lines(t,a);d=len(H[0]);J,E,R=peel(H,q,3*w)
        assert not J and E==0
        ell=w+1;T=[x+zero(d-ell) for x in canonical_lines(q,ell)]
        W=witnesses(H,w,q);assert any(x not in W for x in T)
        nullity=len(H)-rank(H,q);assert nullity>0
        print(f'affine-line design t={t}, a={a}, q={q}, w={w}, d={d}, K={len(H)}, full nullity={nullity}, E={E}, probes={len(T)}, absent={sum(x not in W for x in T)}')


def graph(L,Delta):
    return [(u,L+(u+j)%L) for j in range(Delta) for u in range(L)]


def footprints(L,E):return [(u,) for u in range(2*L)]+list(E)


def conflict_masks(L,E):
    A=footprints(L,E);N=[{u} for u in range(2*L)]
    for u,v in E:N[u].add(v);N[v].add(u)
    result=[]
    for i,a in enumerate(A):
        near=set.union(*(N[u] for u in a))
        result.append({j for j,b in enumerate(A) if j!=i and any(v in near for v in b)})
    return result


def signature(c,L,E,q):
    D={}
    for i,a in c:
        if i<2*L:D[i]=(D.get(i,0)+a)%q
        else:
            for u in E[i-2*L]:D[u]=(D.get(u,0)-a)%q
    return tuple(sorted((u,a) for u,a in D.items() if a))


def good(c,conf):
    I={i for i,a in c}
    return all(not (conf[i]&I) for i in I)


def build_dict(L,Delta,q,d,H0,seed):
    E=graph(L,Delta);assert len(H0)<=len(E);rng=Random(seed)
    X=[tuple(rng.randrange(q) for _ in range(d)) for u in range(2*L)]
    C=list(H0)+[zero(d)]*(len(E)-len(H0))
    H=X+[sub(c,add(X[u],X[v],q),q) for (u,v),c in zip(E,C)]
    return H,E,X,C


def check_planted():
    print('\nPLANTED HIGH-GROWTH SPARSE-SPAN DICTIONARIES')
    total=0;ngood=0;collisionchecks=0
    for q,d,L,Delta,w in [(2,16,8,3,3),(2,21,16,4,3),(3,12,7,3,3),(5,12,4,3,3)]:
        E=graph(L,Delta);K=(Delta+2)*L;Q=q**d;assert ball(K,w,q)<Q
        # Arbitrary prescribed targets: a coordinate subspace plus random points.
        rng=Random(811);T=[a+zero(d-2) for a in product(range(q),repeat=2)]
        T=T[:len(E)];T+= [tuple(rng.randrange(q) for _ in range(d)) for _ in range(len(E)-len(T))]
        for seed in range(100,150):
            H,E,X,C=build_dict(L,Delta,q,d,T,seed)
            lines={scale(a,pow(next(v for v in a if v),-1,q),q) for a in H if any(a)}
            if len(lines)==K:break
        else:raise AssertionError('no projectively distinct finite fixture')
        conf=conflict_masks(L,E);Dmax=2*(Delta+1)**2;assert max(map(len,conf))<=Dmax
        for j,((u,v),c) in enumerate(zip(E,C)):
            assert eval_sparse(H,tuple(sorted(((u,1),(v,1),(2*L+j,1)))),q)==c
        allsets=[set() for j in range(w+1)];G=[Counter() for j in range(w+1)];sig={};cnt=[0]*(w+1)
        for c in sparse_patterns(K,w,q):
            j=len(c);y=eval_sparse(H,c,q);allsets[j].add(y);total+=1
            if good(c,conf):
                z=signature(c,L,E,q);assert z not in sig;sig[z]=c;G[j][y]+=1;cnt[j]+=1;ngood+=1
        S=set();F=Counter();rows=[];cumgood=0
        for j in range(w+1):
            S|=allsets[j];F.update(G[j]);cumgood+=cnt[j]
            b=ball(K,j,q);coll=sum(comb(a,2) for a in F.values());collisionchecks+=1
            assert len(S)>=cumgood-coll
            if j:
                assert cumgood*(K-1)>=b*((K-1)-comb(j,2)*Dmax)
                rows.append((j,len(S),b,cumgood,coll))
        assert set(T)<=S
        print(f'q={q}, d={d}, K={K}, Delta={Delta}, L={L}, planted={len(set(T))}, all radii={rows}')
    print(f'{total} sparse patterns enumerated, {ngood} injective formal signatures, {collisionchecks} exact collision inequalities')


def check_marginals():
    print('\nUNIFORM MATCHING MARGINALS AND POLYNOMIAL/RANK CONTROLS')
    # Exhaustively check joint independence of labels in a disjoint matching.
    checks=0
    for q,d in [(2,1),(3,1),(2,2)]:
        V=list(product(range(q),repeat=d));Q=len(V);C=Counter()
        c0=V[-1];c1=V[1]
        for x0,x1,y0,y1 in product(V,repeat=4):
            C[(sub(c0,add(x0,y0,q),q),sub(c1,add(x1,y1,q),q))]+=1;checks+=1
        assert len(C)==Q**2 and set(C.values())=={Q**2}
    matchecks=0
    for q,a,b in [(2,2,2),(3,2,2),(5,2,2),(3,2,3)]:
        low=0;total=q**(a*b)
        for x in product(range(q),repeat=a*b):
            low+=rank([x[i*b:(i+1)*b] for i in range(a)],q)<min(a,b);matchecks+=1
        assert low*(q-1)<=total
        print(f'q={q}, shape={a}x{b}, singular fraction={low}/{total}')
    polychecks=0
    for q,d,t in [(3,2,1),(3,2,2),(5,2,1),(5,1,3)]:
        X=list(product(range(q),repeat=d));E=[e for e in product(range(t+1),repeat=d) if sum(e)<=t]
        from math import prod
        table=[[prod(pow(x,k,q) for x,k in zip(a,e))%q for e in E] for a in X]
        for c in product(range(q),repeat=len(E)):
            if not any(c):continue
            z=sum(sum(a*b for a,b in zip(c,row))%q==0 for row in table)
            assert z*q<=len(X)*t;polychecks+=1
    print(f'{checks} exact matching assignments, {matchecks} matrix rank comparisons, {polychecks} nonzero polynomial zero-count checks')

if __name__=='__main__':
    check_defects();check_invalid_witnesses();check_sharpness();check_peel();check_planted();check_marginals();print('\nAll checks passed.')
