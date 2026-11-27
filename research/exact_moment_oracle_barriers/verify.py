#!/usr/bin/env python3
"""Exact checks of higher-order moment and distinct-tag obstructions.
No cryptographic curve computations. NetworkX is used only by graph_base.py.
"""
from collections import Counter
from itertools import combinations, product
from math import ceil
import random
from graph_base import (word_for, graph_from_word, incidence, blocks,
                        null_vector, gram, det, stats)


def rank_mod(rows, p):
    piv = {}
    for row in rows:
        a = [v % p for v in row]
        for i, u in sorted(piv.items()):
            if a[i]:
                c = a[i]
                a = [(v-c*w) % p for v,w in zip(a,u)]
        j = next((j for j,v in enumerate(a) if v),None)
        if j is not None:
            c = pow(a[j],-1,p)
            piv[j] = [(v*c) % p for v in a]
            if len(piv)==len(a):
                return len(piv)
    return len(piv)


def coset(p,k,s):
    b=(p-1).bit_length()
    bits=[tuple((x>>j)&1 for j in range(b)) for x in range(p)]
    return [(1,)+sum((bits[x] for x in a),())+bits[(s-sum(a))%p]
            for a in product(range(p),repeat=k-1)]


def binary_moments(A,t):
    out=Counter()
    for a in A:
        support=[j for j,v in enumerate(a) if v]
        for k in range(min(t,len(support))+1):
            out.update(combinations(support,k))
    return out


def dot(a,z,p):
    return sum(x*y for x,y in zip(a,z))%p


def check_cosets():
    print('HIGHER-ORDER BINARY COSETS')
    print('r | t | d | rows/coset | #nonzero moment keys | rank H0 | rank H0+H1 | rational rank certificate')
    fixtures=[(2,2),(3,2),(3,4),(3,6),(5,2),(5,3),(5,4),(5,5),(7,2),(7,3),(7,4),(11,2),(11,3),(13,2),(13,3),(17,2),(31,2)]
    for p,t in fixtures:
        k=t+1;b=(p-1).bit_length();D=1+k*b
        H=coset(p,k,0);J=coset(p,k,1)
        assert len(H)==len(J)==p**t
        assert len(set(H))==len(H) and len(set(J))==len(J)
        assert not (set(H)&set(J))
        h=binary_moments(H,t);j=binary_moments(J,t)
        assert h==j
        z=(0,)+tuple(1<<i for i in range(b))*k
        assert all(dot(a,z,p)==0 for a in H) and all(dot(a,z,p)==1 for a in J)
        assert rank_mod(H,p)==rank_mod(J,p)==D-1 and rank_mod(H+J,p)==D
        # A full rank over another prime certifies rational full rank, exactly.
        q=1009
        assert p!=q and rank_mod(H,q)==rank_mod(J,q)==D
        print(f'{p} | {t} | {D} | {len(H)} | {len(h)} | {D-1} | {D} | mod {q}')
    print(f'{len(fixtures)} exact coset fixtures passed, through degree 6.')


def dense_graph_blocks(p):
    ell=p.bit_length(); G,e,(t,u)=graph_from_word(word_for(p,ell)); a,b=e
    G.remove_edge(a,b); path=[a]
    for _ in range(ell-1):
        v=G.number_of_nodes();G.add_node(v);path.append(v)
    path.append(b);G.add_edges_from(zip(path,path[1:]))
    U,V,E=incidence(G)
    assert det(gram(U))==t+ell*u and (t+ell*u)%p==0 and u%p
    A=[a+b for a,b in zip(U,V)]; B=[a+[-v for v in b] for a,b in zip(U,V)]
    H=sorted(tuple(a) for a in blocks(A));J=sorted(tuple(a) for a in blocks(B))
    assert stats(H,len(H[0]))==stats(J,len(J[0]))
    z=tuple(null_vector(H,p))
    assert len(set(H))==len(H) and len(set(J))==len(J)
    assert sum(dot(a,z,p)!=0 for a in J)>=2*ell
    return H,J,z


def tag(a,q):
    return tuple(int(x==j) for x in a for j in range(1,q))


def max_radix(p,L,h):
    q=1
    while L*(q+1)**h<p:
        q+=1
    return q


def fixed_tag(a):
    return sum(x!=0 for x in a)<=1


def check_tags():
    print('\nDISTINCT SPARSE TAGS (full exact finite ranks)')
    print('r | h | q | blocks | rows | d | max norm^2 | rank H | rank mixed | violations')
    cases=[(257,2),(1009,2),(1009,3),(1009,4),(4099,2),(4099,3)]
    rng=random.Random(20261127); interval_checks=0
    for p,h in cases:
        H,J,z=dense_graph_blocks(p);L=len(H);q=max_radix(p,L,h)
        tags=sorted((tag(a,q),a) for a in product(range(q),repeat=h))
        free=[i for i,(_,a) in enumerate(tags) if not fixed_tag(a)]
        assert q>=2 and free
        js=set(free[:max(1,len(free)//2)])
        base=[];mixed=[];known=[]
        for i,(v,a) in enumerate(tags):
            rows=[v+w for w in H]
            base+=rows
            mixed += [v+w for w in (J if i in js else H)]
            if fixed_tag(a):known+=rows
        assert len(base)<p and len(base)==len(mixed)
        assert len(set(base))==len(base) and len(set(mixed))==len(mixed)
        assert base==sorted(base) and mixed==sorted(mixed)
        d=len(base[0]); Z=(0,)*(h*(q-1))+z
        assert rank_mod(known,p)==rank_mod(base,p)==d-1 and rank_mod(mixed,p)==d
        assert rank_mod(base,65537)==rank_mod(mixed,65537)==d
        assert all(dot(a,Z,p)==0 for a in base)
        violations=sum(dot(a,Z,p)!=0 for a in mixed)
        assert violations>=len(js)*sum(dot(a,z,p)!=0 for a in J)
        norm=max(sum(v*v for v in a) for a in base+mixed)
        assert norm<=5+h
        # Exact full-block invariance with tags, for each changed block.
        for i in js:
            assert stats(base[i*L:(i+1)*L],d)==stats(mixed[i*L:(i+1)*L],d)
        for _ in range(20):
            i=rng.choice(free);l,u=sorted(rng.sample(range(len(base)+1),2))
            if l<=i*L and (i+1)*L<=u:
                v=tags[i][0]; changed=[v+w for w in J]
                assert stats(base[i*L:(i+1)*L],d)==stats(changed,d)
                interval_checks+=1
            elif u<=i*L or (i+1)*L<=l:
                interval_checks+=1
        print(f'{p} | {h} | {q} | {len(tags)} | {len(base)} | {d} | {norm} | {d-1} | {d} | {violations}')
    print(f'{len(cases)} full tagged rank checks; {interval_checks} interval checks.')


def check_score_support():
    # Exhaust all multisets of size < p via counts for p <= 7.
    print('\nNO-WRAP SCORE-SUPPORT BARRIER')
    checked=0; equal_examples=0
    for p in (3,5,7):
        for ns in product(range(p),repeat=p-1):
            if not 0<sum(ns)<p:continue
            supp=sum(n>0 for n in ns)
            sums=[sum(n*pow(s,j,p) for s,n in enumerate(ns,1))%p
                  for j in range(1,p)]
            for t in range(1,p):
                checked+=1
                if all(x==0 for x in sums[:t]):
                    assert supp>=t+1
                    equal_examples+=1
    print(f'{checked} exact score/moment cases; {equal_examples} zero-moment cases; all passed.')
    # Sharpness for powers of roots of unity: each score appears once.
    ex=0
    for p in (5,7,11,13,17,19,31,41):
        for k in range(2,p):
            if (p-1)%k:continue
            roots=[x for x in range(1,p) if pow(x,k,p)==1]
            assert len(roots)==k
            assert all(sum(pow(x,j,p) for x in roots)%p==0 for j in range(1,k))
            assert sum(pow(x,k,p) for x in roots)%p==k
            ex+=1
    print(f'{ex} sharp finite-field score distributions.')

if __name__=='__main__':
    check_cosets();check_tags();check_score_support()
    print('\nAll checks passed.')
