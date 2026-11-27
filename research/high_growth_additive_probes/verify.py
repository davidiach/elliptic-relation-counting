#!/usr/bin/env python3
"""Exact finite checks for higher-order probes and high-growth additive covers.
All examples are abstract finite additive planes. No curve implementation.
"""
from collections import Counter, deque
from itertools import combinations, combinations_with_replacement, product
from math import comb, isqrt, prod
from random import Random

if not __debug__:
    raise RuntimeError('Run with assertions enabled.')


def isprime(p):
    if p < 2: return False
    if p % 2 == 0: return p == 2
    return all(p % d for d in range(3, isqrt(p)+1, 2))


def nextprime(p):
    while not isprime(p): p += 1
    return p


def ceilroot(a,h):
    lo,hi=0,1
    while hi**h<a: hi*=2
    while hi-lo>1:
        z=(lo+hi)//2
        if z**h<a:lo=z
        else:hi=z
    return hi


def encode_powers(K,h,q,r,guard=True):
    assert q>=K and q>h and isprime(q)
    base=h*(q-1)+1; cut=(h+1)//2
    if guard: assert r>=base**cut
    H=[]
    for t in range(K):
        d=[pow(t,j,q) for j in range(1,h+1)]
        H.append((sum(x*base**j for j,x in enumerate(d[:cut]))%r,
                  sum(x*base**j for j,x in enumerate(d[cut:]))%r))
    return H


def sumkey(xs,r):
    return sum(x for x,y in xs)%r,sum(y for x,y in xs)%r


def add(a,b,r):return ((a[0]+b[0])%r,(a[1]+b[1])%r)
def sub(a,b,r):return ((a[0]-b[0])%r,(a[1]-b[1])%r)


def sumsets(A,r):
    S={(0,0)}
    for B in A:S={add(a,b,r) for a in S for b in B}
    return S


def check_bh():
    print('PACKED HIGHER-ORDER SIDON PROBES')
    print('K | order | auxiliary prime | plane prime | multiset sums checked')
    total=0
    for K,h in [(7,2),(8,4),(8,6),(7,8),(9,5),(37,4)]:
        q=nextprime(max(K,h+1,3));base=h*(q-1)+1;r=nextprime(base**((h+1)//2))
        H=encode_powers(K,h,q,r);assert len(set(H))==K
        count=0
        for j in range(1,h+1):
            sums=set()
            for xs in combinations_with_replacement(H,j):
                key=sumkey(xs,r);assert key not in sums;sums.add(key)
            assert len(sums)==comb(K+j-1,j);count+=len(sums)
        total+=count;print(f'{K} | {h} | {q} | {r} | {count}')
    H=encode_powers(8,4,11,17,guard=False)
    assert len({sumkey(xs,17) for xs in combinations_with_replacement(H,4)})<comb(11,4)
    print(f'{total} exact distinct multiset sums; failed-guard control exhibits a collision.')


def graph_girth(adj):
    best=len(adj)+1
    for s in range(len(adj)):
        dist={s:0};parent={s:-1};todo=deque([s])
        while todo:
            u=todo.popleft()
            if 2*dist[u]>=best:continue
            for v in adj[u]:
                if v not in dist:
                    dist[v]=dist[u]+1;parent[v]=u;todo.append(v)
                elif parent[u]!=v and parent[v]!=u:
                    best=min(best,dist[u]+dist[v]+1)
    return best


def adj_from_edges(l,edges):
    A=[set() for _ in range(2*l)]
    for u,v in edges:
        assert 0<=u<l and 0<=v<l and v+l not in A[u]
        A[u].add(v+l);A[v+l].add(u)
    return A


def check_moore():
    print('\nEXHAUSTIVE SMALL BIPARTITE GRAPH CONTROLS')
    for n in (3,4):
        E=list(product(range(n),repeat=2)); counts={2:0,3:0};best={2:0,3:0}
        for mask in range(1<<len(E)):
            edges=[e for j,e in enumerate(E) if (mask>>j)&1]
            girth=graph_girth(adj_from_edges(n,edges))
            for h in (2,3):
                if girth>2*h:
                    # The real Moore inequality implies this rounded integer cap.
                    assert len(edges)<=n*(1+ceilroot(n,h))
                    counts[h]+=1;best[h]=max(best[h],len(edges))
        print(f'n={n} graphs={1<<len(E)}: no short cycles h=2 {counts[2]} max edges {best[2]}; h=3 {counts[3]} max edges {best[3]}')


def check_two_groups():
    print('\nTWO-FACTOR PROBE CONSTRUCTIONS')
    rng=Random(251127);total=0
    for n,h in [(12,4),(20,4),(12,6),(32,6)]:
        K=1+n*(1+ceilroot(n,h));q=nextprime(max(K,h+1,3));B=h*(q-1)+1;r=nextprime(B**((h+1)//2))
        H=encode_powers(K,h,q,r);assert len(set(H))==K
        A={H[0],*H[1:n]};C={(0,0),*(sub(v,H[0],r) for v in H[n:2*n-1])}
        assert len(A)==len(C)==n and set(H[:2*n-1])<=sumsets([A,C],r)
        maxhit=len(set(H)&sumsets([A,C],r));assert maxhit<K
        for _ in range(30):
            A={divmod(x,r) for x in rng.sample(range(r*r),n)}
            C={divmod(x,r) for x in rng.sample(range(r*r),n)}
            hit=len(set(H)&sumsets([A,C],r));assert hit<K;maxhit=max(maxhit,hit)
        total+=31
        print(f'n={n} order={h} probes={K} plane prime={r} max hit={maxhit}')
    print(f'{total} finite cover exclusions (one constructed large intersection per fixture).')


def build_graph(m,l):
    # Exact proof supplies a guard for the greedy case. Smaller fixtures can
    # still succeed; verify their girth independently rather than assume it.
    if m==3:edges=[(u,u) for u in range(l)]
    elif m==4:
        assert l>m
        edges=[(u,u) for u in range(l)]+[(u,(u+1)%l) for u in range(l)]
    else:
        edges=[];adj=[set() for _ in range(2*l)];D=2*m
        for _ in range((m-2)*l):
            u=min(range(l),key=lambda u:(len(adj[u]),u));assert len(adj[u])<D
            seen={u};front={u}
            for j in range(max(4,2*(m-2))-1):
                front={w for v in front for w in adj[v]}-seen
                seen|=front
                if not front:break
            cand=[v for v in range(l,2*l) if v not in seen and len(adj[v])<D]
            assert cand,'finite greedy graph failed'
            v=min(cand,key=lambda v:(len(adj[v]),v));adj[u].add(v);adj[v].add(u);edges.append((u,v-l))
    assert len(edges)==(m-2)*l
    A=adj_from_edges(l,edges);assert graph_girth(A)>max(4,2*(m-2))
    assert max(map(len,A))<=2*m
    return [edges[i*l:(i+1)*l] for i in range(m-2)],graph_girth(A),max(map(len,A))


def tuple_signature(I,js,l,edges):
    choice=dict(zip(I,js));z=Counter();left=set();right=set()
    for i,j in choice.items():
        if i<2:continue
        u,v=edges[i-2][j];left.add(u);right.add(v)
        z[u]-=1;z[l+v]-=1
    if 0 in choice:
        if choice[0] in left:return None
        z[choice[0]]+=1
    if 1 in choice:
        if choice[1] in right:return None
        z[l+choice[1]]+=1
    return tuple(sorted((j,c) for j,c in z.items() if c))


def double_signature(i,js,l,edges):
    C=Counter()
    for j in js:
        if i==0:C[j]+=1
        elif i==1:C[l+j]+=1
        else:
            u,v=edges[i-2][j];C[u]-=1;C[l+v]-=1
    return tuple(sorted(C.items()))


def subsets(m):
    for k in range(1,m+1):yield from combinations(range(m),k)


def matchings(edges):
    colors=[];L=[];R=[]
    for u,v in edges:
        c=next((c for c in range(len(colors)) if u not in L[c] and v not in R[c]),None)
        if c is None:c=len(colors);colors.append([]);L.append(set());R.append(set())
        colors[c].append((u,v));L[c].add(u);R[c].add(v)
    assert sum(map(len,colors))==len(edges)
    return colors


def check_signatures():
    print('\nHIGH-GIRTH GRAPH SIGNATURE CONTROLS')
    rng=Random(251127);total=doubletotal=0
    for m,l,exhaust in [(3,6,True),(4,8,True),(5,300,False),(6,1000,False)]:
        E,g,D=build_graph(m,l);n=0;nd=0
        for I in subsets(m):
            if exhaust:rows=product(range(l),repeat=len(I))
            else:rows=(tuple(rng.randrange(l) for _ in I) for _ in range(300))
            seen={};valid=0
            for js in rows:
                z=tuple_signature(I,js,l,E)
                if z is None:continue
                assert z not in seen or seen[z]==js;seen[z]=js;valid+=1
            if exhaust:
                assert len(seen)==valid
                assert valid*l>=l**len(I)*max(0,l-2*(m-2))
            n+=valid
        for i in range(m):
            if l<=300:rows=combinations_with_replacement(range(l),2)
            else:rows=(tuple(sorted((rng.randrange(l),rng.randrange(l)))) for _ in range(1500))
            seen={}
            for js in rows:
                z=double_signature(i,js,l,E);assert z not in seen or seen[z]==js;seen[z]=js;nd+=1
        C=matchings([e for A in E for e in A]);assert len(C)<=2*D-1
        total+=n;doubletotal+=nd
        print(f'm={m} random slots={l} edges={(m-2)*l} girth={g if g<=2*l else "forest"} degree={D} matching colors={len(C)}: {n} good signatures, {nd} double signatures, exhaustive={exhaust}')
    print(f'{total} good tuple signatures and {doubletotal} double signatures checked.')


def peval(a,x,p):
    y=0
    for v in reversed(a):y=(y*x+v)%p
    return y


def check_independence():
    print('\nLIMITED-INDEPENDENCE AND CHARACTER-MOMENT CONTROLS')
    p,k=5,3;arrays=[tuple(peval(a,x,p) for x in range(p)) for a in product(range(p),repeat=k)]
    tests=0
    for j in range(1,k+1):
        for I in combinations(range(p),j):
            C=Counter(tuple(a[i] for i in I) for a in arrays)
            assert len(C)==p**j and set(C.values())=={p**(k-j)};tests+=1
    linear=0
    for j in range(1,k+1):
        for I in combinations(range(p),j):
            for cs in product(range(1,p),repeat=j):
                C=Counter(sum(c*a[i] for c,i in zip(cs,I))%p for a in arrays)
                assert len(C)==p and set(C.values())=={p**(k-1)};linear+=1
    moments=0
    for L in range(1,5):
        for t in range(1,5):
            C=Counter(tuple(sorted(a)) for a in product(range(L),repeat=t))
            moment=sum(c*c for c in C.values())
            assert moment<=(2*t)**(2*t)*L**t;moments+=1
    print(f'{tests} uniform marginals, {linear} nonzero linear forms, {moments} exact character-moment combinatorial bounds.')


def make_lists(m,r,l,E,h,seed,axes):
    rng=Random(seed);k=64*m
    ax=[rng.randrange(r) for _ in range(k)];ay=[rng.randrange(r) for _ in range(k)]
    V=[(peval(ax,j,r),peval(ay,j,r)) for j in range(2*l)]
    vals=[V[:l],V[l:]]
    for i in range(m-2):
        vals.append([sub(h[i][j],add(V[u],V[l+v],r),r) for j,(u,v) in enumerate(E[i])])
    A=[set(B)|set(D) for B,D in zip(vals,axes)]
    return A,vals


def axes_sets(m,r,full):
    if not full:return [{(0,0)} for _ in range(m)]
    q=ceilroot(r,m)
    return [{(a*pow(q,i,r)%r,0) for a in range(q)}|{(0,a*pow(q,i,r)%r) for a in range(q)} for i in range(m)]


def check_covers():
    print('\nFINITE HIGH-GROWTH COVER FIXTURES')
    print('m | r | n | prescribed points | full axes | |S| / n^m | min subset fraction | min double fraction | good-tuple collisions')
    for m,r,l,full in [(3,101,6,True),(4,1009,10,True),(3,10007,39,False),(4,100003,15,False)]:
        assert isprime(r) and r>2*l
        E,g,D=build_graph(m,l);F=axes_sets(m,r,full);d=len(F[0]);assert all(len(B)==d for B in F)
        n=l+d;M=n**m;assert M<r*r
        # Deliberately structured prescribed points, plus a common translation.
        H=[(j,j*j%r) for j in range(1,(m-2)*l+1)]
        hh=[H[i*l:(i+1)*l] for i in range(m-2)]
        for seed in range(251127,251227):
            A,V=make_lists(m,r,l,E,hh,seed,F)
            if all(len(B)==n for B in A):break
        else:raise AssertionError('no collision-free individual fixture')
        for i in range(m-2):
            for j,(u,v) in enumerate(E[i]):assert sumkey([V[0][u],V[1][v],V[i+2][j]],r)==hh[i][j]
        H=set(H)|{(0,0)};shift=(3,5);A[0]={add(x,shift,r) for x in A[0]};V[0]=[add(x,shift,r) for x in V[0]]
        minratio=(1,1);mindouble=(1,1);collisions=0;fullsize=0
        for I in subsets(m):
            S=sumsets([A[i] for i in I],r)
            count=Counter()
            for js in product(range(l),repeat=len(I)):
                z=tuple_signature(I,js,l,E)
                if z is not None:count[sumkey([V[i][j] for i,j in zip(I,js)],r)]+=1
            good=sum(count.values());C=sum(comb(c,2) for c in count.values());collisions+=C
            assert len(S)>=good-C and good*l>=l**len(I)*max(0,l-2*(m-2))
            den=n**len(I)
            if len(S)*minratio[1]<minratio[0]*den:minratio=(len(S),den)
            if len(I)==m:
                fullsize=len(S);assert {add(x,shift,r) for x in H}<=S
                if full:assert all(add((j,0),shift,r) in S and add((0,j),shift,r) in S for j in range(r))
        for i in range(m):
            S=sumsets([A[i],A[i]],r);den=comb(n+1,2)
            C=Counter(sumkey([V[i][u],V[i][v]],r) for u,v in combinations_with_replacement(range(l),2))
            assert len(S)>=comb(l+1,2)-sum(comb(c,2) for c in C.values())
            if len(S)*mindouble[1]<mindouble[0]*den:mindouble=(len(S),den)
        print(f'{m} | {r} | {n} | {len(H)} | {full} | {fullsize}/{M} | {minratio[0]}/{minratio[1]} | {mindouble[0]}/{mindouble[1]} | {collisions}')
    print('4 full finite covers; all subset and double-sum bounds verified. Axis-free controls are not full-projection claims.')


if __name__=='__main__':
    check_bh();check_moore();check_two_groups();check_signatures();check_independence();check_covers()
    print('\nAll checks passed.')
