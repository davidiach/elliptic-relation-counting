#!/usr/bin/env python3
"""Exact toy checks of the planar cut/cycle Gram indistinguishability lemma.
Requires networkx (planar embedding only); all arithmetic is integer or mod p.
No elliptic-curve implementation. The small graph search is for fixtures only.
"""
from collections import deque
from itertools import combinations_with_replacement
import random
from math import prod
import networkx as nx


def trans(A): return list(map(list, zip(*A)))

def gram(A):
    T=trans(A)
    return tuple(tuple(sum(x*y for x,y in zip(a,b)) for b in T) for a in T)

def moment(A,k):
    out={}
    for row in A:
        support=[i for i,v in enumerate(row) if v]
        for inds in combinations_with_replacement(support,k):
            out[inds]=out.get(inds,0)+prod(row[i] for i in inds)
    return {inds:v for inds,v in out.items() if v}

def stats(A, d):
    if not A: return (0, (0,)*d, tuple((0,)*d for _ in range(d)), {})
    return len(A), tuple(map(sum,trans(A))), gram(A), moment(A,3)

def rank(A,p):
    if not A: return 0
    A=[[v%p for v in row] for row in A]; k=0
    for j in range(len(A[0])):
        i=next((i for i in range(k,len(A)) if A[i][j]),None)
        if i is None: continue
        A[k],A[i]=A[i],A[k]
        z=pow(A[k][j],-1,p); A[k]=[(x*z)%p for x in A[k]]
        for i in range(k+1,len(A)):
            if A[i][j]:
                z=A[i][j]; A[i]=[(x-z*y)%p for x,y in zip(A[i],A[k])]
        k+=1
        if k==len(A): break
    return k

def det(A):
    A=[list(row) for row in A]; n=len(A); sign=1; prev=1
    if n==0:return 1
    for k in range(n-1):
        i=next((i for i in range(k,n) if A[i][k]),None)
        if i is None: return 0
        if i!=k: A[k],A[i]=A[i],A[k]; sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=A[i][j]*pivot-A[i][k]*A[k][j]
                q,rem=divmod(z,prev); assert rem==0
                A[i][j]=q
        for i in range(k+1,n): A[i][k]=0
        prev=pivot
    return sign*A[-1][-1]

def word_for(p, L=1):
    # Graph operations A: subdivide witness edge; C: attach a triangle, new witness edge.
    # The pair tracks (tau(G/e), tau(G-e)) modulo p. BFS is only a finite test generator.
    q=deque([(1,0)]); parent={(1,0):None}
    while q:
        x,y=q.popleft()
        if x and y and (x+L*y)%p==0:
            out=[]; state=(x,y)
            while parent[state] is not None:
                state,c=parent[state]; out.append(c)
            return ''.join(reversed(out))
        for c,z in [('A',((x+y)%p,y)),('C',((2*x+y)%p,(x+y)%p))]:
            if z not in parent:parent[z]=((x,y),c);q.append(z)
    raise AssertionError('No witness found')

def graph_from_word(word):
    G=nx.Graph(); G.add_edge(0,1); a,b=0,1; t,u=1,0
    for c in word:
        v=G.number_of_nodes()
        if c=='A': G.remove_edge(a,b); t,u=t+u,u
        else: t,u=2*t+u,t+u
        G.add_edges_from([(a,v),(v,b)]); b=v
    return G,(a,b),(t,u)

def incidence(G):
    ok,emb=nx.check_planarity(G); assert ok and nx.is_connected(G)
    faces={}; marked=set(); nf=0
    for a,b in emb.edges():
        if (a,b) in marked:continue
        walk=emb.traverse_face(a,b,marked)
        for x,y in zip(walk,walk[1:]+walk[:1]): faces[x,y]=nf
        nf+=1
    E=sorted((min(a,b),max(a,b)) for a,b in G.edges()); N=G.number_of_nodes()
    assert N-len(E)+nf==2
    U=[];V=[]
    for a,b in E:
        u=[0]*N;u[a]-=1;u[b]+=1;U.append(u[:-1])
        v=[0]*nf;v[faces[a,b]]+=1;v[faces[b,a]]-=1;V.append(v[:-1])
    return U,V,E

def null_vector(A,p):
    A=[[v%p for v in row] for row in A]; n=len(A[0]); k=0; piv=[]
    for j in range(n):
        i=next((i for i in range(k,len(A)) if A[i][j]),None)
        if i is None: continue
        A[k],A[i]=A[i],A[k]; z=pow(A[k][j],-1,p); A[k]=[(x*z)%p for x in A[k]]
        for i in range(len(A)):
            if i!=k and A[i][j]:
                z=A[i][j];A[i]=[(x-z*y)%p for x,y in zip(A[i],A[k])]
        piv.append(j);k+=1
        if k==len(A):break
    free=[j for j in range(n) if j not in piv]; assert len(free)==1
    u=[0]*n;u[free[0]]=1
    for i,j in enumerate(piv):u[j]=(-A[i][free[0]])%p
    assert all(sum(x*y for x,y in zip(row,u))%p==0 for row in A)
    return u

def blocks(A):
    return [sign_row+[1] for row in A for sign_row in (row,[-v for v in row])]

def dense_fixture(p):
    L=p.bit_length(); G,e,(t,u)=graph_from_word(word_for(p,L)); a,b=e
    G.remove_edge(a,b); path=[a]
    for _ in range(L-1):
        v=G.number_of_nodes();G.add_node(v);path.append(v)
    path.append(b); G.add_edges_from(zip(path,path[1:]))
    path_edges={tuple(sorted(x)) for x in zip(path,path[1:])}
    U,V,E=incidence(G); D=len(E)
    assert (t+L*u)%p==0 and u%p
    assert det(gram(U))==t+L*u
    assert all(sum(x*y for x,y in zip(a,b))==0 for a in trans(U) for b in trans(V))
    for ee in path_edges:
        # Each deletion leaves precisely the original deleted-edge spanning tree count.
        Gdel=[row for row,key in zip(U,E) if key!=ee]
        assert det(gram(Gdel))==u
    A=[a+b for a,b in zip(U,V)]; B=[a+[-v for v in b] for a,b in zip(U,V)]
    H,J=blocks(A),blocks(B); d=D+1
    assert stats(H,d)==stats(J,d) and rank(H,p)==d-1 and rank(H+J,p)==d
    z=null_vector(H,p)
    violated=[sum(x*y for x,y in zip(row,z))%p!=0 for row in J]
    assert all(violated[2*i] and violated[2*i+1] for i,ee in enumerate(E) if ee in path_edges)
    assert sum(violated)>=2*L
    assert max(sum(x*x for x in row) for row in H+J)<=5
    # Use an even number of hidden positions with total rows < p.
    M=max(0, ((p-1)//len(J)-1)//2*2)
    if M>=2: assert (M+1)*len(J)<p
    return D,L,sum(violated),len(J),M

def main():
    ps=[3,5,7,11,13,17,19,31,101,257,1009,65537,104729]
    rng=random.Random(20261127); checks=0
    print('r | V | E | rank(+) | rank(mixed) | max squared norm | tau/r')
    for p in ps:
        word=word_for(p); G,e,(t,u)=graph_from_word(word)
        U,V,E=incidence(G); D=len(E)
        assert all(sum(x*y for x,y in zip(a,b))==0 for a in trans(U) for b in trans(V))
        assert det(gram(U))==t+u and (t+u)%p==0 and u%p
        Udel=[row for row,ee in zip(U,E) if ee!=tuple(sorted(e))]
        assert det(gram(Udel))==u
        A=[a+b for a,b in zip(U,V)]; B=[a+[-v for v in b] for a,b in zip(U,V)]
        assert det(A)!=0 and abs(det(A))==t+u
        assert gram(A)==gram(B)
        assert rank(A,p)==rank(B,p)==D-1 and rank(A+B,p)==D
        H,J=blocks(A),blocks(B); d=D+1
        assert stats(H,d)==stats(J,d)
        assert moment(H,4)!=moment(J,4)
        assert all(row[-1]==1 for row in H+J)
        assert rank(H,p)==d-1 and rank(H+J,p)==d
        assert max(sum(x*x for x in row) for row in H+J)<=5
        z=null_vector(H,p); assert z[-1]==0
        assert any(sum(x*y for x,y in zip(row,z))%p for row in J)
        # For each interval, hidden replacement cannot be observed when its block
        # lies fully inside or outside the interval. Sample exact oracle comparisons.
        M=9; L=len(H); base=H*(M+1)
        for _ in range(100):
            h=rng.randrange(1,M+1); a,b=sorted(rng.sample(range(len(base)+1),2))
            if (a<=h*L and (h+1)*L<=b) or b<=h*L or (h+1)*L<=a:
                other=base[:h*L]+J+base[(h+1)*L:]
                # additivity, compute only overlap with replaced block
                lo=max(a,h*L);hi=min(b,(h+1)*L)
                if hi>lo:
                    assert stats(base[lo:hi],d)==stats(other[lo:hi],d)
                else: assert base[a:b]==other[a:b]
                checks+=1
        print(f'{p} | {G.number_of_nodes()} | {D} | {d-1} | {d} | {max(sum(x*x for x in row) for row in H+J)} | {(t+u)//p}')
    print(f'{len(ps)} graph fixtures; {checks} exact interval-invisibility checks; all passed.')
    print('All 13 paired blocks agree through moment degree 3 and differ at degree 4.')
    print('Dense variant: r | E | path length | violated rows / total | no-wrap M')
    for p in ps:
        D,L,viol,total,M=dense_fixture(p)
        print(f'{p} | {D} | {L} | {viol}/{total} | {M}')
    print(f'{len(ps)} dense fixtures; all path-edge deletions and nonzero scores checked exactly.')

if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions required')
    main()
