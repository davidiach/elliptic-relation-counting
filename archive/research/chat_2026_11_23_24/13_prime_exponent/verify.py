from itertools import combinations, product
from random import Random
from math import comb, isqrt

rng = Random(20261124)
p = 65521
a = 44
c = 1
O = None

def prime(n):
    return n >= 2 and all(n % j for j in range(2,isqrt(n)+1))

def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x,y=P; u,v=Q
    if x==u and (y+v)%p==0: return None
    s=((3*x*x+a)*pow(2*y,-1,p) if x==u else (v-y)*pow(u-x,-1,p))%p
    z=(s*s-x-u)%p
    return z,(s*(x-z)-y)%p

def neg(P): return None if P is None else (P[0],-P[1]%p)
def total(Ps):
    S=None
    for P in Ps:S=add(S,P)
    return S

# Uniform sampling from the actual complete group (toy field only).
roots={}
for y in range(p): roots.setdefault(y*y%p,[]).append(y)
G=[None]+[(x,y) for x in range(p) for y in roots.get((x*x*x+a*x+c)%p,[])]
r=len(G)
assert prime(p) and prime(r), (p,r)
T=(0,1); assert T in G and add(T,T) is not None

def basis(m):
    d=m+1
    return [(j,0) for j in range(1,d//2+1)]+[(j,1) for j in range((d-3)//2+1)]

def row(P,m):
    x,y=P; u,v=neg(T)
    return [(pow(x,j,p)*pow(y,k,p)-pow(u,j,p)*pow(v,k,p))%p for j,k in basis(m)]

def det(A):
    A=[list(row) for row in A]; n=len(A); z=1
    for i in range(n):
        j=next((j for j in range(i,n) if A[j][i]%p),None)
        if j is None:return 0
        if j!=i:A[i],A[j]=A[j],A[i];z=-z
        u=A[i][i]%p;z=z*u%p; inv=pow(u,-1,p)
        for j in range(i+1,n):
            f=A[j][i]*inv%p
            for k in range(i+1,n): A[j][k]=(A[j][k]-f*A[i][k])%p
            A[j][i]=0
    return z%p

def planted(m,w):
    n=m*w
    idx=[i*w+rng.randrange(w) for i in range(m)]
    s=[rng.randrange(2) for _ in range(m)]
    pts=[rng.choice(G) for _ in range(n)]
    marked=[neg(pts[j]) if e else pts[j] for j,e in zip(idx[:-1],s[:-1])]
    P=add(T,neg(total(marked)));pts[idx[-1]]=neg(P) if s[-1] else P
    points=[Q for P in pts for Q in (P,neg(P))]
    hit=tuple(sorted(2*j+e for j,e in zip(idx,s)))
    assert total(points[j] for j in hit)==T
    return points,hit

print(f'Curve p={p}, a={a}, c={c}, prime group order={r}.')
for m,w,trials in [(3,2,120),(4,2,80)]:
    N=2*m*w; reg=unique=tests=other=cb=0
    for _ in range(trials):
        pts,hit=planted(m,w)
        if None in pts or T in pts or neg(T) in pts or len(set(pts))!=N: continue
        reg+=1;rows=[row(P,m) for P in pts]; hits=[]
        assert len(rows[0])==m
        # Compute both the elliptic sum and the determinant for every m-subset.
        for I in combinations(range(N),m):
            expected=total(pts[j] for j in I)==T
            actual=det([rows[j] for j in I])==0
            assert actual==expected
            tests+=1
            if expected:hits.append(I)
        assert hit in hits
        other+=len(hits)-1
        unique+=hits==[hit]
        # Check Cauchy-Binet for dense public linear-combination queries.
        for _ in range(2):
            C=[[rng.randrange(p) for j in range(N)] for i in range(m)]
            CR=[[sum(C[i][j]*rows[j][k] for j in range(N))%p for k in range(m)] for i in range(m)]
            rhs=sum(det([[C[i][j] for j in I] for i in range(m)])*det([rows[j] for j in I]) for I in combinations(range(N),m))%p
            assert det(CR)==rhs
            cb+=1
    print(f'm={m}, w={w}: trials={trials}, regular={reg}, unique={unique}, extra circuits={other}, determinant/group comparisons={tests}, dense Cauchy-Binet identities={cb}, union bound for extra circuits={(comb(N,m)-1)/r:.8f}')
print('All assertions passed. Finite diagnostics; no faster relation evaluator implemented.')
