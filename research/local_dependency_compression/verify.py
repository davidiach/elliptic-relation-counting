#!/usr/bin/env python3
"""Exact tests for abstract local-dependency compression and prefix probes."""
from itertools import combinations, product
from math import comb
from random import Random
from collections import defaultdict

if not __debug__: raise RuntimeError('Assertions required')

def add(x,y,p): return tuple((a+b)%p for a,b in zip(x,y))
def sub(x,y,p): return tuple((a-b)%p for a,b in zip(x,y))
def scale(x,c,p): return tuple(a*c%p for a in x)
def dot(x,y,p): return sum(a*b for a,b in zip(x,y))%p
def norm(x,p):
    j=next((j for j,v in enumerate(x) if v),None)
    return scale(x,pow(x[j],-1,p),p) if j is not None else tuple(x)
def rank(rows,p):
    B={}
    for row in rows:
        a=list(row)
        for i,b in sorted(B.items()):
            if a[i]:
                z=a[i];a=[(x-z*y)%p for x,y in zip(a,b)]
        j=next((j for j,x in enumerate(a) if x),None)
        if j is not None:B[j]=scale(a,pow(a[j],-1,p),p)
    return len(B)
def basis(rows,p):
    B=[]
    for x in rows:
        if rank(B+[x],p)>len(B):B.append(x)
    return B
def contains(B,x,p):return rank(B+[x],p)==len(B)
def zero(d):return (0,)*d
def coord(i,d):return tuple(int(j==i) for j in range(d))
def linear(c,B,p):
    assert len(c)==len(B)
    return tuple(sum(a*b for a,b in zip(c,z))%p for z in zip(*B))

def divrem(f,g,p):
    f=list(f);g=list(g)
    while f and not f[-1]:f.pop()
    assert g and g[-1]==1
    while len(f)>=len(g):
        c=f[-1];j=len(f)-len(g)
        for i,v in enumerate(g):f[i+j]=(f[i+j]-c*v)%p
        while f and not f[-1]:f.pop()
    return f
class Extension:
    def __init__(self,p,k):
        self.p=p;self.k=k;self.Q=p**k;self.mod=()
        if k>1:
            for c in product(range(p),repeat=k):
                if not c[0]:continue
                f=c+(1,)
                if all(divrem(f,g+(1,),p) for j in range(1,k//2+1) for g in product(range(p),repeat=j)):
                    self.mod=c;break
            assert self.mod
    def digits(self,x):return tuple((x//self.p**j)%self.p for j in range(self.k))
    def mul(self,x,y):
        p,k=self.p,self.k
        if k==1:return x*y%p
        a,b=self.digits(x),self.digits(y);z=[0]*(2*k-1)
        for i,v in enumerate(a):
            for j,u in enumerate(b):z[i+j]=(z[i+j]+v*u)%p
        for i in range(2*k-2,k-1,-1):
            c=z[i]
            for j,v in enumerate(self.mod):z[i-k+j]=(z[i-k+j]-c*v)%p
        return sum(z[j]*p**j for j in range(k))
def projective(p,d):
    for j in range(d):
        for tail in product(range(p),repeat=d-j-1):yield (0,)*j+(1,)+tail

def prefix_probes(p,l,b):
    assert b>=1 and l>=1
    n=min(b+1,l)
    if p>b:m=n*(l-1)+1;mode='hyperplane'
    else:m=((p**n-1)//(p-1))*(l-1)+1;mode='independent'
    k=1
    while p**k<m:k+=1
    F=Extension(p,k);G=[]
    for t in range(m):
        a=[1]
        for i in range(1,l):a.append(F.mul(a[-1],t))
        digits=[F.digits(x) for x in a]
        for f in projective(p,k):G.append(tuple(dot(f,c,p) for c in digits))
    H={coord(i,l) for i in range(l)}
    for x in G:
        for j in range(1,l+1):
            y=x[:j]+(0,)*(l-j)
            if any(y):H.add(norm(y,p))
    G={norm(x,p) for x in G if any(x)}
    for x in G:
        assert all(not any(y) or norm(y,p) in H for j in range(1,l+1) for y in [x[:j]+(0,)*(l-j)])
    assert len(H)<=l*(1+2*m*m)
    return sorted(H),sorted(G),m,k,mode

def check_hit():
    rng=Random(271126);total=0
    print('PREFIX-COMPLETE HYPERPLANE TESTS')
    for p,l,b,exhaust in [(2,3,1,True),(3,3,2,True),(3,4,2,True),(5,3,2,True),(7,4,3,False),(13,5,2,False),
                           (2,4,2,True),(2,5,3,False),(3,5,3,False)]:
        H,G,m,k,mode=prefix_probes(p,l,b);forms=list(projective(p,l));h=min(b+1,l)
        masks=[sum(1<<i for i,x in enumerate(G) if dot(a,x,p)) for a in forms]
        seq=combinations(range(len(forms)),h) if exhaust else (rng.sample(range(len(forms)),h) for _ in range(1200))
        n=0
        for I in seq:
            if mode=='independent' and rank([forms[i] for i in I],p)<h:continue
            z=(1<<len(G))-1
            for i in I:z&=masks[i]
            assert z;n+=1
        total+=n
        print(f'q={p} dim={l} budget={b} mode={mode} ext={k} nodes={m} terminal={len(G)} probes={len(H)} tests={n}')
    print(f'{total} exact hyperplane-family checks.')

def patterns(K,p,w):
    for j in range(min(w,K)+1):
        for I in combinations(range(K),j):
            for c in product(range(1,p),repeat=j):yield I,c

def syn(H,I,c,p):
    z=zero(len(H[0]))
    for i,a in zip(I,c):z=add(z,scale(H[i],a,p),p)
    return z

def sparse(H,p,w):return {syn(H,I,c,p) for I,c in patterns(len(H),p,w)}
def partitions(K):
    if K==0:yield ();return
    for P in partitions(K-1):
        yield P+((K-1,),)
        for j in range(len(P)):yield P[:j]+(P[j]+(K-1,),)+P[j+1:]
def local_words(H,p,s):return [(I,c) for I,c in patterns(len(H),p,s) if I and not any(syn(H,I,c,p))]
def local_ok(H,P,words,p):
    group={i:j for j,B in enumerate(P) for i in B}
    for I,c in words:
        D=defaultdict(list)
        for i,a in zip(I,c):D[group[i]].append((i,a))
        for v in D.values():
            J,z=zip(*v)
            if any(syn(H,J,z,p)):return False
    return True

def check_partitions():
    print('\nEXACT SMALL LOCALITY AND PREFIX PROBES')
    allp=guard=good=total_words=negative=0
    for p,d,w in [(2,3,1),(3,2,1),(2,3,2)]:
        V=list(projective(p,d));count=passed=trigger=words_n=bad=0
        Hprobe,_,_,_,_=prefix_probes(p,d,w)
        for mask in range(1,1<<len(V)):
            H=[x for i,x in enumerate(V) if mask>>i&1]
            C=local_words(H,p,3*w);words_n+=len(C);S=sparse(H,p,w)
            for P in partitions(len(H)):
                count+=1
                if not local_ok(H,P,C,p):continue
                passed+=1;ranks=sorted((rank([H[i] for i in B],p) for B in P),reverse=True)
                if sum(ranks[:w])<d:assert set(Hprobe)-S;trigger+=1
                elif set(Hprobe)<=S:bad+=1
        print(f'q={p} d={d} w={w}: partitions={count} local={passed} guard={trigger} covered no-guard={bad} short words={words_n}')
        allp+=count;good+=passed;guard+=trigger;negative+=bad;total_words+=words_n
    print(f'{allp} exact partitions, {good} local, {guard} triggered certificates, {negative} covered nontriggering controls; {total_words} short words.')

def info_sets(B,p):
    r=len(B)
    if not r:return []
    R=list(zip(*B));avail=list(range(len(R)));out=[]
    while True:
        I=[];rows=[]
        for j in avail:
            if rank(rows+[R[j]],p)>len(rows):I.append(j);rows.append(R[j])
            if len(rows)==r:break
        if len(rows)<r:break
        out.append(set(I));used=set(I);avail=[j for j in avail if j not in used]
    return out

def support(B):return {j for j,z in enumerate(zip(*B)) if any(z)} if B else set()
def peel(B,p,s):
    R=[info_sets(A,p) for A in B];S=[support(A) for A in B];active=set(range(len(B)));rounds=[]
    touch=[[sum(bool(U & T) for U in R[i]) for T in S] for i in range(len(B))]
    while active:
        remove=[]
        for i in sorted(active):
            if not B[i] or sum(sorted((touch[i][j] for j in active if i!=j),reverse=True)[:s-1])<len(R[i]):remove.append(i)
        if not remove:break
        rounds.append(remove);active.difference_update(remove)
    return rounds,sorted(active),[len(x) for x in R]

def affine_lines(t):
    for a in range(t):
        for b in range(t):yield [x*t+(a*x+b)%t for x in range(t)]
    for b in range(t):yield [b*t+y for y in range(t)]
def line_subspaces(t,p,a):
    assert p>=t;B=[]
    for line in affine_lines(t):
        cols=[]
        for j in range(a):
            v=[0]*(t*t)
            for x,i in enumerate(line):v[i]=pow(x,j,p)
            cols.append(tuple(v))
        B.append(cols)
    return B

def complete_group(B,p):return [norm(linear(x,B,p),p) for x in projective(p,len(B))]
def check_info():
    print('\nDISJOINT-INFORMATION-SET LOCALITY CERTIFICATES')
    for t,p,a,w in [(7,7,2,1),(13,13,2,2)]:
        B=line_subspaces(t,p,a);rounds,J,ns=peel(B,p,3*w)
        assert not J and len(rounds)==1 and min(ns)>=3*w
        H=[h for A in B for h in complete_group(A,p)];l=w*a+1;T=[]
        for x in H:
            if len(T)>=l:break
            if not contains(T,x,p):T.append(x)
        assert len(T)==l
        P,_,_,_,_=prefix_probes(p,l,w);Q=[linear(x,T,p) for x in P]
        if w==1:S={norm(h,p) for h in H};miss=next(x for x in Q if norm(x,p) not in S)
        else:
            S={scale(h,c,p) for h in H for c in range(1,p)}|{zero(len(H[0]))}
            miss=next(x for x in Q if not any(sub(x,y,p) in S for y in S))
        assert any(miss);e_lower=len(H)-len(B)*a
        r=rank([v for A in B for v in A],p)
        print(f'q={p} lines={len(B)} d={t*t} span-rank={r} columns={len(H)} local short-rank lower={e_lower} info sets={min(ns)} probe-dim={l} probes={len(P)} absent checked')
    p=3;d=9
    B=[[coord(0,d),coord(1,d)],[coord(1,d),coord(2,d)],
       [tuple(int(i in (3,4,5)) for i in range(d))],[tuple(int(i in (5,6,7)) for i in range(d))]]
    R,J,_=peel(B,p,3);assert J==[0,1] and set(sum(R,[]))=={2,3}
    H=[h for A in B for h in complete_group(A,p)];ids=[];j=0
    for A in B:
        n=len(complete_group(A,p));ids.append(tuple(range(j,j+n)));j+=n
    P=[ids[i] for i in sum(R,[])]+[tuple(x for i in J for x in ids[i])]
    C=local_words(H,p,3);assert local_ok(H,P,C,p)
    print(f'overlapping-core fixture: {len(H)} columns, {len(C)} short words, removed groups={sum(R,[])}, remaining merged groups={J}')

def check_general_lifts():
    print('\nCOARSENED OVERLAPPING-LIFT CERTIFICATES')
    rng=Random(1927);cases=deps=success=merged=0
    for p in (2,3,5):
        d=6;w=1;t=2;b=t*w
        E,_,_,_,_=prefix_probes(p,d,b)
        for rep in range(20):
            g=6
            if rep%3==0:B=[[coord(i,d)] for i in range(g)]
            elif rep%3==1:B=[[coord(0,d)],[coord(0,d)]]+[[coord(i,d)] for i in range(1,5)]
            else:
                B=[basis([tuple(rng.randrange(p) if rng.randrange(3)==0 else 0 for _ in range(d))
                          for j in range(rng.randint(1,2))],p) for i in range(g)]
            # Build arbitrary one- or two-block atoms and retain their block lifts.
            P=[];H=[]
            for j in range(8):
                I=rng.sample(range(g),rng.randint(1,t));v=[]
                for i,A in enumerate(B):
                    v.append(linear(tuple(rng.randrange(p) for _ in A),A,p) if A and i in I else zero(d))
                P.append(v);H.append(tuple(sum(x[k] for x in v)%p for k in range(d)))
            R,J,_=peel(B,p,3*t*w);groups=[[i] for i in sum(R,[])]+([J] if J else [])
            ranks=[rank([x for i in A for x in B[i]],p) for A in groups]
            C=local_words(H,p,3*w)
            for I,c in C:
                for G in groups:
                    assert all(sum(a*P[j][i][k] for j,a in zip(I,c) for i in G)%p==0 for k in range(d))
                deps+=1
            if sum(sorted(ranks,reverse=True)[:b])<d:
                assert set(E)-sparse(H,p,w);success+=1
            merged+=len(J)>1;cases+=1
    print(f'{cases} overlapping lifts, {deps} exact short-word checks, {merged} nontrivial merges, {success} triggered probe certificates.')

def linekey(x,y,p):
    a,b=list(x),list(y);i=next(i for i,v in enumerate(a) if v);j=next(i for i,v in enumerate(b) if v)
    if j<i:a,b=b,a;i=j
    a=list(norm(a,p));b=list(sub(b,scale(a,b[i],p),p));j=next((j for j,v in enumerate(b) if v),None)
    assert j is not None and j>i
    b=list(norm(b,p));a=list(sub(a,scale(b,a[j],p),p))
    return tuple(a),tuple(b)
class UF:
    def __init__(self,n):self.p=list(range(n))
    def root(self,x):
        while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
        return x
    def join(self,x,y):self.p[self.root(y)]=self.root(x)
    def groups(self):
        D=defaultdict(list)
        for i in range(len(self.p)):D[self.root(i)].append(i)
        return list(D.values())
def dense_line_groups(H,p,minimum):
    T=defaultdict(set)
    for i,j in combinations(range(len(H)),2):T[linekey(H[i],H[j],p)].update((i,j))
    U=UF(len(H))
    for S in T.values():
        if len(S)>=minimum:
            a=min(S)
            for b in S:U.join(a,b)
    return sorted((sorted(v) for v in U.groups()),key=lambda v:v[0]),len(T),sum(len(S)>=minimum for S in T.values())
def vand_blocks(p,a,s,g):
    assert g<=p
    return [[tuple(pow(i,k,p)*int(j==z)%p for k in range(s) for z in range(a)) for j in range(a)] for i in range(g)]
def check_bridges():
    print('\nOVERLAPPING BLOCK LIFTS AND AUTOMATIC DENSE-LINE RECOGNITION')
    for p,a,w,g in [(11,2,1,8),(17,2,2,14)]:
        s=6*w;B=vand_blocks(p,a,s,g);d=s*a;checks=0
        for I in combinations(range(g),s):
            assert rank([x for i in I for x in B[i]],p)==s*a;checks+=1
        H=[h for A in B for h in complete_group(A,p)];m=len(H)//g;core=len(H)
        H += [norm(add(B[i][0],B[i+1][0],p),p) for i in range(g-1)]
        assert len(set(H))==len(H) and all(any(x) for x in H)
        P,nlines,ndense=dense_line_groups(H,p,4)
        expected=[list(range(i*m,(i+1)*m)) for i in range(g)]+[[i] for i in range(core,len(H))]
        assert P==expected
        C,_,_=dense_line_groups(H,p,3);assert len(C)==1
        A=[basis([H[i] for i in I],p) for I in P if len(I)>1]
        assert len(A)==g and all(len(x)==a for x in A);count=0
        for I in P:
            if len(I)>1:continue
            h=H[I[0]];pairs=[(i,j) for i,j in combinations(range(g),2) if contains(A[i]+A[j],h,p)]
            assert len(pairs)==1;count+=1
        assert count==g-1
        b=2*w;l=b*a+1;assert l<=d and p>b
        E,_,_,_,_=prefix_probes(p,l,b);Q=[x+zero(d-l) for x in E]
        if w==1:S={norm(x,p) for x in H};assert any(norm(x,p) not in S for x in Q)
        else:
            S={scale(x,c,p) for x in H for c in range(1,p)}|{zero(d)}
            assert any(not any(sub(x,y,p) in S for y in S) for x in Q)
        print(f'q={p} w={w} d={d} groups={g} columns={len(H)} line spans={nlines} dense lines={ndense} triple components={len(C)} block-independence checks={checks} lifted probes={len(E)}')

def Ball(q,n,w):return sum(comb(n,j)*(q-1)**j for j in range(min(n,w)+1))
def compress_lower(q,K,w,count):
    lo,hi=0,K
    while lo<hi:
        m=(lo+hi)//2
        if Ball(q,m,w)>=count:hi=m
        else:lo=m+1
    return lo

def connected_code_fixture(p,w,L,d,seed):
    rng=Random(seed);a=2*w;D=L*a;assert a<p and d<D
    while True:
        T=[tuple(rng.randrange(p) for _ in range(d)) for _ in range(D)]
        if rank(T,p)==d:break
    B=[T[i*a:(i+1)*a] for i in range(L)];assert all(rank(A,p)==a for A in B)
    H=[norm(linear(tuple(pow(t,j,p) for j in range(a)),A,p),p) for A in B for t in range(p)]
    H += [norm(add(H[i*p],H[(i+1)*p],p),p) for i in range(L-1)]
    assert len(set(H))==len(H) and all(any(x) for x in H)
    circuits=[]
    for i in range(L):
        I=tuple(range(i*p,i*p+a))
        for j in range(i*p+a,(i+1)*p):circuits.append(I+(j,))
    for i in range(L-1):circuits.append((i*p,(i+1)*p,L*p+i))
    U=UF(len(H))
    for I in circuits:
        assert len(I)<=3*w and rank([H[i] for i in I],p)==len(I)-1
        assert all(rank([H[i] for i in I if i!=j],p)==len(I)-1 for j in I)
        for j in I[1:]:U.join(I[0],j)
    assert len(U.groups())==1 and len(circuits)==len(H)-L*a and rank(H,p)==d
    return H,len(circuits)
def check_compression_barrier():
    print('\nCONNECTED HIGH-SHORT-RANK COMPRESSION BARRIER');total=0
    for p,w,L,d,seed in [(7,2,2,7,194),(11,2,3,10,91),(17,2,3,11,731),(7,3,2,11,606),(31,1,4,7,17)]:
        H,e=connected_code_fixture(p,w,L,d,seed);K=len(H);v=Ball(p,K,w);S=sparse(H,p,w);total+=v
        lower=compress_lower(p,K,w,len(S))
        print(f'q={p} w={w} d={d} K={K} independent-short lower={e} span={len(S)}/{v} min radius-preserving replacement={lower} circuit components=1')
    print(f'{total} exact sparse coefficient patterns across five fixtures.')

if __name__=='__main__':
    check_hit();check_partitions();check_info();check_general_lifts();check_bridges();check_compression_barrier()
    print('\nAll checks passed.')
