"""Exact synthetic cyclic-group checks of quotient-sumset certificates.
No cryptographic targets or tensor contractor are used.
"""
from collections import Counter
from itertools import product
from random import Random
from math import prod
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def elimination(rows,r,n,rhs=None):
    piv={}; chosen=[]
    for ix,row in enumerate(rows):
        a=np.asarray([*row, rhs[ix] if rhs is not None else 0],dtype=np.int64)%r
        for j in sorted(piv):
            if a[j]:a=(a-int(a[j])*piv[j])%r
        js=np.flatnonzero(a[:n])
        if len(js):
            j=int(js[0]);piv[j]=a*pow(int(a[j]),-1,r)%r;chosen.append(ix)
        else:assert a[n]==0
    return piv,chosen


def rank_mod(rows,r,n):return len(elimination(rows,r,n)[0])


def nullspace(rows,r,n):
    piv,_=elimination(rows,r,n);out=[]
    for f in sorted(set(range(n))-set(piv)):
        x=np.zeros(n,dtype=np.int64);x[f]=1
        for j in sorted(piv,reverse=True):x[j]=-int(piv[j][j+1:n]@x[j+1:])%r
        out.append([int(a) for a in x])
    return out


def solve(rows,r,n,rhs):
    piv,_=elimination(rows,r,n,rhs);assert len(piv)==n
    x=np.zeros(n,dtype=np.int64)
    for j in sorted(piv,reverse=True):x[j]=(int(piv[j][n])-int(piv[j][j+1:n]@x[j+1:]))%r
    return [int(a) for a in x]


def dot(a,b,r):return sum(int(x)*int(y) for x,y in zip(a,b))%r


def palette_chain(r,random_values):
    vals={0,1};ops=[]
    for n in [r,*random_values]:
        x=1
        for bit in bin(n)[3:]:
            y=2*x%r;ops.append((x,x,y));vals.add(y);x=y
            if bit=='1':
                y=(x+1)%r;ops.append((x,1,y));vals.add(y);x=y
        assert x==n%r
    return sorted(vals),ops


def digits(t,D,m):
    out=[]
    for _ in range(m): t,j=divmod(t,D);out.append(j)
    assert not t
    return out


def make_fixture(r,m,extra,seed,num_random):
    rng=Random(seed)
    D=1<<((r.bit_length()+m-1)//m)
    T,ops=palette_chain(r,rng.sample(range(2,r),num_random))
    ds=[{j*D**i%r for j in range(D+1)} for i in range(m)]
    for zeta in range(1,r):
        cs={(c+zeta)%r for c in T}
        if all(not cs.intersection(s) for s in ds):break
    else:raise ValueError('no zeta')
    used=set.union(cs,*ds)
    tau=next(a for a in range(1,r) if pow(a,-1,r) not in used and -pow(a,-1,r)%r not in used)
    blocks=[];special=D+1+len(T)
    for i in range(m):
        pairs=[(0,j*D**i%r) for j in range(D+1)]+[(c,(c+zeta)%r) for c in T]+[(0,pow(tau,-1,r))]
        F=[((1+tau*t)%r,c,t) for c,t in pairs]
        used={a for a,c,t in F}
        assert 0 not in used and len(used)==len(F)
        for a in rng.sample([a for a in range(1,r) if a not in used],extra):
            F.append((a,rng.choice(T),(a-1)*pow(tau,-1,r)%r))
        blocks.append(F)
    return D,T,ops,zeta,tau,blocks


def genuine_rows(r,m,D,T,ops,zeta,tau,B):
    special=D+1+len(T);idx={c:D+1+j for j,c in enumerate(T)}
    start=[];n=0
    for F in B:start.append(n);n+=len(F)-1
    u=[c for F in B for a,c,t in F[1:]];lam=[(a-1)%r for F in B for a,c,t in F[1:]]
    zero=(0,)*(2*m);rows=[]
    def change(X,updates):
        Y=list(X)
        for i,slot,j in updates:Y[2*i+slot]=j
        return tuple(Y)
    def inc(X):
        v=[0]*n
        for j,k in enumerate(X):
            if k:v[start[j//2]+k-1]+=1
        return v
    def vals(X):return tuple(sum(B[i//2][j][a] for i,j in enumerate(X))%r for a in range(3))
    def add(X,Y):
        g,c,t=vals(X);h,d,v=vals(Y)
        assert (h-g)%r==1 and d==c and (v-t)%r==pow(tau,-1,r)
        row=[b-a for a,b in zip(inc(X),inc(Y))]
        assert sum(a*b for a,b in zip(row,u))%r==0
        assert sum(a*b for a,b in zip(row,lam))%r==1
        rows.append(row)
    def contrast(X,Xp):
        assert vals(X)==vals(Xp)
        j=next(j for j,x in enumerate(X) if not x)
        Y=change(X,[(j//2,j%2,special)])
        add(X,Y);add(Xp,Y)
    def dtup(t):return change(zero,[(i,0,j) for i,j in enumerate(digits(t,D,m))])
    for i,F in enumerate(B):
        add(zero,change(zero,[(i,0,special)]))
        contrast(change(zero,[(i,0,special)]),dtup(pow(tau,-1,r)))
        contrast(change(zero,[(i,0,idx[0])]),dtup(zeta))
        for j in range(2,D+1):
            contrast(change(zero,[(i,0,j)]),change(zero,[(i,0,j-1),(i,1,1)]))
        if i:
            contrast(change(zero,[(i-1,0,D)]),change(zero,[(i,0,1)]))
            for c in T:contrast(change(zero,[(i,0,idx[c])]),change(zero,[(0,0,idx[c])]))
        else:
            for a,b,c in ops:
                contrast(change(zero,[(0,0,idx[a]),(0,1,idx[b])]),change(zero,[(0,0,idx[c]),(0,1,idx[0])]))
        for j,(a,c,t) in enumerate(F[special+1:],special+1):
            X=change(zero,[(i,0,j)])
            Xp=change(dtup((t-c-zeta)%r),[(0,1,idx[c])])
            contrast(X,Xp)
    contrast(zero,dtup(r))
    assert rank_mod(rows,r,n)==n-1
    assert rank_mod(rows,1000003,n)==n # certifies rational rank
    if n<=200:
        H=smith_normal_form(sp.Matrix(rows),domain=sp.ZZ)
        assert [abs(int(H[i,i])) for i in range(n)]==[1]*(n-1)+[r]
    return rows,u,lam,start


def min_range(r,B):
    # Exhaust all nonzero augmented solutions (h,k), minimize lift range over
    # independent circular cuts in each block, including the reference zero.
    ans=r*len(B);per=[r]*len(B);arg=None
    k=np.arange(r,dtype=np.int64)[:,None]
    arrays=[(np.array([c for a,c,t in F],dtype=np.int64)[None,:],np.array([(a-1)%r for a,c,t in F],dtype=np.int64)[None,:]) for F in B]
    for h in range(r):
        scores=np.zeros(r,dtype=np.int64)
        for i,(c,l) in enumerate(arrays):
            x=np.sort((h*l+k*c)%r,axis=1)
            gap=np.maximum(np.diff(x,axis=1).max(axis=1),r+x[:,0]-x[:,-1])
            w=r-gap;scores+=w
            per[i]=min(per[i],int(w[1 if h==0 else 0:].min()))
        if h==0:scores[0]=r*len(B)
        j=int(scores.argmin());v=int(scores[j])
        if v<ans:ans=v;arg=(h,j)
    return ans,arg,per


def sumset_dp(sets,r,cap):
    S={0}
    for A in sets:
        T=set()
        for a in A:
            for x in S:
                T.add((a+x)%r)
                if len(T)>cap:return None
        S=T
    return S


def counts(r,B):
    def conv(Ss):
        x=np.zeros(r,dtype=np.int64);x[0]=1
        for S in Ss:
            y=np.zeros_like(x)
            for a in S:y+=np.roll(x,a)
            x=y
        return x
    F=[[a for a,c,t in Q] for Q in B];W=prod(map(len,F));x=conv(F);assert x.sum()==W and x.min()>0
    P=[conv([[-a%r for a in S] for S in F]+[S for j,S in enumerate(F) if j!=i]) for i in range(len(F))]
    return (int(x.min()),int(x.max())),(min(int(a.min()) for a in P),max(int(a.max()) for a in P))


def short_tuple(r,B,g):
    for X in product(*(range(len(F)) for F in B)):
        if sum(B[i][a][0] for i,a in enumerate(X))%r==g:return X
    raise AssertionError('coverage failed')


def incidence(B,X,copies):
    start=[];n=0
    for F in B:start.append(n);n+=len(F)-1
    v=[0]*n
    for i,j in enumerate(X):
        if j:v[start[i%len(B)]+j-1]+=1
    return v


def verify_completion(r,m,D,T,zeta,tau,B,rows,u,lam,start,copies):
    n=len(u)
    U=nullspace(rows,r,n);assert len(U)==1
    v=U[0];factor=next(v[j]*pow(u[j],-1,r)%r for j in range(n) if u[j]);assert factor and all(v[j]==factor*u[j]%r for j in range(n))
    labelsets=[]
    for i,F in enumerate(B):labelsets.append({0,*v[start[i]:start[i]+len(F)-1]})
    S=sumset_dp(labelsets*copies,r,r)
    assert S is not None
    # Independent pairwise convolution of support at each step.
    brute={0}
    for C in labelsets*copies:brute={ (x+y)%r for x in brute for y in C }
    assert S==brute
    if len(S)==r:return (len(S),None,None,None)
    t=next(j for j in range(r) if j not in S)
    assert sumset_dp(labelsets*copies,r,len(S)-1) is None
    j=next(j for j,a in enumerate(v) if a)
    lift=[0]*n;lift[j]=t*pow(v[j],-1,r)%r
    g=(copies*m+dot(lift,lam,r))%r
    # Global short coverage suffices; pad a short representation with references.
    X=short_tuple(r,B,(g-(copies-1)*m)%r)+(0,)*((copies-1)*m)
    assert sum(B[i%m][a][0] for i,a in enumerate(X))%r==g
    c=incidence(B,X,copies);new=[(a-b)%r for a,b in zip(lift,c)]
    assert dot(new,lam,r)==0 and dot(new,v,r)!=0
    assert rank_mod([*rows,new],r,n)==n
    assert solve([*rows,new],r,n,[1]*len(rows)+[0])==lam
    return len(S),t,g,dot(new,v,r)


def joint_dp(sets,r,cap,d):
    S={(0,)*d}
    for A in sets:
        T=set()
        for a in A:
            for x in S:
                T.add(tuple((j+k)%r for j,k in zip(a,x)))
                if len(T)>cap:return None
        S=T
    return S


def joint_branch_checks():
    for r,d in [(101,2),(101,3),(257,4)]:
        # A covered cyclic toy domain; the known rows are dense, verified
        # equations, not asserted to be a canonical fixed-shift extraction.
        B=[[(a,0,0) for a in range(1,r)]]+[[ (a,0,0) for a in (1,2)] for _ in range(d)]
        m=len(B);n=r-2+d;lam=[(a-1)%r for F in B for a,_,_ in F[1:]]
        rows=[[int(i==j) for i in range(n)] for j in range(r-2)];rhs=[dot(a,lam,r) for a in rows]
        K=nullspace(rows,r,n);assert len(K)==d
        start=0;sets=[]
        for F in B:
            V={(0,)*d}
            for j in range(len(F)-1):V.add(tuple(u[start+j] for u in K))
            start+=len(F)-1;sets.append(V)
        S=joint_dp(sets,r,r-1,d);assert S is not None and len(S)==2**d
        assert joint_dp(sets,r,len(S)-1,d) is None
        for j in range(d):
            V=nullspace(rows,r,n);assert len(V)==d-j
            v=V[0]
            # Initial K is a unit-coordinate basis in this toy; for general
            # K use field elimination to express v in its row span.
            coeff=v[-d:];assert all(v[i]==sum(coeff[k]*K[k][i] for k in range(d))%r for i in range(n))
            vals={dot(coeff,s,r) for s in S};t=next(x for x in range(r) if x not in vals)
            i=next(i for i,x in enumerate(v) if x);z=[0]*n;z[i]=t*pow(v[i],-1,r)%r
            g=(m+dot(z,lam,r))%r;X=short_tuple(r,B,g);a=incidence(B,X,1)
            new=[(x-y)%r for x,y in zip(z,a)]
            assert dot(new,lam,r)==0 and dot(new,v,r)!=0
            rows.append(new);rhs.append(0)
        assert solve(rows,r,n,rhs)==lam
        print(f'Joint support branch: r={r}, initial corank {d}, support {len(S)}; exactly {d} new verified zero-value rows recover full rank.')


def prior_check(path):
    import importlib.util
    spec=importlib.util.spec_from_file_location('prior',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    for r,m,extra,seed in [(101,2,8,17),(257,3,12,22),(1009,2,10,17)]:
        D,zeta,tau,B=mod.dictionary(r,m,extra,seed)
        rows,u,lam,make,vals=mod.selected_rows(r,m,D,zeta,tau,B)
        S=sumset_dp([{c for a,c,t in F} for F in B],r,r);assert len(S)==r
        # The dictionary's pure d and c encodings prove surjectivity of every
        # nonzero augmented projection; verify both short encodings here.
        for x in range(r):
            ds=mod.digits(x,D,m)
            assert sum(B[i][j][2] for i,j in enumerate(ds))%r==x
            assert sum(B[i][D+1+j][1] for i,j in enumerate(ds))%r==x
        print(f'PR12 control r={r}: short homogeneous score support {len(S)}/{r}; every nonzero augmented scalar projection is surjective by the two dictionary encodings.')


def main():
    for r,m,extra,seed,nrand in [(257,2,8,8,2),(1009,2,20,18,4),(2039,2,4,18,2)]:
        D,T,ops,zeta,tau,B=make_fixture(r,m,extra,seed,nrand)
        allrows,u,lam,start=genuine_rows(r,m,D,T,ops,zeta,tau,B);n=len(u)
        _,ids=elimination(allrows,r,n);basis=[allrows[i] for i in ids]
        assert len(basis)==n-1 and rank_mod(basis,1000003,n)==n-1
        _,extraids=elimination(basis+allrows,1000003,n)
        assert len(extraids)==n
        basis=[(basis+allrows)[i] for i in extraids]
        assert len(basis)==n and rank_mod(basis,r,n)==n-1
        R,arg,per=min_range(r,B);assert 2*R>=r
        fibers,pins=counts(r,B)
        frac=[str(sp.Rational(max(Counter(c for a,c,t in F).values()),len(F))) for F in B]
        print(f'r={r}, m={m}, block sizes={[len(F) for F in B]}, palette={len(T)}: {len(allrows)} genuine rows, selected rational rank {n}, modular rank {n-1}.')
        print(f'  Exact minimum augmented circular width {R} at {arg}; per-block minima {per}; largest quotient fractions {frac}.')
        print(f'  Short fibers {fibers}; pin counts {pins}; no-wrap {0<pins[0]<=pins[1]<r}.')
        for cpb in (1,2):
            S,miss,g,score=verify_completion(r,m,D,T,zeta,tau,B,basis,u,lam,start,cpb)
            expect_recovery = cpb == 1 or r == 2039
            assert (miss is not None) == expect_recovery, (
                f'C{r}, copies={cpb}: expected recovery={expect_recovery}, missing residue={miss}'
            )
            if miss is None:print(f'  Copies={cpb}: projection support {S}/{r}, missing-value branch correctly reports no missing value.')
            else:print(f'  Copies={cpb}: projection support {S}/{r}; missing {miss}; group fiber {g}; verified zero-value score {score}; full recovery.')
        if r==2039:assert 0<pins[0]<=pins[1]<r
    joint_branch_checks()
    import sys
    if len(sys.argv)>1:prior_check(sys.argv[1])


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions required')
    main()
