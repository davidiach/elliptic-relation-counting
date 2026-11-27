#!/usr/bin/env python3
"""Exact toy checks for graded sumset interpolation and Hankel rank certificates."""
from itertools import product, combinations, combinations_with_replacement
from math import comb, prod
from random import Random

if not __debug__:
    raise RuntimeError('Run with assertions enabled')


def exps(p,d):
    return sorted(product(range(p),repeat=d),key=lambda a:(sum(a),a))


def pivots(rows,p):
    B={}
    for row in rows:
        a=[x%p for x in row]
        for i,b in sorted(B.items()):
            if a[i]:
                c=a[i];a=[(x-c*y)%p for x,y in zip(a,b)]
        i=next((i for i,x in enumerate(a) if x),None)
        if i is not None:
            c=pow(a[i],-1,p);B[i]=[x*c%p for x in a]
    return set(B)


def profile(A,p,E):
    return {E[i] for i in pivots(([prod(comb(x,k) for x,k in zip(a,e))%p for e in E] for a in A),p)}


def plus(A,B,p):
    return {tuple((x+y)%p for x,y in zip(a,b)) for a in A for b in B}


def lead_plus(A,B,p):
    return {tuple(x+y for x,y in zip(a,b)) for a in A for b in B if all(x+y<p for x,y in zip(a,b))}


def complete(L,E):
    return min((sum(e)-1 for e in E if e not in L),default=sum(E[-1]))


def check_profiles():
    print('EXHAUSTIVE GRADED PROFILE INCLUSIONS')
    for p,d in [(2,2),(2,3),(3,2)]:
        X=list(product(range(p),repeat=d)); E=exps(p,d);Q=len(X)
        index={x:i for i,x in enumerate(X)}
        trans=[[index[tuple((x+y)%p for x,y in zip(a,b))] for b in X] for a in X]
        subsets=[[]]+[[i for i in range(Q) if z>>i&1] for z in range(1,1<<Q)]
        L=[set()]+[profile([X[i] for i in I],p,E) for I in subsets[1:]]
        C=[-1]+[complete(l,E) for l in L[1:]]
        shifts=[[sum(1<<trans[i][j] for j in I) for I in subsets] for i in range(Q)]
        strict=0;count=0
        for a in range(1,1<<Q):
            assert len(L[a])==len(subsets[a])
            for b in range(a,1<<Q):
                z=0
                for i in subsets[a]:z|=shifts[i][b]
                lower=lead_plus(L[a],L[b],p)
                assert lower<=L[z]
                assert C[z]>=min(d*(p-1),C[a]+C[b])
                strict+=lower!=L[z];count+=1
        print(f'p={p} d={d}: {count} unordered pairs, {strict} strict leading-profile inclusions')
    rng=Random(271126);count=0;strict=0;stronger=0
    for p,d in [(2,4),(3,3),(5,2)]:
        X=list(product(range(p),repeat=d));E=exps(p,d)
        for _ in range(120):
            m=rng.randrange(2,5)
            A=[set(rng.sample(X,rng.randrange(1,min(len(X),9)+1))) for i in range(m)]
            S={tuple(0 for _ in range(d))};L=S.copy();c=0
            for a in A:
                P=profile(a,p,E);L=lead_plus(L,P,p);c+=complete(P,E);S=plus(S,a,p)
            F=profile(S,p,E);assert L<=F
            assert complete(F,E)>=min(d*(p-1),c)
            strict+=L!=F;stronger+=complete(L,E)>min(d*(p-1),c);count+=1
    print(f'{count} seeded multi-factor checks; {strict} strict inclusions; {stronger} refined bounds stronger than scalar completeness')


class Field:
    """Prime fields or small extensions, encoded by base-p coefficient words."""
    def __init__(self,p,mod=()):
        self.p=p;self.mod=tuple(mod);self.s=len(mod) or 1;self.Q=p**self.s
        self.digits=[tuple((a//p**i)%p for i in range(self.s)) for a in range(self.Q)]
        if self.s>1:
            assert all(self.power(a,self.Q-1)==1 for a in range(1,self.Q))
    def add(self,a,b):
        if self.s==1:return (a+b)%self.p
        return sum(((x+y)%self.p)*self.p**i for i,(x,y) in enumerate(zip(self.digits[a],self.digits[b])))
    def neg(self,a):
        return sum((-x%self.p)*self.p**i for i,x in enumerate(self.digits[a]))
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        if self.s==1:return a*b%self.p
        if not a or not b:return 0
        s=self.s;p=self.p;c=[0]*(2*s-1)
        for i,x in enumerate(self.digits[a]):
            for j,y in enumerate(self.digits[b]):c[i+j]=(c[i+j]+x*y)%p
        for i in range(2*s-2,s-1,-1):
            for j,x in enumerate(self.mod):c[i-s+j]=(c[i-s+j]-c[i]*x)%p
        return sum(c[i]*p**i for i in range(s))
    def power(self,a,n):
        y=1
        while n:
            if n&1:y=self.mul(y,a)
            a=self.mul(a,a);n//=2
        return y
    def inv(self,a):
        assert a
        return self.power(a,self.Q-2)
    def div(self,a,b):return self.mul(a,self.inv(b))
    def sum(self,A):
        t=0
        for a in A:t=self.add(t,a)
        return t
    def dot(self,a,b):return self.sum(self.mul(x,y) for x,y in zip(a,b))


def matvec(M,x,F):return tuple(F.dot(a,x) for a in M)
def matmul(A,B,F):return [list(matvec(list(zip(*B)),a,F)) for a in A]
def vadd(a,b,F):return tuple(F.add(x,y) for x,y in zip(a,b))
def vscale(a,z,F):return tuple(F.mul(x,z) for x in a)


def rref(A,F):
    A=[list(a) for a in A];n=len(A);d=len(A[0]);k=0;piv=[]
    for j in range(d):
        i=next((i for i in range(k,n) if A[i][j]),None)
        if i is None:continue
        A[i],A[k]=A[k],A[i];c=F.inv(A[k][j]);A[k]=list(vscale(A[k],c,F))
        for i in range(n):
            if i!=k and A[i][j]:
                c=A[i][j];A[i]=[F.sub(x,F.mul(c,y)) for x,y in zip(A[i],A[k])]
        piv.append(j);k+=1
        if k==n:break
    return A,piv


def inverse(A,F):
    n=len(A);R,piv=rref([list(a)+[int(i==j) for j in range(n)] for i,a in enumerate(A)],F)
    if piv[:n]!=list(range(n)):return None
    return [a[n:] for a in R]


def rank(A,F):return len(rref(A,F)[1])
def curve(t,k,F):return tuple(F.power(t,j) for j in range(k+1))
def hankel(a):
    k=(len(a)-1)//2
    return [a[i:i+k+1] for i in range(k+1)]


def recover(frame,F):
    d=len(frame[0])-1;assert len(frame)==d+3
    V=list(zip(*frame[:d+1]));Vi=inverse(V,F)
    if Vi is None:return None
    a=matvec(Vi,frame[d+1],F);b=matvec(Vi,frame[d+2],F)
    if not all(a) or not all(b):return None
    t=[F.div(y,x) for x,y in zip(a,b)]
    if len(set(t))!=d+1:return None
    C=[]
    for i,u in enumerate(t):
        weight=1
        for j,v in enumerate(t):
            if i!=j:weight=F.mul(weight,F.div(F.neg(v),F.sub(u,v)))
        C.append(vscale(curve(u,d,F),F.div(weight,a[i]),F))
    return matmul(list(zip(*C)),Vi,F)


def select(A,F,k,frame=None):
    # Lists consist of affine points. Select a hidden-curve frame if supplied.
    if frame is None:
        if len(A[0])<2*k+3:return None
        frame=[(1,)+tuple(x) for x in A[0][:2*k+3]]
    T=recover(frame,F)
    if T is None:return None
    Ti=inverse(T,F)
    if Ti is None:return None
    budgets=[max(rank(hankel(matvec(T,(1,)+tuple(x),F)),F) for x in B) for B in A]
    if sum(budgets)>k:return None
    L=Ti[0];v=tuple(int(i==k) for i in range(2*k+1));c=F.dot(L,v)
    if c==0:
        j=next(i for i,x in enumerate(L) if x)
        for a in range(1,F.Q):
            w=list(v);w[j]=F.add(w[j],a)
            if rank(hankel(w),F)==k+1:
                v=tuple(w);c=F.dot(L,v);assert c;break
        else:return None
    target=matvec(Ti,vscale(v,F.div(len(A)%F.p,c),F),F)
    assert target[0]==len(A)%F.p
    assert rank(hankel(matvec(T,target,F)),F)==k+1
    return target[1:],T,budgets


def random_matrix(d,F,rng,affine=False):
    while True:
        M=[[rng.randrange(F.Q) for j in range(d)] for i in range(d)]
        if affine:M[0]=[1]+[0]*(d-1)
        if inverse(M,F) is not None:return M


def make_affine_curve(F,k,M):
    P=[]
    for t in range(F.Q):
        a=matvec(M,curve(t,2*k,F),F)
        if a[0]:P.append(vscale(a,F.inv(a[0]),F)[1:])
    a=matvec(M,(0,)*(2*k)+(1,),F)
    if a[0]:P.append(vscale(a,F.inv(a[0]),F)[1:])
    return P


def check_hidden():
    print('\nHIDDEN HANKEL CERTIFICATES')
    rng=Random(417);checks=0;negative=0;sampled=0;exact=0;adjust=0;remote=0
    for p,k,m in [(17,2,2),(23,3,3),(31,4,3),(37,5,4)]:
        F=Field(p)
        for affine in (True,False):
            for _ in range(8):
                M=random_matrix(2*k+1,F,rng,affine);P=make_affine_curve(F,k,M)
                A=[rng.sample(P,2*k+3)]+[rng.sample(P,3) for i in range(m-1)]
                ans=select(A,F,k);assert ans is not None
                z,T,B=ans;assert B==[1]*m
                Ti=inverse(T,F);adjust+=Ti[0][k]==0
                S={(0,)*(2*k)}
                for a in A:S={vadd(x,y,F) for x in S for y in a}
                assert z not in S;exact+=len(S);checks+=1
                for _ in range(25):
                    xs=[rng.choice(a) for a in A];s=(0,)*(2*k)
                    for x in xs:s=vadd(s,x,F)
                    assert rank(hankel(matvec(T,(m,)+s,F)),F)<=m;sampled+=1
                # The same recovered map handles arbitrary coefficients in a dictionary.
                target=(0,)*k+(1,)+(0,)*k
                for _ in range(10):
                    v=(0,)*(2*k+1)
                    for i in range(k):
                        x=rng.choice(P);v=vadd(v,vscale(matvec(T,(1,)+x,F),rng.randrange(p),F),F)
                    assert rank(hankel(v),F)<=k and v!=target;remote+=1
                if k==m:
                    # Append a point whose recovered Hankel rank exceeds one.
                    bad=next(x for x in (tuple(rng.randrange(p) for i in range(2*k)) for _ in range(100))
                             if rank(hankel(matvec(T,(1,)+x,F)),F)>1)
                    A[-1].append(bad);assert select(A,F,k) is None;negative+=1
    # One rank-two list, with total budget two + one = three.
    F=Field(31);k=3;M=random_matrix(7,F,rng,True);P=make_affine_curve(F,k,M)
    B=[vscale(vadd(P[i],P[i+1],F),F.inv(2),F) for i in range(0,10,2)]
    A=[P[:9],B];ans=select(A,F,k);assert ans and ans[2]==[1,2]
    S={vadd(a,b,F) for a in A[0] for b in A[1]};assert ans[0] not in S
    print(f'{checks} hidden-curve successes, {negative} honest failures, {exact} enumerated sumset points, {sampled} sampled rank bounds, {remote} sparse-span rank bounds')
    print(f'{adjust} affine-slice adjustments; one heterogeneous rank budget [1,2] passed')


def trace(a,F):return F.sum(F.power(a,F.p**i) for i in range(F.s))


def baseword(x,F):return sum((F.digits[a] for a in x),())


def check_norms():
    print('\nNORM-COSET MOMENT CURVES')
    rng=Random(271126);count=0;ds=0;proj=0;char=0;balls=0
    # x^3+x+1 over F_5 and x^3+2 over F_7 are irreducible.
    for F,m in [(Field(5,(1,1,0)),2),(Field(7,(2,0,0)),3)]:
        n=(F.Q-1)//(F.p-1);k=2*m
        A=[[curve(t,k,F)[1:] for t in range(1,F.Q) if F.power(t,n)==c] for c in range(1,m+1)]
        assert all(len(a)==n for a in A)
        # Power-sum uniqueness checked by exact enumeration (small fixtures).
        for I in range(1,1<<m):
            S={(0,)*k};expected=1
            for i,a in enumerate(A):
                if I>>i&1:S={vadd(x,y,F) for x in S for y in a};expected*=n
            assert len(S)==expected;count+=len(S)
        for a in A:
            assert len({vadd(x,y,F) for x,y in combinations_with_replacement(a,2)})==comb(n+1,2);ds+=comb(n+1,2)
            assert len(pivots(([1]+list(baseword(x,F)) for x in a),F.p))==k*F.s+1
            # Exact marginal counts for sampled scalar directions. We use the
            # crude analytic Weil bound, not floating-point Fourier arithmetic.
            for _ in range(20):
                u=tuple(rng.randrange(F.Q) for i in range(k));assert any(u)
                counts=[0]*F.p
                for x in a:
                    y=trace(F.dot(x,u),F);assert y<F.p;counts[y]+=1
                # Orthogonality: squared sum of all nontrivial base characters.
                assert F.p*sum(z*z for z in counts)-n*n <= (F.p-1)*k*k*F.Q
                char+=1
            for _ in range(8):
                u=tuple(rng.randrange(F.Q) for i in range(k));assert any(u)
                got={trace(F.dot(x,u),F) for x in a};assert all(t<F.p for t in got);proj+=1
        ans=select(A,F,m);assert ans is not None and ans[2]==[1]*m
        z,T,_=ans;assert z not in S
        print(f'p={F.p} s={F.s} m={m} n={n} full support={len(S)} hidden rank certificate passed')
    print(f'{count} exact retained-support points, {ds} distinct unordered doubles, {char} sampled character-count checks, {proj} sampled scalar projections')
    for p,d,m,t in [(2,5,2,1),(3,4,2,1),(3,3,2,2),(2,6,3,1)]:
        X=list(product(range(p),repeat=d));B={a for a in X if sum(x!=0 for x in a)<=t}
        S={(0,)*d}
        for _ in range(m):S=plus(S,B,p)
        assert S=={a for a in X if sum(x!=0 for x in a)<=m*t}
        E=exps(p,d);assert complete(profile(S,p,E),E)>=min(m*t,d*(p-1));balls+=1
    print(f'{balls} augmented-ball support and completeness controls')

def row_basis(rows,p):
    B={}
    for row in rows:
        a=[x%p for x in row]
        for i,b in sorted(B.items()):
            if a[i]:
                c=a[i];a=[(x-c*y)%p for x,y in zip(a,b)]
        i=next((i for i,x in enumerate(a) if x),None)
        if i is not None:
            c=pow(a[i],-1,p);B[i]=[x*c%p for x in a]
    return list(B.values())


def evalrow(a,E,p):return [prod(comb(x,k) for x,k in zip(a,e))%p for e in E]


def product_space(A,E,p):
    index={e:i for i,e in enumerate(E)};d=len(E[0]);D=len(E)
    pairs=[(i,j,index[e]) for i,a in enumerate(E) for j,b in enumerate(E)
           if (e:=tuple(x+y for x,y in zip(a,b))) in index]
    def mul(a,b):
        z=[0]*D
        for i,j,k in pairs:z[k]=(z[k]+a[i]*b[j])%p
        return z
    W=[evalrow((0,)*d,E,p)]
    for B in A:
        U=row_basis((evalrow(a,E,p) for a in B),p)
        W=row_basis((mul(x,y) for x in W for y in U),p)
    return W,mul


def kernel_vector(A,p,D):
    F=Field(p);R,P=rref(A,F);free=next((j for j in range(D) if j not in P),None)
    if free is None:return None
    u=[0]*D;u[free]=1
    for i,j in enumerate(P):u[j]=-R[i][free]%p
    assert all(sum(x*y for x,y in zip(a,u))%p==0 for a in A)
    return u


def nonzero_point(E,c,p):
    T={e:a%p for e,a in zip(E,c) if a%p};x=[]
    while len(next(iter(T)))>0:
        deg=max(e[0] for e in T)
        for a in range(deg+1):
            U={}
            for e,v in T.items():U[e[1:]]=(U.get(e[1:],0)+v*comb(a,e[0]))%p
            U={e:v for e,v in U.items() if v}
            if U:T=U;x.append(a);break
        else:raise AssertionError('nonzero polynomial vanished on all trial values')
    assert T[()]!=0
    return tuple(x)


def check_products():
    print('\nEXACT TRUNCATED PRODUCT-SPACE ALGORITHM')
    rng=Random(271126);checks=0;holes=0;nondegree=0;high=0;algebra=0
    for p,d in [(2,4),(3,3),(5,2),(7,2)]:
        X=list(product(range(p),repeat=d));F=exps(p,d)
        for i in range(70):
            m=rng.randrange(2,5)
            A=[rng.sample(X,rng.randrange(1,7)) for j in range(m)]
            if i%2:
                top=rng.sample(F,rng.randrange(1,4));E=[e for e in F if any(all(x<=y for x,y in zip(e,v)) for v in top)]
                nondegree+=1
            else:
                t=rng.randrange(1,d*(p-1)+1);E=[e for e in F if sum(e)<=t]
            high+=max(map(sum,E))>=p
            W,mul=product_space(A,E,p)
            S={(0,)*d}
            for a in A:S=plus(S,a,p)
            U=row_basis((evalrow(x,E,p) for x in S),p)
            assert len(W)==len(U)==len(row_basis(W+U,p))
            for _ in range(3):
                x,y=rng.choice(X),rng.choice(X);z=tuple((a+b)%p for a,b in zip(x,y))
                assert mul(evalrow(x,E,p),evalrow(y,E,p))==evalrow(z,E,p);algebra+=1
            u=kernel_vector(W,p,len(E))
            if u is not None:
                z=nonzero_point(E,u,p);assert z not in S
                assert sum(a*b for a,b in zip(evalrow(z,E,p),u))%p!=0
                holes+=1
            else:assert len(W)==len(E)
            checks+=1
    print(f'{checks} exact product-space comparisons, {nondegree} general lower sets, {high} degree >= characteristic, {algebra} group-algebra identities, {holes} certified absent points')


def shift_rank(values,X,p):
    index={a:i for i,a in enumerate(X)}
    return len(pivots(([values[index[tuple((a+b)%p for a,b in zip(x,y))]] for x in X] for y in X),p))


def check_thick():
    print('\nTHICK-BALL TRANSLATION-COMPLEXITY THRESHOLD')
    rng=Random(417);cases=0;tests=0
    for p,d,L,h in [(2,3,1,1),(3,2,2,1),(3,3,2,1),(3,3,1,1),(5,2,3,0),(2,4,1,2)]:
        X=list(product(range(p),repeat=d));missing=[i for i,x in enumerate(X) if sum(a>=L for a in x)>h]
        threshold=(L+1)**(h+1);vs=[]
        for x in X:vs.append(prod(a-j for a in x[:h+1] for j in range(L))%p)
        assert shift_rank(vs,X,p)==threshold
        if p**len(missing)<=100:
            vectors=product(range(p),repeat=len(missing))
        else:vectors=(tuple(rng.randrange(p) for i in missing) for _ in range(40))
        n=0
        for u in vectors:
            if not any(u):continue
            v=[0]*len(X)
            for i,a in zip(missing,u):v[i]=a
            assert shift_rank(v,X,p)>=threshold;n+=1
        tests+=n;cases+=1
        print(f'p={p} d={d} L={L} h={h}: missing={len(missing)}, threshold={threshold}, {n} separator checks')
    print(f'{cases} sharp threshold fixtures, {tests} exact separator translation ranks')


def packed_core(m,n,p,q):
    base=m*(q-1)+1;cut=(m+1)//2;assert q>=m*n and q>m and base**cut<=p
    out=[]
    for t in range(m*n):
        a=[pow(t,j,q) for j in range(1,m+1)]
        out.append((sum(x*base**j for j,x in enumerate(a[:cut])),sum(x*base**j for j,x in enumerate(a[cut:]))))
    return [out[i*n:(i+1)*n] for i in range(m)]


def check_plane():
    print('\nORIGINAL-PLANE HIGH-GROWTH GRID CONTROLS')
    for m,n,p,q,g in [(2,7,101,17,3),(3,6,4001,19,3),(4,5,10009,23,2)]:
        A=packed_core(m,n,p,q);assert all(len(set(a))==n for a in A)
        base=[{(a*g**i%p,b*g**i%p) for a in range(g) for b in range(g)} for i in range(m)]
        axes=[{(a*g**i%p,0) for a in range(g)}|{(0,a*g**i%p) for a in range(g)} for i in range(m)]
        S={(0,0)}
        for B in A:S=plus(S,B,p)
        assert len(S)==n**m
        assert all(len(plus(a,a,p))==comb(n+1,2) for a in A)
        B=[set(a)|b|c for a,b,c in zip(A,base,axes)];S={(0,0)}
        for C in B:S=plus(S,C,p)
        L=g**m;assert L<p and all((a,b) in S for a in range(L) for b in range(L))
        assert len(S)>=n**m
        print(f'm={m} p={p} core={n} grid={L}x{L} actual support={len(S)} translation dimension >= {L+1}')
    print('3 packed-core and square-grid controls (small fixtures do not impose the asymptotic axis/product guards)')


if __name__=='__main__':
    check_profiles();check_products();check_thick();check_plane();check_hidden();check_norms()
    print('\nAll checks passed.')
