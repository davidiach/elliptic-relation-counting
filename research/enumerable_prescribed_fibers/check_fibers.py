"""Exact finite checks for prescribed-fiber covariance and batched mask evaluation.
Only synthetic cyclic groups and one small curve; no asymptotic contractor.
"""
from itertools import combinations, product
from fractions import Fraction as F
from math import prod, lcm
from random import Random
from collections import defaultdict
import numpy as np


def conv(r, blocks):
    a = np.zeros(r, dtype=np.int64); a[0] = 1
    for B in blocks:
        z = np.zeros(r, dtype=np.int64)
        for x in B: z += np.roll(a, x)
        a = z
    return a


def pin_moment_expectations(r, m, num, den, pin_sets):
    masks = []
    for mask in range(1 << (r-1)):
        block = [value for value in range(1, r) if mask & (1 << (value-1))]
        size = len(block)
        masks.append((block, num**size * (den-num)**(r-1-size)))
    sums = {pins: [0]*r for pins in pin_sets}
    for selection in product(masks, repeat=m):
        blocks = [block for block, _ in selection]
        weight = prod(w for _, w in selection)
        for pins, totals in sums.items():
            free_blocks = [block for i, block in enumerate(blocks) if i not in pins]
            mass = prod(map(len, free_blocks))
            squared_error = (r*conv(r, free_blocks) - mass)**2
            for values in product(*(blocks[i] for i in pins)):
                shifted = np.roll(squared_error, sum(values) % r)
                for target, error in enumerate(shifted):
                    totals[target] += weight * int(error)
    scale = den**((r-1)*m) * r**2
    return {pins: tuple(F(total, scale) for total in totals) for pins, totals in sums.items()}


def pairwise_moment_checks():
    cases = 0; pin_cases = 0; max_ratio = F(0)
    for r,m,num,den in [(3,2,1,2),(3,3,1,3),(5,2,1,3),(5,3,1,2)]:
        q = F((r-1)*num,den); mu = q**m/r
        # The bound requires at least two unpinned blocks.
        pin_sets = (I for s in range(min(2, m-2)+1) for I in combinations(range(m), s))
        for I, moments in pin_moment_expectations(r, m, num, den, pin_sets).items():
            for ex in moments:
                assert ex <= 2*mu
                max_ratio=max(max_ratio,ex/mu)
                if I: pin_cases+=1
                else: cases+=1
    print(f'Pairwise moment bounds: {cases} unpinned and {pin_cases} selected-pin expectations checked exactly; maximum E/mean {max_ratio}.')

    r, m = 3, 4
    pin_pairs = list(combinations(range(m), 2))
    moments = pin_moment_expectations(r, m, 1, 2, pin_pairs)
    bound = 2 * F(r-1, 2)**m / r
    expected = F(1, 8)
    for pins in pin_pairs:
        assert moments[pins] == (expected,)*r, (pins, moments[pins])
        assert all(value <= bound for value in moments[pins])
    print(f'C{r} four-block, two-pin moments: {len(pin_pairs)} pin pairs x {r} targets; all equal {expected}, bound {bound}.')


def curve_affine_checks():
    p=13
    pts=[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-2)%p==0]
    assert len(pts)==18
    checks=0
    for T in (1,3,7):
        for parity in (0,1):
            X=[x for x,y in pts if y%2==parity]; assert len(X)==len(set(X))
            seeds=list(product(range(p),repeat=2))
            for x in X:
                assert sum((a*x+c)%p<T for a,c in seeds)==p*T; checks+=1
            for x,y in combinations(X,2):
                assert sum((a*x+c)%p<T and (a*y+c)%p<T for a,c in seeds)==T*T; checks+=1
    print(f'Actual E/F13, order 19: {checks} exact affine-mask marginal and within-parity pair checks (including zero slopes).')


def exact_pin_data(r,Bs,t):
    m=len(Bs); n=[len(B) for B in Bs]; W=prod(n); d0=conv(r,Bs); d=int(d0[t]); assert d
    cache={I:conv(r,[Bs[j] for j in range(m) if j not in I]) for s in (1,2) for I in combinations(range(m),s)}
    errors={}; masses={}
    for I,a in cache.items():
        w=W//prod(n[i] for i in I)
        z=[]; zcnt=[]
        for P in product(*(Bs[i] for i in I)):
            c=int(a[(t-sum(P))%r]); zcnt.append(c); z.append(r*c-w)
        errors[I]=sum(x*x for x in z); masses[I]=zcnt
    return W,d,errors,masses


def pin_identity_checks():
    rng=Random(513); checks=0
    for r,m,k in [(5,2,3),(7,3,3),(11,3,4)]:
        for _ in range(5):
            Bs=[rng.sample(range(r),k) for i in range(m)]
            X=list(product(*Bs)); vals=[sum(x)%r for x in X]
            for t in [0,1,r-1]:
                if t not in vals: continue
                W,d,E,ms=exact_pin_data(r,Bs,t)
                for I in E:
                    lhs=sum((r*int(a==t)-1)*(r*int(b==t)-1) for x,a in zip(X,vals) for y,b in zip(X,vals) if all(x[i]==y[i] for i in I))
                    assert lhs==E[I]; checks+=1
    print(f'Shared-pin count identity: {checks} exact direct pair comparisons.')


def integer_pd(a):
    # Fraction-free symmetric elimination. All positive leading minors iff PD.
    A=[list(map(int,x)) for x in a]; n=len(A); old=1
    for k in range(n):
        piv=A[k][k]
        if piv<=0:return False
        for i in range(k+1,n):
            for j in range(i,n):
                z=A[i][j]*piv-A[i][k]*A[k][j]
                assert z%old==0
                A[i][j]=A[j][i]=z//old
        for i in range(k+1,n): A[k][i]=A[i][k]=0
        old=piv
    return True


def covariance_check(r,Bs,t,check_pd=True):
    m=len(Bs); ns=[len(B) for B in Bs]; W,d,E,masses=exact_pin_data(r,Bs,t)
    # For equal-size fixtures q is exactly the common size; this is a test
    # of the deterministic conclusion, not the randomized mean.
    assert len(set(ns))==1; q=ns[0]; mu=F(W,r)
    phi=F((r*d-W)**2,r*r)/(mu*mu)+F(q*q,r*r)*sum(E.values())/(mu*mu)
    coords=[(i,j) for i,B in enumerate(Bs) for j in range(1,len(B))]; D=len(coords)
    si=[masses[(i,)][j] for i,j in coords]
    C=[]; V=[]; scale=lcm(*(n*n for n in ns))
    for a,(i,j) in enumerate(coords):
        crow=[];vrow=[]
        for b,(k,l) in enumerate(coords):
            if i==k:
                Q=si[a] if j==l else 0
                vv=scale*((ns[i] if j==l else 0)-1)//(ns[i]**2)
            elif i<k: Q=masses[(i,k)][j*ns[k]+l];vv=0
            else: Q=masses[(k,i)][l*ns[i]+j];vv=0
            crow.append(d*Q-si[a]*si[b]); vrow.append(vv)
        C.append(crow);V.append(vrow)
    if check_pd:
        # C/d^2 and V/scale are the exact conditional and independent covariances.
        lo=[[2*scale*C[i][j]-d*d*V[i][j] for j in range(D)] for i in range(D)]
        hi=[[2*d*d*V[i][j]-scale*C[i][j] for j in range(D)] for i in range(D)]
        assert integer_pd(lo) and integer_pd(hi)
    return d,phi,D


def fiber_rows(r,Bs,t):
    offsets=[]; n=0
    for B in Bs: offsets.append(n);n+=len(B)-1
    last={x:i for i,x in enumerate(Bs[-1])};m=len(Bs)
    for J in product(*(range(len(B)) for B in Bs[:-1])):
        s=(t-sum(Bs[i][j] for i,j in enumerate(J)))%r
        if s not in last: continue
        row=np.zeros(n,dtype=np.int64)
        for i,j in enumerate((*J,last[s])):
            if j: row[offsets[i]+j-1]=1
        yield row


def fiber_rank(r,Bs,t):
    D=sum(len(B)-1 for B in Bs);piv={}; gen=fiber_rows(r,Bs,t); ref=next(gen,None)
    if ref is None:return 0,0
    cnt=1
    for x in gen:
        cnt+=1; a=(x-ref)%r
        for j in sorted(piv):
            if a[j]:a=(a-int(a[j])*piv[j])%r
        inds=np.flatnonzero(a)
        if len(inds):
            j=int(inds[0]);piv[j]=a*pow(int(a[j]),-1,r)%r
        if len(piv)==D-1:
            # Upper bound D-1 follows from true scalar difference direction.
            assert any((x-B[0])%r for B in Bs for x in B)
            return len(piv),cnt
    return len(piv),cnt


def fiber_checks():
    pd_cases=0; guard_cases=0; ranks=0; full=0; minimum=1; rng=Random(319)
    for r,m in [(7,6),(11,6),(13,6)]:
        Bs=[list(range(1,r)) for _ in range(m)]
        for t in (0,1):
            d,phi,D=covariance_check(r,Bs,t)
            assert phi<F(1,64)**2;guard_cases+=1;pd_cases+=2
            rank,cnt=fiber_rank(r,Bs,t);assert rank==D-1;full+=1;ranks+=1
    for r,m,k,reps in [(19,4,9,6),(31,4,12,8),(101,4,12,8),(101,3,6,8)]:
        successes=0;tested=0; deficits=[]
        for _ in range(reps):
            Bs=[rng.sample(range(1,r),k) for i in range(m)];d=conv(r,Bs);D=m*(k-1)
            for t in (0,1,2):
                if d[t]==0:continue
                rank,cnt=fiber_rank(r,Bs,t);assert rank<=D-1;tested+=1;ranks+=1;successes+=rank==D-1;full+=rank==D-1;deficits.append(D-1-rank)
        print(f'Explicit cyclic fibers r={r}, m={m}, n={k}: {successes}/{tested} attain maximum modular contrast rank; largest deficiency from maximum {max(deficits)}.')
    print(f'Covariance guard: {guard_cases} qualifying dense fixtures, {pd_cases} exact positive-definiteness inequalities; modular maximum rank in {full}/{ranks} total nonempty fibers.')


def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a or [0]


def add(a,b,p):return trim([( (a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%p for i in range(max(len(a),len(b)))])
def neg(a,p):return [(-x)%p for x in a]
def mul_naive(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
    return trim(out)


def inverse_matrix(A,p):
    n=len(A); M=[[x%p for x in row]+[int(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k=next(i for i in range(j,n) if M[i][j]);M[j],M[k]=M[k],M[j]; inv=pow(M[j][j],-1,p);M[j]=[x*inv%p for x in M[j]]
        for i in range(n):
            if i!=j:
                c=M[i][j];M[i]=[(x-c*y)%p for x,y in zip(M[i],M[j])]
    return [row[n:] for row in M]


class Toom:
    def __init__(self,p,k):
        self.p=p;self.k=k;assert p>2*k-2
        self.inv=inverse_matrix([[pow(x,j,p) for j in range(2*k-1)] for x in range(2*k-1)],p);self.calls=0
    def mul(self,a,b):
        p=self.p;k=self.k; self.calls+=1; n=max(len(a),len(b))
        if n<=2*k:return mul_naive(a,b,p)
        s=(n+k-1)//k
        A=[list(a[i*s:(i+1)*s])+[0]*max(0,s-len(a[i*s:(i+1)*s])) for i in range(k)]
        B=[list(b[i*s:(i+1)*s])+[0]*max(0,s-len(b[i*s:(i+1)*s])) for i in range(k)]
        values=[]
        for x in range(2*k-1):
            a0=[sum(A[j][i]*pow(x,j,p) for j in range(k))%p for i in range(s)]
            b0=[sum(B[j][i]*pow(x,j,p) for j in range(k))%p for i in range(s)]
            c=self.mul(a0,b0);values.append(c+[0]*(2*s-1-len(c)))
        out=[0]*(2*k*s-1)
        for j in range(2*k-1):
            for i in range(2*s-1):out[j*s+i]=(out[j*s+i]+sum(self.inv[j][l]*values[l][i] for l in range(2*k-1)))%p
        return trim(out)
    def series_inv(self,f,n):
        p=self.p;assert f[0]==1; g=[1];k=1
        while k<n:
            k=min(2*k,n);t=self.mul(f[:k],g)[:k];t=t+[0]*(k-len(t));t=[(-x)%p for x in t];t[0]=(t[0]+2)%p
            g=self.mul(g,t)[:k];g+= [0]*(k-len(g))
        assert self.mul(f[:n],g)[:n]==[1] or (self.mul(f[:n],g)+[0]*n)[:n]==[1]+[0]*(n-1)
        return g
    def remainder(self,a,b):
        p=self.p;a=trim(a);b=trim(b);assert b[-1]==1
        if len(a)<len(b):return a
        n=len(a)-len(b)+1;g=self.series_inv(b[::-1],n);q=self.mul(a[::-1][:n],g)[:n];q+= [0]*(n-len(q));q=q[::-1]
        r=add(a,neg(self.mul(q,b),p),p);assert len(r)<len(b) or r==[0]
        return r
    def eval_many(self,f,X):
        p=self.p;assert X
        def tree(xs):
            if len(xs)==1:return [(-xs[0])%p,1],None,None
            k=len(xs)//2;A=tree(xs[:k]);B=tree(xs[k:]);return self.mul(A[0],B[0]),A,B
        def descend(a,T):
            a=self.remainder(a,T[0]);
            if T[1] is None:return [a[0]]
            return descend(a,T[1])+descend(a,T[2])
        return descend(f,tree(X))


def eval_horner(f,x,p):
    u=0
    for a in reversed(f):u=(u*x+a)%p
    return u


def multipoint_checks():
    rng=Random(329);mul=0;ev=0;reps=0;calls=0
    for p in (101,1009,65537):
        for k in (2,3,5,8):
            T=Toom(p,k)
            for n in (1,3,7,16,33,64,129):
                a=[rng.randrange(p) for _ in range(n)];b=[rng.randrange(p) for _ in range(n+1)];assert T.mul(a,b)==mul_naive(a,b,p);mul+=1
                X=[rng.randrange(p) for _ in range(n)];X[-1]=X[0];reps+=n>1
                f=[rng.randrange(p) for _ in range(n+1)];f[-1]=1
                assert T.eval_many(f,X)==[eval_horner(f,x,p) for x in X];ev+=1
            calls+=T.calls
    print(f'Multipoint compiler arithmetic: {mul} Toom-Cook multiplication comparisons, {ev} product/remainder-tree evaluations ({reps} with repeated points), {calls} recursive multiplication calls.')


if __name__=='__main__':
    pairwise_moment_checks();curve_affine_checks();pin_identity_checks();fiber_checks();multipoint_checks()
