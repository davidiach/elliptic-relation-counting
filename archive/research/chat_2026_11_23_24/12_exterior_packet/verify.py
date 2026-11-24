from fractions import Fraction as F
from itertools import combinations
from math import prod, log2
from collections import Counter
from random import Random
rng=Random(20261124)
def v2(x):
    assert x; x=abs(x);return (x&-x).bit_length()-1
def det(A):
    # Exact rational elimination, or integer fraction-free elimination.
    A=[list(r) for r in A];n=len(A);sign=1
    if not n:return 1
    if any(isinstance(x,F) for r in A for x in r):
        A=[[F(x) for x in r] for r in A];z=F(1)
        for i in range(n):
            j=next((j for j in range(i,n) if A[j][i]),None)
            if j is None:return 0
            if j!=i:A[i],A[j]=A[j],A[i];sign=-sign
            p=A[i][i];z*=p
            for j in range(i+1,n):
                q=A[j][i]/p
                for k in range(i+1,n):A[j][k]-=q*A[i][k]
                A[j][i]=0
        return sign*z
    prev=1
    for i in range(n-1):
        j=next((j for j in range(i,n) if A[j][i]),None)
        if j is None:return 0
        if j!=i:A[i],A[j]=A[j],A[i];sign=-sign
        p=A[i][i]
        for j in range(i+1,n):
            for k in range(i+1,n):
                x=A[j][k]*p-A[j][i]*A[i][k];assert x%prev==0;A[j][k]=x//prev
            A[j][i]=0
        prev=p
    return sign*A[-1][-1]
def mf(x,Q):
    x=F(x);assert x.denominator%2;return x.numerator*pow(x.denominator,-1,Q)%Q
def c(L,T):return F(prod(T-s for s in range(-L,L)),prod(-s for s in range(-L,L) if s))
def ev(L,k):return {s:max(0,k-L+abs(s)) for s in range(-L,L)}
def ss(L,k):return list(combinations(range(L),k))
def check_common():
    n=tight=0
    for L in (2,4,8,16):
        cases=[(X,Y) for k in range(L+1) for X in ss(L,k) for Y in ss(L,k)] if L<=4 else []
        if L>4:
            for _ in range(120):
                k=rng.randrange(L+1);cases.append((tuple(sorted(rng.sample(range(L),k))),tuple(sorted(rng.sample(range(L),k)))))
        for X,Y in cases:
            k=len(X);e=ev(L,k);ns=Counter(y-x for x in Y for y in X);ex={s:k-ns[s]-e[s] for s in e}
            assert min(ex.values())>=0 and sum(ex.values())==2*k*(L-k)
            T=L+rng.randrange(1,L+5);C=c(L,T);d=det([[C/F(T+x-y) for y in X] for x in Y])
            vd=prod(x2-x1 for x1,x2 in combinations(X,2))*prod(y2-y1 for y1,y2 in combinations(Y,2))
            cp=prod(-s for s in range(-L,L) if s);g0=prod((-s)**e[s] for s in e if e[s]);gt=prod((T-s)**e[s] for s in e)
            stripped=F(vd*g0,cp**k)*prod((T-s)**ex[s] for s in ex)
            assert (d*F(g0,gt))**2==stripped**2
            at0=F(vd*g0,cp**k)*prod((-s)**ex[s] for s in ex)
            assert at0**2==int(X==Y);n+=1
        for k in range(L+1):
            for t in range(-L,L):
                pairs=[(x,x+t) for x in range(L) if 0<=x+t<L][:k];Y={x for x,y in pairs};X={y for x,y in pairs}
                Y.update([x for x in range(L) if x not in Y][:k-len(Y)]);X.update([y for y in range(L) if y not in X][:k-len(X)])
                assert k-sum(x+t in X for x in Y)==ev(L,k)[t];tight+=1
    print('Exact rational minor, reduced degree and origin checks:',n)
    print('Tight common-root multiplicity checks:',tight)
def check_val():
    levels=nodes=direct=0;worst=F(0)
    for L in (2,4,8,16,32,64):
        a=L.bit_length()-1
        for k in range(L+1):
            e=ev(L,k);vals=[0]*(2*L);assert sum(e.values())==k*k
            for j in range(1,a+2):
                q=2**j;cnt=[0]*q
                for s,w in e.items():cnt[s%q]+=w
                bins=[0]*q
                for x in range(k):bins[x%q]+=1
                for d in range(q):
                    assert cnt[d]==sum(bins[x]*bins[(x+d+L)%q] for x in range(q))
                    if j<=a:assert cnt[d]<=cnt[0];levels+=1
                for t in range(-L,L):vals[t+L]+=cnt[t%q]
            for t in range(-L,L):
                delta=vals[t+L]-vals[L];assert delta<=k;nodes+=1
                if k:worst=max(worst,F(delta,k))
                if L<=16:
                    base=sum(w*v2(-s) for s,w in e.items() if w)
                    for z in (1,3):
                        assert sum(w*v2(t+2*L*z-s) for s,w in e.items() if w)-base==delta;direct+=1
    print('Residue autocorrelation / maximum comparisons:',levels)
    print('Normalization valuation inequalities:',nodes)
    print('Direct valuation comparisons:',direct,'maximum positive loss/k:',worst)
def modes_at(L,i,t,P):
    Q=2**P;T=t+2*L;C=c(L,T);q=2**(i+1);row=lambda x:(x>>i)&1;col=lambda y:1-((y>>i)&1)
    W=[[mf(C/F(T+x-y),Q) for y in range(L)] for x in range(L)];A=[[W[x][y] if row(x)==col(y) else 0 for y in range(L)] for x in range(L)];ms=[]
    for y0 in range(q):
        sg=col(y0);tg=1-sg
        for j in range(P-1):
            w=j+1;u=[mf(C/F(T-y0+x)*F(q,T-y0+x)**j/F(2**w),Q) if row(x)==tg else 0 for x in range(L)]
            v=[pow((y-y0)//q,j,Q) if (y-y0)%q==0 else 0 for y in range(L)]
            if any(u) and any(v):ms.append((w,sg,tg,u,v))
    for x in range(L):
        for y in range(L):assert (A[x][y]+sum(2**w*u[x]*v[y] for w,sg,tg,u,v in ms))%Q==W[x][y]
    return W,A,ms,row,col
def sm(ms,P,k):
    out=[]
    def rec(j,left,bits):
        out.append((P-left,tuple(bits)))
        if len(bits)==k:return
        for x in range(j,len(ms)):
            w=ms[x][0]
            if w<left:rec(x+1,left-w,bits+[x])
    rec(0,P,[]);return out
def par(ls):return sum(x>y for i,x in enumerate(ls) for y in ls[i+1:])%2
def local(A,ms,I,Y,X,g):
    Is=[j for j in I if ms[j][1]==g];It=[j for j in I if ms[j][2]==g]
    if len(Y)+len(Is)!=len(X)+len(It):return 0
    return det([[A[y][x] for x in X]+[ms[j][3][y] for j in It] for y in Y]+[[ms[j][4][x] for x in X]+[0]*len(It) for j in Is])
def sign(ms,I,ny,nx):return (-1)**(len(I)+par([0]*ny[0]+[1]*ny[1]+[ms[j][1] for j in I])+par([0]*nx[0]+[1]*nx[1]+[ms[j][2] for j in I]))
def border(A,ms,I,Y,X):return (-1)**len(I)*det([[A[y][x] for x in X]+[ms[j][3][y] for j in I] for y in Y]+[[ms[j][4][x] for x in X]+[0]*len(I) for j in I])
def check_modes():
    entries=terms=paired=n=0
    for L,i,P in [(2,0,5),(4,0,6),(4,1,6),(8,0,5),(8,1,5),(8,2,5)]:
        for t in (2**i,-2**i):
            W,A,ms,row,col=modes_at(L,i,t,P);Q=2**P;entries+=L*L
            for _ in range(8):
                k=rng.randrange(L+1);X=tuple(sorted(rng.sample(range(L),k)));Y=tuple(sorted(rng.sample(range(L),k)))
                Xg=[tuple(x for x in X if col(x)==g) for g in (0,1)];Yg=[tuple(y for y in Y if row(y)==g) for g in (0,1)];nx=list(map(len,Xg));ny=list(map(len,Yg));bas=(-1)**(par([row(y) for y in Y])+par([col(x) for x in X]));S=sm(ms,P,k);total=0;cache={}
                for w,I in S:
                    aa=border(A,ms,I,Y,X)%Q;bb=sign(ms,I,ny,nx)*prod(local(A,ms,I,Yg[g],Xg[g],g) for g in (0,1))%Q
                    assert aa==bas*bb%Q;cache[I]=(aa,bb);total=(total+2**w*aa)%Q;terms+=1
                assert total==det([[W[y][x] for x in X] for y in Y])%Q;n+=1
                for _ in range(12):
                    _,I=rng.choice(S);_,J=rng.choice(S);aa,bb=cache[I];cc,dd=cache[J];assert aa*cc%Q==bb*dd%Q;paired+=1
    print('Actual Cauchy crossing entries:',entries)
    print('Bordered determinant / grouped-sign terms:',terms)
    print('Complete minor expansions:',n)
    print('Two-copy sign cancellation checks:',paired)
def check_joins():
    n=terms=0
    for L,i,P in ((2,0,5),(4,0,5),(4,1,5)):
        W,A,ms,row,col=modes_at(L,i,2**i,P);Q=2**P
        ys=[tuple(y for y in range(L) if row(y)==g) for g in (0,1)];xs=[tuple(x for x in range(L) if col(x)==g) for g in (0,1)]
        yo=[sum((list(combinations(q,k)) for k in range(len(q)+1)),[]) for q in ys];xo=[sum((list(combinations(q,k)) for k in range(len(q)+1)),[]) for q in xs]
        for _ in range(2):
            weights=[{(Y,X):rng.randrange(Q) for Y in yo[g] for X in xo[g]} for g in (0,1)]
            for k in range(L+1):
                direct=sum(prod(weights[g][tuple(y for y in Y if row(y)==g),tuple(x for x in X if col(x)==g)] for g in (0,1))*det([[W[y][x] for x in X] for y in Y])**2 for X in ss(L,k) for Y in ss(L,k))%Q
                actual=0;S=sm(ms,P,k)
                for ny0 in range(len(ys[0])+1):
                    ny=[ny0,k-ny0]
                    if not 0<=ny[1]<=len(ys[1]):continue
                    for nx0 in range(len(xs[0])+1):
                        nx=[nx0,k-nx0]
                        if not 0<=nx[1]<=len(xs[1]):continue
                        for w,I in S:
                            for v,J in S:
                                if w+v>=P:continue
                                f=sign(ms,I,ny,nx)*sign(ms,J,ny,nx)*2**(w+v)
                                for g in (0,1):f=f*sum(weights[g][Y,X]*local(A,ms,I,Y,X,g)*local(A,ms,J,Y,X,g) for Y in yo[g] if len(Y)==ny[g] for X in xo[g] if len(X)==nx[g])%Q
                                actual=(actual+f)%Q;terms+=1
                assert actual==direct;n+=1
    print('Full squared-minor group contractions:',n)
    print('Full grouped contraction summands:',terms)
def catalogue(c,P,distinct=False):
    dp=[0]*P;dp[0]=1
    for w in range(1,P):
        for _ in range(2*c if distinct else c):
            for j in (range(P-1,w-1,-1) if distinct else range(w,P)):dp[j]+=dp[j-w]
    return sum(dp)
if __name__=='__main__':
    check_common();check_val();check_modes();check_joins();
    print('Illustrative catalogue (c=2, half precision for exterior, not runtime):')
    for P in (20,40,80,160):print(f' P={P:3}: unrestricted log2={log2(catalogue(2,P)):.6f}, exterior log2={log2(catalogue(2,P//2,True)):.6f}')
    print('ALL ASSERTIONS PASSED')
