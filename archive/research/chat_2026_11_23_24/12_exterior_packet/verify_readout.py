from verify import F,det,v2,c,ev,ss,mf,rng
from math import prod,comb,ceil
from collections import Counter

def mul(a,b,cap=None):
    d=len(a)+len(b)-1 if cap is None else min(len(a)+len(b)-1,cap)
    out=[F(0)]*d
    for i,x in enumerate(a):
        for j,y in enumerate(b[:max(0,d-i)]):out[i+j]+=x*y
    return out

def at(a,x):
    y=F(0)
    for v in reversed(a):y=y*x+v
    return y

def jet(a,t,alpha,q=None):
    q=len(a) if q is None else q
    return [sum((a[k]*comb(k,j)*t**(k-j)*alpha**j for k in range(j,len(a))),F(0)) for j in range(q)]

def recip(a,q):
    assert a[0]==1
    b=[F(1)]
    for i in range(1,q):b.append(-sum((a[j]*b[i-j] for j in range(1,min(len(a),i+1))),F(0)))
    return b

def check():
    cases=jets=values=0
    for L in (2,4,8):
        a=L.bit_length()-1;alpha=2*L;u=a
        nodes=[t for t in range(-L,L) if t and v2(t)<u]
        for k in range(L+1):
            for trial in range(2):
                d=4*k*(L-k);e=ev(L,k);g0=prod((-s)**e[s] for s in e if e[s]);cp=prod(-s for s in range(-L,L) if s)
                terms=[];G=[F(0)]*(d+1)
                opts=ss(L,k)
                for j in range(12):
                    X=rng.choice(opts);Y=X if j%3==0 else rng.choice(opts);w=rng.randrange(16);terms.append((X,Y,w))
                    ns=Counter(y-x for x in Y for y in X);ex={s:k-ns[s]-e[s] for s in e}
                    vd=prod(x2-x1 for i,x1 in enumerate(X) for x2 in X[i+1:])*prod(y2-y1 for i,y1 in enumerate(Y) for y2 in Y[i+1:])
                    a0=F(vd*g0,cp**k);f=[F(w*2**(2*k))*a0*a0]
                    for s in ex:
                        for _ in range(2*ex[s]):f=mul(f,[-s,1])
                    assert len(f)==d+1
                    for i,x in enumerate(f):G[i]+=x
                Z=sum(w for X,Y,w in terms if X==Y);assert G[0]==2**(2*k)*Z
                # Whole local polynomial integrality, not just integral evaluations.
                for t in nodes:
                    J=jet(G,t,alpha)
                    assert all(x.denominator%2 for x in J);jets+=len(J)
                    for z in (1,3):
                        T=t+alpha*z;C=c(L,T);gt=prod((T-s)**e[s] for s in e)
                        literal=sum((w*det([[C/F(T+x-y) for y in X] for x in Y])**2 for X,Y,w in terms),F(0))
                        assert at(G,T)==2**(2*k)*F(g0,gt)**2*literal;values+=1
                h=ceil((d+1)/sum(u-v2(t) for t in nodes));qs={t:h*(u-v2(t)) for t in nodes};assert sum(qs.values())>d
                P0=8;target=P0+2*k;R=F(0)
                for t,q in qs.items():
                    D=L//2**v2(t);p=target+8*h*D+q*(a+1);Q=2**p
                    E0=F(1);E=[F(1)]
                    for s,m in qs.items():
                        if s==t:continue
                        E0*=F(-s,t-s)**m
                        lin=[F(1),F(alpha,t-s)]
                        for _ in range(m):E=mul(E,lin,q)
                    I=recip(E,q);J=jet(G,t,alpha,q);obs=[F(mf(x,Q)) for x in J]
                    B=mul(obs,I,q);R+=E0*at(B,F(-t,alpha))
                assert R.denominator%2
                assert mf(R,2**target)==2**(2*k)*Z%(2**target)
                assert mf(R,2**target)//(2**(2*k))==Z%(2**P0)
                cases+=1
    print('Complete fixed-cardinality Hermite readouts:',cases)
    print('Locally integral Taylor coefficients:',jets)
    print('Independent normalized determinant evaluation comparisons:',values)
    print('ALL ASSERTIONS PASSED')
if __name__=='__main__':check()
