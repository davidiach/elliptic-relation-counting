#!/usr/bin/env python3
"""Exact finite-lane tests. Not a collective trace algorithm or large-field benchmark."""
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import product
from math import comb
from random import Random
from previous_fiber_verifier import Ring,roots_poly,pe,interp,res,s3
rng=Random(2026112421)
class Scalar:
    def __init__(self,p):self.p=p;self.zero=0;self.one=1
    def c(self,a):return a%self.p
    def add(self,a,b):return (a+b)%self.p
    def sub(self,a,b):return (a-b)%self.p
    def mul(self,a,b):return a*b%self.p
    def pw(self,a,k):return pow(a,k,self.p)
class Quotient:
    def __init__(self,R,f):
        self.R=R;self.f=f;self.d=len(f)-1;assert self.d and f[-1]==R.o
        self.zero=self.c(0);self.one=self.c(1);self.x=tuple(R.rem([R.z,R.o],f))
    def c(self,a):return (self.R.c(a),)+(self.R.z,)*(self.d-1)
    def add(self,a,b):return tuple(self.R.add(x,y) for x,y in zip(a,b))
    def sub(self,a,b):return tuple(self.R.sub(x,y) for x,y in zip(a,b))
    def mul(self,a,b):return tuple(self.R.rem(self.R.pm(a,b),self.f))
    def pw(self,a,k):
        z=self.one
        while k:
            if k&1:z=self.mul(z,a)
            k//=2
            if k:a=self.mul(a,a)
        return z

def division(n,x,a,c,R):
    """phi_n=psi_n for odd n, psi_n/(2y) for even n. No division."""
    def sq(v):return R.mul(v,v)
    def cube(v):return R.mul(sq(v),v)
    def sm(*v):
        out=R.zero
        for u in v:out=R.add(out,u)
        return out
    x2=sq(x);x3=R.mul(x2,x);x4=sq(x2);x6=sq(x3)
    D=R.mul(R.c(4),sm(x3,R.mul(R.c(a),x),R.c(c)));D2=sq(D)
    base={0:R.zero,1:R.one,2:R.one,
          3:sm(R.mul(R.c(3),x4),R.mul(R.c(6*a),x2),R.mul(R.c(12*c),x),R.c(-a*a)),
          4:R.mul(R.c(2),sm(x6,R.mul(R.c(5*a),x4),R.mul(R.c(20*c),x3),R.mul(R.c(-5*a*a),x2),R.mul(R.c(-4*a*c),x),R.c(-8*c*c-a**3)))}
    @lru_cache(None)
    def phi(k):
        if k<=4:return base[k]
        j=k//2
        if k%2==0:return R.mul(phi(j),R.sub(R.mul(phi(j+2),sq(phi(j-1))),R.mul(phi(j-2),sq(phi(j+1)))))
        L=R.mul(phi(j+2),cube(phi(j)));R0=R.mul(phi(j-1),cube(phi(j+1)))
        if j%2==0:L=R.mul(D2,L)
        else:R0=R.mul(D2,R0)
        return R.sub(L,R0)
    val=phi(n);states=phi.cache_info().currsize
    assert states<=8*n.bit_length()+10
    return val,states

def points(p,a,c):
    assert (4*a**3+27*c*c)%p
    ys=defaultdict(list)
    for y in range(p):ys[y*y%p].append(y)
    return [(x,y) for x in range(p) for y in ys[(x**3+a*x+c)%p]]
def add(P,Q,p,a):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if (y+v)%p==0:return None
        t=(3*x*x+a)*pow(2*y,-1,p)%p
    else:t=(v-y)*pow(u-x,-1,p)%p
    z=(t*t-x-u)%p;return z,(t*(x-z)-y)%p
def mul(P,n,p,a):
    Q=None
    while n:
        if n&1:Q=add(Q,P,p,a)
        P=add(P,P,p,a);n//=2
    return Q
@lru_cache(None)
def aux(p,ell):
    assert p%ell not in (0,1,ell-1)
    tries=0
    while True:
        tries+=1;a,c=rng.randrange(p),rng.randrange(p)
        if (4*a**3+27*c*c)%p==0:continue
        ps=points(p,a,c);N=len(ps)+1
        if N%ell:continue
        H=next((Q for P in ps if (Q:=mul(P,N//ell,p,a)) is not None),None)
        assert H is not None and mul(H,ell,p,a) is None
        S=tuple(sorted({mul(H,k,p,a)[0] for k in range(1,ell)}));assert len(S)==(ell-1)//2
        return a,c,N,H,S,tries
FIX=[(101,7),(127,11),(257,13),(509,13),(1009,19),(2003,29),(4093,29)]

def support_checks():
    C=Counter();largest=0
    for p,ell in FIX:
        a,c,N,H,S,tr=aux(p,ell);R=Scalar(p);ss=set(S)
        for x in range(p):
            q,k=division(ell,x,a,c,R);assert (q==0)==(x in ss);C['abscissas']+=1
        ps=points(p,a,c)
        for bits in (10,32,64,128,256,512):
            for j in range(12):
                P=rng.choice(ps);n=(1<<(bits-1))|rng.getrandbits(bits-1)|1
                if n%p==0:n+=2
                q,k=division(n,P[0],a,c,R);assert (q==0)==(mul(P,n,p,a) is None)
                C['multiplication_comparisons']+=1;largest=max(largest,k)
        C['curve_proposals']+=tr;C['auxiliary_curves']+=1
        print(f'aux p={p}, ell={ell}, #Eprime={N}, width={len(S)}, trials={tr}',flush=True)
    print('Auxiliary supports and recurrence:',dict(C),'max memoized states:',largest,flush=True)

def fiber_gcd(R,F,ell,a,c,square_filter=None):
    d=len(F)-1;masks,zero=R.mask(F);G=[masks[0]]+[R.z]*d
    for k in range(1,d+1):
        f=R.monic(F,k);Q=Quotient(R,f)
        psi,_=division(ell,Q.x,a,c,Q);xp=Q.sub(Q.pw(Q.x,R.p),Q.x)
        g=R.gcd(f,list(psi));g=R.gcd(g,list(xp))
        if square_filter is not None:
            A,t,aa,cc=square_filter;x=Q.add(Q.c(A),Q.mul(Q.c(t),Q.x))
            fx=Q.add(Q.add(Q.mul(Q.mul(x,x),x),Q.mul(Q.c(aa),x)),Q.c(cc))
            qr=Q.sub(Q.pw(fx,(R.p-1)//2),Q.one)
            g=R.gcd(g,list(qr))
        for j,x in enumerate(g):G[j]=R.add(G[j],R.mul(masks[k],x))
    return G,zero

def moments(R,G,zero,S,K):
    d=len(G)-1;mask,_=R.mask(G);H=R.sm(R.mul(R.c(j),e) for j,e in enumerate(mask))
    C=[R.o]+[R.sm(R.mul(mask[j],G[j-i]) for j in range(i,d+1)) for i in range(1,d+1)];P=[H]
    for j in range(1,K+1):
        v=R.mul(R.c(-j),C[j]) if j<=d else R.z
        for i in range(1,min(d,j-1)+1):v=R.sub(v,R.mul(C[i],P[j-i]))
        P.append(v)
    return [R.add(v,R.mul(zero,R.c(sum(pow(s,j,R.p) for s in S)))) for j,v in enumerate(P)]

def check_lanes(p,ell,a,c,S,fs,d,K=5,square_filter=None):
    R=Ring(p,len(fs));F=list(zip(*[f+[0]*(d+1-len(f)) for f in fs]));G,zero=fiber_gcd(R,F,ell,a,c,square_filter)
    if square_filter is not None:
        A,t,aa,cc=square_filter
        S=tuple(u for u in S if pow(((A+t*u)**3+aa*(A+t*u)+cc)%p,(p-1)//2,p)==1)
    P=moments(R,G,zero,S,K)
    out=[]
    for j,f in enumerate(fs):
        true=[s for s in S if pe(f,s,p)==0];out.append(true)
        assert bool(zero[j])==(not any(f))
        if not zero[j]:assert [x[j] for x in G]==roots_poly(p,true)+[0]*(d-len(true))
        assert [x[j] for x in P]==[sum(pow(s,k,p) for s in true)%p for k in range(K+1)]
    return out

def mixed_lanes():
    C=Counter()
    for _ in range(100):
        p,ell=rng.choice(FIX);a,c,N,H,S,tr=aux(p,ell);d=rng.randrange(1,8);n=rng.randrange(2,7);fs=[]
        for i in range(n):
            if rng.randrange(6)==0:f=[0];C['zero']+=1
            else:
                k=rng.randrange(d+1);rs=[rng.choice(S) if rng.randrange(2) else rng.randrange(p) for _ in range(k)];f=roots_poly(p,rs);v=rng.randrange(1,p);f=[x*v%p for x in f]
                C['repeated']+=len(rs)!=len(set(rs));C['degree_drop']+=k<d
            fs.append(f)
        L=check_lanes(p,ell,a,c,S,fs,d)
        A,t=rng.randrange(p),rng.randrange(1,p);aa,cc=a,c
        L2=check_lanes(p,ell,a,c,S,fs,d,5,(A,t,aa,cc))
        for L0,L1 in zip(L,L2):assert set(L1)<=set(L0)
        C['lanes']+=n;C['roots']+=sum(map(len,L));C['moments']+=6*n;C['fixtures']+=1
        C['square_filtered_lanes']+=n;C['square_filtered_roots']+=sum(map(len,L2));C['square_filtered_moments']+=6*n
    print('Masked torsion gcd:',dict(C),flush=True)

class Curve2:
    """Independent group truth over F_p^2, including mixed twists and branch points."""
    def __init__(self,p,a,c):
        self.p=p;self.a=a;self.c=c;self.nr=next(t for t in range(2,p) if pow(t,(p-1)//2,p)==p-1);self.zero=(0,0)
        sq=defaultdict(list)
        for y in range(p):sq[y*y%p].append(y)
        self.ys={}
        for x in range(p):
            f=(x**3+a*x+c)%p
            self.ys[x]=tuple((y,0) for y in sq[f]) if f in sq else tuple((0,y) for y in sq[f*pow(self.nr,-1,p)%p])
    def A(self,x,y):return ((x[0]+y[0])%self.p,(x[1]+y[1])%self.p)
    def S(self,x,y):return ((x[0]-y[0])%self.p,(x[1]-y[1])%self.p)
    def M(self,x,y):return ((x[0]*y[0]+self.nr*x[1]*y[1])%self.p,(x[0]*y[1]+x[1]*y[0])%self.p)
    def I(self,x):
        v=pow((x[0]*x[0]-self.nr*x[1]*x[1])%self.p,-1,self.p);return (x[0]*v%self.p,-x[1]*v%self.p)
    def plus(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u:
            if self.A(y,v)==self.zero:return None
            t=self.M(self.A(self.M((3,0),self.M(x,x)),(self.a,0)),self.I(self.M((2,0),y)))
        else:t=self.M(self.S(v,y),self.I(self.S(u,x)))
        z=self.S(self.S(self.M(t,t),x),u);return z,self.S(self.M(t,self.S(x,z)),y)
    def truth(self,xs,T,h):
        for ys in product(*(self.ys[x] for x in xs)):
            Q=None
            for x,y in zip(xs,ys):Q=self.plus(Q,((x,0),y))
            P=None
            for _ in range(h):P=self.plus(P,Q)
            if P is not None and P[0]==(T,0):return True
        return False

def sum_poly(xs,p,a,c):
    if len(xs)==1:return [(-xs[0])%p,1]
    if len(xs)==2:return s3(*xs,a,c,p)
    prev=sum_poly(xs[:-1],p,a,c);d=2**(len(xs)-1)
    return interp([res(prev,s3(xs[-1],v,a,c,p),p) for v in range(d+1)],p)
def projection_poly(p,a,c,h,T):
    if h==1:return [(-T)%p,1]
    R=Scalar(p);vv=[]
    for x in range(h*h+1):
        z=[division(j,x,a,c,R)[0] for j in (h-1,h,h+1)];D=4*(x**3+a*x+c)%p
        V=z[1]**2*(D if h%2==0 else 1)%p;W=z[0]*z[2]*(D if h%2 else 1)%p;vv.append(((x-T)*V-W)%p)
    return interp(vv,p)
def elliptic_grids():
    C=Counter()
    for p,ell in FIX[:4]:
        aa,cc,N,H,S,tr=aux(p,ell)
        for _ in range(8):
            while True:
                a,c=rng.randrange(p),rng.randrange(p)
                if (4*a**3+27*c*c)%p:break
            E=Curve2(p,a,c);ps=[P for P in points(p,a,c) if P[1]];m=rng.choice((2,3));h=rng.choice((1,1,2));d=h*h*2**(m-1)
            t=rng.randrange(1,p);A=rng.randrange(p)
            rest_lists=[rng.sample(range(p),2) for _ in range(m-1)]
            T=rng.choice(ps)[0]
            if _%2==0:
                # Force genuine hits, without changing the torsion support.
                while True:
                    leaf=rng.choice(ps);u=rng.choice(S);A=(leaf[0]-t*u)%p
                    others=[rng.choice(ps) for _ in range(m-1)];P=leaf
                    for Q in others:P=add(P,Q,p,a)
                    P=mul(P,h,p,a)
                    if P is not None and P[1]:break
                T=P[0];rest_lists=[[Q[0],rng.randrange(p)] for Q in others];C['planted']+=1
            if _==1:
                # Degree-drop lane with the known target abscissa as the other leaf.
                m=2;h=1;d=2;rest_lists=[[T,(T+1)%p]]
            pr=projection_poly(p,a,c,h,T);rest=list(product(*rest_lists));fs=[]
            def val(z,x):
                if h==1:return res(sum_poly(z,p,a,c),s3(x,T,a,c,p),p)
                return res(sum_poly(z+(x,),p,a,c),pr,p)
            for z in rest:
                f=interp([val(z,(A+t*u)%p) for u in range(d+1)],p);assert any(f);fs.append(f);C['degree_drop']+=f[-1]==0
            L=check_lanes(p,ell,aa,cc,S,fs,d,4)
            Lr=check_lanes(p,ell,aa,cc,S,fs,d,4,(A,t,a,c))
            for z,f,rs,rsr in zip(rest,fs,L,Lr):
                for u in S:
                    x=(A+t*u)%p;assert (pe(f,u,p)==0)==(val(z,x)==0)==E.truth(z+(x,),T,h);C['addresses']+=1;C['positive']+=u in rs
                    lift=any(y[0] and not y[1] for y in E.ys[x]);assert (u in rsr)==(u in rs and lift)
                    C['rational_filter_comparisons']+=1;C['retained_positive']+=u in rsr
                C['fibers']+=1
            C['grids']+=1;C['projected']+=h!=1
    assert C['positive']>=C['planted'] and C['degree_drop']>0 and C['projected']>0
    print('Actual elliptic relation grids:',dict(C),flush=True)

def dense_delta():
    C=Counter()
    for p,ell in FIX[:4]:
        a,c,N,H,S,tr=aux(p,ell);w=len(S)
        for s in S:
            q=roots_poly(p,[x for x in S if x!=s]);q=[x*pow(pe(q,s,p),-1,p)%p for x in q]
            for j in range(w):
                vals=[sum(q[k]*comb(k,j)*pow(-A,k-j,p) for k in range(j,w))%p for A in range(p)]
                assert vals.count(0)<=w-1-j;C['coefficient_root_bounds']+=1;C['translation_evaluations']+=p
        bad=sum(not all(all(roots_poly(p,[(x+A)%p for x in S if x!=s])) for s in S) for A in range(p))
        assert bad<=w*w*(w-1)//2;C['union_bounds']+=1
    print('Dense singleton coefficient tests:',dict(C),flush=True)
if __name__=='__main__':
    assert __debug__
    support_checks();mixed_lanes();elliptic_grids();dense_delta()
    print('All assertions passed. Explicit finite-lane component checks only.')
