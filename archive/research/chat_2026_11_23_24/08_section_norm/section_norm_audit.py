"""Finite checks for the section norm compiler. Not a runtime benchmark.
Requires SymPy; prior toy principal-ideal construction is reused.
"""
from itertools import combinations, product
from random import Random
from math import gcd
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from polynomial_ring_verifier import add_ec, fold, ev, add, scale, lll_exact
from short_projector_audit import toy_construct


def section_degs(m):
    d=m+1
    return d//2,(d-3)//2


def coefficients(m):
    da,db=section_degs(m)
    slots=[(0,i) for i in range(da+1)]+[(1,i) for i in range(db+1)]
    leading=(0,da) if (m+1)%2==0 else (1,db)
    return [v for v in slots if v!=leading],leading


def lists_from(m,cs):
    da,db=section_degs(m); A=[0]*(da+1); B=[0]*(db+1)
    names,lead=coefficients(m)
    assert len(names)==m
    (A,B)[lead[0]][lead[1]]=1
    for (t,j),c in zip(names,cs):(A,B)[t][j]=c
    return A,B


def peval(A,x,p):
    v=0
    for a in reversed(A):v=(v*x+a)%p
    return v


def sum_ec(Ps,p,a):
    S=None
    for P in Ps:S=add_ec(S,P,p,a)
    return S


def model(p,a,c):
    u=c*pow(a,-1,p)%p if a and c else 1
    f=[c,a*u%p,0,pow(u,3,p)]
    return u,f


def finite_sections():
    checked=norm_good=direct_good=target_checks=0
    for p,a,c in [(7,0,3),(11,1,0),(13,4,1),(17,2,2),(19,1,1)]:
        assert (4*a**3+27*c*c)%p
        G=[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-a*x-c)%p==0]
        u,f=model(p,a,c); ui=pow(u,-1,p); fs=[peval(f,z,p) for z in range(p)]
        h=2%p;s=3%p; D=h*s*s%p
        for m in (2,3,4):
            da,db=section_degs(m); d=m+1
            direct=set()
            for Ps in combinations(G,d):
                if len({P[0] for P in Ps})<d:continue
                if sum_ec(Ps,p,a) is None:direct.add(tuple(sorted(Ps)))
            got=set()
            powers=[[pow(x,j,p) for j in range(da+1)] for x in range(p)]
            for cs in product(range(p),repeat=m):
                A,B=lists_from(m,cs)
                roots=[]
                for z in range(p):
                    az=sum(x*y for x,y in zip(A,powers[z]))%p
                    bz=sum(x*y for x,y in zip(B,powers[z]))%p
                    if (D*az*az-h*fs[z]*bz*bz)%p==0:roots.append((z,az,bz))
                checked+=1
                if len(roots)!=d:continue
                Ps=[]
                for z,az,bz in roots:
                    assert bz!=0
                    y=-s*az*pow(bz,-1,p)%p
                    x=u*z%p
                    assert (y*y-x**3-a*x-c)%p==0
                    Ps.append((x,y))
                assert sum_ec(Ps,p,a) is None
                for i,P in enumerate(Ps):
                    assert sum_ec(Ps[:i]+Ps[i+1:],p,a)==(P[0],-P[1]%p)
                    target_checks+=1
                got.add(tuple(sorted(Ps)));norm_good+=1
            assert got==direct
            direct_good+=len(direct)
            print(f'p={p}, m={m}: {len(direct)} distinct reduced divisors, exact support agreement',flush=True)
    print(f'Enumerated {checked:,} normalized sections; {norm_good:,} norm hits; {direct_good:,} direct divisors; {target_checks:,} target sign checks.',flush=True)


def pair_short(p,gamma,n,d,r):
    # Lattice { (U,V): phi(U)-r phi(V)=0 }. Basis rows.
    coeff=[pow(gamma,i,p) for i in range(d)]+[(-r*pow(gamma,i,p))%p for i in range(d)]
    L=sp.zeros(2*d,2*d);L[0,0]=p
    for i in range(1,2*d):L[i,0]=-coeff[i];L[i,i]=1
    v=lll_exact(L.tolist())[0]
    U=v[:d]+[0]*(n-d);V=v[d:]+[0]*(n-d)
    assert ev(U,gamma,p)==r*ev(V,gamma,p)%p and ev(V,gamma,p)!=0
    # Finite LLL fixture check: a loose bound on the asymptotic SVP radius.
    cap=4*d*2**((p.bit_length()+2*d-1)//(2*d))
    assert max(map(abs,U+V))<=cap
    return U,V


def short_data(p,a,c,gamma,lam):
    n=5;d=n-1 if lam==1 else n
    u,f=model(p,a,c)
    # At most two distinct nonzero coefficients, after scaling.
    vals=list(dict.fromkeys(v for v in f if v));assert 1<=len(vals)<=2
    c0=vals[0];c1=vals[-1]
    U,V=pair_short(p,gamma,n,d,c1*pow(c0,-1,p)%p)
    h=ev(V,gamma,p)*pow(c0,-1,p)%p
    assert h and ev(U,gamma,p)==h*c1%p
    weights=[([0]*n if v==0 else V if v==c0 else U) for v in f]
    if pow(h,(p-1)//2,p)==1:
        D=[1]+[0]*(n-1)
    else:
        nonsq=next(x for x in range(2,10000) if pow(x,(p-1)//2,p)==p-1)
        J,K=pair_short(p,gamma,n,d,nonsq)
        D=next(x for x in (J,K) if pow(ev(x,gamma,p),(p-1)//2,p)==p-1)
    s=int(sp.sqrt_mod(ev(D,gamma,p)*pow(h,-1,p)%p,p))
    assert s and h*s*s%p==ev(D,gamma,p)
    assert all(ev(W,gamma,p)==h*v%p for W,v in zip(weights,f))
    return u,f,h,s,D,weights


def solve_mod(rows,rhs,p):
    R=[list(v)+[r%p] for v,r in zip(rows,rhs)]; k=len(rows[0]);i=0
    for j in range(k):
        a=next(a for a in range(i,len(R)) if R[a][j]%p);R[i],R[a]=R[a],R[i]
        t=pow(R[i][j]%p,-1,p);R[i]=[x*t%p for x in R[i]]
        for a in range(len(R)):
            if a==i:continue
            t=R[a][j]%p
            R[a]=[(x-t*y)%p for x,y in zip(R[a],R[i])]
        i+=1
    assert all(not any(v) for v in R[k:])
    return [R[i][-1] for i in range(k)]


def mulpoly(A,B,mul,zero,plus):
    C=[zero() for _ in range(len(A)+len(B)-1)]
    for i,a in enumerate(A):
        for j,b in enumerate(B):C[i+j]=plus(C[i+j],mul(a,b))
    return C


def lift_checks():
    rng=Random(20261124)
    for bits,cyc in [(96,False),(128,True)]:
        p=int(sp.nextprime(2**(bits-1)+314159))
        while ((p-1)%5==0)!=cyc:p=int(sp.nextprime(p))
        gamma,lam,M,det,lift,quotient=toy_construct(p)
        n=len(M); d=n-1 if cyc else n; c=det//p;e=c*pow(c,-1,p)
        unit=[1]+[0]*(n-1);zero=lambda:[0]*n
        mat=sp.Matrix.hstack(*[sp.Matrix(fold(M,[int(i==j) for i in range(n)],lam)) for j in range(n)])
        J=hermite_normal_form(mat.row_join(e*sp.eye(n)));E=lll_exact(J.T.tolist())[0]
        assert ev(E,gamma,p)!=0 and quotient(scale(p,E)) is not None
        # Check thin support injectivity on a small exact cube.
        d0=min(d,3); cube=[list(v)+[0]*(n-d0) for v in product(range(8),repeat=d0)]
        assert len({ev(v,gamma,p) for v in cube})==len(cube)
        print(f'{bits}-bit field, lambda={lam}: {len(cube)} injective thin addresses.',flush=True)
        # Draw actual random curve points, with all x distinct, and complete a target.
        ok=bad=eqs=norm_checks=target_lifts=0;wmax=0;ratio_checks=0
        for a,c0 in [(0,3),(1,0),(7,11),(rng.randrange(1,p),rng.randrange(1,p))]:
            assert (4*a**3+27*c0*c0)%p
            u,f,h,s,D,weights=short_data(p,a,c0,gamma,lam);ui=pow(u,-1,p);ratio_checks+=1
            wmax=max(wmax,max(abs(x).bit_length() for P in [D]+weights for x in P))
            def point():
                while True:
                    x=rng.randrange(p);g=(x**3+a*x+c0)%p
                    if pow(g,(p-1)//2,p) in (0,1):return x,int(sp.sqrt_mod(g,p))
            for m in (2,3,4,5,6):
                for trial in range(2):
                    while True:
                        Ps=[point() for _ in range(m)];T=sum_ec(Ps,p,a)
                        if T is not None and len({P[0] for P in Ps+[T]})==m+1:break
                    roots=[(P[0]*ui%p,P[1]) for P in Ps]+[(T[0]*ui%p,-T[1]%p)]
                    names,lead=coefficients(m)
                    def elt(t,j,x,y):return (s if t==0 else y)*pow(x,j,p)%p
                    cs=solve_mod([[elt(t,j,x,y) for t,j in names] for x,y in roots],[-elt(*lead,x,y) for x,y in roots],p)
                    A,B=lists_from(m,cs)
                    assert all((s*peval(A,x,p)+y*peval(B,x,p))%p==0 for x,y in roots)
                    Ar=list(map(lift,A));Br=list(map(lift,B));(Ar,Br)[lead[0]][lead[1]]=unit
                    mult=lambda x,y:fold(x,y,lam)
                    N1=[mult(D,v) for v in mulpoly(Ar,Ar,mult,zero,add)]
                    N2=mulpoly(weights,mulpoly(Br,Br,mult,zero,add),mult,zero,add)
                    N=[add(N1[i] if i<len(N1) else zero(),scale(-1,N2[i] if i<len(N2) else zero())) for i in range(m+2)]
                    # A separate linear lift fixes the actual target sign.
                    xT,yneg=roots[-1]
                    Ftar=zero()
                    for j,v in enumerate(Ar): Ftar=add(Ftar,mult(lift(s*pow(xT,j,p)%p),v))
                    for j,v in enumerate(Br): Ftar=add(Ftar,mult(lift(yneg*pow(xT,j,p)%p),v))
                    assert quotient(mult(E,Ftar)) is not None
                    target_lifts+=1
                    Lc=ev(N[-1],gamma,p); assert Lc!=0
                    us=[];aa=[];tt=[]
                    for x,y in roots[:-1]:
                        z=[rng.randrange(8) if j<d else 0 for j in range(n)];t=rng.randrange(1,p)
                        us.append(z);aa.append((x-t*ev(z,gamma,p))%p);tt.append(t)
                    for delta in (0,1):
                        av=aa.copy();av[0]=(av[0]+delta)%p
                        # Public coefficients compressed before multiplying private digit polynomials.
                        R=[zero() for _ in range(m+2)]
                        for mask in range(1<<m):
                            coeff=[-roots[-1][0]%p,1]; Uprod=unit
                            for j in range(m):
                                if mask>>j&1:
                                    coeff=[-tt[j]*v%p for v in coeff];Uprod=mult(Uprod,us[j])
                                else:coeff=mulpoly(coeff,[-av[j]%p,1],lambda x,y:x*y%p,lambda:0,lambda x,y:(x+y)%p)
                            for j,v in enumerate(coeff):
                                if v%p==0:continue
                                C=lift(Lc*v%p)
                                if j==m+1:C=N[-1]
                                R[j]=add(R[j],mult(C,Uprod))
                        exact=0
                        for v,w in zip(N,R):
                            F=add(v,scale(-1,w));EF=mult(E,F);Q=quotient(EF)
                            assert (Q is not None)==(ev(F,gamma,p)==0)
                            if Q is not None:
                                assert mult(M,Q)==EF;exact+=1;eqs+=1
                        assert (exact==m+2)==(delta==0)
                        if delta:bad+=1
                        else:ok+=1
                    norm_checks+=1
        print(f'{bits}-bit field: {ratio_checks} short multiplier systems; max public weight bits {wmax}; {ok} exact lifted norm systems, {bad} mutated systems rejected; {eqs} exact projected quotient identities; {target_lifts} separate target sign lifts.',flush=True)

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    finite_sections();lift_checks()
