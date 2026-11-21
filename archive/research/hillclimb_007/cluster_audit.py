"""Exact cluster-grid formulas, Hermite readout losses and candidate budgets."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, isqrt, prod
from pathlib import Path
import json
import sys
import time

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'hillclimb_006'))
from nonuniform_grid import coordinates, differences, derivative_values, v2, qv, multiply, w_jet


@lru_cache(None)
def coarse(n):
    D=differences(coordinates(n))
    E=tuple(sorted(set(D)|{d-1 for d in D}))
    J=derivative_values(E)
    V=derivative_values(D)
    previous=differences(coordinates(n-1))
    assert set(E)=={2*x+r for x in previous for r in (-1,0)}
    c=(n-3)*(n-4)+1
    e=2*(n*n-5*n+9)
    assert len(E)==e and len(E)-len(D)==c
    assert J[0]==V[0]
    assert min(J.values())==n*n-2*n-3
    for x in E:
        j=(x+1)//2
        assert J[x]==len(previous)-1+sum(v2(j-y) for y in previous if y!=j)
    return D,E,V,J,c,e


def admissible_depth(n,u):
    m=n-2
    c=(n-3)*(n-4)+1
    A=(m**3+5*m)//3-(n*n-2*n-3)
    gap=max(1,(A+c-1)//c)
    return max(0,u-gap+1)


def A(h,z):
    return 3*(1<<h)-(1<<(h-z))*(z+3)-z*(h+1)


def central_loss(n,u,v):
    m=n-2
    origin=(m**3+5*m)//3
    e=2*(n*n-5*n+9)
    return e*A(u,v)+v*((e-1)*u+origin)


def nodes_for(G,v,h0=1):
    return {x:h0*(v-v2(x)) for x in G if x and v2(x)<v}


def highest_loss(nodes,tau,ell):
    value=sum(q*(v2(tau-s)-v2(s)) for s,q in nodes.items() if s!=tau)
    return value+(nodes[tau]-1)*(ell-v2(tau))


def hermite_weights(nodes,ell):
    """Full Hermite evaluation functional at zero, in scaled jet coordinates."""
    alpha=1<<ell
    out={}
    for tau,q in nodes.items():
        E0=Q(1)
        inverse=[Q(1)]
        for s,mult in nodes.items():
            if s==tau:
                continue
            E0*=Q(-s,tau-s)**mult
            ratio=Q(alpha,tau-s)
            factor=[(-1)**k*comb(mult+k-1,k)*ratio**k for k in range(q)]
            inverse=multiply(inverse,factor,q)
        inverse += [Q(0)]*(q-len(inverse))
        w=Q(-tau,alpha)
        out[tau]=[E0*w**j*sum((inverse[k]*w**k for k in range(q-j)),Q(0)) for j in range(q)]
    return out


def gaussian(A):
    n=len(A)
    R=[row[:] for row in A]
    L=[[Q(int(i==j)) for j in range(n)] for i in range(n)]
    U=[[Q(int(i==j)) for j in range(n)] for i in range(n)]
    diagonal=[]
    for k in range(n):
        pivot=R[k][k]
        assert pivot
        diagonal.append(pivot)
        for i in range(k+1,n):L[i][k]=R[i][k]/pivot
        for j in range(k+1,n):U[k][j]=R[k][j]/pivot
        for i in range(k+1,n):
            for j in range(k+1,n):R[i][j]-=L[i][k]*pivot*U[k][j]
    for i in range(n):
        for j in range(n):
            assert sum((L[i][k]*diagonal[k]*U[k][j] for k in range(n)),Q(0))==A[i][j]
    assert all(qv(L[i][k])>=0 for k in range(n) for i in range(k,n))
    assert all(qv(U[k][j])>=0 for k in range(n) for j in range(k,n))
    return [qv(d) for d in diagonal]


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    start=time.monotonic()
    cases=[]
    valuation_checks=mode_checks=hermite_loss_checks=0
    gaussian_checks=0
    for n,u in product((4,5,6,8,10,12),(1,2,3,4,5)):
        B=1<<u
        D,E,V,J,c,e=coarse(n)
        X=tuple(B*x+r for x in coordinates(n) for r in range(B))
        G=differences(X)
        expected={B*t for t in D}|{B*t+r for t in E for r in range(1,B)}
        assert set(G)==expected and len(G)==B*e-c
        origin=u*(len(D)-1)+V[0]+e*(B-1-u)
        assert origin==sum(v2(-s) for s in G if s)
        selected=G if len(G)<=600 else tuple(sorted({0,1,-1,G[0],G[-1]}|
                    {G[j*len(G)//23] for j in range(23)}))
        for tau in selected:
            t,r=divmod(tau,B)
            predicted=(V[t]-V[0]) if r==0 else c*(u-v2(r))+J[t]-V[0]
            actual=sum(v2(tau-s) for s in G if s!=tau)-origin
            assert actual==predicted
            valuation_checks+=1
        v=admissible_depth(n,u)
        if v:
            worst=c*(u-v+1)+min(J.values())-V[0]
            assert worst>=0 and n<=3*u+3
            if v<u:
                assert c*(u-v)+min(J.values())-V[0]<0
            assert all(sum(x%(1<<v)==r for x in X)==len(X)//(1<<v) for r in range(1<<v))
            profile=nodes_for(G,v)
            for tau in (1,1<<max(0,v-1)):
                i=v2(tau);h=u-i;z=v-i
                direct=sum(q*(v2(tau-s)-v2(s)) for s,q in profile.items() if s!=tau)
                predicted=e*A(h,z)+z*((e-1)*h+J[tau//B])
                assert direct==predicted
                hermite_loss_checks+=1
            degree=len(X)*(len(G)-1)
            denominator=e*(B*(v-1)+(B>>v))
            h0=(degree+1+denominator-1)//denominator
            assert h0*denominator>degree
            assert h0*v>=Q(len(X),2)
            L0=h0*central_loss(n,u,v)
        else:
            h0=L0=None
        block=tuple(x for x in X if x%2==0)
        prefix=[0]*len(block)
        spectrum=[]
        for j,x in enumerate(block):
            # Independent greedy-minimum certificate for integral normalized
            # triangular factors, not merely the diagonal formula.
            assert prefix[j]==min(prefix[j:])
            a,t=divmod(j,B//2)
            predicted=4*j-2*t.bit_count()+a*(a-3)
            assert 2*prefix[j]==predicted
            spectrum.append(predicted)
            mode_checks+=1
            for k in range(j+1,len(block)):prefix[k]+=v2(block[k]-x)
        if len(block)<=12:
            matrix=[[Q(1,1+x-y) for y in block] for x in block]
            assert gaussian(matrix)==spectrum
            gaussian_checks+=1
        cases.append(dict(n=n,u=u,ports=len(X),grid_nodes=len(G),coarse_runs=c,
                          admissible_retained_levels=v,mode_count_one_orientation=len(block),
                          maximum_diagonal_increment=max(spectrum),
                          minimum_triangular_height=h0,central_E0_loss=L0))

    # Direct negative control: central values alone miss bad coarse nodes.
    n,u=8,2
    B=1<<u;D,E,V,J,c,e=coarse(n)
    X=tuple(B*x+r for x in coordinates(n) for r in range(B));G=differences(X)
    t=min(E,key=lambda x:J[x]);tau=B*t+B//2
    bad_node=tau
    raw=w_jet(G,tau,tau,u+n-1,2)
    assert qv(raw[0])==-16

    # Exact interpolation identities on every monomial, including the
    # highest coefficient whose readout loss has no cancellation.
    interpolation=[]
    for n,u,v,policy in ((4,1,1,'all_coarse'),(4,2,2,'all_coarse'),
                         (4,2,2,'central'),(6,3,2,'central')):
        B=1<<u;D,E,V,J,c,e=coarse(n)
        X=tuple(B*x+r for x in coordinates(n) for r in range(B))
        G=differences(X)
        domain=G if policy=='all_coarse' else range(1,B)
        nodes=nodes_for(domain,v)
        ell=u+n-1
        weights=hermite_weights(nodes,ell)
        N=sum(nodes.values())
        for degree in range(N):
            got=sum((w*comb(degree,j)*tau**(degree-j)*(1<<ell)**j
                     for tau,ws in weights.items() for j,w in enumerate(ws) if j<=degree),Q(0))
            assert got==int(degree==0)
        for tau,ws in weights.items():
            assert -qv(ws[-1])==highest_loss(nodes,tau,ell)
        interpolation.append(dict(n=n,u=u,v=v,policy=policy,monomials=N,
                                  central_highest_coefficient_loss=-qv(weights[1][-1])))

    # Large-parameter certificate uses exact integers only; no exponentially
    # large packet, matrix, or candidate table is materialized.
    asymptotic=[]
    for u in (32,40,64,96,128):
        B=1<<u
        assert B>=144*u*(u+1)**3+16*u
        for n in (4,8,16,3*u):
            c=(n-3)*(n-4)+1
            e=2*(n*n-5*n+9)
            vmax=admissible_depth(n,u)
            assert vmax>=1
            for v in sorted({1,vmax}):
                N=n*B;M=16*N
                denominator=e*(B*(v-1)+(B>>v))
                h0=(M*(B*e-c-1)+1+denominator-1)//denominator
                q=h0*v
                L0=h0*central_loss(n,u,v)
                assert L0>=Q(M*n*n*B,8*v)
                k=min(n*B//16,isqrt(n**3*B*B//(256*v)))
                K=16*k
                assert c*u+n*n<=k
                assert K<q and K<=M//16
                structural=16*(k*(c*u+n*n)+2*k*(k-1))
                assert structural<L0
                central_denominator=B*(v-1)+(B>>v)
                central_h0=(M*(B*e-c-1)+1+central_denominator-1)//central_denominator
                central_L0=central_h0*A(u,v)
                assert structural<central_L0 and K<central_h0*v
                if n>=v:assert K>=M//32
                else:assert (32*K)**2*v>=M*M*n
                asymptotic.append(dict(u=u,n=n,v=v,M=str(M),
                                       log2_feasible_binary_candidates_at_least=str(K),
                                       central_only_policy_also_certified=True,
                                       method='closed-form integer certificate; no packet allocation'))

    central=[]
    for u in (1,2,3,4,8,16,32):
        B=1<<u
        for n in (max(4,u),max(8,4*u),max(16,B+1)):
            # Closed forms avoid constructing huge geometric coordinate sets.
            c=(n-3)*(n-4)+1;e=2*(n*n-5*n+9)
            for v in sorted({1,u}):
                assert n>=v
                M=128*n*B
                denominator=B*(v-1)+(B>>v)
                h0=(M*(B*e-c-1)+1+denominator-1)//denominator
                q=h0*v
                L0=h0*A(u,v)
                max_structural=c*u+2*n*B+n*n
                if u>=2:
                    assert A(u,v)>=B//4
                    assert Q(M,128)*max_structural<L0
                else:
                    # Only node 1: q-1 is the full deformation degree.
                    # ell=0 fails on the leading W coefficient; for ell>=1,
                    # the ell-dependent readout budget already pays for all.
                    assert q-1==M*(B*e-c-1)
                    assert Q(M,128)*(1+max_structural)<q-1
                assert M//128<q
                central.append(dict(u=u,n=n,v=v,log2_candidate_bound=str(M//128)))

    result=dict(status='passed',finite_cases=cases,exact_derivative_checks=valuation_checks,
                exact_mode_weights=mode_checks,independent_Gaussian_checks=gaussian_checks,
                exact_Hermite_E0_checks=hermite_loss_checks,interpolation_checks=interpolation,
                inadmissible_coarse_node_negative_control=dict(n=8,u=2,node=bad_node,constant_valuation=-16),
                all_coarse_asymptotic_certificates=asymptotic,
                central_large_cluster_certificates=central,
                scope='Specific balanced geometric-cluster grids with triangular multiplicities and the standard full Hermite readout, followed by unrestricted feasible-mode enumeration. No lower bound for arbitrary interpolation functionals, coefficient-aware precision, support pruning, or ECDLP.',
                seconds=round(time.monotonic()-start,3))
    (ROOT/'cluster_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',finite_cases=len(cases),derivative_checks=valuation_checks,
                         mode_weights=mode_checks,Gaussian_checks=gaussian_checks,
                         interpolation_monomials=sum(x['monomials'] for x in interpolation),
                         asymptotic_certificates=len(asymptotic),central_certificates=len(central),
                         seconds=result['seconds']),indent=2))


if __name__=='__main__':
    main()
