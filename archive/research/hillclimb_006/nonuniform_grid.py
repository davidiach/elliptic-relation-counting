"""Exact dyadic grid, repair, mode-spectrum and readout-cost audit.

Independent standard-library implementation. No asymptotic contractor is
assumed or executed. Assertions are mandatory.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import prod
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parent


def v2(n):
    if not n:
        raise ValueError('valuation of zero is not finite')
    n=abs(n)
    return (n&-n).bit_length()-1


def qv(x):
    return v2(x.numerator)-v2(x.denominator)


def coordinates(n):
    return (0,)+tuple(1<<j for j in range(n-1))


def differences(X):
    return tuple(sorted({x-y for x in X for y in X}))


def derivative_values(grid):
    return {t:sum(v2(t-s) for s in grid if s!=t) for t in grid}


def repair(grid, values):
    # Smallest shift exponent in this paired construction, subject to
    # exceeding every old pair-difference valuation.
    separation=max(v2(x-y) for x,y in combinations(grid,2))+1
    required=max(2*values[0]-2*values[d]+v2(d) for d in grid if d)
    E=max(separation,required)
    added=tuple(d+(1<<E) for d in grid if d)
    out=tuple(sorted(set(grid)|set(added)))
    assert len(out)==2*len(grid)-1
    return out,E


def multiply(a,b,q):
    out=[Q(0)]*min(q,len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:q-i]):
            out[i+j]+=x*y
    return out


def w_jet(grid,delta,tau,ell,q):
    cp=prod(-s for s in grid if s)
    polynomial=[Q(1,cp)]
    for s in grid:
        if s!=delta:
            polynomial=multiply(polynomial,[Q(tau-s),Q(1<<ell)],q)
    return polynomial+[Q(0)]*max(0,q-len(polynomial))


def all_one_packet_jet(X,grid,tau,ell,q):
    """A full packet with every old Y and X pinned to one; F(0)=1."""
    W={d:w_jet(grid,d,tau,ell,q) for d in differences(X)}
    polynomial=[Q(1)]
    for row in X:
        row_form=[sum((W[col-row][j] for col in X),Q(0)) for j in range(q)]
        polynomial=multiply(polynomial,row_form,q)
    return polynomial+[Q(0)]*max(0,q-len(polynomial))


def gaussian_modes(B,tau):
    """Independent exact elimination on a scalar Cauchy matrix."""
    size=len(B)
    A=[[Q(1,tau+x-y) for y in B] for x in B]
    R=[row[:] for row in A]
    L=[[Q(int(i==j)) for j in range(size)] for i in range(size)]
    U=[[Q(int(i==j)) for j in range(size)] for i in range(size)]
    diagonal=[]
    for k in range(size):
        pivot=R[k][k]
        assert pivot
        diagonal.append(pivot)
        for i in range(k+1,size):
            L[i][k]=R[i][k]/pivot
        for j in range(k+1,size):
            U[k][j]=R[k][j]/pivot
        for i in range(k+1,size):
            for j in range(k+1,size):
                R[i][j]-=L[i][k]*pivot*U[k][j]
    for i in range(size):
        for j in range(size):
            assert sum((L[i][k]*diagonal[k]*U[k][j] for k in range(size)),Q(0))==A[i][j]
    assert all(qv(L[i][k])==(i-k if k else 0) for k in range(size) for i in range(k,size))
    assert all(qv(U[k][j])==(j-k if k else 0) for k in range(size) for j in range(k,size))
    return [qv(d) for d in diagonal]


def value(grid,delta,T):
    return prod((T-s for s in grid if s!=delta),start=Q(1))/prod(-s for s in grid if s)


def main():
    if not __debug__:
        raise RuntimeError('Do not disable assertions')
    start=time.monotonic()
    records=[]
    derivative_checks=mode_checks=0
    for n in range(3,19):
        X=coordinates(n)
        D=differences(X)
        V=derivative_values(D)
        m=n-2
        assert len(D)==m*m+m+3
        assert V[0]==(m**3+5*m)//3
        assert V[1]==(n*n-n-4)//2
        assert V[1]-V[0]<0
        assert min(V[t] for t in D if t%2)==3*n-8
        assert min(2*V[t]-v2(t) for t in D if t)==2*(3*n-8)
        G,E=repair(D,V)
        U=derivative_values(G)
        assert E==max(n,2*V[0]-6*n+16)
        assert U[0]==2*V[0]
        for d in D:
            if d:
                predicted=E+2*V[d]-v2(d)
                assert U[d]==U[d+(1<<E)]==predicted
        assert min(U.values())==U[0]
        assert max(v2(x-y) for x,y in combinations(G,2))==E
        ell=E
        delta={t:ell+U[t]-U[0] for t in G if t%2}
        assert min(delta.values())==(E if n>=4 else 4)
        # First, the stronger arbitrary-polynomial Hermite witness.
        whole_binary_budget=n*(n-1)*(n-2)//3
        if n>=7:
            assert min(delta.values())>=whole_binary_budget
            # Stronger relevance: a realizable all-one packet deformation.
            # At most two rows have an old difference equal to the odd node.
            # All other row forms are coefficientwise divisible by 2^(delta-n).
            for d in delta.values():
                precision_lower=Q((n-2)*(d-n),n) # per original edge
                max_weight=d+(n-2)*(n-1)
                assert precision_lower>Q(max_weight,4)
        if n<=10:
            B=tuple(x for x in X if x%2==0)
            for tau in (1,1+(1<<E)):
                actual=gaussian_modes(B,tau)
                assert actual==[k*(k+1) for k in range(n-1)]
                mode_checks+=1
        derivative_checks+=len(D)+len(G)
        records.append(dict(n=n,original_nodes=len(D),derivative_at_zero=V[0],
                            derivative_at_one=V[1],unrepaired_constant_valuation=V[1]-V[0],
                            repaired_nodes=len(G),minimal_pair_shift=E,local_scale=ell,
                            odd_base_weight_min=min(delta.values()),odd_base_weight_max=max(delta.values()),
                            quadratic_mode_increments=[k*(k+1) for k in range(n-1)],
                            exponential_candidate_bound_applies=n>=7))

    # Direct products independently confirm the derivative-valuation sums.
    direct_derivatives=0
    jets=[]
    for n in (3,4,5,7):
        X=coordinates(n);D=differences(X);V=derivative_values(D)
        G,E=repair(D,V);U=derivative_values(G);ell=E
        if n>=4:
            # A minimally separated pair at d=3 has no constant margin.
            # Scaling by one fewer power of two makes its linear term
            # nonintegral, so the chosen common local scale is necessary.
            assert U[3]==U[0]
            deficient=w_jet(G,3,3,E-1,2)
            assert qv(deficient[1])==-1
        for grid,vals in ((D,V),(G,U)):
            for tau in grid:
                assert v2(prod(tau-s for s in grid if s!=tau))==vals[tau]
                direct_derivatives+=1
        # A negative constant term is unchanged by increasing jet scale.
        for extra in (0,5,11):
            raw=w_jet(D,1,1,ell+extra,3)
            assert qv(raw[0])==V[1]-V[0]
        selected=sorted({1,-1,1+(1<<E),-1+(1<<E)}|({3,3+(1<<E)} if 3 in D else set()))
        for tau in selected:
            for d in D:
                local=w_jet(G,d,tau,ell,3)
                assert all(not x or qv(x)>=0 for x in local)
            q=n+2
            F=all_one_packet_jet(X,G,tau,ell,q)
            base=ell+U[tau]-U[0]
            bound=(n-2)*(base-n)
            assert all(not x or qv(x)>=bound for x in F)
            first=next((j for j,x in enumerate(F) if x),None)
            assert first is not None
            jets.append(dict(n=n,node=str(tau),jet_length=q,proved_divisibility=bound,
                             first_nonzero_degree=first,observed_minimum_valuation=min(qv(x) for x in F if x)))
        assert all(value(G,d,Q(0))==int(d==0) for d in D)
        # All Y and X pins are one: W(0)=I, so the original count is one.
        assert prod(sum(value(G,col-row,Q(0)) for col in X) for row in X)==1

    affine_checks=0
    for n in (3,5,8):
        X=coordinates(n);D=differences(X)
        for shift in (1,3,8):
            scale=1<<shift
            expanded=tuple(scale*d for d in D)
            for T in (Q(1,3),Q(7,5)):
                for d in (0,1,max(D)):
                    assert value(expanded,scale*d,scale*T)==value(D,d,T)
                    affine_checks+=1
    clustered=[]
    for u in (1,2,3):
        B=1<<u
        for n in (3,4,5,6,8):
            X=tuple(B*x+r for x in coordinates(n) for r in range(B))
            D=differences(X);V=derivative_values(D)
            retained=[t for t in D if t and v2(t)<u]
            bad=min(retained,key=lambda t:V[t]-V[0])
            clustered.append(dict(u=u,clusters=n,ports=len(X),difference_nodes=len(D),
                                  central_constant_valuation=V[1]-V[0],
                                  worst_retained_node=bad,worst_retained_constant_valuation=V[bad]-V[0]))
    assert all(r['worst_retained_constant_valuation']<0 for r in clustered if r['clusters']==8)
    result=dict(status='passed',families=records,derivative_valuation_checks=derivative_checks,
                independent_direct_derivative_products=direct_derivatives,
                exact_Cauchy_elimination_checks=mode_checks,all_one_packet_jets=jets,
                uniform_dilation_identities=affine_checks,
                balanced_cluster_pilot=clustered,
                conclusions=[
                    'The natural geometric difference grid has negative local constant valuations of cubic magnitude.',
                    'A minimal separated paired-node repair restores local integrality with polynomially many nodes and polynomial bit lengths.',
                    'For n>=7, odd-node Hermite reconstruction plus unrestricted normalized-mode enumeration still has a 2^(M/4) feasible-candidate lower bound at some required node.',
                    'Uniform dilation of all coordinates and interpolation nodes is exactly invariant.'
                ],
                scope='A specific geometric coordinate family, a specific minimal paired-node repair, and reconstruction using only odd-node jet residues with unrestricted candidate enumeration. This is not a lower bound for all nonuniform grids, support-aware methods, or ECDLP algorithms.',
                seconds=round(time.monotonic()-start,3))
    (ROOT/'grid_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',families=len(records),derivative_checks=derivative_checks,
                         Cauchy_eliminations=mode_checks,packet_jet_checks=len(jets),
                         dilation_checks=affine_checks,seconds=result['seconds']),indent=2))


if __name__=='__main__':
    main()
