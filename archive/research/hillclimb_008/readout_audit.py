"""Exact operator optima, independent certificates and asymptotic witnesses."""
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import json
import random
import sys
import time

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'hillclimb_007'))
from cluster_audit import hermite_weights, coarse
from optimal_readout import (optimize,optimize_geometric_profile,jet_matrix,
                             row_module,v2,qv)
from polynomial_obstruction import (consecutive,cluster,verify_obstruction,
                                    large_parameter_checks)
import verify_certificate as independent


def bareiss(A):
    A=[list(row) for row in A];n=len(A);sign=1;old=1
    for k in range(n-1):
        pivot=next(i for i in range(k,n) if A[i][k])
        if pivot!=k:A[k],A[pivot]=A[pivot],A[k];sign=-sign
        current=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*current-A[i][k]*A[k][j]
                assert numerator%old==0
                A[i][j]=numerator//old
            A[i][k]=0
        old=current
    return sign*A[-1][-1]


def full_cost(nodes,ell):
    weights=hermite_weights(nodes,ell)
    common=max([0]+[-qv(w) for ws in weights.values() for w in ws if w])
    weighted=max([0]+[-(1<<v2(t))*qv(w) for t,ws in weights.items() for w in ws if w])
    return common,weighted,weights


def inverse_check(matrix,bound):
    n=len(matrix)
    A=[[Q(x) for x in row]+[Q(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if A[i][k])
        A[k],A[pivot]=A[pivot],A[k]
        scale=A[k][k];A[k]=[x/scale for x in A[k]]
        for i in range(n):
            if i!=k:
                scale=A[i][k]
                A[i]=[x-scale*y for x,y in zip(A[i],A[k])]
    assert all(A[i][j]==int(i==j) for i in range(n) for j in range(n))
    assert all(not x or qv(x)>=-bound for row in A for x in row[n:])
    return n*n


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    start=time.monotonic();certificates={};results=[]
    counts=dict(common_certificates=0,geometric_certificates=0,minor_determinants=0,
                inverse_entry_checks=0,module_index_checks=0,unique_hermite_controls=0,negative_controls=0)

    def common_case(name,nodes,degree,scales,full=None):
        cert=optimize(nodes,degree,scales)
        independent.common(json.loads(json.dumps(cert)))
        A,_=jet_matrix(nodes,degree,scales)
        if degree<=25:
            minor=[A[i] for i in cert['minor_rows']]
            assert v2(bareiss(minor))==cert['determinant_valuation_bound']
            counts['minor_determinants']+=1
            if degree<=11:
                counts['inverse_entry_checks']+=inverse_check(minor,cert['inverse_entry_bound'])
        augmented=A+[[int(k==0) for k in range(degree+1)]]
        index=sum(row_module(augmented,cert['working_precision'],track_rows=False)['valuations'])
        assert sum(cert['invariant_valuations'])-index==cert['optimal_denominator_exponent']
        counts['module_index_checks']+=1;counts['common_certificates']+=1
        certificates[name]=dict(common=cert)
        row=dict(name=name,degree=degree,jets=sum(nodes.values()),
                 optimal_common=cert['optimal_denominator_exponent'])
        if full is not None:
            row['full_common'],row['full_geometric'],weights=full
            geometric=optimize_geometric_profile(nodes,degree,scales,row['full_geometric'])
            independent.geometric(json.loads(json.dumps(geometric)))
            counts['geometric_certificates']+=1
            row['optimal_geometric']=geometric['optimal_scale']
            certificates[name]['geometric']=geometric
            assert row['optimal_common']<=row['full_common']
            assert row['optimal_geometric']<=row['full_geometric']
            if degree==sum(nodes.values())-1:
                expected=[w for t in sorted(nodes) for w in weights[t]]
                assert [Q(w) for w in cert['weights']]==expected
                assert row['optimal_common']==row['full_common']
                assert row['optimal_geometric']==row['full_geometric']
                counts['unique_hermite_controls']+=1
        results.append(row)
        print(json.dumps(row),flush=True)
        return cert

    # Unique interpolation, overdetermined interpolation, a zero center,
    # and unequal local scales are separate controls.
    for name,nodes,d,ell in [
            ('constant',{-1:1,1:1},0,2),('two_values',{-1:1,1:1},1,2),
            ('single_odd',{1:5},4,3),('single_even',{2:5},4,3),
            ('overdetermined',{-1:2,1:2},2,3),('unique',{-1:2,1:2},3,3)]:
        common_case(name,nodes,d,{t:ell for t in nodes},full_cost(nodes,ell))
    common_case('zero_node',{-2:2,0:2,3:1},4,{-2:2,0:4,3:1})
    common_case('unequal_scales',{-3:2,2:3,5:2},5,{-3:1,2:4,5:2})
    rng=random.Random(80631)
    for i in range(24):
        centers=sorted(rng.sample([x for x in range(-8,9) if x],rng.randrange(1,5)))
        nodes={t:rng.randrange(1,4) for t in centers}
        d=rng.randrange(sum(nodes.values()))
        scales={t:rng.randrange(6) for t in centers}
        common_case('random_'+str(i),nodes,d,scales)

    # M is the deformation mass. Only mass % packet_size == 0 has a
    # complete packet interpretation; all rows remain polynomial tests.
    for a,u,mass in [(2,2,1),(2,2,2),(2,2,4),(3,2,1),(3,2,2),
                     (3,3,1),(3,3,2),(3,3,4),(3,3,8)]:
        L=1<<a;h=(mass+u-2)//(u-1)
        nodes={t:h*(u-v2(t)) for t in range(-L,L) if t and v2(t)<u}
        common_case(f'consecutive_{a}_{u}_{mass}',nodes,mass*(2*L-1),
                    {t:a+1 for t in nodes},full_cost(nodes,a+1))
        results[-1]['mass_divisible_by_packet_size']=mass%L==0

    # Small balanced-cluster degree fixtures, not full network executions.
    for n,u,v,mass,central in [(4,1,1,1,False),(4,2,2,1,False),
                              (4,2,2,1,True),(5,1,1,1,False)]:
        c=cluster(n,u,v,mass,central)
        h=(c['degree']+c['J'])//c['J']
        nodes={t:h*q for t,q in c['nodes'].items()}
        common_case(f'cluster_{n}_{u}_{v}_{mass}_{central}',nodes,c['degree'],
                    {t:c['ell'] for t in nodes},full_cost(nodes,c['ell']))

    # Corrupt three distinct certificate assertions; all must be rejected.
    base=certificates['overdetermined']['common']
    for mutation in ('weight','denominator','witness'):
        bad=deepcopy(base)
        if mutation=='weight':bad['weights'][0]=str(Q(bad['weights'][0])+1)
        elif mutation=='denominator':bad['optimal_denominator_exponent']-=1
        else:bad['lower_certificate']['polynomial']=[1]+[0]*bad['degree']
        try:independent.common(bad)
        except AssertionError:counts['negative_controls']+=1
        else:raise AssertionError('A corrupted certificate was accepted: '+mutation)

    obstructions=[verify_obstruction(c,expanded=True) for c in [
        consecutive(2,2,4),consecutive(3,3,4),cluster(4,1,1,1),
        cluster(4,2,2,2),cluster(5,2,2,2),cluster(4,2,2,1,True)]]
    large=large_parameter_checks()
    result=dict(status='passed',seconds=round(time.monotonic()-start,3),counts=counts,
                cases=results,polynomial_obstructions=obstructions,
                large_formula_certificates=large,certificates=certificates,
                scope='Exact linear readout on arbitrary rational polynomials; no complete network contraction or ECDLP runtime improvement claimed')
    (ROOT/'readout_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],seconds=result['seconds'],counts=counts,
                         expanded_jet_checks=sum(x['expanded_jet_checks'] for x in obstructions),
                         large_formula_certificates=len(large))),flush=True)


if __name__=='__main__':main()
