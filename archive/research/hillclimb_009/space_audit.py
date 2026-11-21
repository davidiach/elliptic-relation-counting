"""Reconstruct legal packet families and certify rank and parity readout."""
from copy import deepcopy
from fractions import Fraction as Q
from itertools import product
from math import comb,isqrt
from pathlib import Path
import json
import random
import sys
import time

ROOT=Path(__file__).resolve().parent
for previous in ('hillclimb_003','hillclimb_007','hillclimb_008'):
    sys.path.insert(0,str(ROOT.parent/previous))
from c38_slices import copy_restriction,sc,common
from cluster_audit import hermite_weights,qv
from optimal_readout import optimize_geometric_profile
import verify_certificate as exact_verify
from packet_space import kernels,family,complete_rows,copy_polynomial,multiply,add,rational_basis,modular_basis
from parity_readout import jet_data,optimize,verify,v2


def independent_kernel(a):
    L=1<<a;C=[1]
    for t in range(-L,L):
        new=[0]*(len(C)+1)
        for k,x in enumerate(C):new[k]-=t*x;new[k+1]+=x
        C=new
    cp=C[1];out={}
    for delta in range(-L+1,L):
        quotient=[0]*(len(C)-1);quotient[-1]=C[-1]
        for k in range(len(quotient)-1,0,-1):quotient[k-1]=C[k]+delta*quotient[k]
        assert C[0]+delta*quotient[0]==0
        out[delta]=quotient
    return out,cp


def indicator(a,X,Y,W,cp):
    out=[1]
    for x,bit in enumerate(Y):
        form=[sum(W[y-x][k]*X[y] for y in range(len(X))) for k in range(2*(1<<a))]
        if not bit:form=[cp-form[0]]+[-z for z in form[1:]]
        out=multiply(out,form)
    return out


def determinant_mod(matrix,prime):
    n=len(matrix);assert all(len(row)==n for row in matrix)
    A=[[x%prime for x in row] for row in matrix];answer=1
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:A[j],A[pivot]=A[pivot],A[j];answer=-answer
        value=A[j][j];answer=answer*value%prime;inverse=pow(value,-1,prime)
        for i in range(j+1,n):
            factor=A[i][j]*inverse%prime
            for k in range(j+1,n):A[i][k]=(A[i][k]-factor*A[j][k])%prime
            A[i][j]=0
    return answer%prime


def copy_fixture(a):
    L=1<<a;net,_,pairs,accepted=copy_restriction(L//2)
    labels={t.uid:(0,)*a for t in net.tensors}
    packets=sc.geom.packetize(net,labels,a,a)
    assert len(packets)==1 and set(packets[0]['positions'])==set(range(L))
    position={e:x for x,e in packets[0]['positions'].items()}
    paired=[(position[x],position[y]) for x,y in pairs]
    return net,packets,position,paired,accepted


def native_checks(a,fixture):
    L=1<<a;W,cp=kernels(a);net,packets,position,pairs,accepted=fixture
    q=5;ring=sc.Ring(96,(('check',q),));axis=ring.variable('check')
    centers=[-1,1]+([2] if a>=2 else [])+([4] if a>=3 else [])
    checks=kernel_checks=0
    sample_X=[(0,)*L,(1,)*L,tuple(int(j%2==0) for j in range(L)),tuple(int(j==0) for j in range(L))]
    for tau in centers:
        local,_=sc.geom.local_polynomials(a,tau,q)
        for delta,poly in W.items():
            shifted=tuple(Q(sum(comb(k,j)*poly[k]*tau**(k-j) for k in range(j,len(poly)))<<(a+1)*j,cp)
                          for j in range(q))
            assert shifted==local[delta]+(Q(0),)*(q-len(local[delta]));kernel_checks+=q
        values={delta:ring.polynomial(coefficients,axis) for delta,coefficients in local.items()}
        for X in sample_X:
            forms={}
            for x in range(L):
                value=ring.zero
                for y,bit in enumerate(X):
                    if bit:value=ring.add(value,values[y-x])
                forms[x]=value
            actual=ring.zero
            for Y in accepted:
                assert net.rows[0].entry(tuple(Y[e] for e in net.rows[0].ports),ring)==ring.one
                term=ring.one
                for e,x in position.items():term=ring.mul(term,net.alphabet.selector(Y[e],(forms[x],),ring))
                actual=ring.add(actual,term)
            f=copy_polynomial(a,X,pairs,W,cp)
            rows,_,_,odd=jet_data([f],cp**L,{tau:q},a+1)
            expected=tuple(row[0]*pow(odd,-1,ring.N)%ring.N for row in rows)
            assert actual==expected;checks+=q
    # Full original-assignment sums on small generic indicator networks.
    full_checks=0
    if a<=2:
        ports=tuple(range(L));alphabet=sc.Alphabet(1,1)
        prototype=sc.network('packet-prototype',[(ports,lambda s:1)],[(ports,lambda s:1)],alphabet,p=96)
        packed=sc.geom.packetize(prototype,{t.uid:(0,)*a for t in prototype.tensors},a,a)
        pos={e:x for x,e in packed[0]['positions'].items()}
        for X,Y in [(sample_X[0],sample_X[0]),(sample_X[2],sample_X[2]),(sample_X[3],sample_X[1])]:
            network=sc.network('literal-indicators',
                [(ports,lambda s,Y=Y:int(all(s[e]==Y[pos[e]] for e in ports)))],
                [(ports,lambda s,X=X:int(all(s[e]==X[pos[e]] for e in ports)))],alphabet,p=96)
            for tau in centers:
                local,_=sc.geom.local_polynomials(a,tau,q)
                values={delta:ring.polynomial(c,axis) for delta,c in local.items()}
                actual=sc.direct_deformation(network,packed,values,ring)
                f=indicator(a,X,Y,W,cp)
                rows,_,_,odd=jet_data([f],cp**L,{tau:q},a+1)
                assert actual==tuple(row[0]*pow(odd,-1,ring.N)%ring.N for row in rows)
                full_checks+=q
    return dict(kernel_coefficients=kernel_checks,actual_C38_selector_coefficients=checks,
                complete_indicator_deformation_coefficients=full_checks)


def observation_trials(rows,target,specs,odd,result,seed):
    rng=random.Random(seed);scale=result['optimal_scale'];upper=result['upper_certificate']
    budgets={int(t):b for t,b in upper['budgets'].items()};B=max(budgets.values())
    weights=upper['weight_numerators'];checks=0
    for trial in range(20):
        coefficients=[rng.randrange(-1000,1001) for _ in target]
        answer=sum(x*y for x,y in zip(coefficients,target))%2
        observed=[]
        for row,(tau,j) in zip(rows,specs):
            modulus=1<<(budgets[tau]+1)
            # Actual jet data includes the global odd denominator.
            residue=sum(x*y for x,y in zip(row,coefficients))*pow(odd,-1,modulus)%modulus
            residue+=rng.randrange(-9,10)*modulus
            observed.append(residue)
        numerator=sum(w*odd*(1<<(B-budgets[t]))*x for w,x,(t,j) in zip(weights,observed,specs))
        assert numerator%(1<<(B+1))==(answer<<B)
        checks+=1
    return checks


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    started=time.monotonic();common.verify_sources()
    prime=1000000007
    assert prime%2 and all(prime%d for d in range(3,isqrt(prime)+1,2))
    labels=json.loads((ROOT/'basis_labels.json').read_text())
    ranks=[];native=[];copies={};families={};negative=0
    for a in (1,2,3):
        L=1<<a;width=L*(2*L-1)+1;W,cp=kernels(a);independent,other_cp=independent_kernel(a)
        assert independent==W and cp==other_cp
        generators=labels[str(L)];rows=[]
        for X,Y in generators:
            f=indicator(a,X,Y,independent,cp)
            assert f[0]==cp**L*int(X==Y)
            rows.append(f+[0]*(width-len(f)))
        assert len(rows)==width
        det=determinant_mod(rows,prime);assert det
        ranks.append(dict(family='generic',L=L,dimension=width,rank=width,
                          prime=prime,determinant_mod_prime=det,generator_labels=generators))
        if a<=2:
            complete=list(complete_rows(a));exact=rational_basis(complete,width)
            assert exact['rank']==width
            for (X,Y),f in complete:assert f==indicator(a,X,Y,independent,cp)
            families[a,'generic']=complete
        fixture=copy_fixture(a);copies[a]=fixture
        pairs=fixture[3]
        entries=[(X,copy_polynomial(a,X,pairs,W,cp)) for X in product((0,1),repeat=L)]
        families[a,'copy']=entries
        rank=modular_basis(entries,width,prime)
        selected=[dict(entries)[tuple(X)] for X in rank['selected']]
        columns=rank['pivots'];minor=[[f[j] if j<len(f) else 0 for j in columns] for f in selected]
        det_copy=determinant_mod(minor,prime);assert det_copy
        assert rank['rank']==min(len(entries),width)
        ranks.append(dict(family='actual_C38_copy_predicate',L=L,dimension=width,
                          rank=rank['rank'],generators=len(entries),prime=prime,
                          determinant_mod_prime=det_copy,selected=rank['selected'],
                          coefficient_columns=columns,pairs=pairs))
        if a<=2:
            for X,f in entries:
                literal=family(a,X,W,cp);total=[0]
                for Y in fixture[4]:
                    key=tuple(Y[fixture[1][0]['positions'][x]] for x in range(L))
                    total=add(total,literal[key])
                assert f==total
        native.append(dict(L=L,**native_checks(a,fixture)))
        print(json.dumps(dict(L=L,generic_rank=width,C38_copy_rank=rank['rank'])),flush=True)
    damaged=[row[:] for row in rows];damaged[0]=[0]*len(damaged[0])
    assert determinant_mod(damaged,prime)==0;negative+=1

    previous=json.loads((ROOT.parent/'hillclimb_008/readout_results.json').read_text())
    readouts=[];certificates={};observations=0
    for a,kind in ((1,'generic'),(1,'copy'),(2,'generic'),(2,'copy'),(3,'copy')):
        start=time.monotonic();L=1<<a;degree=L*(2*L-1);u=a
        h=(L+u-2)//(u-1) if u>1 else (degree+2)//2
        nodes={t:h*(u-v2(t)) for t in range(-L,L) if t and v2(t)<u}
        full=hermite_weights(nodes,a+1)
        bound=max([0]+[-(1<<v2(t))*qv(w) for t,ws in full.items() for w in ws if w])
        if a==1:
            exact=optimize_geometric_profile(nodes,degree,{t:a+1 for t in nodes},bound)
        else:
            exact=previous['certificates'][f'consecutive_{a}_{a}_{L}']['geometric']
            cert=exact['upper_certificate']
            assert {int(t):q for t,q in cert['nodes'].items()}==nodes and cert['degree']==degree
            assert set(cert['scales'].values())=={a+1}
        exact_verify.geometric(exact)
        entries=families[a,kind];_,cp=kernels(a)
        data,target,specs,odd=jet_data([f for X,f in entries],cp**L,nodes,a+1)
        result=optimize(data,target,specs,bound)
        verify(data,target,specs,json.loads(json.dumps(result)))
        observations+=observation_trials(data,target,specs,odd,result,900+a)
        key=f'{kind}_L{L}'
        certificates[key]=dict(parity=result,full_space_exact=exact,
                               odd_denominator=odd,specs=specs,
                               generator_labels=[X for X,f in entries])
        row=dict(family=kind,L=L,degree=degree,generators=len(entries),jets=len(data),
                 standard_full_Hermite_scale=bound,optimal_full_space_exact_scale=exact['optimal_scale'],
                 optimal_integer_family_parity_scale=result['optimal_scale'],seconds=round(time.monotonic()-start,3))
        readouts.append(row);print(json.dumps(row),flush=True)
        for mutation in ('weights','lower'):
            bad=deepcopy(result)
            if mutation=='weights':bad['upper_certificate']['weight_numerators']=[0]*len(data)
            else:bad['lower_certificate']['indistinguishable_odd_combination']=[0]*len(target)
            try:verify(data,target,specs,bad)
            except AssertionError:negative+=1
            else:raise AssertionError('A corrupted parity certificate was accepted')

    # Same rational span, different coefficient lattice: the solver must
    # distinguish genuinely modular recovery from exact rational readout.
    controls=[]
    for multiplier,expected in ((1,1),(2,0)):
        data,target,specs,odd=jet_data([[1],[0,multiplier]],1,{1:2},1)
        result=optimize(data,target,specs,2)
        verify(data,target,specs,result);assert result['optimal_scale']==expected
        controls.append(dict(generators=[[1],[0,multiplier]],parity_scale=expected,certificate=result))
    output=dict(status='passed',seconds=round(time.monotonic()-started,3),rank_certificates=ranks,
                native_checks=native,readouts=readouts,certificates=certificates,
                lattice_controls=controls,random_observation_trials=observations,
                negative_controls=negative,
                scope='Finite complete packets and actual C38 predicate restrictions, plus the proved fixed-packet product extension; no full C38 contraction or ECDLP runtime improvement')
    (ROOT/'space_results.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',seconds=output['seconds'],readout_certificates=len(readouts),
                         observation_trials=observations,negative_controls=negative)),flush=True)


if __name__=='__main__':main()
