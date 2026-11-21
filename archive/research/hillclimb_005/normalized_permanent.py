"""Test a Cauchy permanent shortcut on the actual normalized mode factors."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'hillclimb_001'))
import common
import simplex_contraction as sc


def residual(X,k,x,y,T):
    value=F(1,T+x-y)
    for earlier in X[:k]:
        value *= F((x-earlier)*(earlier-y),(T+x-earlier)*(T+earlier-y))
    return value


def rational_detper(A):
    determinant=permanent=F(0)
    for p in permutations(range(len(A))):
        value=F(1)
        for i,j in enumerate(p):
            value *= A[i][j]
        permanent += value
        determinant += (-1)**sum(p[i]>p[j] for i in range(len(A)) for j in range(i+1,len(A)))*value
    return determinant,permanent


def ring_detper(A,ring):
    determinant=permanent=ring.zero
    for p in permutations(range(len(A))):
        value=ring.one
        for i,j in enumerate(p):
            value=ring.mul(value,A[i][j])
        permanent=ring.add(permanent,value)
        sign=(-1)**sum(p[i]>p[j] for i in range(len(A)) for j in range(i+1,len(A)))
        determinant=ring.add(determinant,ring.scale(value,sign))
    return determinant,permanent


def main():
    common.verify_sources()
    X=list(range(0,16,2)); rows=(6,10,14)
    A=[[residual(X,k,x,X[k],1)/residual(X,k,X[k],X[k],1) for k in range(3)] for x in rows]
    d,p=rational_detper(A)
    squared=[[v*v for v in row] for row in A]
    defect=d*p-rational_detper(squared)[0]
    assert defect == F(-8,495990495)
    inverse_entry_det=rational_detper([[1/v for v in row] for row in A])[0]
    assert inverse_entry_det == F(-22,135)
    C=[[F(1,1+x-X[k]) for k in range(3)] for x in rows]
    cd,cp=rational_detper(C)
    assert cd*cp == rational_detper([[v*v for v in row] for row in C])[0]

    # Fully occupied L=16 packet. The selected three rows and first three
    # actual normalized modes yield a nonzero identity defect modulo 2^64.
    ring=sc.Ring(64,(('normalization',7),))
    packet={'key':((0,1),(0,0)),'positions':{x:x for x in range(16)}}
    W,modes,forms=sc.modes_and_forms([packet],4,1,ring,'normalization',1)
    ids=[j for j,m in enumerate(modes) if m['target']==0][:3]
    assert [modes[j]['weight'] for j in ids] == [5,7,11]
    matrix=[[forms[0,x]['modes'][j] for j in ids] for x in rows]
    det,per=ring_detper(matrix,ring)
    square_det=ring_detper([[ring.mul(v,v) for v in row] for row in matrix],ring)[0]
    finite_defect=ring.sub(ring.mul(det,per),square_det)
    assert finite_defect != ring.zero
    assert all(c==0 for c in finite_defect[:6]) and sc.geom.v2(finite_defect[6])==49

    beta=[0]*len(modes)
    for j in ids: beta[j]=1
    specs=[('Y',x) for x in range(0,16,2)]+[('X',x) for x in range(1,16,2)]
    symbols=tuple(int(side=='Y' and x in rows) for side,x in specs)
    tensor=sc.Packet('permanent-probe',specs,specs,0,0,beta,modes,forms,ring,Counter(),sc.Alphabet(1,1))
    actual=tensor.entry(symbols,ring)
    assert actual==per
    dense_rows=[]
    for a in (3,4,5,6,7):
        L=1<<a
        local=sc.Ring(2*L+1,(('dense',2),))
        full={'key':((0,1),(0,0)),'positions':{x:x for x in range(L)}}
        _,all_modes,all_forms=sc.modes_and_forms([full],a,1,local,'dense',1)
        incoming=[j for j,m in enumerate(all_modes) if m['target']==0]
        valuations=[sc.geom.v2(all_forms[0,L-2]['modes'][j][1]) for j in incoming]
        expected=[a+1+4*k-2*k.bit_count() for k in range(L//2)]
        assert valuations==expected and len(incoming)==L//2
        assert max(valuations)==2*L-a-1
        dense_rows.append(dict(L=L, nonzero_modes=len(incoming), precision=local.p,
                               largest_leading_valuation=max(valuations)))
    result=dict(status='passed',normalized_scalar_matrix=[[str(v) for v in row] for row in A],
                scalar_identity_defect=str(defect), reciprocal_matrix_determinant=str(inverse_entry_det),
                raw_Cauchy_positive_control=True, actual_packet_entry_equals_permanent=True,
                actual_mode_weights=[modes[j]['weight'] for j in ids],
                finite_jet_defect=list(finite_defect), protected_valuation=49, precision=64, jet_length=7,
                dense_last_row_checks=dense_rows,
                scope='Rejects applying the raw Cauchy Borchardt identity directly to these normalized mode coefficients, including nonzero diagonal row/column rescalings over Q. Does not exclude different structured-permanent algorithms.',
                primary_reference='https://arxiv.org/abs/2309.05886')
    (ROOT/'normalized_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
