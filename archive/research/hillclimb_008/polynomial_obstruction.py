"""Explicit dual polynomials for degree-aware linear readout.

These are linear-operator lower bounds, not lower bounds for ECDLP or
for nonlinear readout on a structured class of deformation polynomials.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sys

from optimal_readout import v2

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'hillclimb_007'))
from cluster_audit import coarse, A as central_A


def consecutive(a,u,mass):
    assert 1<=u<=a
    L=1<<a
    nodes={t:u-v2(t) for t in range(-L,L) if t and v2(t)<u}
    degree=mass*(2*L-1)
    J=2*L*(u-1)+(2*L>>u)
    assert J==sum(nodes.values())
    h=degree//J
    losses={}
    for tau in nodes:
        i=v2(tau)
        H=Q(L,1<<i)*sum((Q(j,1<<j)*min(j,u-i) for j in range(1,a-i+1)),Q(0))
        assert H.denominator==1
        losses[tau]=h*(int(H)+(a+1-i)*(u-i))
    return dict(kind='consecutive',parameters=dict(a=a,u=u,mass=mass),degree=degree,
                nodes=nodes,ell=a+1,J=J,h=h,losses=losses,
                lower_scale=min((1<<v2(t))*losses[t] for t in nodes),
                simple_lower=Q(h*L,2))


def cluster(n,u,v,mass,central=False):
    assert n>=4 and 1<=v<=u
    B=1<<u
    D,E,_,VE,c,e=coarse(n)
    grid=range(1,B) if central else sorted({B*t for t in D}|{B*t+r for t in E for r in range(1,B)})
    nodes={t:v-v2(t) for t in grid if t and v2(t)<v}
    degree=mass*(B*e-c-1)
    J=(B*(v-1)+(B>>v))*(1 if central else e)
    assert J==sum(nodes.values())
    h=degree//J
    ell=u+n-1
    losses={}
    for tau in nodes:
        i=v2(tau);z=v-i;low=u-i
        H=central_A(low,z)
        if not central:
            t,_=divmod(tau,B)
            H=e*H+z*((e-1)*low+VE[t])
        losses[tau]=h*(H+(ell-i)*z)
    return dict(kind='central_cluster' if central else 'full_cluster',
                parameters=dict(n=n,u=u,v=v,mass=mass),degree=degree,
                nodes=nodes,ell=ell,J=J,h=h,losses=losses,
                lower_scale=min((1<<v2(t))*losses[t] for t in nodes),
                simple_lower=Q(h*B,2) if central else Q(h*e*B,4))


def polynomial(orders):
    out=[1]
    for tau,multiplicity in sorted(orders.items()):
        for _ in range(multiplicity):
            new=[0]*(len(out)+1)
            for k,c in enumerate(out):
                new[k]-=tau*c
                new[k+1]+=c
            out=new
    return out


def verify_obstruction(case,expanded=False):
    """Check the node-difference proof, optionally every polynomial jet."""
    nodes=case['nodes'];ell=case['ell'];h=case['h']
    assert 0 not in nodes and h>=1
    orders={tau:h*q for tau,q in nodes.items()}
    assert sum(orders.values())<=case['degree']
    assert all(v2(t-s)<=ell for t in nodes for s in nodes if t!=s)
    root_v=sum(q*v2(t) for t,q in orders.items())
    for tau,a in orders.items():
        direct=ell*a+sum(q*v2(tau-s) for s,q in orders.items() if s!=tau)-root_v
        assert direct==case['losses'][tau]
    assert case['lower_scale']>=case['simple_lower']
    checks=0
    if expanded:
        f=polynomial(orders)
        assert len(f)-1==sum(orders.values()) and v2(f[0])==root_v
        for tau,a in orders.items():
            jets=[sum(comb(k,j)*f[k]*tau**(k-j) for k in range(j,len(f)))<<(ell*j)
                  for j in range(len(f))]
            assert all(x==0 for x in jets[:a]) and jets[a]
            target=root_v+case['losses'][tau]
            assert v2(jets[a])==target
            assert all(not x or v2(x)>=target for x in jets)
            checks+=len(jets)
    return dict(kind=case['kind'],parameters=case['parameters'],degree=case['degree'],
                witness_degree=sum(orders.values()),node_count=len(nodes),
                base_jet_count=case['J'],h=h,lower_scale=case['lower_scale'],
                simple_lower=str(case['simple_lower']),expanded_jet_checks=checks)


def large_parameter_checks():
    """Closed-form inequalities only; never construct exponential grids."""
    records=[]
    # Both balanced-cluster policies, with physical mass at least one packet.
    for u in (32,48,64,128,256):
        B=1<<u
        for n in (4,u,3*u+3):
            c=(n-3)*(n-4)+1;e=2*(n*n-5*n+9)
            mass=n*B;degree=mass*(B*e-c-1)
            for v in (1,max(1,u//2),u):
                base=B*(v-1)+(B>>v)
                for central in (False,True):
                    J=base*(1 if central else e)
                    h=degree//J
                    lower=Q(h*B,2) if central else Q(h*e*B,4)
                    target=Q(mass*e*B,8*v if central else 16*v)
                    assert lower>=target
                    assert e>=Q(n*n,2) and c+1<=e
                    records.append(dict(u=u,n=n,v=v,central=central,
                                        lower=str(lower),target=str(target)))
    # Uniform H_i estimates used in the proof, independent of node count.
    for low in range(1,129):
        for z in range(1,low+1):
            A=central_A(low,z)
            assert A>= (1<<low)-low-1
            assert A+Q(low,2)*z>=Q(1<<low,4)
            assert A+low*z>=Q(1<<low,2)
    return records
