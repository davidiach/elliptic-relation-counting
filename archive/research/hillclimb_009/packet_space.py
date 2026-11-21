"""Literal complete binary-packet deformation polynomials.

All polynomial families here are generated from legal indicator tensors.
No readout or lower-bound witness is assumed to belong to their span.
"""
from fractions import Fraction as Q
from itertools import product
from math import gcd,prod


def trim(poly):
    while len(poly)>1 and not poly[-1]:poly.pop()
    return poly


def add(a,b):
    out=list(a)+[0]*max(0,len(b)-len(a))
    for i,x in enumerate(b):out[i]+=x
    return trim(out)


def multiply(a,b):
    if not any(a) or not any(b):return [0]
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:out[i+j]+=x*y
    return trim(out)


def evaluate(poly,x):
    out=0
    for c in reversed(poly):out=out*x+c
    return out


def kernels(a):
    L=1<<a;grid=range(-L,L);cp=prod(-t for t in grid if t)
    W={}
    for delta in range(-L+1,L):
        polynomial=[1]
        for t in grid:
            if t!=delta:polynomial=multiply(polynomial,[-t,1])
        assert polynomial[0]==(cp if delta==0 else 0)
        W[delta]=polynomial
    return W,cp


def family(a,X,W=None,cp=None):
    """Return all Y-indicator numerators with common denominator cp^L."""
    L=1<<a
    assert len(X)==L and all(x in (0,1) for x in X)
    if W is None:W,cp=kernels(a)
    states={():[1]}
    for x in range(L):
        one=[0]
        for y,value in enumerate(X):
            if value:one=add(one,W[y-x])
        zero=add([cp],[-c for c in one])
        states={bits+(bit,):multiply(poly,factor)
                for bits,poly in states.items() for bit,factor in ((0,zero),(1,one))}
    denominator=cp**L
    assert all(f[0]==denominator*int(tuple(X)==Y) for Y,f in states.items())
    total=[0]
    for f in states.values():total=add(total,f)
    assert total==[denominator]
    return states


def modular_basis(rows,width,prime=1000000007,stop_at_full=False):
    """A rank lower bound over Q, certified by the recorded prime pivots."""
    basis={};selected=[]
    for label,poly in rows:
        row=[x%prime for x in poly]+[0]*(width-len(poly))
        assert len(row)==width
        for p,b in sorted(basis.items()):
            value=row[p]
            if value:row=[(x-value*y)%prime for x,y in zip(row,b)]
        pivot=next((i for i,x in enumerate(row) if x),None)
        if pivot is not None:
            inverse=pow(row[pivot],-1,prime)
            basis[pivot]=[(x*inverse)%prime for x in row]
            selected.append(label)
            if stop_at_full and len(basis)==width:break
    return dict(rank=len(basis),pivots=sorted(basis),selected=selected,prime=prime)


def primitive(row):
    divisor=0
    for x in row:divisor=gcd(divisor,x)
    if divisor:row=[x//divisor for x in row]
    first=next((x for x in row if x),0)
    return [-x for x in row] if first<0 else row


def rational_basis(rows,width):
    """Exact primitive integer elimination, with all rows actually checked."""
    basis={};selected=[];count=0
    for label,poly in rows:
        row=list(poly)+[0]*(width-len(poly));assert len(row)==width
        for p,b in sorted(basis.items()):
            if row[p]:
                left,right=b[p],row[p]
                row=primitive([left*x-right*y for x,y in zip(row,b)])
        pivot=next((i for i,x in enumerate(row) if x),None)
        if pivot is not None:
            basis[pivot]=primitive(row);selected.append(label)
        count+=1
    return dict(rank=len(basis),pivots=sorted(basis),selected=selected,
                echelon=[basis[p] for p in sorted(basis)],rows_checked=count)


def complete_rows(a):
    W,cp=kernels(a)
    for X in product((0,1),repeat=1<<a):
        for Y,f in family(a,X,W,cp).items():yield (X,Y),f


def copy_polynomial(a,X,pairs,W=None,cp=None):
    """Sum the actual equality-predicate selectors, retaining every factor."""
    L=1<<a
    assert sorted(x for pair in pairs for x in pair)==list(range(L))
    if W is None:W,cp=kernels(a)
    forms=[]
    for x in range(L):
        form=[0]
        for y,bit in enumerate(X):
            if bit:form=add(form,W[y-x])
        forms.append(form)
    out=[1]
    for x,y in pairs:
        factor=add([cp*cp],[-cp*c for c in add(forms[x],forms[y])])
        factor=add(factor,[2*c for c in multiply(forms[x],forms[y])])
        out=multiply(out,factor)
    assert out[0]==cp**L*int(all(X[x]==X[y] for x,y in pairs))
    return out
