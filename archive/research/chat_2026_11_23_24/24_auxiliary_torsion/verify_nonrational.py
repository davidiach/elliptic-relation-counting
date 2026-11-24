"""Negative control: auxiliary division polynomial alone admits nonrational roots."""
from verify_chart import *
p,ell=13,5
a,c,N,H,S,tr=aux(p,ell);K=K2(p)
zs=[z for z in product(range(p),range(1,p)) if division(ell,z,a,c,K)[0]==K.zero]
assert zs
z=zs[0];z2=K.pw(z,p);s=K.add(z,z2);v=K.mul(z,z2);assert s[1]==v[1]==0
f=[v[0],-s[0]%p,1];alpha=next(x for x in product(range(p),range(1,p)) if x!=z and x!=z2)
assert all(pe(f,x,p)!=0 for x in range(p))
K,Ms=chart_matrices(p,ell,a,c,[f],2,alpha,(0,1,a,c))
assert det(K,Ms[0][0])==K.zero
bad=0
for eta in product(range(p),repeat=2):
    e2=K.mul(eta,eta)
    M=[[K.add(Ms[0][0][i][j],K.add(K.mul(eta,Ms[1][0][i][j]),K.mul(e2,Ms[2][0][i][j]))) for j in range(2)] for i in range(2)]
    bad+=det(K,M)==K.zero
assert bad<=4
print(f'p={p}, ell={ell}, curve coefficients=({a},{c}), order={N}, rational support={S}')
print(f'Nonrational roots: {len(zs)}. Division polynomial alone gives a false hit on irreducible quadratic {f}.')
print(f'Filtered determinant: {bad} false hits out of {p*p} eta values; rigorous bound 4. All assertions passed.')
