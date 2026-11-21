"""Certified parity recovery on an integer-generated polynomial family.

This differs from an exact identity on its rational span. Upper bounds
apply to every integer combination of the listed generators; lower
bounds concern that additive closure and the specified jet observations.
"""
from math import comb


def v2(x):
    assert x
    x=abs(x);return (x&-x).bit_length()-1


def jet_data(polynomials,denominator,nodes,ell):
    """Return odd-unit-normalized integral jets and exact root parities."""
    loss=v2(denominator);odd=denominator//(1<<loss)
    target=[]
    for f in polynomials:
        assert f[0]%denominator==0
        target.append((f[0]//denominator)%2)
    rows=[];specs=[]
    for tau,q in sorted(nodes.items()):
        for j in range(q):
            values=[]
            for f in polynomials:
                value=sum(comb(k,j)*f[k]*tau**(k-j) for k in range(j,len(f)))<<(ell*j)
                assert value%(1<<loss)==0
                values.append(value//(1<<loss))
            rows.append(values);specs.append((tau,j))
    return rows,target,specs,odd


def diagonalize(A,precision,track_rows):
    """Capped Smith elimination, including rectangular/rank-deficient cases."""
    N=1<<precision;m=len(A);n=len(A[0])
    R=[[x%N for x in row] for row in A]
    U=[[int(i==j) for j in range(m)] for i in range(m)] if track_rows else None
    V=[[int(i==j) for j in range(n)] for i in range(n)]
    vals=[]
    for k in range(min(m,n)):
        candidate=min(((v2(R[i][j]),i,j) for i in range(k,m) for j in range(k,n) if R[i][j]),default=None)
        if candidate is None:break
        value,r,c=candidate
        R[k],R[r]=R[r],R[k]
        if U is not None:U[k],U[r]=U[r],U[k]
        if c!=k:
            for row in R:row[k],row[c]=row[c],row[k]
            for row in V:row[k],row[c]=row[c],row[k]
        inverse=pow(R[k][k]>>value,-1,N)
        R[k]=[(x*inverse)%N for x in R[k]]
        if U is not None:U[k]=[(x*inverse)%N for x in U[k]]
        for i in range(k+1,m):
            assert R[i][k]%(1<<value)==0
            factor=R[i][k]>>value
            if factor:
                R[i]=[(x-factor*y)%N for x,y in zip(R[i],R[k])]
                if U is not None:U[i]=[(x-factor*y)%N for x,y in zip(U[i],U[k])]
        for j in range(k+1,n):
            assert R[k][j]%(1<<value)==0
            factor=R[k][j]>>value
            if factor:
                for i in range(m):R[i][j]=(R[i][j]-factor*R[i][k])%N
                for i in range(n):V[i][j]=(V[i][j]-factor*V[i][k])%N
        vals.append(value)
    assert vals==sorted(vals)
    assert all(R[i][j]==((1<<vals[i]) if i==j and i<len(vals) else 0)
               for i in range(m) for j in range(n))
    return vals,U,V


def feasible(rows,target,specs,budgets,construct=False):
    """Each jet has 1+budget_tau bits; recover target modulo two."""
    assert all(b>=0 for b in budgets.values())
    B=max(budgets.values());K=B+1;N=1<<K
    powers=[B-budgets[t] for t,j in specs]
    scaled=[[(x<<power)%N for x in row] for row,power in zip(rows,powers)]
    vals,U,V=diagonalize(scaled,K,construct)
    transformed=[sum(x*row[j] for x,row in zip(target,V))%2 for j in range(len(target))]
    missing=next((j for j in range(len(vals),len(target)) if transformed[j]),None)
    result=dict(feasible=missing is None,budgets=budgets,working_precision=K,
                modular_invariant_valuations=vals)
    if missing is not None:
        combination=[row[missing] for row in V]
        assert sum(c*x for c,x in zip(target,combination))%2==1
        assert all(sum(a*x for a,x in zip(row,combination))%(1<<(budgets[t]+1))==0
                   for row,(t,j) in zip(rows,specs))
        result['indistinguishable_odd_combination']=combination
    elif construct:
        x=[(transformed[i]<<B)>>value for i,value in enumerate(vals)]
        coefficients=[sum(c*row[j] for c,row in zip(x,U))%(1<<(budgets[t]+1))
                      for j,(t,order) in enumerate(specs)]
        # The weights on odd-normalized data are coefficients / 2^budget.
        for k,c in enumerate(target):
            value=sum((coefficient<<(B-budgets[t]))*row[k]
                      for coefficient,row,(t,j) in zip(coefficients,rows,specs))-(c<<B)
            assert value%N==0
        result['weight_numerators']=coefficients
    return result


def optimize(rows,target,specs,upper,geometric=True):
    centers={t for t,j in specs}
    assert 0 not in centers
    def budgets(scale):return {t:scale//(1<<v2(t)) if geometric else scale for t in centers}
    assert feasible(rows,target,specs,budgets(upper))['feasible']
    lo=-1;hi=upper;calls=1
    while hi-lo>1:
        mid=(lo+hi)//2
        if feasible(rows,target,specs,budgets(mid))['feasible']:hi=mid
        else:lo=mid
        calls+=1
    good=feasible(rows,target,specs,budgets(hi),True)
    bad=None if hi==0 else feasible(rows,target,specs,budgets(hi-1))
    assert good['feasible'] and (bad is None or not bad['feasible'])
    return dict(optimal_scale=hi,profile='geometric' if geometric else 'common',
                binary_search_calls=calls,upper_certificate=good,lower_certificate=bad)


def verify(rows,target,specs,result):
    """Direct certificate checks without using the module elimination."""
    scale=result['optimal_scale'];geometric=result['profile']=='geometric'
    def budget(t,amount):return amount//(1<<v2(t)) if geometric else amount
    upper=result['upper_certificate'];coefficients=upper['weight_numerators']
    expected={t:budget(t,scale) for t,j in specs}
    assert expected=={int(t):b for t,b in upper['budgets'].items()}
    assert len(coefficients)==len(rows)
    B=max(expected.values())
    for k,c in enumerate(target):
        residual=sum((w<<(B-expected[t]))*row[k] for w,row,(t,j) in zip(coefficients,rows,specs))-(c<<B)
        assert residual%(1<<(B+1))==0
    if scale:
        lower=result['lower_certificate'];f=lower['indistinguishable_odd_combination']
        assert len(f)==len(target) and all(isinstance(x,int) for x in f)
        assert sum(c*x for c,x in zip(target,f))%2==1
        expected={t:budget(t,scale-1) for t,j in specs}
        assert expected=={int(t):b for t,b in lower['budgets'].items()}
        assert all(sum(a*x for a,x in zip(row,f))%(1<<(expected[t]+1))==0
                   for row,(t,j) in zip(rows,specs))
    else:assert result['lower_certificate'] is None
    return True
