"""Certified optimal dyadic denominators for degree-limited linear readout.

The model is an exact linear identity on arbitrary rational polynomials.
No claim is made about nonlinear readout with additional input promises.
"""
from fractions import Fraction as Q
from math import comb


def v2(x):
    if not x:
        raise ValueError('zero has no finite valuation')
    x=abs(x)
    return (x&-x).bit_length()-1


def qv(x):
    return v2(x.numerator)-v2(x.denominator)


def jet_matrix(nodes,degree,scales):
    """Nodes and scales are dictionaries indexed by distinct integer centers."""
    assert degree>=0 and nodes and all(q>=1 for q in nodes.values())
    assert sum(nodes.values())>=degree+1
    assert set(nodes)==set(scales) and all(ell>=0 for ell in scales.values())
    rows=[];specs=[]
    for tau,q in sorted(nodes.items()):
        for j in range(q):
            rows.append([comb(k,j)*tau**(k-j)*(1<<(scales[tau]*j)) if j<=k else 0
                         for k in range(degree+1)])
            specs.append((tau,j))
    return rows,specs


def minor_certificate(nodes,degree,scales,specs):
    """Greedily select initial jet prefixes; optimality of this choice is unneeded."""
    used={tau:0 for tau in nodes}
    C=0
    for _ in range(degree+1):
        def increment(tau):
            return scales[tau]*used[tau]+sum(q*v2(tau-s) for s,q in used.items() if s!=tau)
        tau=min((t for t in nodes if used[t]<nodes[t]),key=lambda t:(increment(t),t))
        C+=increment(tau)
        used[tau]+=1
    selected=[i for i,(tau,j) in enumerate(specs) if j<used[tau]]
    assert len(selected)==degree+1
    return selected,C,used


def inverse_entry_bound(prefixes,scales,specs,powers,selected):
    """Bound every inverse-minor coefficient using cardinal polynomials.

    This avoids paying the sum of all invariant factors as working
    precision when only their maximum is needed. The proof is in
    readout.tex; no numerical estimate or guessed cutoff is used.
    """
    used={t:q for t,q in prefixes.items() if q}
    separation=max([0]+[v2(t-s) for t in used for s in used if t!=s])
    bounds=[]
    for r in selected:
        tau,j=specs[r]
        cross=sum(q*v2(tau-s) for s,q in used.items() if s!=tau)
        bounds.append(cross+(used[tau]-1-j)*separation+scales[tau]*j+powers[r])
    return max(bounds)


def row_module(A,precision,track_rows=True):
    """Diagonalize modulo 2^precision using invertible row/column operations.

    The caller supplies a bound on the negative valuations of inverse
    entries of a full-rank minor, below precision. Thus all invariant
    factors relevant to membership occur before the modular cutoff.
    """
    modulus=1<<precision
    m,n=len(A),len(A[0])
    assert m>=n
    matrix=[[x%modulus for x in row] for row in A]
    U=[[int(i==j) for j in range(m)] for i in range(m)] if track_rows else None
    V=[[int(i==j) for j in range(n)] for i in range(n)]
    valuations=[]
    for k in range(n):
        candidates=((v2(matrix[i][j]),i,j) for i in range(k,m) for j in range(k,n) if matrix[i][j])
        try:
            valuation,pivot_row,pivot_col=min(candidates)
        except ValueError as error:
            raise AssertionError('Working precision did not expose full column rank') from error
        matrix[k],matrix[pivot_row]=matrix[pivot_row],matrix[k]
        if U is not None:U[k],U[pivot_row]=U[pivot_row],U[k]
        if pivot_col!=k:
            for row in matrix:row[k],row[pivot_col]=row[pivot_col],row[k]
            for row in V:row[k],row[pivot_col]=row[pivot_col],row[k]
        unit=matrix[k][k]>>valuation
        inverse=pow(unit,-1,modulus)
        matrix[k]=[(x*inverse)%modulus for x in matrix[k]]
        if U is not None:U[k]=[(x*inverse)%modulus for x in U[k]]
        assert matrix[k][k]==1<<valuation
        for i in range(k+1,m):
            assert matrix[i][k]%(1<<valuation)==0
            factor=matrix[i][k]>>valuation
            if factor:
                matrix[i]=[(x-factor*y)%modulus for x,y in zip(matrix[i],matrix[k])]
                if U is not None:U[i]=[(x-factor*y)%modulus for x,y in zip(U[i],U[k])]
        for j in range(k+1,n):
            assert matrix[k][j]%(1<<valuation)==0
            factor=matrix[k][j]>>valuation
            if factor:
                for i in range(m):matrix[i][j]=(matrix[i][j]-factor*matrix[i][k])%modulus
                for i in range(n):V[i][j]=(V[i][j]-factor*V[i][k])%modulus
        valuations.append(valuation)
    assert valuations==sorted(valuations)
    assert all(matrix[i][j]==((1<<valuations[i]) if i==j else 0) for i in range(m) for j in range(n))
    target=V[0] # e_0 times the column transform
    losses=[max(0,v-(v2(x) if x else precision)) for v,x in zip(valuations,target)]
    return dict(valuations=valuations,target=target,losses=losses,
                optimum=max(losses),U=U,V=V,modulus=modulus)


def solve_exact(matrix,rhs):
    """Fraction-free elimination, used to lift a modular upper certificate."""
    n=len(matrix)
    assert all(isinstance(x,int) for row in matrix for x in row)
    assert all(isinstance(y,int) for y in rhs)
    A=[list(row)+[y] for row,y in zip(matrix,rhs)]
    previous=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if A[i][k]),None)
        if pivot is None:raise AssertionError('Singular certified minor')
        A[k],A[pivot]=A[pivot],A[k]
        current=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n+1):
                numerator=A[i][j]*current-A[i][k]*A[k][j]
                assert numerator%previous==0
                A[i][j]=numerator//previous
            A[i][k]=0
        previous=current
    answer=[Q(0)]*n
    for i in range(n-1,-1,-1):
        answer[i]=(Q(A[i][-1])-sum(A[i][j]*answer[j] for j in range(i+1,n)))/A[i][i]
    return answer


def verify_upper(A,weights,loss):
    assert len(A)==len(weights)
    for k in range(len(A[0])):
        assert sum((w*row[k] for w,row in zip(weights,A)),Q(0))==int(k==0)
    assert all(not w or qv(w)>=-loss for w in weights)


def lower_certificate(A,module):
    L=module['optimum']
    if L==0:
        return dict(reason='The allowed denominator exponent is nonnegative by definition',loss=0)
    i=module['losses'].index(L)
    polynomial=[row[i] for row in module['V']]
    values=[sum(x*y for x,y in zip(row,polynomial)) for row in A]
    constant_valuation=v2(polynomial[0])
    data_valuation=min(v2(x) for x in values if x)
    assert data_valuation-constant_valuation==L
    return dict(loss=L,polynomial=polynomial,constant_valuation=constant_valuation,
                minimum_jet_valuation=data_valuation,
                argument='Weights of denominator exponent at most L-1 would make the recovered constant divisible by one extra power of two')


def optimize(nodes,degree,scales,extra_row_powers=None,construct=True):
    """Optimize common denominator on optional dyadically scaled jet rows."""
    A,specs=jet_matrix(nodes,degree,scales)
    powers=[0]*len(A) if extra_row_powers is None else list(extra_row_powers)
    assert len(powers)==len(A) and all(x>=0 for x in powers)
    scaled=[[x<<power for x in row] for row,power in zip(A,powers)]
    selected,C,prefixes=minor_certificate(nodes,degree,scales,specs)
    C+=sum(powers[i] for i in selected)
    cardinal_bound=inverse_entry_bound(prefixes,scales,specs,powers,selected)
    inverse_bound=min(C,cardinal_bound)
    precision=inverse_bound+2
    module=row_module(scaled,precision,track_rows=construct)
    assert sum(module['valuations'])<=C
    L=module['optimum']
    result=dict(degree=degree,nodes=nodes,scales=scales,specs=specs,
                row_powers=powers,minor_rows=selected,minor_prefixes=prefixes,
                determinant_valuation_bound=C,working_precision=precision,
                inverse_entry_bound=inverse_bound,cardinal_inverse_bound=cardinal_bound,
                invariant_valuations=module['valuations'],optimal_denominator_exponent=L,
                lower_certificate=lower_certificate(scaled,module))
    if construct:
        modulus=module['modulus']
        x=[0]*len(A)
        for i,(value,v) in enumerate(zip(module['target'],module['valuations'])):
            numerator=(value<<L)%modulus
            assert numerator%(1<<v)==0
            x[i]=numerator>>v
        modular=[sum(x[i]*module['U'][i][j] for i in range(len(A)))%modulus for j in range(len(A))]
        modular=[value if value<=modulus//2 else value-modulus for value in modular]
        residual=[(1<<L)*int(k==0)-sum(c*row[k] for c,row in zip(modular,scaled)) for k in range(degree+1)]
        assert all(x%modulus==0 for x in residual)
        minor_transpose=[[scaled[i][k] for i in selected] for k in range(degree+1)]
        correction=solve_exact(minor_transpose,[x//modulus for x in residual])
        lifted=[Q(x) for x in modular]
        for i,value in zip(selected,correction):
            assert not value or qv(value)>=-inverse_bound
            lifted[i]+=modulus*value
        assert all(not x or qv(x)>=0 for x in lifted)
        weights=[x/(1<<L) for x in lifted]
        verify_upper(scaled,weights,L)
        assert max([0]+[-qv(w) for w in weights if w])==L
        result['weights']=[str(w) for w in weights]
    return result


def profile_feasible(nodes,degree,scales,budgets,construct=False):
    """Check whether some exact weights satisfy v2(w_tau,j)>=-budget_tau."""
    assert set(nodes)==set(budgets) and all(x>=0 for x in budgets.values())
    ceiling=max(budgets.values())
    _,specs=jet_matrix(nodes,degree,scales)
    powers=[ceiling-budgets[tau] for tau,j in specs]
    result=optimize(nodes,degree,scales,powers,construct=construct)
    feasible=result['optimal_denominator_exponent']<=ceiling
    result['budgets']=budgets
    result['feasible']=feasible
    if feasible and construct:
        weights=[Q(w)*(1<<power) for w,power in zip(result['weights'],powers)]
        A,_=jet_matrix(nodes,degree,scales)
        for k in range(degree+1):
            assert sum((w*row[k] for w,row in zip(weights,A)),Q(0))==int(k==0)
        assert all(not w or qv(w)>=-budgets[tau] for w,(tau,j) in zip(weights,specs))
        result['original_weights']=[str(w) for w in weights]
    if not feasible:
        polynomial=result['lower_certificate']['polynomial']
        A,_=jet_matrix(nodes,degree,scales)
        values=[sum(x*y for x,y in zip(row,polynomial)) for row in A]
        root_v=v2(polynomial[0])
        assert min(v2(value)-budgets[tau] for value,(tau,j) in zip(values,specs) if value)>root_v
    return result


def optimize_geometric_profile(nodes,degree,scales,upper):
    """Minimize lambda with node budget floor(lambda/2^v2(tau))."""
    assert 0 not in nodes
    def budget(value):return {tau:value//(1<<v2(tau)) for tau in nodes}
    assert profile_feasible(nodes,degree,scales,budget(upper))['feasible']
    lower=-1;high=upper;steps=1
    while high-lower>1:
        mid=(lower+high)//2
        feasible=profile_feasible(nodes,degree,scales,budget(mid))['feasible']
        steps+=1
        if feasible:high=mid
        else:lower=mid
    upper_certificate=profile_feasible(nodes,degree,scales,budget(high),construct=True)
    lower_certificate_=None if high==0 else profile_feasible(nodes,degree,scales,budget(high-1))
    assert upper_certificate['feasible']
    assert lower_certificate_ is None or not lower_certificate_['feasible']
    return dict(optimal_scale=high,binary_search_calls=steps,
                upper_certificate=upper_certificate,lower_certificate=lower_certificate_)
