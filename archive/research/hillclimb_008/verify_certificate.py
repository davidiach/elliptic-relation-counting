"""Direct certificate checker; does not use the module diagonalization.

Rows are rebuilt by multiplying linear polynomials, independently of the
optimizer's binomial matrix construction. Both exact upper identities
and integer-polynomial lower witnesses are checked.
"""
from fractions import Fraction as Q


def valuation(x):
    x=Q(x)
    assert x
    a,b=abs(x.numerator),x.denominator
    return (a&-a).bit_length()-(b&-b).bit_length()


def shifted_monomials(tau,ell,degree):
    out=[[1]]
    for _ in range(degree):
        previous=out[-1]
        following=[0]*(len(previous)+1)
        for j,x in enumerate(previous):
            following[j]+=tau*x
            following[j+1]+=(1<<ell)*x
        out.append(following)
    return out


def matrix(certificate):
    nodes={int(t):q for t,q in certificate['nodes'].items()}
    scales={int(t):q for t,q in certificate['scales'].items()}
    d=certificate['degree'];rows=[];specs=[]
    for tau,q in sorted(nodes.items()):
        columns=shifted_monomials(tau,scales[tau],d)
        for j in range(q):
            rows.append([col[j] if j<len(col) else 0 for col in columns])
            specs.append((tau,j))
    assert [list(x) for x in specs]==[list(x) for x in certificate['specs']]
    return rows,specs


def identity(rows,weights):
    assert len(weights)==len(rows)
    for k in range(len(rows[0])):
        assert sum((Q(w)*row[k] for row,w in zip(rows,weights)),Q(0))==int(k==0)


def common(certificate,require_upper=True):
    A,_=matrix(certificate)
    powers=certificate['row_powers']
    assert len(powers)==len(A) and all(isinstance(p,int) and p>=0 for p in powers)
    B=[[x<<p for x in row] for row,p in zip(A,powers)]
    L=certificate['optimal_denominator_exponent']
    assert isinstance(L,int) and L>=0
    if require_upper:
        w=[Q(x) for x in certificate['weights']]
        identity(B,w)
        assert all(not x or valuation(x)>=-L for x in w)
    if L:
        lower=certificate['lower_certificate']
        f=lower['polynomial']
        assert len(f)==certificate['degree']+1 and all(isinstance(x,int) for x in f) and f[0]
        jets=[sum(x*y for x,y in zip(row,f)) for row in B]
        assert any(jets)
        t=valuation(f[0]);s=min(valuation(x) for x in jets if x)
        assert s-t==L
        assert lower['constant_valuation']==t and lower['minimum_jet_valuation']==s
    return True


def profile(certificate):
    A,specs=matrix(certificate)
    budgets={int(t):b for t,b in certificate['budgets'].items()}
    assert all(isinstance(x,int) and x>=0 for x in budgets.values())
    if certificate['feasible']:
        w=[Q(x) for x in certificate['original_weights']]
        identity(A,w)
        assert all(not x or valuation(x)>=-budgets[t] for x,(t,j) in zip(w,specs))
    else:
        f=certificate['lower_certificate']['polynomial']
        assert len(f)==certificate['degree']+1 and f[0]
        jets=[sum(x*y for x,y in zip(row,f)) for row in A]
        assert min(valuation(x)-budgets[t] for x,(t,j) in zip(jets,specs) if x)>valuation(f[0])
    return True


def geometric(result):
    scale=result['optimal_scale']
    upper=result['upper_certificate']
    assert upper['feasible']
    expected={int(t):scale//(1<<valuation(int(t))) for t in upper['nodes']}
    assert expected=={int(t):b for t,b in upper['budgets'].items()}
    profile(upper)
    lower=result['lower_certificate']
    if scale:
        assert lower and not lower['feasible']
        assert lower['degree']==upper['degree'] and lower['nodes']==upper['nodes'] and lower['scales']==upper['scales']
        expected={int(t):(scale-1)//(1<<valuation(int(t))) for t in lower['nodes']}
        assert expected=={int(t):b for t,b in lower['budgets'].items()}
        profile(lower)
    else:
        assert lower is None
    return True
