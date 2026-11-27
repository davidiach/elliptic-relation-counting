"""Exact synthetic checks of CRT lifting and modular tensor-train counting.
No deployed systems, implicit contractor, or claimed general fast avoidance.
"""
from itertools import product
from collections import Counter
from math import prod
from random import Random
import numpy as np


def isprime(p):
    return p >= 2 and all(p % j for j in range(2, int(p**.5)+1))


def crt(vals, mods):
    P=prod(mods)
    return sum(x*(P//p)*pow(P//p,-1,p) for x,p in zip(vals,mods))%P


def enc(z, mods):
    return tuple(x % p for x in z for p in mods)


def dec(z, mods):
    s=len(mods)
    return crt(z[:s],mods),crt(z[s:],mods)


def add(x,y,dims):
    return tuple((a+b)%p for a,b,p in zip(x,y,dims))


def conv_dict(A,B,dims):
    C=Counter()
    for x,a in A.items():
        for y,b in B.items():
            C[add(x,y,dims)]+=a*b
    return C


def support(lists,dims):
    S={tuple(0 for p in dims)}
    for A in lists:S={add(x,y,dims) for x in S for y in A}
    return S


def wrap_lift(lists, r, mods):
    assert prod(mods)<=r
    return [sorted(set(enc(z,mods) for z in A)) for A in lists]+[
        sorted(set(enc((-a*r,-b*r),mods) for a in range(len(lists)) for b in range(len(lists))))]


def factor(A,p):
    """Exact rank factorization by pivot columns; integers mod small p."""
    A=np.asarray(A,dtype=np.int64)%p
    if A.shape[0]==1:
        js=np.flatnonzero(A[0])
        if not len(js):return A[:,:0],A[:0,:]
        j=int(js[0]);return A[:,j:j+1],(A*pow(int(A[0,j]),-1,p))%p
    if A.shape[1]==1:
        return (A,np.ones((1,1),dtype=np.int64)) if np.any(A) else (A[:,:0],A[:0,:])
    R=A.copy(); k=0; piv=[]
    for j in range(A.shape[1]):
        ix=np.flatnonzero(R[k:,j])
        if not len(ix):continue
        i=k+int(ix[0]); R[[i,k]]=R[[k,i]]
        R[k]=(R[k]*pow(int(R[k,j]),-1,p))%p
        a=R[:,j].copy(); a[k]=0
        R=(R-np.outer(a,R[k]))%p
        piv.append(j);k+=1
        if k==A.shape[0]:break
    return A[:,piv],R[:k,:]


def zero(dims):
    return [np.zeros((1,n,1),dtype=np.int64) for n in dims]


def compress(G,p):
    """Two exact sweeps; each internal dimension equals the unfolding rank
    for a nonzero tensor. A zero tensor uses rank-one zero cores."""
    G=[a.copy()%p for a in G]; dims=[a.shape[1] for a in G]
    for j in range(len(G)-1):
        a,n,b=G[j].shape
        X,Y=factor(G[j].reshape(a*n,b),p)
        k=X.shape[1]
        if k==0:return zero(dims)
        G[j]=X.reshape(a,n,k)
        _,n2,b2=G[j+1].shape
        G[j+1]=(Y @ G[j+1].reshape(b,n2*b2)%p).reshape(k,n2,b2)
    for j in range(len(G)-1,0,-1):
        a,n,b=G[j].shape
        X,Y=factor(G[j].reshape(a,n*b),p)
        k=X.shape[1]
        if k==0:return zero(dims)
        G[j]=Y.reshape(k,n,b)
        a2,n2,_=G[j-1].shape
        G[j-1]=(G[j-1].reshape(a2*n2,a) @ X%p).reshape(a2,n2,k)
    return G


def sparse_tt(table,dims,p):
    """Right-to-left prefix trie, using a rank factorization at each level.
    The table may have duplicate words represented by integer multiplicities.
    """
    vals={tuple(x):np.array([a%p],dtype=np.int64) for x,a in table.items() if a%p}
    if not vals:return zero(dims)
    G=[None]*len(dims); R=1
    for j in range(len(dims)-1,-1,-1):
        prefixes=sorted({x[:-1] for x in vals}); index={x:i for i,x in enumerate(prefixes)}
        A=np.zeros((len(prefixes),dims[j]*R),dtype=np.int64)
        for x,v in vals.items():A[index[x[:-1]],x[-1]*R:(x[-1]+1)*R]=v
        if j==0:
            assert len(prefixes)==1
            G[0]=A.reshape(1,dims[0],R)
        else:
            X,Y=factor(A,p); k=X.shape[1]
            if not k:return zero(dims)
            G[j]=Y.reshape(k,dims[j],R)
            vals={x:X[i,:] for i,x in enumerate(prefixes)};R=k
    return compress(G,p)


def conv_tt(A,B,p):
    C=[]
    for a,b in zip(A,B):
        r,n,s=a.shape;t,n2,u=b.shape;assert n==n2
        if r==s==1:
            ix=np.flatnonzero(a[0,:,0])
            if len(ix)==1:
                x=int(ix[0]);C.append(np.roll(b,x,axis=1)*int(a[0,x,0])%p);continue
        if t==u==1:
            ix=np.flatnonzero(b[0,:,0])
            if len(ix)==1:
                x=int(ix[0]);C.append(np.roll(a,x,axis=1)*int(b[0,x,0])%p);continue
        c=np.zeros((r*t,n,s*u),dtype=np.int64)
        # skip identically zero slices; every step is exact modulo p
        ax=[x for x in range(n) if np.any(a[:,x,:])]
        by=[y for y in range(n) if np.any(b[:,y,:])]
        for x in ax:
            for y in by:
                z=(x+y)%n
                c[:,z,:]=(c[:,z,:]+np.kron(a[:,x,:],b[:,y,:]))%p
        C.append(c)
    return compress(C,p)


def value(G,x,p):
    u=np.array([1],dtype=np.int64)
    for a,v in zip(G,x):u=(u@a[:,v,:])%p
    return int(u[0])


def suffixes(G,p):
    S=[None]*(len(G)+1);S[-1]=np.array([1],dtype=np.int64)
    for j in range(len(G)-1,-1,-1):S[j]=(G[j].sum(axis=1)%p @ S[j+1])%p
    return S


def small_primes_bound(M):
    ps=[];P=1;p=2
    while P<=M:
        if isprime(p):ps.append(p);P*=p
        p+=1
    return ps


def prepare(tables,dims,cap):
    M=prod(sum(A.values()) for A in tables); ps=small_primes_bound(M)
    GG=[]; maxrank=1
    for p in ps:
        G=None
        for A in tables:
            T=sparse_tt(A,dims,p)
            rank=max(a.shape[-1] for a in T)
            maxrank = max(maxrank, rank)
            if rank>cap:return None,rank
            G=T if G is None else conv_tt(G,T,p)
            rank=max(a.shape[-1] for a in G)
            maxrank = max(maxrank, rank)
            if rank>cap:return None,rank
        GG.append(G)
    return (M,ps,GG),maxrank


def absent(tables,dims,cap):
    prep,width=prepare(tables,dims,cap)
    if prep is None:return None,width,0,0
    M,ps,GG=prep; assert M<prod(dims)
    SS=[suffixes(G,p) for G,p in zip(GG,ps)]
    P=[np.array([1],dtype=np.int64) for p in ps];out=[]; queries=0
    count=M; remain=prod(dims)
    for j,n in enumerate(dims):
        remain//=n
        children=[];next_states=[]
        for a in range(n):
            states=[(v@G[j][:,a,:])%p for G,v,p in zip(GG,P,ps)]
            residues=[int((v@S[j+1])%p) for v,S,p in zip(states,SS,ps)]
            c=crt(residues,ps);assert c<=M
            children.append(c);next_states.append(states);queries+=1
        assert sum(children)==count
        a=next(a for a,c in enumerate(children) if c<remain)
        out.append(a);count=children[a];P=next_states[a]
    assert count==0
    return tuple(out),width,len(ps),queries


def rank_of(G,x,p):
    return value(G,x,p)


def lift_tests():
    rng=Random(817); cases=0; inclusions=0; absent_checks=0; reverse=0
    for mods,r in [([2,3],7),([3,5],17),([2,5],11),([2,3,5],31)]:
        dims=tuple(mods)*2;R=prod(mods)
        for _ in range(9):
            A=[rng.sample(list(product(range(r),repeat=2)),rng.randrange(2,5)) for _ in range(3)]
            direct=support(A,(r,r)); lifts=wrap_lift(A,r,mods); S=support(lifts,dims)
            assert all(enc(x,mods) in S for x in direct); inclusions+=len(direct)
            for z in product(*(range(n) for n in dims)):
                if z not in S:
                    assert dec(z,mods) not in direct;absent_checks+=1
            cases+=1
            # Reverse reduction with a supplied prime r>=R: every z in
            # the prime plane whose reduction mod R is present is in the
            # corrected prime sumset.
            B=[rng.sample(list(product(range(R),repeat=2)),rng.randrange(2,5)) for _ in range(3)]
            SB=support(B,(R,R)); h=(r-1)//R
            C=list(product(range(-(len(B)-1),h+1),repeat=2))
            T=support(B+[[((a*R)%r,(b*R)%r) for a,b in C]],(r,r))
            for x in product(range(r),repeat=2):
                if tuple(v%R for v in x) in SB:assert x in T
                reverse+=1
    print(f'CRT lift: {cases} exact plane fixtures, {inclusions} represented-point inclusions, {absent_checks} safe absent decodings, {reverse} reverse-lift point checks.')


def tt_tests():
    rng=Random(73);cases=0; values=0; ranks=0; convs=0; gaps=0
    for dims in [(2,3),(2,2,3),(3,3,2),(2,2,2,3)]:
        pts=list(product(*(range(x) for x in dims)))
        for p in [2,3,5,11]:
            for _ in range(5):
                A=Counter({x:rng.randrange(1,7) for x in rng.sample(pts,rng.randrange(1,len(pts)+1))})
                B=Counter({x:rng.randrange(1,4) for x in rng.sample(pts,rng.randrange(1,len(pts)+1))})
                for table in [A,B,conv_dict(A,B,dims)]:
                    if table is A or table is B:G=sparse_tt(table,dims,p)
                    else:G=conv_tt(sparse_tt(A,dims,p),sparse_tt(B,dims,p),p);convs+=1
                    arr=np.array([table[x]%p for x in pts],dtype=np.int64)
                    assert all(value(G,x,p)==table[x]%p for x in pts);values+=len(pts)
                    for j in range(1,len(dims)):
                        u=arr.reshape(prod(dims[:j]),-1);rnk=factor(u,p)[0].shape[1]
                        assert G[j-1].shape[-1]==max(1,rnk);ranks+=1
                    cases+=1
    for dims in [(2,3,5),(3,3,3,2),(2,3,2,3,2)]:
        pts=list(product(*(range(x) for x in dims)))
        for _ in range(10):
            tables=[Counter(rng.sample(pts,2)) for j in range(3)]
            # nontrivial multiplicities in two cases where product permits
            if len(pts)>=54:tables[0][next(iter(tables[0]))]+=2
            exact=Counter({tuple(0 for n in dims):1})
            for A in tables:exact=conv_dict(exact,A,dims)
            z,rank,primes,queries=absent(tables,dims,200)
            assert z is not None and exact[z]==0;gaps+=1
    # Need >1 coefficient prime: tensor total 32, the first primes alone are insufficient.
    dims=(3,3,5); A=Counter({(0,0,0):4,(1,2,0):4});B=Counter({(0,0,0):3,(0,1,1):1})
    z,rank,ps,qs=absent([A,B],dims,20);assert ps>=3 and conv_dict(A,B,dims)[z]==0
    # Capped rank branch must honestly fail on a diagonal table.
    G,rank=prepare([Counter({(i,i):1 for i in range(5)})],(5,5),2)
    assert G is None and rank==5
    print(f'Tensor train: {cases} exact tensors, {values} point evaluations, {ranks} unfolding-rank comparisons, {convs} convolutions, {gaps+1} certified absent points; one rank-cap failure retained.')
    input_rank_test()
    intermediate_rank_cap_test()


def input_rank_test():
    dims = (2, 2, 3)
    first = Counter([(0, 0, 0), (1, 0, 0)])
    second = Counter([(0, 0, 0), (1, 1, 0)])
    tables = [first, second]
    final = conv_dict(first, second, dims)
    # The third coordinate is fixed, so only the first cut can have rank above one.
    for p in small_primes_bound(sum(final.values())):
        ranks = [sparse_matrix_rank(t, (0,), (1, 2), p) for t in (first, second, final)]
        assert ranks == [1, 2, 1]
    rejected, width = prepare(tables, dims, cap=1)
    assert rejected is None and width == 2
    point, peak, _, _ = absent(tables, dims, cap=2)
    assert point is not None and final[point] == 0
    assert peak == 2, 'reported peak must include the rank-two input'
    print(f'Input ranks 1 and 2 convolved to rank 1. Reported peak {peak}, cap 1 rejected.')


def intermediate_rank_cap_test():
    dims = (4, 4, 2)
    first = Counter([(0, 0, 0), (1, 1, 0)])
    second = Counter([(0, 0, 0), (2, 2, 0)])
    third = Counter((i, 0, 0) for i in range(4))
    tables = [first, second, third]
    intermediate = conv_dict(first, second, dims)
    final = conv_dict(intermediate, third, dims)
    for p in small_primes_bound(sum(final.values())):
        ranks = [sparse_matrix_rank(t, (0,), (1, 2), p)
                 for t in (*tables, intermediate, final)]
        assert ranks == [2, 2, 1, 4, 1]
    rejected, width = prepare(tables, dims, cap=2)
    assert rejected is None and width == 4, 'cap 2 must reject the rank-four intermediate'
    point, peak, _, _ = absent(tables, dims, cap=4)
    assert point is not None and final[point] == 0 and peak == 4
    print(f'Intermediate rank {width} exceeded cap 2. Cap {peak} passed with final rank 1.')


def block_codes(q,L,k,blocks):
    assert isprime(q) and L<=q
    dims=(q,)*(L*blocks); local=[]
    for c in product(range(q),repeat=k):
        local.append(tuple(sum(a*pow(x,j,q) for j,a in enumerate(c))%q for x in range(L)))
    assert len(set(local))==q**k
    # All k-coordinate projections full, checked on one repeated code block.
    from itertools import combinations
    for I in combinations(range(L),k):assert len({tuple(x[i] for i in I) for x in local})==q**k
    tables=[]
    for a in range(blocks):
        tables.append(Counter({(0,)*(a*L)+x+(0,)*((blocks-a-1)*L):1 for x in local}))
    z,rank,ps,queries=absent(tables,dims,q**k)
    assert z is not None and any(z[i*L:(i+1)*L] not in set(local) for i in range(blocks))
    print(f'Product-code fixture: q={q}, {len(dims)} coordinates, {blocks} explicit lists of {len(local)} entries; ambient={q}^{len(dims)}, choices={q**k}^{blocks}; all projections onto <= {k} coordinates full; exact peak bond {rank}, {ps} coefficient primes, {queries} child counts.')


def mixed_code_fixture():
    # Reorder pair coordinates consecutively by modulus for the tensor order.
    mods=[2,3,5,7]; R=prod(mods); r=R+1
    while not isprime(r):r+=1
    dims=tuple(p for p in mods for _ in range(2)); Q=53
    assert Q>=max(mods)**2 and isprime(Q)
    def enc_interleaved(z):return tuple(x%p for p in mods for x in z)
    def dec_interleaved(z):return crt(z[0::2],mods),crt(z[1::2],mods)
    lists=[]
    for i in range(0,len(mods),2):
        F=set()
        for t in range(Q):
            x=[0]*len(dims)
            for j in (i,i+1):
                p=mods[j];a=t%(p*p);x[2*j]=a%p;x[2*j+1]=a//p
            F.add(tuple(x))
        for j in (i,i+1):assert len({x[2*j:2*j+2] for x in F})==mods[j]**2
        lists.append(sorted(F))
    original=[[dec_interleaved(x) for x in A] for A in lists]
    assert all(set(enc_interleaved(x) for x in A)==set(B) for A,B in zip(original,lists))
    C=sorted({enc_interleaved((-a*r,-b*r)) for a in range(len(lists)) for b in range(len(lists))})
    tables=[Counter(A) for A in lists+[C]]
    M=prod(map(len,tables));assert M<R*R
    z,rank,ps,queries=absent(tables,dims,1000)
    assert z is not None
    # Direct verification using block product: z absent iff for each wrap
    # correction some block slice is not in that block's local code.
    for c in C:
        y=tuple((x-a)%p for x,a,p in zip(z,c,dims))
        assert any(y[4*i:4*i+4] not in {x[4*i:4*i+4] for x in A} for i,A in enumerate(lists))
    print(f'Mixed CRT fixture: r={r}, R={R}, distinct primes={mods}, list sizes={[len(A) for A in tables]}, {M} tuple choices < {R*R} ambient; every single-prime pair projection full; peak bond {rank}, {ps} coefficient primes, {queries} child counts.')


def shadow(lists,dims,corrections,Q,cap):
    """Cross-characteristic affine-space overapproximation, exact rational guard.
    Returns an absent point, or honest failure if the guard/cap does not pass.
    """
    from fractions import Fraction
    d=len(dims);L=len(corrections)
    assert all(lists) and corrections and isprime(Q)
    refs=[A[0] for A in lists]
    a0=np.array([sum(a[j] for a in refs)%Q for j in range(d)],dtype=np.int64)
    columns=np.array([[(x-y)%Q for x,y in zip(a,ref)] for A,ref in zip(lists,refs) for a in A],dtype=np.int64).T
    V,_=factor(columns,Q);k=V.shape[1]
    # Pivot columns of V^T pick an information coordinate set.
    if k:
        I=[]; selected=[]
        for j in range(d):
            test=selected+[V[j,:].tolist()]
            if factor(np.array(test),Q)[0].shape[1]>len(selected):
                selected=test;I.append(j)
                if len(I)==k:break
        assert len(I)==k
        J=V[I,:]; Aug=np.concatenate([J,np.eye(k,dtype=np.int64)],axis=1)
        X,Y=factor(Aug,Q);inv=Y[:,k:];assert np.array_equal(J@inv%Q,np.eye(k,dtype=np.int64))
        T=(V@inv)%Q
    else:I=[];T=np.zeros((d,0),dtype=np.int64)
    U=[sum(max(a[j] for a in A) for A in lists) for j in range(d)]
    widths=[1+u//p for u,p in zip(U,dims)]
    n0=L*prod(widths[j] for j in I)
    beta=[min(Fraction(1),Fraction(w*((p+Q-1)//Q),p)) for w,p in zip(widths,dims)]
    guard=Fraction(n0)*prod(beta[j] for j in range(d) if j not in I)
    if n0>cap or guard>=1:return None,k,n0,guard
    def choices(j,z):
        return {(z+t*dims[j])%Q for t in range(widths[j]) if z+t*dims[j]<=U[j]}
    z=[0]*d;candidates=[]
    for c in corrections:
        for vals in product(*(sorted(choices(j,(-c[j])%dims[j])) for j in I)):
            v=(a0+T@(np.array(vals,dtype=np.int64)-a0[I]))%Q
            candidates.append((c,tuple(map(int,v))))
    actual=len(candidates)
    for j,p in enumerate(dims):
        if j in I:continue
        options=[]
        for a in range(p):
            new=[(c,v) for c,v in candidates if v[j] in choices(j,(a-c[j])%p)]
            options.append(new)
        a=min(range(p),key=lambda a:len(options[a]));z[j]=a;candidates=options[a]
    assert not candidates
    return tuple(z),k,actual,guard


def mds_lists(dims,Q,k,m):
    """Split consecutive coefficients of a common-field RS evaluation.
    Lists are then regarded as mixed-modulus residues, not a homomorphism.
    """
    d=len(dims);assert d<Q<=min(dims) and m<=k<=d
    parts=[];start=0
    for i in range(m):
        a=k//m+int(i<k%m) # equal-as-possible consecutive groups
        F=[]
        for cs in product(range(Q),repeat=a):
            F.append(tuple(sum(c*pow(t,start+j,Q) for j,c in enumerate(cs))%Q for t in range(1,d+1)))
        assert len(set(F))==Q**a
        parts.append(F);start+=a
    assert start==k
    return parts


def shadow_tests():
    from itertools import combinations
    # Full low-coordinate projections in mixed moduli, true distinct-prime lifts.
    # Small sizes deliberately do not impose the asymptotic guards.
    cases=0;cuts=0;projection_checks=0;certs=0
    for dims,Q,k,m in [((5,5,7,7),5,3,2),((11,11,13,13,17,17),11,3,2),((11,11,13,13,17,17),11,4,2),((7,7,11,11),7,2,2)]:
        if Q<=len(dims):continue
        A=mds_lists(dims,Q,k,m);M=prod(map(len,A));d=len(dims);s=d//2
        exact=Counter({(0,)*d:1})
        for B in A:exact=conv_dict(exact,Counter(B),dims)
        mods=dims[0::2];R=prod(mods);r=R+1
        while not isprime(r):r+=1
        C=sorted({tuple((-a*r if j%2==0 else -b*r)%p for j,p in enumerate(dims)) for a in range(m) for b in range(m)})
        shifted=conv_dict(exact,Counter(C),dims)
        assert sum(exact.values())==M and sum(shifted.values())==M*len(C)
        C0=1+ m*(Q-1)//min(dims);Bnd=C0**s*Q**max(0,k-s)
        for I in combinations(range(d),s):
            J=tuple(j for j in range(d) if j not in I)
            row=Counter();col=Counter();r2=Counter();c2=Counter()
            for x,v in exact.items():row[tuple(x[j] for j in I)]+=v;col[tuple(x[j] for j in J)]+=v
            assert max(row.values())<=Bnd and max(col.values())<=Bnd
            for x,v in shifted.items():r2[tuple(x[j] for j in I)]+=1;c2[tuple(x[j] for j in J)]+=1
            assert max(r2.values())<=len(C)*Bnd and max(c2.values())<=len(C)*Bnd
            assert len(c2)*Bnd>=M;cuts+=1
        if m*(Q-1)>=max(dims)-1:
            for I in combinations(range(d),k//m):
                assert len({tuple(x[j] for j in I) for x in exact})==prod(dims[j] for j in I);projection_checks+=1
        z,rank,n,guard=shadow(A,dims,C,Q,200000)
        assert rank==k
        if z is not None:assert shifted[z]==0;certs+=1
        cases+=1
        print(f'Common-field shadow fixture: Q={Q}, d={d}, rank={rank}, choices={M}, correction size={len(C)}, shadow guard {"passes" if z is not None else "fails"}.')
    # A direct positive guard case with all single-coordinate marginals full:
    # each one-list evaluation is a constant word and every alphabet <= Q.
    # Use the large-alphabet mixed case d=12; Q=17, p_j=17..29, m=2.
    dims=(17,17,19,19,23,23,29,29,31,31,37,37);Q=19
    # choose Q<=min dims: use Q=17, m=3, constant/multi-coefficient lists
    Q=17;k=3;m=3;A=mds_lists(dims,Q,k,m)
    C=[(0,)*len(dims)]
    z,rank,n,guard=shadow(A,dims,C,Q,50000)
    assert z is not None and rank==k and z not in support(A,dims)
    certs+=1
    for j,p in enumerate(dims):assert len({sum(x)%p for x in product(*(sorted({a[j] for a in B}) for B in A))})==p
    print(f'Shadow rescue: d={len(dims)}, all {len(dims)} single-coordinate projections full, rank={rank}, {n} initial affine candidates, exact guard={guard}; absent point verified.')
    # full-rank guard failure, despite enormous volume slack
    A=[[(0,)*5]+[tuple(int(i==j) for i in range(5)) for j in range(5)]]
    z,k,n,g=shadow(A,(11,)*5,[(0,)*5],7,1000)
    assert z is None and k==5
    print(f'All-order width/carry tests: {cases} mixed MDS fixtures, {cuts} balanced cuts, {projection_checks} full low-coordinate projections; {certs} certified shadow successes; full-rank negative control retained.')
    translated_shadow_tests()


def translated_shadow_tests():
    from fractions import Fraction
    dims = (5, 5)
    lists = [[(0, 1), (1, 2)]]
    corrections = [(0, 0)]
    ambient = set(product(*(range(p) for p in dims)))
    point, rank, count, guard = shadow(lists, dims, corrections, 5, cap=1)
    assert (rank, count, guard) == (1, 1, Fraction(1, 5))
    assert point in ambient - support(lists + [corrections], dims), 'translated shadow returned a represented or invalid point'
    print(f'Translated affine line over F5 returned {point}, verified absent with guard {guard}.')

    corrections = [(0, 0), (0, 1)]
    rejected, rank, bound, guard = shadow(lists, dims, corrections, 5, cap=1)
    assert (rank, bound, guard) == (1, 2, Fraction(2, 5))
    assert rejected is None, 'candidate cap 1 must reject a bound of two'
    point, rank, count, same_guard = shadow(lists, dims, corrections, 5, cap=2)
    assert (rank, count, same_guard) == (1, bound, guard)
    assert point in ambient - support(lists + [corrections], dims)
    print(f'Candidate bound {bound} with guard {guard} rejected at cap 1 and certified {point} absent at cap 2.')


def sparse_matrix_rank(table,I,J,p):
    """Independent sparse elimination on the *explicit tiny* unfolding."""
    from collections import defaultdict
    rows=defaultdict(dict)
    for x,a in table.items():
        if a%p:rows[tuple(x[j] for j in I)][tuple(x[j] for j in J)]=a%p
    basis={}
    for row in rows.values():
        row=row.copy()
        while row:
            j=min(row);a=row[j]
            if j not in basis:
                inv=pow(a,-1,p);basis[j]={h:(v*inv)%p for h,v in row.items()};break
            for h,v in basis[j].items():
                x=(row.get(h,0)-a*v)%p
                if x:row[h]=x
                else:row.pop(h,None)
    return len(basis)


def rank_barrier_checks():
    from itertools import combinations
    from fractions import Fraction
    dims=(11,11,13,13,17,17);Q=11;k=3;m=2;d=len(dims);s=d//2
    A=mds_lists(dims,Q,k,m);M=prod(map(len,A));R=11*13*17;r=R+1
    while not isprime(r):r+=1
    C=sorted({tuple((-a*r if j%2==0 else -b*r)%p for j,p in enumerate(dims)) for a in range(m) for b in range(m)})
    S=Counter({(0,)*d:1})
    for B in A:S=conv_dict(S,Counter(B),dims)
    S=conv_dict(S,Counter(C),dims)
    C0=1+m*(Q-1)//min(dims);Bnd=C0**s
    assert max(S.values())<=len(C)*Bnd<101
    ranks=[]
    for I in combinations(range(d),s):
        if 0 not in I:continue
        J=tuple(j for j in range(d) if j not in I)
        h=sparse_matrix_rank(S,I,J,101)
        assert h*len(C)*Bnd**2>=M
        ranks.append(h)
    z,k2,n,g=shadow(A,dims,C,Q,20000);assert z is not None and S[z]==0
    print(f'Exact high-width rescue: {len(ranks)} distinct balanced unfoldings have ranks {min(ranks)}..{max(ranks)} over F101; common-field affine rank {k2}, {n} initial candidates, shadow certifies absence.')
    # Check the asymptotic integer inequalities without constructing a huge
    # prime or enumerating the astronomical dictionary. 2**b-1 is only a
    # numeric upper-endpoint surrogate, not claimed to be prime.
    def sieve(n):
        a=bytearray(b'\x01')*(n+1);a[0:2]=b'\x00\x00'
        for p in range(2,int(n**.5)+1):
            if a[p]:a[p*p:n+1:p]=b'\x00'*(((n-p*p)//p)+1)
        return a
    for b,m in [(1024,4),(8192,5),(65536,6),(262144,7)]:
        r0=2**b-1;R=1;ps=[];flags=sieve(4*b)
        for p in range(b,4*b+1):
            if not flags[p]:continue
            if R*p>r0:break
            R*=p;ps.append(p)
        assert ps and R<=r0<4*b*R
        s=len(ps);d=2*s;Q=(max(ps)-1+m-1)//m+1
        while not flags[Q]:Q+=1
        k=0;M=1
        while M<r0:M*=Q;k+=1
        C0=1+m*(Q-1)//min(ps)
        # Conservative uniform bound on the shadow's failure mass:
        # ceil(p_j/Q)/p_j <= 1/Q+1/p_min.
        beta=Fraction(C0,Q)+Fraction(C0,min(ps));L=m*m
        guard=Fraction(L*C0**k)*beta**(d-k) if k<=d else Fraction(10**100)
        ok=(d<Q<=min(ps) and m*(Q-1)>=max(ps)-1 and s<=k<d and L*M<R*R and guard<1)
        # Rational (over Q) rank lower bound, and modular bound with <= b*2+10 primes.
        num=M;den=L*(C0**s*Q**max(0,k-s))**2*(2*b+10)
        print(f'Numeric guard b={b}, m={m}: d={d}, Q={Q}, k={k}, C={C0}; {"PASS" if ok else "FAIL retained"}; conservative log2(rank lower bound) >= {num.bit_length()-den.bit_length()-1}.')
        if b>=8192:assert ok


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    lift_tests();tt_tests();mixed_code_fixture();block_codes(5,4,2,12);shadow_tests();rank_barrier_checks()
