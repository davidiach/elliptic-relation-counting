"""Recheck PR #12 with ordinary complete short-fiber enumeration.
Imports the existing finite fixture generator from a supplied directory.
"""
import sys
from itertools import product
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1])))
from check_dual_lift import dictionary, selected_rows, rank_mod

for r,m,extra,seed in [(101,2,8,17),(257,3,12,22),(1009,2,10,17)]:
    D,zeta,tau,B=dictionary(r,m,extra,seed)
    A,u,lam,_,_=selected_rows(r,m,D,zeta,tau,B)
    n=len(u);starts=[];j=0
    for F in B:starts.append(j);j+=len(F)-1
    # The selected basis in PR12 has exactly the same row space as A.
    assert rank_mod(A,r,n)==n-1
    Xs=[]
    for X in product(*(range(len(F)) for F in B)):
        if sum(B[i][j][0] for i,j in enumerate(X))%r:continue
        c=[0]*n
        for i,j in enumerate(X):
            if j:c[starts[i]+j-1]+=1
        Xs.append(c)
    scores={sum(x*y for x,y in zip(c,u))%r for c in Xs}
    assert len(scores)>1
    contrast=next([x-y for x,y in zip(c,Xs[0])] for c in Xs if sum(x*y for x,y in zip(c,u))%r !=sum(x*y for x,y in zip(Xs[0],u))%r)
    assert sum(x*y for x,y in zip(contrast,lam))%r==0
    assert rank_mod(A+[contrast],r,n)==n
    print(f'PR12 C{r}: identity short fiber {len(Xs)} tuples, {len(scores)} distinct quotient scores; one genuine zero-shift contrast raises rank {n-1} to {n}.')
