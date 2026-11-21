"""Polynomial certificate over Z for the target-dependent second moment.

The specification uses a literal determinant and an exhaustive three-sign sum.
The proposed readout uses only squared Vandermonde polynomials. Neither side
imports the project solver, experiment code, or curve fixtures.
"""
from itertools import permutations, product

from sympy.polys.domains import ZZ
from sympy.polys.rings import ring


def determinant(matrix, algebra):
    total = algebra.zero
    for perm in permutations(range(4)):
        inv = sum(perm[i] > perm[j] for i in range(4) for j in range(i+1,4))
        term = algebra.one
        for i, j in enumerate(perm):
            term *= matrix[i][j]
        total += (-1)**inv * term
    return total


def check():
    algebra, x1,x2,x3,X,y1,y2,y3,Y = ring("x1,x2,x3,X,y1,y2,y3,Y",ZZ)
    xs = [x1,x2,x3,X]
    ys = [y1,y2,y3,Y]
    specification = algebra.zero
    for signs in product((-1,1),repeat=3):
        signed = [e*y for e,y in zip(signs,ys[:3])] + [-Y]
        rows = [[algebra.one,x,x*x,v] for x,v in zip(xs,signed)]
        specification += determinant(rows,algebra)**2
    readout = algebra.zero
    for i in range(4):
        u,v,w = xs[:i]+xs[i+1:]
        readout += 8*ys[i]**2*((v-u)*(w-u)*(w-v))**2
    if specification != readout:
        raise AssertionError("Nonzero polynomial discrepancy")
    print(f"Literal sign-summed square and four squared minors agree over Z "
          f"({len(readout)} nonzero monomials).")


if __name__ == "__main__":
    check()
