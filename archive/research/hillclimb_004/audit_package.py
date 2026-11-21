"""Independent exact checks of selected additions in the supplied ZIP.

This is a lemma audit, not an implementation of its assumed contractor.
Only the selector implementation is imported; probability and join checks
are independent. Run with Python assertions enabled.
"""
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parent
WORK = ROOT.parent.parent
PACKAGE = WORK/'_review/compressed_relation_search_package/compressed_relation_search'
RERUN = WORK/'_review/reruns/compressed_relation_search_20261031'
ARCHIVE_HASH = '95f2078edf80cc3ce1090be7a4216fa1f407e91044d376d035c99ebf7b9204c5'


def spectral():
    edges = list(combinations(range(9), 2))
    parts = []
    for a in combinations(range(1, 9), 2):
        A = {0, *a}
        rest = set(range(9))-A
        first = min(rest)
        for b in combinations(sorted(rest-{first}), 2):
            B = {first, *b}
            parts.append((A, B, rest-B))
    assert len(parts) == 280
    counts = [[0]*36 for _ in edges]
    for part in parts:
        label = {v:i for i, group in enumerate(part) for v in group}
        bundles = defaultdict(list)
        for e, (v, w) in enumerate(edges):
            if label[v] != label[w]:
                bundles[tuple(sorted((label[v], label[w])))].append(e)
        for bundle in bundles.values():
            for e in bundle:
                for f in bundle:
                    counts[e][f] += 1
    for e, edge in enumerate(edges):
        for f, other in enumerate(edges):
            expected = F(1, 7)+F(len(set(edge)&set(other)), 14)+F(13*(e == f), 28)
            assert F(counts[e][f], 280) == expected

    # Derive the low-t averaged forms from their actual 36 edge coefficients,
    # not from the closed expressions for m_k and R_k in the package.
    beta, delta = F(39, 68), F(347, 425)
    B2, residual = F(34641, 2800), F(1751, 5600)
    triple = {0, 1, 2}
    external = [i for i, e in enumerate(edges) if len(set(e)&triple) == 1]
    assert len(external) == 18
    bounds = [F(n, 1000) for n in (3545, 3598, 3642, 3673)]
    means = []
    for types in product((0, 1), repeat=3):
        c = [F(0) for _ in edges]
        for last in triple:
            alpha, gamma = (beta, delta) if types[last] else (delta, beta)
            pair = triple-{last}
            for i, edge in enumerate(edges):
                e = set(edge)
                if last in e and e & pair:
                    c[i] += alpha/3
                elif len(e & triple) == 1:
                    c[i] += (gamma if last in e else F(1))/3
        mean = sum(c[i] for i in external)/18
        orthogonal2 = sum(x*x for x in c)-18*mean*mean
        k = sum(types)
        expected_mean = (6+k*delta+(3-k)*beta)/9
        expected_orthogonal2 = ((k*beta+(3-k)*delta)**2+k*beta**2+
                               (3-k)*delta**2+2*k*(3-k)*(delta-beta)**2)/9
        assert mean == expected_mean and orthogonal2 == expected_orthogonal2
        assert mean*mean >= orthogonal2/6
        u, v, cap = mean*mean*B2, orthogonal2*residual, bounds[k]
        gap = cap*cap-u-v
        assert gap > 0 and gap*gap > 4*u*v and cap*cap < 14
        means.append(dict(types=types, mean=str(mean), residual_squared=str(orthogonal2)))
    return dict(partitions=280, exact_averaging_matrix_entries=1296,
                independently_derived_averaged_forms=means)


def extend(span, x):
    return span | sum(1 << (v ^ x) for v in range(8) if span >> v & 1)


def inverse_moments():
    rng = random.Random(421759)
    weights = [[1]*8, [10000, 1, 1, 1, 1, 1, 1, 1],
               [0, 10000, 1, 1, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0, 0, 0]]
    weights += [[rng.randrange(0, 30) for _ in range(8)] for _ in range(12)]
    checks = 0
    for row in weights:
        law = [F(n, sum(row)) for n in row]
        support = 1
        for x, p in enumerate(law):
            if p:
                support = extend(support, x)
        d = support.bit_count().bit_length()-1
        spans = {1: F(1)}
        for t in range(1, 8):
            following = defaultdict(F)
            for span, mass in spans.items():
                for x, p in enumerate(law):
                    if p:
                        following[extend(span, x)] += mass*p
            spans = following
            assert sum(spans.values()) == 1
            for j in range(1, 4):
                if t <= d+j-1:
                    continue
                moment = sum(prob/sum(law[x] for x in range(8) if span >> x & 1)**j
                             for span, prob in spans.items())
                bound = F(1)
                for i in range(j):
                    bound *= F(t-i, t-d-i)
                assert moment <= bound
                checks += 1
    return dict(distributions=len(weights), exact_inverse_moment_inequalities=checks,
                field='F_2', ambient_dimension=3, maximum_training_length=7)


def selectors():
    source = PACKAGE/'checks/unit_simplex.py'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == '615af6e85bc05a511c273ddbe72dc4d2704b22f64e057303ec6fe58567b41cbb'
    spec = importlib.util.spec_from_file_location('reviewed_unit_simplex', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # New degree-four extension fixture; x^4+x+1 has all 15 nonzero residues
    # invertible, checked here rather than assumed by the Ring constructor.
    ring = module.Ring(11, (1, 1, 0, 0, 1), q=3)
    for v in range(1, 16):
        assert ring.mul(ring.symbol(v), ring.inv(ring.symbol(v))) == ring.one
    rng = random.Random(87112)
    cardinal_checks = off_grid_checks = 0
    for r, k in [(2, 4), (3, 4)]:
        alphabet = module.Alphabet(r, k, ring)
        polys = [alphabet.expanded(i) for i in range(alphabet.Q)]
        total = defaultdict(lambda: ring.zero)
        for poly in polys:
            for e, v in poly.items():
                assert sum(e) <= k
                total[e] = ring.add(total[e], v)
        assert {e:v for e, v in total.items() if v != ring.zero} == {(0,)*r:ring.one}

        def evaluate(poly, point):
            out = ring.zero
            for powers, coefficient in poly.items():
                term = coefficient
                for x, power in zip(point, powers):
                    term = ring.mul(term, ring.power(x, power))
                out = ring.add(out, term)
            return out

        for i, poly in enumerate(polys):
            for j, point in enumerate(alphabet.points):
                assert evaluate(poly, point) == (ring.one if i == j else ring.zero)
                cardinal_checks += 1
        for _ in range(6):
            point = [tuple(rng.randrange(ring.N) for _ in range(ring.dim)) for _ in range(r)]
            actual = [evaluate(poly, point) for poly in polys]
            assert actual == alphabet.selectors(point)
            off_grid_checks += len(actual)
    return dict(nonzero_residue_inverses=15, monomial_cardinal_checks=cardinal_checks,
                off_grid_nilpotent_comparisons=off_grid_checks,
                fixture='(Z/2048)[u,e]/(u^4+u+1,e^3), (r,k)=(2,4),(3,4)')


def main():
    if not __debug__:
        raise RuntimeError('Do not disable assertions')
    archive = WORK/'compressed_relation_search_package.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == ARCHIVE_HASH
    manifest = json.loads((PACKAGE.parent/'review_manifest.json').read_text())
    for name, expected in manifest['files'].items():
        assert hashlib.sha256((PACKAGE.parent/name).read_bytes()).hexdigest() == expected
    rerun = json.loads((RERUN/'checks/run_report.json').read_text())
    assert rerun['status'] == 'passed' and len(rerun['checks']) == 8
    for record in rerun['checks']:
        result = RERUN/'checks'/record['record']
        assert hashlib.sha256(result.read_bytes()).hexdigest() == record['record_sha256']
    for name, expected in rerun['source_sha256'].items():
        assert hashlib.sha256((RERUN/name).read_bytes()).hexdigest() == expected
    result = dict(status='passed', archive_sha256=ARCHIVE_HASH,
                  package_suites_passed=8, spectral=spectral(),
                  inverse_moments=inverse_moments(), selectors=selectors(),
                  scope='Selected finite algebra and probability checks. No proof or execution of the assumed asymptotic packet contractor.',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  python=sys.version.split()[0])
    run = ROOT/'runs'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    run.mkdir(parents=True)
    (run/'independent_audit.json').write_text(json.dumps(result, indent=2)+'\n')
    (ROOT/'latest_run.json').write_text(json.dumps(dict(path=str(run.relative_to(ROOT)), status='passed'), indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
