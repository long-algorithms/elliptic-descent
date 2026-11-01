"""Independent release checks using eclib, rational torsion and division points.

No class-group computation, normal-field cocycle, or expected addition table is
used for the reference. Four small fixtures also compare EVERY named covering
and EVERY sum with E(Q)/4E(Q), using exact rational points on the actual models.
Point-search failure is an audit failure, never evidence of insolubility.
"""
import argparse
from functools import lru_cache
from itertools import product
import json
from math import gcd, isqrt, lcm
from pathlib import Path
from time import perf_counter

from foundation import metadata

# None denotes the origin. The first point is used only for the rank-one lower
# bound; rank-zero bounds come independently from eclib and rational torsion.
FIXTURES = {
    'rank_one': ((0, 0, 0, -1, 1), (0, 1), True),
    'split_two_torsion': ((0, 0, 0, -1, 0), None, True),
    'cyclic_two_torsion': ((0, 0, 0, 1, 0), None, True),
    'cyclic_four_torsion': ((0, 0, 0, 4, 0), None, True),
    'full_two_and_four': ((0, -5, 0, -16, 80), None, False),
    'zero_selmer': ((0, -1, 1, -10, -20), None, False),
}


def integer_roots(a, b, c):
    if a:
        disc = b*b-4*a*c
        if disc < 0:
            return ()
        root = isqrt(disc)
        if root*root != disc:
            return ()
        return tuple(set(v//(2*a) for v in (-b+root, -b-root) if v % (2*a) == 0))
    if b:
        return (-c//b,) if -c % b == 0 else ()
    return None if c == 0 else ()


def small_point(quadrics, *, bound=32):
    """Positive witness discovery only; exact integer quadratic elimination.

    This independent audit path does not use the package's local solver,
    transformations, torsion or cocycle code. It is deliberately bounded.
    """
    import sympy as sp
    symbols = sorted(set().union(*(q.free_symbols for q in quadrics)), key=str)
    if len(symbols) != 4:
        raise ValueError('the audit expects four quadric coordinates')
    pairs = tuple((i, j) for i in range(4) for j in range(i, 4))
    rows = []
    for q in quadrics:
        p = sp.Poly(q, *symbols, domain=sp.QQ)
        cs = [p.coeff_monomial(symbols[i]*symbols[j]) for i, j in pairs]
        den = lcm(*(int(v.q) for v in cs))
        rows.append(tuple(int(v*den) for v in cs))
    for i in range(4):
        p = tuple(int(i == j) for j in range(4))
        if all(sum(c*p[i]*p[j] for c, (i, j) in zip(row, pairs)) == 0 for row in rows):
            return p
    heights = [1]
    while heights[-1] < bound:
        heights.append(min(bound, 2*heights[-1]))
    for height in heights:
        for a, b, c in product(range(-height, height+1), repeat=3):
            if not (a or b or c) or next(v for v in (a, b, c) if v) < 0:
                continue
            equations = [(r[9], r[3]*a+r[6]*b+r[8]*c,
                          r[0]*a*a+r[1]*a*b+r[2]*a*c+r[4]*b*b+r[5]*b*c+r[7]*c*c)
                         for r in rows]
            roots, other = integer_roots(*equations[0]), equations[1]
            if roots is None:
                roots, other = integer_roots(*equations[1]), equations[0]
            for d in (0,) if roots is None else roots:
                if (abs(d) <= height and gcd(a, b, c, d) == 1
                        and other[0]*d*d+other[1]*d+other[2] == 0):
                    return a, b, c, d
    raise ArithmeticError('no bounded rational witness found; no insolubility claim follows')


def evaluate_map(covering, point, E):
    """Evaluate the exported polynomial maps with Sage rational arithmetic.

    This does not call MarkedFourCover.evaluate or the package group law.
    Its result is checked by Sage's elliptic-curve constructor.
    """
    from sage.all import QQ
    import sympy as sp
    from elliptic_descent.geometry import U
    from elliptic_descent.two_cover import X, Y, Z

    def evaluate(expr, variables, values):
        result = QQ(0)
        for powers, coeff in sp.Poly(expr, *variables, domain=sp.QQ).terms():
            term = QQ(str(coeff))
            for value, power in zip(values, powers):
                term *= value**power
            result += term
        return result

    p = tuple(QQ(v) for v in point)
    if not any(p) or any(evaluate(q, U, p) for q in covering.model.quadrics):
        raise ArithmeticError('invalid rational covering point')
    q = tuple(evaluate(v, U, p) for v in covering.model.map_xyz)
    a, b, c = (evaluate(v, (X, Y, Z), q) for v in covering.marking.source.map_weighted)
    J = covering.marking.change.source
    a1, a2, a3, a4, a6 = (QQ(str(v)) for v in J.ainvs)
    if b*b+a1*a*b*c+a3*b*c**3 != a**3+a2*a*a*c*c+a4*a*c**4+a6*c**6:
        raise ArithmeticError('exported map failed its weighted target equation')
    if covering.marking.translation is not None:
        raise ValueError('this audit expects the assembled origin-preserving marking')
    if not c:
        if not (a or b):
            raise ArithmeticError('exported covering map has an indeterminacy')
        return E(0)
    u, r, s, t = (QQ(str(v)) for v in covering.marking.change.urst)
    xx = (a-r*c*c)/(u*u*c*c)
    yy = (b-s*(a-r*c*c)*c-t*c**3)/(u**3*c**3)
    return E(xx, yy)


def audit_fixture(name, *, require_installed=False):
    from sage.all import EllipticCurve
    from sage.env import SAGE_VERSION
    import elliptic_descent as ed
    if require_installed and Path(ed.__file__).resolve().is_relative_to(Path(__file__).resolve().parents[1]):
        raise ArithmeticError('the release audit imported the checkout, not the installed wheel')
    ainvs, free_point, labelled = FIXTURES[name]
    E = EllipticCurve(list(ainvs))
    torsion = tuple(E.torsion_points())
    torsion_two = tuple(P for P in torsion if (2*P).is_zero())
    t = len(torsion_two).bit_length()-1
    if 2**t != len(torsion_two):
        raise ArithmeticError('rational two-torsion is not an elementary two-group')
    eclib_dimension = int(E.global_minimal_model().mwrank_curve(verbose=False).selmer_rank())
    rank = int(free_point is not None)
    if free_point is not None:
        P = E(*free_point)
        if P in torsion or P.division_points(2):
            raise ArithmeticError('the supplied point does not prove the required mod-two lower bound')
    # E(Q)/2E(Q) injects into Sel2. Independent torsion plus the checked
    # infinite point saturate the eclib upper bound, so Sha[2]=Sha[4]=0.
    if eclib_dimension != t+rank:
        raise ArithmeticError('independent bounds do not prove Sha[2]=0 for this fixture')
    four_torsion = {4*P for P in torsion}
    representatives = []
    for P in torsion:
        if not any(P-Q in four_torsion for Q in representatives):
            representatives.append(P)
    twice = []
    for P in representatives:
        if not any(2*P-Q in four_torsion for Q in twice):
            twice.append(2*P)
    n = len(representatives).bit_length()-1
    r = len(twice).bit_length()-1
    if 2**n != len(representatives) or 2**r != len(twice):
        raise ArithmeticError('torsion quotient is not a two-group')
    expected = (2,)*(n-2*r)+(4,)*(r+rank)
    start = perf_counter()
    result = ed.four_descent(ainvs, proof=True)
    G = result.require_group(certified=True)
    if not result.verify_evidence() or G.invariant_factors != expected:
        raise ArithmeticError('computed Sel4 disagrees with the independently determined rational quotient')
    point_records, sums, separations = [], 0, 0
    if labelled:
        images = []
        for element in G.elements:
            point = small_point(element.covering.model.quadrics)
            image = evaluate_map(element.covering, point, E)
            images.append(image)
            point_records.append({'key': element.key, 'point': point,
                                  'image': None if image.is_zero() else [str(v) for v in image.xy()]})

        @lru_cache(maxsize=None)
        def divisible_by_four(P):
            halves = P.division_points(4)  # All rational division points, not a bounded search.
            if any(4*Q != P for Q in halves):
                raise ArithmeticError('Sage division-point witness failed exact multiplication')
            return bool(halves)

        for i, P in enumerate(images):
            for j, Q in enumerate(images):
                if divisible_by_four(P-Q) != (i == j):
                    raise ArithmeticError('distinct named coverings have equal rational Kummer classes')
                separations += 1
                matches = [k for k, R in enumerate(images) if divisible_by_four(P+Q-R)]
                if matches != [G.law.table[i][j]]:
                    raise ArithmeticError('covering addition disagrees with independent rational point division')
                sums += 1
    return {**metadata(), 'version': ed.__version__, 'sage': SAGE_VERSION, 'fixture': name,
            'scope': 'eclib/torsion bound; rational Kummer-4 checks where labelled=true',
            'independence': 'reference uses no package square classes, pairing, or cocycle law',
            'verified': True, 'certified': result.certified, 'assumptions': result.assumptions,
            'complete': result.complete, 'eclib_selmer_dimension': eclib_dimension,
            'proved_rank': rank, 'sha_two_and_four_trivial': True,
            'expected_invariant_factors': expected, 'invariant_factors': G.invariant_factors,
            'classes': G.cardinality, 'labelled': labelled, 'rational_point_witnesses': point_records,
            'ordered_sums_checked': sums, 'pairwise_separations_checked': separations,
            'seconds': perf_counter()-start, 'installed_module': ed.__file__}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=['all', *FIXTURES], default='all')
    parser.add_argument('--require-installed', action='store_true')
    args = parser.parse_args()
    for name in FIXTURES if args.fixture == 'all' else (args.fixture,):
        print(json.dumps(audit_fixture(name, require_installed=args.require_installed), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
