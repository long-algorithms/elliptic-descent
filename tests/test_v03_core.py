from dataclasses import replace
from itertools import product
from random import Random
from types import SimpleNamespace

import pytest
import sympy as sp

from elliptic_descent.arithmetic import affine_f2
from elliptic_descent.geometry import U, construct_four_cover
from elliptic_descent.local import QuadraticSystem, padic_solubility
from elliptic_descent.local_complete import (
    _Budget, _initial_solutions, _linear_solutions, padic_solubility_complete,
    complete_local_screen,
)
from elliptic_descent.quotient import row_basis, normal_form, quotient_fibre
from elliptic_descent.real import real_squareclass_solubility
from elliptic_descent.status import LocalStatus
from elliptic_descent.two_cover import QuarticTwoCover


def test_quotient_against_exhaustive_cosets():
    rng = Random(8814)
    for n in range(1, 8):
        for _ in range(20):
            columns = [rng.randrange(8) for j in range(n)]
            fibre = affine_f2(columns, rng.randrange(8))
            if fibre.particular is None:
                continue
            scalars = tuple(b for b in fibre.basis if rng.randrange(2))
            quotient, rows = quotient_fibre(fibre, scalars)
            expected = {normal_form(v, rows) for v in fibre.elements()}
            assert set(quotient.elements()) == expected
            assert fibre.cardinality == quotient.cardinality*(1 << len(rows))
            assert len(list(quotient.elements())) == len(expected)
    empty, _ = quotient_fibre(affine_f2([0], 1), (1,))
    assert empty.cardinality == 0
    with pytest.raises(ValueError):
        quotient_fibre(affine_f2([1], 0), (1,))


@pytest.mark.parametrize('p', [2, 3, 5, 7])
def test_linear_lifts_exhaustively(p):
    rng = Random(p)
    for _ in range(25):
        matrix = [[rng.randrange(p) for j in range(3)] for i in range(2)]
        rhs = [rng.randrange(p) for i in range(2)]
        actual = set(_linear_solutions(matrix, rhs, p))
        expected = {x for x in product(range(p), repeat=3)
                    if all(sum(a*b for a,b in zip(row,x)) % p == r for row,r in zip(matrix,rhs))}
        assert actual == expected


@pytest.mark.parametrize('p', [2, 3, 5, 7])
def test_finite_field_elimination_exhaustively(p):
    rng = Random(100+p)
    for _ in range(20):
        equations = tuple(sum(rng.randrange(-3,4)*U[i]*U[j]
                              for i in range(4) for j in range(i,4)) for k in range(2))
        if any(q == 0 for q in equations):
            continue
        system = QuadraticSystem.from_quadrics(equations)
        actual = set(_initial_solutions(system, p, _Budget(None)))
        expected = set()
        for point in product(range(p), repeat=4):
            if not any(point):
                continue
            pivot = next(i for i,v in enumerate(point) if v)
            inverse = pow(point[pivot], -1, p)
            canonical = tuple(v*inverse % p for v in point)
            if all(v % p == 0 for v in system.values(point)):
                expected.add((pivot,canonical))
        assert actual == expected


def test_complete_local_and_budgets():
    u,v,w,z=U
    Q=(u*u+v*v+w*w+z*z,u*u+2*v*v+3*w*w+4*z*z)
    for p in (2,3,5,7):
        old = padic_solubility(Q,p)
        new = padic_solubility_complete(Q,p)
        assert new.status == old.status and new.verify(Q)
    for kwargs in ({'max_nodes':0}, {'max_depth':0}):
        assert padic_solubility_complete(Q,2,**kwargs).status == LocalStatus.UNKNOWN
    singular=(u*u-17*v*v,w*w-17*z*z)
    with pytest.raises(NotImplementedError):
        padic_solubility_complete(singular,2)
    r=padic_solubility_complete(singular,2,max_depth=6)
    assert r.status==LocalStatus.SOLUBLE and r.verify(singular)


def test_real_sign_certificates_and_tampering():
    # Root signs +,-,+,- cannot arise as a threshold, even up to common sign.
    f=(4,0,-5,0,1)
    t=sp.Symbol('t')
    polynomial=sp.interpolate([(-2,1),(-1,-1),(1,1),(2,-1)],t)
    xi=tuple(sp.Poly(polynomial,t).nth(i) for i in range(4))
    D=construct_four_cover(f,xi)
    r=real_squareclass_solubility(D)
    assert r.status==LocalStatus.INSOLUBLE and r.verify(D)
    assert not replace(r,status=LocalStatus.SOLUBLE).verify(D)
    assert not replace(r,signs=(1,1,1,1)).verify(D)
    for f,xi in [((1,-1,0,0,1),(0,-1,0,0)), ((-1,-1,0,0,1),(0,-1,0,0)),
                 ((4,0,-5,0,1),(1,0,0,0))]:
        model=SimpleNamespace(f=f,xi=xi)
        r=real_squareclass_solubility(model)
        assert r.status==LocalStatus.SOLUBLE and r.verify(model)


def test_two_cover_and_orientations():
    f=(1,-1,0,0,1)
    C=QuarticTwoCover(f)
    assert C.ainvs==(0,0,0,-4,1) and C.check_map()
    D=construct_four_cover(f,(0,-1,0,0))
    opposite=replace(D,norm_witness=-D.norm_witness)
    assert D.quadrics==opposite.quadrics
    assert opposite.check_map_identity()
    image=D.evaluate((1,0,0,0))[1]
    opposite_image=opposite.evaluate((1,0,0,0))[1]
    P,minusP=C.evaluate(image),C.evaluate(opposite_image)
    assert P[0]==minusP[0] and P[1]==-minusP[1] and P[1]!=0
    screen=complete_local_screen(D)
    assert screen.status==LocalStatus.SOLUBLE and screen.verify(D)


@pytest.mark.parametrize('f', [(1,-1,0,0,1), (172,-120,20,12,-7),
                               (-64,164,-52,-68,-11), (-7,-68,-142,112,-19)])
def test_canonical_map_identities(f):
    C=QuarticTwoCover(f)
    assert C.check_map()
