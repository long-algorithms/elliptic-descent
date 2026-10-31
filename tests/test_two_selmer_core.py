"""Exact cubic algebra, Kummer, local evidence and coordinate tests without Sage."""
from dataclasses import replace
from itertools import product
from random import Random

import pytest
import sympy as sp

from elliptic_descent.cubic_arithmetic import (CubicEtaleAlgebra, T, cubic_norm, kummer_algebra,
    kummer_value, rational_two_torsion, rational_halves_of_two_torsion)
from elliptic_descent.cubic_local import CubicCover, cubic_local_solubility, cubic_real_solubility
from elliptic_descent.geometry import U
from elliptic_descent.status import LocalStatus
from elliptic_descent.two_selmer import F2Labels
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange


@pytest.mark.parametrize('a', [(-1,1),(-1,0),(4,0),(0,-5,0,-16,80),
                              (0,0,1,-7,6),(1,2,3,4,5),(0,-1,0,5,-5)])
def test_cubic_algebra_norms_and_unit_cover(a):
    E = WeierstrassModel(a)
    A = kummer_algebra(E)
    assert A.verify() and sum(A.degrees) == 3
    rng = Random(20261031)
    for _ in range(10):
        u = tuple(sp.Rational(rng.randrange(-4,5),rng.randrange(1,5)) for j in range(3))
        v = tuple(sp.Rational(rng.randrange(-4,5),rng.randrange(1,5)) for j in range(3))
        assert A.from_components(A.to_components(u)) == u
        assert A.norm(u) == cubic_norm(u,A.f)
        assert A.norm(A.multiply(u,v)) == A.norm(u)*A.norm(v)
    D = CubicCover(E,A.one)
    assert D.verify() and D.evaluate((1,0,0,0)) is None
    local = cubic_local_solubility(D)
    assert local.status == LocalStatus.SOLUBLE and local.verify(D)


def test_kummer_at_every_rational_two_torsion_point():
    E = WeierstrassModel((-1,0)); A=kummer_algebra(E)
    assert len(rational_two_torsion(E)) == 4
    expected = {(1,0):(2,1,2), (0,0):(-1,-1,1), (-1,0):(-2,-1,2)}
    # Factors ordered theta-1, theta, theta+1.
    for P, values in expected.items():
        d = kummer_value(E,P)
        assert A.to_components(d) == tuple((c,) for c in values)
        assert CubicCover(E,d).verify()
        assert not rational_halves_of_two_torsion(E,P)
    assert kummer_value(E,None) == A.one
    E=WeierstrassModel((4,0))
    assert set(rational_halves_of_two_torsion(E,(0,0))) == {(2,4),(2,-4)}
    E=WeierstrassModel((0,-5,0,-16,80))
    halves = rational_halves_of_two_torsion(E,(5,0))
    assert set(halves) == {(2,6),(2,-6),(8,12),(8,-12)}
    assert all(E.double(P) == (5,0) for P in halves)


def test_kummer_homomorphism_by_independent_chord_identity():
    base=WeierstrassModel((-1,1))
    for change in (WeierstrassChange(base),WeierstrassChange(base,(-2,1,2,-3))):
        E=change.target; A=kummer_algebra(E)
        P=change.apply((0,1)); Q=change.apply(base.double((0,1))); R=E.add(P,Q)
        a1,_,a3,_,_=E.ainvs
        vP=P[1]+(a1*P[0]+a3)/2
        vQ=Q[1]+(a1*Q[0]+a3)/2
        slope=(vQ-vP)/(Q[0]-P[0]); intercept=vP-slope*P[0]
        line=(intercept,slope,0)
        triple=A.multiply(A.multiply(kummer_value(E,P),kummer_value(E,Q)),kummer_value(E,R))
        assert triple == A.multiply(line,line)


def test_cubic_real_sign_has_no_spurious_global_sign_quotient():
    E=WeierstrassModel((-1,0)); A=kummer_algebra(E)
    # Factor order is +1,0,-1; real-root order is the reverse.
    for signs in product((-1,1),repeat=3):
        if sp.prod(signs) != 1:
            continue
        delta=A.from_components(tuple((s,) for s in signs))
        D=CubicCover(E,delta); report=cubic_real_solubility(D)
        assert report.verify(D)
        expected = signs in ((1,1,1),(-1,-1,1))
        assert (report.status == LocalStatus.SOLUBLE) == expected
        assert not replace(report,status=LocalStatus.INSOLUBLE if expected else LocalStatus.SOLUBLE).verify(D)
    # Replacing all signs by their negatives is not an equivalence in the cubic norm kernel.
    with pytest.raises(ValueError): CubicCover(E,(-1,0,0))


def test_local_budget_failure_is_not_insolubility():
    D=CubicCover((-1,0),(1,0,0))
    for limits in ({'max_nodes':0},{'max_depth':0}):
        report=cubic_local_solubility(D,**limits)
        assert report.status == LocalStatus.UNKNOWN and not report.verify(D)
    A=D.algebra
    bad=CubicCover(D.curve,A.from_components(((1,),(-1,),(-1,))))
    assert cubic_real_solubility(bad).status == LocalStatus.INSOLUBLE


def test_coordinate_roundtrips_and_alternate_bases():
    rng=Random(5032)
    for n in range(8):
        basis=tuple(1<<i for i in range(n))
        for _ in range(20):
            B=list(basis); rng.shuffle(B)
            if n>1:
                i,j=rng.sample(range(n),2); B[i]^=B[j]
            labels=F2Labels(n,tuple(B))
            for c in range(1<<n):
                assert labels.coordinates(labels.ambient(c)) == c
    L=F2Labels(4,(3,12))
    assert L.coordinates(15)==3
    for v in (-1,1,16):
        with pytest.raises(ValueError): L.coordinates(v)
    for B in ((0,), (1,1), (1,2,3), (16,)):
        with pytest.raises(ValueError): F2Labels(4,B)
    with pytest.raises(TypeError): L.ambient(1.0)
    with pytest.raises(ValueError): L.ambient(4)


def test_input_and_corrupted_model_rejections():
    for f in ((1,0,0,0),(1,0,0,2),(0,0,0,1)):
        with pytest.raises(ValueError): CubicEtaleAlgebra(f)
    with pytest.raises(TypeError): CubicEtaleAlgebra((1.0,-1,0,1))
    A=kummer_algebra(WeierstrassModel((-1,0)))
    with pytest.raises(ValueError): A.from_components(((1,),))
    with pytest.raises(ValueError): A.embed_unit(1,(0,))
    with pytest.raises(ValueError): CubicCover((-1,0),(0,0,0))
    with pytest.raises(ValueError): kummer_value((-1,1),(1,2))
    with pytest.raises(ValueError): rational_halves_of_two_torsion((-1,1),(0,1))
    D=CubicCover((-1,1),(1,0,0))
    with pytest.raises(ValueError): D.evaluate((0,0,0,0))
    with pytest.raises(ValueError): D.evaluate((1,1,1,1))
    local=cubic_local_solubility(D)
    assert not local.verify(CubicCover((-2,1),(1,0,0)))


def test_two_torsion_and_halves_transport_under_nonminimal_changes():
    E=WeierstrassModel((0,-5,0,-16,80))
    C=WeierstrassChange(E,(sp.Rational(1,3),2,-1,4))
    assert set(rational_two_torsion(C.target)) == {C.apply(P) for P in rational_two_torsion(E)}
    for P in rational_two_torsion(E):
        assert set(rational_halves_of_two_torsion(C.target,C.apply(P))) == {
            C.apply(Q) for Q in rational_halves_of_two_torsion(E,P)}


def test_kummer_values_lift_their_rational_points_on_the_auxiliary_models():
    for a, points in [((-1,1), ((0,1),(0,-1),(1,1))),
                      ((-1,0), ((1,0),(0,0),(-1,0))),
                      ((0,-5,0,-16,80), ((8,12),(5,0),(4,0)))]:
        E=WeierstrassModel(a); A=kummer_algebra(E)
        for P in points:
            delta=kummer_value(E,P)
            D=CubicCover(E,delta)
            parts=[(sp.S.One,)+(sp.S.Zero,)*(degree-1) for degree in A.degrees]
            for i,g in enumerate(A.factors):
                if len(g)==2 and P[0]==-g[0]:
                    parts[i]=(sp.S.Zero,)
            u=A.from_components(parts)
            assert D.evaluate((*u,1)) in (E.point(P), E.negate(P))


def test_split_selmer_classes_by_direct_rational_enumeration():
    from elliptic_descent.arithmetic import square_class
    E=WeierstrassModel((-1,0)); A=kummer_algebra(E)
    reps=tuple(map(sp.Rational,(-2,-1,1,2)))
    byclass={square_class(a,(2,)):a for a in reps}
    key=lambda d: tuple(square_class(c[0],(2,)) for c in A.to_components(d))
    accepted=set()
    for a,b in product(reps,repeat=2):
        c=byclass[square_class(a*b,(2,))]
        d=A.from_components(((a,),(b,),(c,)))
        model=CubicCover(E,d); evidence=cubic_local_solubility(model)
        assert evidence.status != LocalStatus.UNKNOWN and evidence.verify(model)
        if evidence.status == LocalStatus.SOLUBLE: accepted.add(key(d))
    assert accepted == {key(kummer_value(E,P)) for P in rational_two_torsion(E)}
    assert len(accepted)==4
