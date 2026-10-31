"""Every vector, not only a basis: labelled quartics with checked original markings."""
import os
from dataclasses import replace
from functools import lru_cache
from itertools import product

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage=pytest.importorskip('sage.all',reason='requires SageMath runtime')

from elliptic_descent import two_selmer, WeierstrassModel, WeierstrassChange
from elliptic_descent.backends.sage_conics import rational_conic_point
from elliptic_descent.selmer_quartic import labelled_quartic
from elliptic_descent.geometry import U
from elliptic_descent.two_cover import X, Z

FIXTURES=(
    (0,0,0,-1,1), (0,0,0,-1,0), (0,0,0,1,0), (0,0,0,4,0),
    (0,-5,0,-16,80), (0,-1,1,-10,-20), (0,-1,1,-929,-10595),
)


@lru_cache(maxsize=16)
def selmer(a):
    return two_selmer(a)


@lru_cache(maxsize=32)
def cover(a,c):
    return selmer(a).quartic_cover(c)


@pytest.mark.parametrize('a',FIXTURES)
def test_all_classes_including_zero_basis_and_sums(a):
    S=selmer(a)
    assert S.certified
    E=sage.EllipticCurve(list(a)).global_minimal_model()
    assert S.cardinality==2**int(E.mwrank_curve(verbose=False).selmer_rank())
    labels=set()
    for c in range(S.cardinality):
        C=cover(a,c)
        assert C.certified and C.verify()
        assert C.target==S.curves.original and C.marking.target==S.curves.original
        assert C.coordinates==c and C.ambient_label==S.ambient_label(c)
        assert S.coordinates(C.delta)==c
        assert C.quartic.f[4] and C.quartic.check_map()
        assert C.marking.translation is None
        assert C.ambient_label not in labels
        labels.add(C.ambient_label)
    assert labels==set(S.known_soluble_labels)
    for left,right in product(range(S.cardinality),repeat=2):
        C=cover(a,S.add(left,right))
        delta=S.arithmetic.algebra.multiply(cover(a,left).delta,cover(a,right).delta)
        assert S.coordinates(delta)==C.coordinates


def test_sha_classes_do_not_require_a_rational_point_on_the_genus_one_curve():
    a=(0,-1,1,-929,-10595)
    S=selmer(a)
    assert S.dimension==2
    for c in (1,2,3):
        C=cover(a,c)
        assert C.certified and C.verify()
        assert C.presentation.parametrization.verify(C.presentation.auxiliary.q[2])
        assert C.target==S.curves.original
    # In particular the sum of the basis vectors has its own marked quartic.
    assert cover(a,3).ambient_label==S.basis_labels[0]^S.basis_labels[1]


def test_alternate_selmer_bases_and_conic_parameters_preserve_actual_labels():
    a=(0,-5,0,-16,80)
    S=selmer(a)
    b0,b1=S.basis_labels
    other=S.with_basis((b0^b1,b0))
    for c in range(S.cardinality):
        base=cover(a,c)
        d=other.coordinates(base.delta)
        C=other.quartic_cover(d,parameter_change=((1,1),(1,2)))
        assert C.certified and C.ambient_label==base.ambient_label and C.delta==base.delta
        assert C.target==base.target
        assert S.transport_to(other,c)==d
        point=tuple(p.subs({X:1,Z:2},simultaneous=True) for p in base.presentation.parametrization.forms)
        alternative=S.quartic_cover(c,conic_point=point)
        assert alternative.verify() and alternative.ambient_label==base.ambient_label


def test_original_nonminimal_model_is_preserved_under_class_transport():
    E=WeierstrassModel((-1,1))
    change=WeierstrassChange(E,(sp.Rational(1,2),1,-1,3))
    S=selmer(E.ainvs)
    moved=two_selmer(change.target)
    for c in range(S.cardinality):
        d=S.transport_to(moved,c,change)
        C=moved.quartic_cover(d)
        assert C.certified and C.target==change.target
        assert C.marking.change.target==change.target
        assert moved.transport_to(S,C.coordinates,change.inverse())==c
    # A different working presentation also produces the same arithmetic class.
    integral_change=WeierstrassChange(E,(1,1,-1,2))
    retained=two_selmer(integral_change.target,minimal=False)
    C=retained.quartic_cover(S.transport_to(retained,1,integral_change))
    assert C.certified and C.target==integral_change.target


def test_tampered_labels_models_and_partial_inputs_are_rejected():
    a=(0,0,0,-1,1)
    S=selmer(a); C=cover(a,0)
    assert not replace(C,coordinates=1).verify()
    assert not replace(C,ambient_label=1<<S.arithmetic.dimension).verify()
    assert not replace(C,presentation=cover(a,1).presentation).verify()
    for c in (-1,S.cardinality):
        with pytest.raises(ValueError):S.quartic_cover(c)
    with pytest.raises(TypeError):S.quartic_cover(1.0)
    partial=two_selmer(a,max_candidates=0)
    with pytest.raises(ValueError):partial.quartic_cover(0)
    conditional=two_selmer(a,proof=False)
    result=conditional.quartic_cover(0)
    assert result.verify() and not result.certified and result.assumptions
    rational=two_selmer((0,0,0,-1,0),proof=False).quartic_cover(0)
    assert rational.certified and not rational.assumptions


def test_relative_integration_keeps_the_label_marking_and_assumptions():
    for c in (0,1):
        C=cover((0,0,0,-1,1),c)
        result=C.descend(proof=True)
        assert result.complete and result.certified and result.verify_evidence()
        assert result.target==C.target
        assert result.representative.ambient_label==C.ambient_label
        assert len(result.coverings)==selmer((0,0,0,-1,1)).cardinality
        assert all(D.target==C.target and D.verify() for D in result.coverings)
    partial=C.descend(max_candidates=0)
    assert not partial.complete and not partial.certified
    assert not partial.verify_evidence()
    conditional=two_selmer((0,0,0,-1,1),proof=False).quartic_cover(0).descend(proof=True)
    assert conditional.complete and conditional.verify_evidence()
    assert not conditional.certified and conditional.assumptions


def test_sage_conic_solver_returns_checkable_points_and_fails_explicitly():
    u,v,w=U[:3]
    q=u*u+v*v-2*w*w
    point=rational_conic_point(q)
    assert any(point) and q.subs(dict(zip(U[:3],point)))==0
    with pytest.raises(ArithmeticError):rational_conic_point(u*u+v*v+w*w)
