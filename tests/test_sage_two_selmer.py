"""Full labelled 2-Selmer data: genuine Sage arithmetic plus independent eclib checks."""
import os
from dataclasses import replace
from functools import lru_cache
from itertools import product

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent.two_selmer import two_selmer, F2Labels
from elliptic_descent.cubic_arithmetic import kummer_algebra, kummer_value
from elliptic_descent.backends.sage_two_selmer import cubic_square_classes
from elliptic_descent.status import ProofStatus
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange

FIXTURES = (
    ('rank_one', (0,0,0,-1,1)),
    ('split_two_torsion', (0,0,0,-1,0)),
    ('cyclic_two_torsion', (0,0,0,1,0)),
    ('cyclic_four_torsion', (0,0,0,4,0)),
    ('full_two_and_four', (0,-5,0,-16,80)),
    ('zero_selmer', (0,-1,1,-10,-20)),
    ('nontrivial_sha', (0,-1,1,-929,-10595)),
)


@lru_cache(maxsize=16)
def compute(a):
    return two_selmer(a)


@pytest.mark.parametrize('name,a', FIXTURES)
def test_complete_labelled_space_matches_independent_eclib(name,a):
    r = compute(a)
    assert r.enumeration_complete and r.complete and r.certified, name
    assert not r.unresolved and not r.assumptions
    assert r.verify_evidence(), name
    E = sage.EllipticCurve(list(a)).global_minimal_model()
    independent_dimension = int(E.mwrank_curve(verbose=False).selmer_rank())
    assert r.dimension == independent_dimension, (name,r.dimension,independent_dimension)
    assert r.cardinality == len(r.known_soluble_labels) == 2**independent_dimension
    assert r.arithmetic.verify_basis()
    assert {d.label for d in r.decisions} == set(r.candidate_space.elements())
    for c in range(r.cardinality):
        delta = r.representative(c)
        assert r.coordinates(delta) == c
        assert r.arithmetic.norm_label(r.ambient_label(c)) == 0
        assert r.add(c,c) == 0 and r.add(c,0) == c
    for a,b in product(range(r.cardinality),repeat=2):
        actual = r.arithmetic.algebra.multiply(r.representative(a),r.representative(b))
        assert r.coordinates(actual) == r.add(a,b)


def test_rank_zero_nontrivial_sha_is_not_a_point_generated_basis():
    r=compute((0,-1,1,-929,-10595))
    assert r.dimension == 2 and r.cardinality == 4
    assert r.torsion_data().image_rank == 0
    # Published 571a1 fixture: rank zero, no rational 2-torsion, nontrivial Selmer.
    # No point-search failure is used to prove either Selmer membership or rank.
    assert r.torsion_data().kernel == (None,)


@pytest.mark.parametrize('a', [(0,0,0,-1,0),(0,0,0,1,0),(0,0,0,4,0),(0,-5,0,-16,80)])
def test_torsion_kummer_images_and_halves_against_sage(a):
    r=compute(a); data=r.torsion_data()
    assert data.verify(r)
    E=sage.EllipticCurve(list(a))
    for entry in data.entries:
        P=E(0) if entry.point is None else E(*[sage.QQ(str(c)) for c in entry.point])
        halves = P.division_points(2)
        actual={None if Q.is_zero() else tuple(sp.Rational(str(c)) for c in Q.xy()) for Q in halves}
        assert set(entry.rational_halves) == actual
        assert (entry.coordinates == 0) == bool(actual)
    assert len(data.kernel)*2**data.image_rank == len(data.entries)
    if a == (0,0,0,4,0):
        assert r.kummer((0,0)) == 0 and r.kummer((2,4)) != 0
        assert data.image_rank == 0 and len(data.four_torsion) == 4
    if a == (0,-5,0,-16,80):
        assert len(data.entries) == 4 and len(data.kernel) == 2 and len(data.four_torsion) == 8


def test_kummer_group_homomorphism_on_nonzero_a1_a3_and_nonminimal_models():
    base=WeierstrassModel((-1,1))
    change=WeierstrassChange(base,(sp.Rational(1,2),1,-1,3))
    r=two_selmer(change.target)
    P=change.apply((0,1)); E=change.target
    points=[None]
    for _ in range(7): points.append(E.add(points[-1],P))
    for A,B in product(points,repeat=2):
        assert r.kummer(E.add(A,B)) == r.kummer(A)^r.kummer(B)
    assert r.kummer(E.double(P)) == 0
    for P in points:
        assert r.kummer(E.negate(P)) == r.kummer(P)
    old=compute(base.ainvs)
    for c in range(old.cardinality):
        d=old.transport_to(r,c,change)
        assert r.transport_to(old,d,change.inverse()) == c
    assert old.transport_to(r,old.kummer((0,1)),change) == r.kummer(change.apply((0,1)))


def test_basis_changes_preserve_actual_classes_and_torsion_images():
    r=compute((0,-5,0,-16,80))
    labels=list(r.basis_labels)
    labels.reverse(); labels[0]^=labels[1]
    alternate=r.with_basis(labels)
    assert alternate.certified and alternate.verify_evidence()
    assert alternate.dimension == r.dimension
    for c in range(r.cardinality):
        d=alternate.coordinates(r.representative(c))
        assert r.coordinates(alternate.representative(d)) == c
        assert r.transport_to(alternate,c) == d
    for entry in r.torsion_data().entries:
        assert alternate.kummer(entry.point) == alternate.coordinates_of_label(r.ambient_label(entry.coordinates))
    for invalid in ((0,), tuple(labels+labels[:1]), ()):
        with pytest.raises(ValueError): r.with_basis(invalid)


def test_proof_modes_and_budget_limits_do_not_claim_a_complete_basis():
    a=(0,0,0,-1,1)
    for limits in ({'max_candidates':0},{'max_nodes':0},{'max_depth':0}):
        r=two_selmer(a,**limits)
        assert not r.complete and not r.certified and not r.verify_evidence()
        with pytest.raises(ValueError): _=r.dimension
        with pytest.raises(ValueError): _=r.basis_labels
        assert r.kummer_ambient(None) == 0
    conditional=two_selmer(a,proof=False)
    assert conditional.complete and not conditional.certified and conditional.assumptions
    assert conditional.verify_evidence()
    rational=two_selmer((0,0,0,-1,0),proof=False)
    assert rational.complete and rational.certified and not rational.assumptions
    for invalid in (0,1,None,'false'):
        with pytest.raises(TypeError): two_selmer(a,proof=invalid)


def test_corrupt_labels_and_membership_fail_closed():
    r=compute((0,0,0,-1,1))
    duplicate=replace(r,decisions=(r.decisions[0],)*len(r.decisions))
    assert not duplicate.enumeration_complete and not duplicate.certified
    assert not duplicate.verify_evidence()
    first=replace(r.decisions[0],label=1<<r.arithmetic.dimension)
    bad=replace(r,decisions=(first,*r.decisions[1:]))
    assert not bad.complete
    for d in r.decisions:
        if d.label not in r.known_soluble_labels:
            with pytest.raises(ValueError): r.coordinates(d.model.delta)
            break
    with pytest.raises(ValueError): r.kummer((0,0))
    with pytest.raises(ValueError): r.add(r.cardinality,0)
    with pytest.raises(ValueError): r.transport_to(compute((0,0,0,-1,0)),0)


def test_factorwise_norm_kernel_is_not_quotiented_by_rational_scalars():
    for a in ((-1,1),(1,0),(-1,0)):
        arithmetic=cubic_square_classes(kummer_algebra(WeierstrassModel(a)))
        for c in (-1,2):
            label=arithmetic.coordinates((c,0,0))
            assert arithmetic.norm_label(label) != 0
        assert arithmetic.coordinates(arithmetic.representative(0)) == 0
        for i,g in enumerate(arithmetic.generators):
            assert arithmetic.coordinates(g) == 1<<i


def test_explicit_nonminimal_working_model_and_natural_algebra_transport():
    E=WeierstrassModel((-1,1))
    change=WeierstrassChange(E,(1,1,-1,2))
    normalized=two_selmer(change.target)
    retained=two_selmer(change.target,minimal=False)
    assert retained.curves.original == retained.curves.working == change.target
    assert retained.certified and retained.verify_evidence()
    assert retained.dimension == normalized.dimension
    for c in range(normalized.cardinality):
        d=normalized.transport_to(retained,c)
        assert retained.transport_to(normalized,d) == c
    P=change.apply((0,1))
    assert normalized.transport_to(retained,normalized.kummer(P)) == retained.kummer(P)
