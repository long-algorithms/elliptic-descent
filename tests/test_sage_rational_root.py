"""All rational-root patterns: complete descents vs independent eclib 2-Selmer counts."""
import os
from dataclasses import replace
from functools import lru_cache

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import relative_four_descent
from elliptic_descent.rational_root import RationalRootIdentification

FIXTURES = (
    ('cubic_root_zero', (0,1,-1,0,1), (3,1)),
    ('cubic_root_three', (-3,4,-1,-3,1), (3,1)),
    ('cubic_all_real', (0,1,-3,0,1), (3,1)),
    ('imaginary_quadratic', (-1,0,0,0,1), (2,1,1)),
    ('quadratic_two_roots', (2,0,-3,0,1), (2,1,1)),
    ('split_four_torsion', (4,0,-5,0,1), (1,1,1,1)),
    ('negative_cubic', (0,-1,1,0,-1), (3,1)),
    ('rational_split', (sp.Rational(1,4),0,sp.Rational(-5,4),0,1), (1,1,1,1)),
)


@lru_cache(maxsize=16)
def descend(f):
    return relative_four_descent(f, proof=True)


@pytest.mark.parametrize('name,f,kind', FIXTURES)
def test_complete_descent_and_independent_second_descent(name,f,kind):
    result = descend(f)
    assert result.factorization_type == kind
    assert result.complete and result.certified and not result.unresolved
    assert result.coverings  # A rational branch point rules out an empty Selmer set.
    assert result.orientation.multiplicity == 1
    assert len(result.coverings) == len(result.unoriented_coverings)
    assert result.verify_evidence(), name
    # Use the derivative-based cubic, not the canonical invariant formula,
    # as input to eclib. The identification checks the marked map as well.
    identification = RationalRootIdentification(result.source)
    assert identification.check_marked_map()
    E = sage.EllipticCurve([sage.QQ(str(c)) for c in identification.branch_cubic_ainvs])
    E = E.global_minimal_model()
    dimension = int(E.mwrank_curve(verbose=False).selmer_rank())
    assert len(result.coverings) == 2**dimension, (name, len(result.coverings), dimension)
    for D in result.coverings:
        witness = result.orientation.identification(D)
        assert witness.verify(D, replace(D, norm_witness=-D.norm_witness))
        assert D.check_construction()


def test_partial_modes_and_genuine_empty_result_are_distinct():
    f = (0,1,-1,0,1)
    partial = relative_four_descent(f, max_candidates=0)
    assert not partial.coverings and partial.total_model_classes > 0
    assert not partial.enumeration_complete and not partial.certified
    assert not partial.verify_evidence()
    for limits in ({'max_nodes':0}, {'max_depth':0}):
        bounded = relative_four_descent(f, **limits)
        assert bounded.enumeration_complete and bounded.unresolved
        assert not bounded.complete and not bounded.certified
        assert not bounded.verify_evidence()


def test_proof_mode_depends_on_arithmetic_not_just_user_flag():
    mixed = relative_four_descent((0,1,-1,0,1), proof=False)
    assert mixed.complete and mixed.assumptions and not mixed.certified
    assert mixed.verify_evidence()
    assert len(mixed.coverings) == len(descend((0,1,-1,0,1)).coverings)
    split = relative_four_descent((4,0,-5,0,1), proof=False)
    assert split.complete and split.certified and not split.assumptions
    assert len(split.coverings) == len(descend((4,0,-5,0,1)).coverings)


def test_orientation_tampering_cannot_certify_or_double_the_output():
    result = descend((2,0,-3,0,1))
    for multiplier in (None, (1,0,0,0), (0,0,0,0)):
        bad = replace(result, orientation=replace(result.orientation,reversing_multiplier=multiplier))
        assert not bad.certified and not bad.verify_evidence()
        with pytest.raises(ValueError):
            _ = bad.coverings


def test_rational_four_torsion_does_not_trigger_an_extra_quotient():
    result = descend((4,0,-5,0,1))
    E = sage.EllipticCurve([sage.QQ(str(c)) for c in result.source.ainvs])
    P = E(8,12)
    assert P.order() == 4 and 2*P == E(5,0)
    dimension = int(E.mwrank_curve(verbose=False).selmer_rank())
    assert len(result.coverings) == 2**dimension
    assert result.orientation.multiplicity == 1


def test_translation_and_nonmonic_square_rescaling_preserve_counts():
    f = (2,0,-3,0,1)
    t = sp.Symbol('t')
    shifted = sp.Poly(sum(c*(t+1)**i for i,c in enumerate(f)),t)
    transformed = (tuple(shifted.nth(i) for i in range(5)), tuple(sp.Rational(c,4) for c in f))
    for changed in transformed:
        result = relative_four_descent(changed)
        assert result.certified and result.verify_evidence()
        assert len(result.coverings) == len(descend(f).coverings)
        assert RationalRootIdentification(changed).check_marked_map()


def test_singular_and_degree_three_presentations_stay_rejected():
    for f in ((0,0,-1,0,1), (1,-1,0,1,0)):
        with pytest.raises(ValueError):
            relative_four_descent(f)
