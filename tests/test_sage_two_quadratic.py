"""Complete 2+2 descents and independent eclib cardinality checks.

The known rational points identify C with its Jacobian, so relative 2-coverings
are counted by the 2-Selmer group. eclib is a different arithmetic algorithm,
not another invocation of our factorwise Sage norm-fibre computation.
"""
import os
from dataclasses import replace
from functools import lru_cache
from itertools import product

import pytest
import sympy as sp

if os.environ.get("ELLIPTIC_DESCENT_REQUIRE_SAGE") == "1":
    import sage.all as sage
else:
    sage = pytest.importorskip("sage.all", reason="requires SageMath runtime")

from elliptic_descent import EtaleQuarticAlgebra, QuarticTwoCover, relative_four_descent
from elliptic_descent.backends.sage_etale import sage_etale_square_classes
from elliptic_descent.equivalence import OrientationCertificate
from elliptic_descent.quotient import quotient_fibre
from elliptic_descent.status import LocalStatus

# Counts below are local replay regressions, independently compared to eclib
# whenever a rational point is supplied. No live Magma run is claimed.
FIXTURES = (
    ("two_real_fields", (6, 0, -5, 0, 1), (1, 1, 0), 1),
    ("isomorphic_fields", (16, 0, -10, 0, 1), (1, 1, 0), 2),
    ("nonmonic_positive", (12, 0, -10, 0, 2), (1, 2, 1), 2),
    ("two_imaginary_fields", (2, 0, 3, 0, 1), (1, 1, 0), 1),
    ("mixed_signature", (-6, 0, 1, 0, 1), (1, 1, 0), 1),
    ("affine_positive", (6, 0, 8, 0, 2), (1, 4, 1), 2),
    ("real_obstruction", (-2, 0, -3, 0, -1), None, 0),
    ("two_adic_obstruction", (6, 0, 9, 0, 3), None, 0),
)


@lru_cache(maxsize=16)
def descend(f):
    return relative_four_descent(f, proof=True)


@pytest.mark.parametrize("name,f,point,expected", FIXTURES)
def test_complete_two_quadratic_descent(name, f, point, expected):
    result = descend(f)
    assert result.factorization_type == (2, 2)
    assert result.enumeration_complete and result.complete and result.certified
    assert not result.unresolved and not result.assumptions
    assert len(result.unoriented_coverings) == expected, name
    assert len(result.coverings) == 2*expected
    assert result.verify_evidence(), name
    for D in result.coverings:
        assert D.check_construction()
    if point is not None:
        # evaluate() verifies the point on the exact input quartic, including
        # the nonmonic model; merely knowing the Jacobian rank is insufficient.
        result.source.evaluate(point)
        E = sage.EllipticCurve([sage.QQ(str(a)) for a in result.source.ainvs])
        independent_dimension = int(E.mwrank_curve().selmer_rank())
        assert len(result.coverings) == 2**independent_dimension, (name, independent_dimension)


def test_two_adic_empty_source_has_an_independent_mod_eight_obstruction():
    f = (6, 0, 9, 0, 3)  # 3*(t^2+1)*(t^2+2)
    squares = {y*y % 8 for y in range(8)}
    assert not any((x % 2 or z % 2) and
                   sum(c*x**i*z**(4-i) for i, c in enumerate(f)) % 8 in squares
                   for x, z in product(range(8), repeat=2))
    result = descend(f)
    # Exercise actual local rejection, not just an empty arithmetic norm fibre.
    assert result.total_model_classes > 0
    assert result.certified and not result.coverings
    assert all(c.local.status == LocalStatus.INSOLUBLE for c in result.candidates)


def test_diagonal_scalar_image_and_isomorphic_factor_dependencies():
    for f in ((6, 0, -5, 0, 1), (16, 0, -10, 0, 1)):
        context = sage_etale_square_classes(f, proof=True)
        scalar_labels = [context.coordinates(context.algebra.scalar(c))
                         for c in (-1, *context.space.support)]
        quotient, rows = quotient_fibre(context.norm_fibre(), scalar_labels)
        result = descend(f)
        assert result.scalar_basis == rows
        assert result.total_model_classes == quotient.cardinality
        assert context.norm_fibre().cardinality == quotient.cardinality*2**len(rows)
        # Both 2 and its inverse are squares in both copies of Q(sqrt(2)).
        if f[0] == 16:
            assert context.coordinates(context.algebra.scalar(2)) == 0
            assert len(rows) < len(context.space.support)+1


def test_two_quadratic_partial_and_conditional_modes():
    f = (6, 0, -5, 0, 1)
    partial = relative_four_descent(f, max_candidates=0)
    assert partial.total_model_classes > 0 and not partial.enumeration_complete
    assert not partial.complete and not partial.certified and not partial.verify_evidence()
    for limits in ({"max_depth": 0}, {"max_nodes": 0}):
        bounded = relative_four_descent(f, **limits)
        assert bounded.enumeration_complete and bounded.unresolved
        assert not bounded.complete and not bounded.certified
        assert not bounded.verify_evidence()
    conditional = relative_four_descent(f, proof=False)
    assert conditional.complete and conditional.assumptions and not conditional.certified
    assert len(conditional.coverings) == len(descend(f).coverings)
    assert conditional.verify_evidence()  # Local/geometry audit does not certify class groups.


def test_affine_change_and_square_scaling_preserve_counts():
    f = (6, 0, -5, 0, 1)
    t = sp.Symbol("t")
    shifted = sp.Poly(sum(c*(t+1)**i for i, c in enumerate(f)), t)
    inputs = (tuple(shifted.nth(i) for i in range(5)), tuple(4*c for c in f))
    expected = len(descend(f).coverings)
    for changed in inputs:
        result = relative_four_descent(QuarticTwoCover(changed))
        assert result.certified and len(result.coverings) == expected
        assert result.verify_evidence()


def test_orientation_evidence_is_bound_to_the_actual_factorization():
    result = descend((6, 0, -5, 0, 1))
    bad = replace(result, orientation=OrientationCertificate(((1, -1, 0, 0, 1),)))
    assert not bad.certified and not bad.verify_evidence()
    with pytest.raises(ValueError):
        _ = bad.coverings
    # Even an isomorphic abstract product is not the same labelled polynomial.
    other = replace(result, orientation=descend((16, 0, -10, 0, 1)).orientation)
    assert not other.certified and not other.verify_evidence()


def test_positive_two_quadratic_models_pass_generic_geometry_audits():
    result = descend((6, 0, -5, 0, 1))
    for D in result.coverings:
        assert D.check_map_identity() and D.check_smoothness() and D.check_map_basepoint_free()
