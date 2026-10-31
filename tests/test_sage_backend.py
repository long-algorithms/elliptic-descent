"""Real Sage integration tests; a dedicated CI job requires these, without mocks."""
import os
import pytest
import sympy as sp

if os.environ.get("ELLIPTIC_DESCENT_REQUIRE_SAGE") == "1":
    import sage.all as sage
else:
    sage = pytest.importorskip("sage.all", reason="requires SageMath runtime")

from elliptic_descent.backends.sage import (
    number_field_square_classes, pari_two_covers, sage_square_classes,
)
from elliptic_descent.foundation import screen_norm_candidates
from elliptic_descent.geometry import _norm, construct_four_cover
from elliptic_descent.status import ProofStatus

F = (1, -1, 0, 0, 1)


def test_quartic_basis_roundtrips_and_independent_norms():
    context = sage_square_classes(F, proof=True)
    space = context.space
    assert space.support == (2, 229)
    assert space.arithmetic_proof == ProofStatus.CERTIFIED
    K = context.field
    for i, g in enumerate(space.generators):
        assert context.coordinates(g) == 1 << i
        element = sum(K(str(c))*K.gen()**j for j, c in enumerate(g))
        assert sp.Rational(str(element.norm())) == _norm(g, F)
    batch = space.enumerate(limit=1024)
    assert batch.norm_fibre_exhausted
    assert len({context.coordinates(c.xi) for c in batch.candidates}) == batch.total_count
    for c in batch.candidates:
        assert c.label == context.coordinates(c.xi)
        assert _norm(c.xi, F) == c.norm_witness**2


def test_class_group_contribution_is_not_lost():
    # Sage's documented class-number-two example: 2 is not a unit but contributes.
    K = sage.QuadraticField(-5, "a")
    assert K.class_number(proof=True) == 2
    V, gens, from_V, to_V = number_field_square_classes(K, [], proof=True)
    assert V.dimension() == 2
    assert to_V(K(2)) != V.zero()
    assert any(abs(g.norm()) != 1 for g in gens)


def test_conditional_arithmetic_is_labelled():
    context = sage_square_classes(F, (2,), proof=False)
    assert context.space.arithmetic_proof == ProofStatus.CONDITIONAL
    assert "GRH" in context.space.enumerate(limit=1).assumptions[0]
    with pytest.raises(TypeError):
        number_field_square_classes(context.field, (2,), proof="false")


def test_nonmonic_field_does_not_change_geometric_quartic():
    f = (1, -1, 0, 0, 4)
    context = sage_square_classes(f, (2,), proof=True)
    assert context.space.f[4] == 4
    assert context.field.defining_polynomial().is_monic()
    for c in context.space.enumerate(limit=4).candidates:
        assert _norm(c.xi, f) == 4*c.norm_witness**2
        assert construct_four_cover(f, c.xi).check_map_identity()


def test_reducible_arithmetic_is_explicitly_unsupported():
    with pytest.raises(NotImplementedError):
        sage_square_classes((4, 0, -5, 0, 1))


def test_pari_basis_with_independently_verified_maps():
    covers = pari_two_covers((-25, 4))
    assert len(covers) == 2
    assert all(c.check_map() for c in covers)
    assert all(c.ainvs == (0, 0, 0, -25, 4) for c in covers)


def test_arithmetic_to_local_screening_pipeline():
    context = sage_square_classes(F, (2,), proof=True)
    result = screen_norm_candidates(context.space, limit=32, max_nodes=5000)
    assert result.arithmetic.norm_fibre_exhausted
    assert not result.relative_descent_complete
    assert len(result.candidates) == result.arithmetic.total_count
    for candidate in result.candidates:
        assert candidate.covering.check_map_identity()
        for report in candidate.local.reports:
            if report.certificate:
                assert report.verify(candidate.covering)
