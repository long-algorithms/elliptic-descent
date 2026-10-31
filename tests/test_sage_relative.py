"""End-to-end relative descent, including published genuinely empty fibres."""
import os
from dataclasses import replace

import pytest

if os.environ.get("ELLIPTIC_DESCENT_REQUIRE_SAGE") == "1":
    import sage.all as sage
else:
    sage = pytest.importorskip("sage.all", reason="requires SageMath runtime")

from elliptic_descent.relative import relative_four_descent
from elliptic_descent.status import LocalStatus

# Magma handbook H133E34--36. Coefficients are ascending, not display order.
# Positive example has one model modulo target Y -> -Y, hence two fixed-target
# oriented coverings. 571a1 has no locally soluble lifts above any of these.
FIXTURES = [
    ("handbook-positive", (172, -120, 20, 12, -7), 1),
    ("571a1-first", (-64, 164, -52, -68, -11), 0),
    ("571a1-second", (-7, -68, -142, 112, -19), 0),
    ("571a1-third", (-3, 52, -232, 60, -4), 0),
]


def description(result):
    return {"total": result.total_model_classes, "support": result.support,
            "scalar_rank": len(result.scalar_basis), "norm_dimension": result.norm_dimension,
            "models": len(result.unoriented_coverings), "complete": result.complete,
            "timings": result.timings,
            "local": [(c.label, c.local.status.value, c.local.real.signs,
                       [(r.place, r.status.value, r.nodes, r.depth) for r in c.local.finite])
                      for c in result.candidates]}


def test_unbounded_relative_positive_and_truncation():
    result = relative_four_descent((1, -1, 0, 0, 1))
    assert result.certified, description(result)
    assert result.coverings and len(result.coverings) == 2*len(result.unoriented_coverings)
    assert result.verify_evidence(), description(result)
    partial = relative_four_descent((1,-1,0,0,1), max_candidates=0)
    assert not partial.enumeration_complete and not partial.complete and not partial.certified
    bounded = relative_four_descent((1,-1,0,0,1), max_depth=0, max_nodes=0)
    assert bounded.enumeration_complete and bounded.unresolved and not bounded.certified


@pytest.mark.parametrize("name,f,expected_models", FIXTURES)
def test_published_relative_descents(name, f, expected_models):
    result = relative_four_descent(f, max_nodes=200_000, max_depth=24)
    print(name, description(result), flush=True)
    assert result.certified, description(result)
    assert len(result.unoriented_coverings) == expected_models, description(result)
    assert len(result.coverings) == 2*expected_models
    for c in result.candidates:
        assert c.local.real.verify(c.model)
        for report in c.local.finite:
            # Cheap positive checks are independent of the search.
            if report.status == LocalStatus.SOLUBLE:
                assert report.verify(c.model)


def test_relative_conditional_status():
    result = relative_four_descent((1,-1,0,0,1), proof=False)
    assert result.complete and not result.certified and result.assumptions


def test_relative_reducible_is_not_an_empty_descent():
    with pytest.raises(NotImplementedError):
        relative_four_descent((4,0,-5,0,1))
