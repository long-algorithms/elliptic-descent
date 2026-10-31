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
    # Includes independent negative-certificate checks, not just point witnesses.
    assert result.verify_evidence(max_nodes=1_000_000), description(result)


def test_relative_conditional_status():
    result = relative_four_descent((1,-1,0,0,1), proof=False)
    assert result.complete and not result.certified and result.assumptions


def test_relative_singular_is_not_an_empty_descent():
    with pytest.raises(ValueError):
        relative_four_descent((1,0,-2,0,1))


def test_large_prime_certification_has_no_64_bit_cutoff():
    from elliptic_descent.arithmetic import prime, prime_support
    p = 2**89 - 1
    assert prime(p) == p
    assert prime_support((p*p,)) == (p,)
    with pytest.raises(ValueError):
        prime(3*p)


def test_irreducible_quartic_with_rational_two_torsion_jacobian():
    # Q(theta), theta^4=-1, contains several quadratic subfields. Both -1 and
    # 2 are squares in the field, so blindly subtracting |S|+1 dimensions is wrong.
    result=relative_four_descent((1,0,0,0,1))
    assert result.certified and len(result.scalar_basis)==0
    assert result.total_model_classes==4
    assert len(result.unoriented_coverings)==2 and len(result.coverings)==4
    assert result.source.ainvs==(0,0,0,-4,0)
    assert result.verify_evidence()


def test_square_rescaling_preserves_the_relative_answer():
    base=relative_four_descent((1,-1,0,0,1))
    scaled=relative_four_descent((4,-4,0,0,4))
    assert base.certified and scaled.certified
    assert len(base.coverings)==len(scaled.coverings)
    assert scaled.verify_evidence()
