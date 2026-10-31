from elliptic_descent.arithmetic import FiniteNormSpace
from elliptic_descent.foundation import screen_norm_candidates
from elliptic_descent.status import ProofStatus


def test_pipeline_keeps_scope_and_budget_explicit():
    space = FiniteNormSpace((1, -1, 0, 0, 1), ((0, -1, 0, 0),), (2, 229))
    result = screen_norm_candidates(space, limit=1, max_nodes=0, max_depth=0)
    assert not result.arithmetic.norm_fibre_exhausted
    assert not result.relative_descent_complete
    assert result.arithmetic.arithmetic_proof == ProofStatus.UNVERIFIED
    assert len(result.candidates) == 1
    full = screen_norm_candidates(space, limit=2, verify_geometry=True)
    assert full.arithmetic.norm_fibre_exhausted
    assert not full.relative_descent_complete
    for candidate in full.candidates:
        for report in candidate.local.reports:
            if report.certificate:
                assert report.verify(candidate.covering)
