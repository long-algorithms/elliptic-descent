"""v0.2 arithmetic-to-local-screening pipeline, deliberately not four_descent()."""
from dataclasses import dataclass

from .arithmetic import CandidateBatch, FiniteNormSpace
from .geometry import construct_four_cover
from .local import LocalScreen, everywhere_local_solubility
from .status import LocalStatus


@dataclass(frozen=True)
class ScreenedCandidate:
    label: int
    covering: object
    local: LocalScreen


@dataclass(frozen=True)
class FoundationResult:
    arithmetic: CandidateBatch
    candidates: tuple[ScreenedCandidate, ...]

    @property
    def relative_descent_complete(self):
        # No rational-scalar/orientation quotient or covering equivalence here.
        return False

    @property
    def unresolved(self):
        return tuple(c for c in self.candidates if c.local.status == LocalStatus.UNKNOWN)

    @property
    def retained(self):
        # Unknown candidates MUST NOT be discarded.
        return tuple(c for c in self.candidates if c.local.status != LocalStatus.INSOLUBLE)


def screen_norm_candidates(space: FiniteNormSpace, *, limit=128,
                           max_depth=6, max_nodes=100_000, verify_geometry=False):
    """Enumerate a bounded norm fibre and certify local decisions where possible.

    Unexpected arithmetic exceptions propagate; they never become empty output.
    Each local budget is per candidate and finite place, not for the whole run.
    """
    batch = space.enumerate(limit=limit)
    records = []
    for candidate in batch.candidates:
        covering = construct_four_cover(space.f, candidate.xi,
                                        norm_witness=candidate.norm_witness)
        if verify_geometry and not all((covering.check_algebra(), covering.check_map_identity(),
                                        covering.check_smoothness(), covering.check_map_basepoint_free())):
            raise ArithmeticError("candidate geometry verification failed")
        local = everywhere_local_solubility(covering, max_depth=max_depth, max_nodes=max_nodes)
        records.append(ScreenedCandidate(candidate.label, covering, local))
    return FoundationResult(batch, tuple(records))
