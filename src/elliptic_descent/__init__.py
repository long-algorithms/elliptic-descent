"""Exact arithmetic and local foundations for elliptic-curve 4-descent."""
from .geometry import FourCoverCandidate, construct_four_cover
from .status import (CandidateDecision, CompletionStatus, DescentResult,
                     EnumerationStatus, EquivalenceStatus, LocalStatus, ProofStatus)
from .arithmetic import FiniteNormSpace, quartic_support
from .foundation import screen_norm_candidates
from .local import padic_solubility, real_solubility, everywhere_local_solubility

__all__ = [
    "FourCoverCandidate", "construct_four_cover", "CandidateDecision",
    "CompletionStatus", "DescentResult", "EnumerationStatus",
    "EquivalenceStatus", "LocalStatus", "ProofStatus", "FiniteNormSpace",
    "quartic_support", "screen_norm_candidates", "padic_solubility",
    "real_solubility", "everywhere_local_solubility",
]
__version__ = "0.2.0"
