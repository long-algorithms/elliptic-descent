"""Reference geometry for explicit 4-descent on elliptic curves."""
from .geometry import FourCoverCandidate, construct_four_cover
from .status import (CandidateDecision, CompletionStatus, DescentResult,
                     EnumerationStatus, EquivalenceStatus, LocalStatus, ProofStatus)

__all__ = [
    "FourCoverCandidate", "construct_four_cover", "CandidateDecision",
    "CompletionStatus", "DescentResult", "EnumerationStatus",
    "EquivalenceStatus", "LocalStatus", "ProofStatus",
]
__version__ = "0.1.0"
