"""Exact relative 4-descent over Q and independently testable arithmetic tools."""
from .geometry import FourCoverCandidate, construct_four_cover
from .etale import EtaleQuarticAlgebra
from .status import (CandidateDecision, CompletionStatus, DescentResult,
                     EnumerationStatus, EquivalenceStatus, LocalStatus, ProofStatus)
from .arithmetic import FiniteNormSpace, quartic_support
from .foundation import screen_norm_candidates
from .local import padic_solubility, real_solubility, everywhere_local_solubility
from .relative import relative_four_descent, RelativeDescentResult
from .two_cover import QuarticTwoCover
from .local_complete import padic_solubility_complete

__all__ = [
    "EtaleQuarticAlgebra", "FourCoverCandidate", "construct_four_cover", "CandidateDecision",
    "CompletionStatus", "DescentResult", "EnumerationStatus",
    "EquivalenceStatus", "LocalStatus", "ProofStatus", "FiniteNormSpace",
    "quartic_support", "screen_norm_candidates", "padic_solubility",
    "real_solubility", "everywhere_local_solubility",
    "relative_four_descent", "RelativeDescentResult", "QuarticTwoCover",
    "padic_solubility_complete",
]
__version__ = "0.3.0"
