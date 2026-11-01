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
from .rational_root import RationalRootIdentification
from .local_complete import padic_solubility_complete
from .weierstrass import WeierstrassModel, WeierstrassChange
from .cubic_conversion import CubicToQuartic, IdentityTwoCover
from .marked_cover import MarkedQuarticCover, MarkedFourCover, MarkedRelativeResult
from .two_selmer import two_selmer, TwoSelmerResult
from .selmer_quartic import LabelledQuarticCover, QuarticPresentation
from .cassels_tate import cassels_tate_pairing, pairing_entry, CasselsTatePairing
from .four_selmer_fibres import assemble_four_selmer_fibres, FourSelmerFibres

from .four_selmer_group import four_selmer_group, FourSelmerGroup, FourSelmerElement

__all__ = [
    "EtaleQuarticAlgebra", "FourCoverCandidate", "construct_four_cover", "CandidateDecision",
    "CompletionStatus", "DescentResult", "EnumerationStatus",
    "EquivalenceStatus", "LocalStatus", "ProofStatus", "FiniteNormSpace",
    "quartic_support", "screen_norm_candidates", "padic_solubility",
    "real_solubility", "everywhere_local_solubility",
    "relative_four_descent", "RelativeDescentResult", "QuarticTwoCover",
    "padic_solubility_complete", "RationalRootIdentification",
    "WeierstrassModel", "WeierstrassChange", "CubicToQuartic", "IdentityTwoCover",
    "MarkedQuarticCover", "MarkedFourCover", "MarkedRelativeResult",
    "two_selmer", "TwoSelmerResult", "LabelledQuarticCover", "QuarticPresentation",
    "cassels_tate_pairing", "pairing_entry", "CasselsTatePairing",
    "assemble_four_selmer_fibres", "FourSelmerFibres",
]
__all__ += ["four_selmer_group", "FourSelmerGroup", "FourSelmerElement"]
__version__ = "0.4.0"

from .descent import four_descent, DescentSession, DescentLimits, DescentStage, FourDescentResult
from .limits import FieldDegreeLimit

__all__ += ["four_descent", "DescentSession", "DescentLimits", "DescentStage",
            "FourDescentResult", "FieldDegreeLimit"]
