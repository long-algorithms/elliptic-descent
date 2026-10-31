"""Auditable *status bookkeeping* for future descent routines.

This module does not calculate Selmer groups or validate mathematical proofs.
ProofStatus values are claims to be assigned by an implementing backend only
when it has established the corresponding mathematical obligations.

The distinction between an empty exhaustive list and an empty partial list is
intentional.  A complete conditional computation is not a certified one.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Tuple


class ProofStatus(str, Enum):
    UNVERIFIED = "unverified"
    CONDITIONAL = "conditional"
    CERTIFIED = "certified"


class EnumerationStatus(str, Enum):
    NOT_STARTED = "not_started"
    PARTIAL = "partial"
    EXHAUSTIVE = "exhaustive"


class LocalStatus(str, Enum):
    UNKNOWN = "unknown"
    SOLUBLE = "soluble"
    INSOLUBLE = "insoluble"


class EquivalenceStatus(str, Enum):
    NOT_CHECKED = "not_checked"
    DEDUPLICATED = "deduplicated"


class CompletionStatus(str, Enum):
    INCOMPLETE = "incomplete"
    CONDITIONAL = "conditional"
    CERTIFIED = "certified"


@dataclass(frozen=True)
class CandidateDecision:
    """Labelled candidate with a local verdict and provenance label.

    The label is an identifier in an *enumerated candidate space*, not proof
    that a covering has been constructed or verified.  'proof' is only an
    annotation; the eventual backend is responsible for producing evidence.
    """

    label: str
    local_status: LocalStatus = LocalStatus.UNKNOWN
    proof: ProofStatus = ProofStatus.UNVERIFIED
    evidence: str = ""

    def __post_init__(self):
        if not self.label:
            raise ValueError("candidate label cannot be empty")
        if self.local_status == LocalStatus.UNKNOWN:
            if self.proof != ProofStatus.UNVERIFIED or self.evidence:
                raise ValueError("an unknown local status cannot have a proof")
        elif self.proof != ProofStatus.UNVERIFIED and not self.evidence:
            raise ValueError("a decided candidate must describe its evidence")


@dataclass(frozen=True)
class DescentResult:
    """A future per-quartic descent's finite decision ledger.

    Every enumerated candidate, including locally insoluble and unresolved
    candidates, must be present in ``candidates``.  An exhaustive enumeration
    requires a proof-labelled enumeration and explicit quotient-equivalence
    check.  This class only checks bookkeeping, not the underlying proofs.
    No descent algorithm produces this object in v0.1.x.
    """

    candidates: Tuple[CandidateDecision, ...] = ()
    enumeration: EnumerationStatus = EnumerationStatus.NOT_STARTED
    enumeration_proof: ProofStatus = ProofStatus.UNVERIFIED
    enumeration_evidence: str = ""
    equivalence: EquivalenceStatus = EquivalenceStatus.NOT_CHECKED
    equivalence_proof: ProofStatus = ProofStatus.UNVERIFIED
    equivalence_evidence: str = ""
    assumptions: Tuple[str, ...] = ()

    def __post_init__(self):
        labels = [c.label for c in self.candidates]
        if len(labels) != len(set(labels)):
            raise ValueError("candidate labels must be unique")
        if self.enumeration != EnumerationStatus.EXHAUSTIVE and self.enumeration_proof != ProofStatus.UNVERIFIED:
            raise ValueError("a partial enumeration cannot claim exhaustive proof")
        if self.enumeration_proof != ProofStatus.UNVERIFIED and not self.enumeration_evidence:
            raise ValueError("enumeration proof must identify its evidence")
        if self.equivalence != EquivalenceStatus.DEDUPLICATED and self.equivalence_proof != ProofStatus.UNVERIFIED:
            raise ValueError("unchecked equivalence cannot claim proof")
        if self.equivalence_proof != ProofStatus.UNVERIFIED and not self.equivalence_evidence:
            raise ValueError("equivalence proof must identify its evidence")
        if any(not s for s in self.assumptions):
            raise ValueError("empty assumption description")

    @property
    def unresolved(self) -> Tuple[CandidateDecision, ...]:
        return tuple(c for c in self.candidates if c.local_status == LocalStatus.UNKNOWN or c.proof == ProofStatus.UNVERIFIED)

    @property
    def locally_soluble(self) -> Tuple[CandidateDecision, ...]:
        return tuple(c for c in self.candidates if c.local_status == LocalStatus.SOLUBLE and c.proof != ProofStatus.UNVERIFIED)

    @property
    def completion_status(self) -> CompletionStatus:
        if (self.enumeration != EnumerationStatus.EXHAUSTIVE
                or self.enumeration_proof == ProofStatus.UNVERIFIED
                or self.equivalence != EquivalenceStatus.DEDUPLICATED
                or self.equivalence_proof == ProofStatus.UNVERIFIED
                or self.unresolved):
            return CompletionStatus.INCOMPLETE
        if (self.assumptions
                or self.enumeration_proof != ProofStatus.CERTIFIED
                or self.equivalence_proof != ProofStatus.CERTIFIED
                or any(c.proof != ProofStatus.CERTIFIED for c in self.candidates)):
            return CompletionStatus.CONDITIONAL
        return CompletionStatus.CERTIFIED

    @property
    def complete(self) -> bool:
        return self.completion_status != CompletionStatus.INCOMPLETE

    @property
    def certified(self) -> bool:
        return self.completion_status == CompletionStatus.CERTIFIED
