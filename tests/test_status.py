"""A result cannot silently upgrade a partial calculation to complete."""
from dataclasses import replace

import pytest

from elliptic_descent.status import (
    CandidateDecision, CompletionStatus, DescentResult, EnumerationStatus,
    EquivalenceStatus, LocalStatus, ProofStatus,
)


def exhaustive(*candidates, proof=ProofStatus.CERTIFIED, assumptions=()):
    return DescentResult(
        candidates=tuple(candidates),
        enumeration=EnumerationStatus.EXHAUSTIVE,
        enumeration_proof=proof,
        enumeration_evidence="finite-support square-class enumeration",
        equivalence=EquivalenceStatus.DEDUPLICATED,
        equivalence_proof=proof,
        equivalence_evidence="verified identification of covering classes",
        assumptions=tuple(assumptions),
    )


def decided(label="class-1", status=LocalStatus.SOLUBLE, proof=ProofStatus.CERTIFIED):
    return CandidateDecision(label, status, proof, "local certificates")


def test_empty_partial_is_not_empty_descent():
    r = DescentResult()
    assert not r.complete and not r.certified
    assert not r.locally_soluble
    assert r.completion_status == CompletionStatus.INCOMPLETE


def test_empty_exhaustive_certified_is_valid_empty_descent():
    r = exhaustive()
    assert r.complete and r.certified
    assert r.locally_soluble == ()
    assert r.completion_status == CompletionStatus.CERTIFIED


def test_unknown_local_decision_blocks_completeness():
    r = exhaustive(decided(), CandidateDecision("undecided"))
    assert not r.complete
    assert [x.label for x in r.unresolved] == ["undecided"]
    assert [x.label for x in r.locally_soluble] == ["class-1"]


def test_conditional_proof_is_complete_but_not_certified():
    r = exhaustive(decided(proof=ProofStatus.CONDITIONAL))
    assert r.complete and not r.certified
    assert r.completion_status == CompletionStatus.CONDITIONAL
    assert exhaustive(proof=ProofStatus.CONDITIONAL).complete
    assert not exhaustive(proof=ProofStatus.CONDITIONAL).certified
    assert not exhaustive(assumptions=("GRH for class groups",)).certified


def test_missing_equivalence_blocks_completeness():
    r = replace(exhaustive(), equivalence=EquivalenceStatus.NOT_CHECKED,
                equivalence_proof=ProofStatus.UNVERIFIED, equivalence_evidence="")
    assert not r.complete


def test_locally_insoluble_decision_is_not_reported_as_soluble():
    r = exhaustive(decided(status=LocalStatus.INSOLUBLE))
    assert r.complete and r.certified
    assert r.locally_soluble == ()


def test_missing_enumeration_proof_blocks_completeness():
    r = replace(exhaustive(), enumeration_proof=ProofStatus.UNVERIFIED,
                enumeration_evidence="")
    assert not r.complete


def test_duplicate_candidate_labels_are_rejected():
    with pytest.raises(ValueError, match="unique"):
        exhaustive(decided("same"), decided("same"))


@pytest.mark.parametrize("kwargs", [
    {"label": ""},
    {"label": "unknown", "proof": ProofStatus.CERTIFIED},
    {"label": "unknown", "evidence": "should be empty"},
    {"label": "decided", "local_status": LocalStatus.SOLUBLE,
     "proof": ProofStatus.CERTIFIED},
])
def test_invalid_candidate_claims(kwargs):
    with pytest.raises(ValueError):
        CandidateDecision(**kwargs)


def test_inconsistent_exhaustiveness_and_equivalence_rejected():
    with pytest.raises(ValueError, match="partial"):
        DescentResult(enumeration=EnumerationStatus.PARTIAL,
                      enumeration_proof=ProofStatus.CERTIFIED)
    with pytest.raises(ValueError, match="equivalence"):
        DescentResult(equivalence_proof=ProofStatus.CERTIFIED)
    with pytest.raises(ValueError, match="evidence"):
        DescentResult(enumeration=EnumerationStatus.EXHAUSTIVE,
                      enumeration_proof=ProofStatus.CERTIFIED)
