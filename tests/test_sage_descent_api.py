"""Genuine-Sage curve-level integration; no simulated number-field results."""
import json
import os
from dataclasses import replace

import pytest

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import (
    DescentSession, DescentLimits, DescentStage, four_descent, two_selmer,
    cassels_tate_pairing, assemble_four_selmer_fibres, WeierstrassChange,
)
from elliptic_descent.fibre_presentations import fibre_representative
from elliptic_descent.status import CompletionStatus


@pytest.mark.parametrize('a,expected', [
    ((0, 0, 0, -1, 1), (4,)), ((0, 0, 0, -1, 0), (2, 2)),
    ((0, -1, 1, -10, -20), ()), ((0, -1, 1, -929, -10595), (2, 2)),
])
def test_curve_api_agrees_with_explicit_pipeline_and_reuses_witnesses(a, expected):
    session = DescentSession(a)
    R = session.run()
    assert R.finished and R.complete and R.certified and R.verify_evidence()
    assert R.completion_status == CompletionStatus.CERTIFIED
    assert R.invariant_factors == expected and R.curve.ainvs == a
    # Same exact arithmetic/basis, but independently reconstruct the lower-stage
    # pairing and fibres; this tests integration, not a second descent algorithm.
    B = cassels_tate_pairing(R.selmer)
    F = assemble_four_selmer_fibres(R.selmer, B)
    G = F.group()
    assert B.matrix == R.pairing.matrix and G.law == R.group.law
    assert R.cardinality == G.cardinality
    for v in R.group.elements:
        assert R.covering(v.key) is v.covering
        assert v.covering.target == R.curve
        assert 2*v == R.group.include_two_selmer(R.group.projection(v))
        for w in R.group.elements:
            assert (v+w).key == (G.element(v.key)+G.element(w.key)).key
    before = session.cache_info()
    again = session.run()
    assert again.group is R.group and again.fibres is R.fibres and again.selmer is R.selmer
    assert session.cache_info() == before and again.verify_evidence()
    assert json.loads(json.dumps(R.summary()))['certified'] is True
    # Cached wrappers hold maps, not a successful verifier's verdict.
    bad = replace(R.group, method='unverified-count-fit')
    assert not replace(R, group=bad).certified
    assert not replace(R, group=bad).verify_evidence()


def test_two_selmer_candidate_and_local_limits_resume_one_arithmetic_context():
    session = DescentSession((-1, 1))
    empty = session.run(limits=DescentLimits(two_selmer_candidates=0))
    assert not empty.finished and not empty.complete and empty.stop_reason == 'two_selmer_incomplete'
    assert empty.selmer.decisions == () and empty.pairing is None
    one = session.run(until='two_selmer', limits=DescentLimits(two_selmer_candidates=1))
    full = session.run(until='two_selmer')
    assert full.finished and not full.complete and not full.certified and full.verify_evidence()
    assert full.selmer.arithmetic is one.selmer.arithmetic is empty.selmer.arithmetic
    assert full.selmer.decisions[0] is one.selmer.decisions[0]
    assert session.cache_info()['arithmetic_contexts'] == 1
    assert session.cache_info()['two_selmer_candidates'] == full.selmer.total_candidates
    independent = two_selmer((-1, 1))
    assert full.selmer.dimension == independent.dimension
    assert full.selmer.known_soluble_labels == independent.known_soluble_labels
    assert all(d.model.delta == e.model.delta for d, e in zip(full.selmer.decisions, independent.decisions))
    for option in ('two_selmer_nodes', 'two_selmer_depth'):
        other = DescentSession((-1, 1))
        partial = other.run(limits=DescentLimits(**{option: 0}))
        assert not partial.complete and not partial.selmer.complete
        assert partial.selmer.unresolved
        done = other.run(until='two_selmer')
        assert done.selmer.certified and done.verify_evidence()
        assert done.selmer.arithmetic is partial.selmer.arithmetic
        assert other.cache_info()['arithmetic_contexts'] == 1


def test_pairing_resume_reuses_all_three_normalized_quartics():
    session = DescentSession((0, -1, 1, -929, -10595))
    partial = session.run(until='pairing', limits=DescentLimits(pairing_trials=0))
    assert partial.selmer.complete and not partial.pairing.complete and partial.fibres is None
    assert partial.stop_reason == 'pairing_incomplete' and not partial.verify_evidence()
    before = session.cache_info()
    full = session.run(until='pairing')
    assert full.finished and not full.complete and full.verify_evidence()
    assert full.selmer is partial.selmer and full.pairing.rank == 2
    assert session.cache_info()['quartic_presentations'] == before['quartic_presentations'] == 3
    assert session.cache_info()['quartic_normalizations'] == before['quartic_normalizations'] == 3
    assert session.cache_info()['pairing_entries'] == before['pairing_entries']+1
    assert session.run(until='pairing').pairing is full.pairing


def test_fibre_limits_keep_earlier_stages_and_reuse_completed_fibres():
    session = DescentSession((-1, 1))
    empty = session.run(limits=DescentLimits(fibres=0))
    assert empty.selmer.complete and empty.pairing.complete and not empty.fibres.complete
    assert empty.fibres.known_classes == () and empty.group is None
    assert not empty.complete and empty.stop_reason == 'fibres_incomplete'
    one = session.run(limits=DescentLimits(fibres=1))
    assert not one.complete and one.fibres.known_classes
    key = one.fibres.known_classes[0].key
    assert one.covering(key).target == one.curve
    done = session.run()
    assert done.certified and done.verify_evidence()
    assert one.fibres.fibres[0] is done.fibres.fibres[0]
    assert session.cache_info()['fibre_assemblies'] == len(done.fibres.fibres)
    assert done.selmer is empty.selmer and done.pairing is empty.pairing
    with pytest.raises(ValueError): empty.require_group()
    other = DescentSession((-1, 1))
    limited = other.run(limits=DescentLimits(relative_candidates=0))
    assert limited.stop_reason == 'fibres_incomplete' and not limited.fibres.complete
    assert other.run().certified


def test_field_degree_limits_retain_fibres_and_do_not_cache_a_false_success():
    session = DescentSession((-1, 0))
    partial = session.run(backend='galois', limits=DescentLimits(field_degree=0))
    assert partial.stop_reason == 'field_degree_limit' and partial.group is None
    assert partial.fibres.complete and partial.verify_evidence() and not partial.complete
    assert dict(partial.limit_detail)['degree'] > 0
    before = session.cache_info()
    done = session.run(backend='galois')
    assert done.certified and done.group.method == 'geometric_galois_cocycles'
    assert done.fibres is partial.fibres and done.selmer is partial.selmer
    after = session.cache_info()
    assert after['group_constructions'] == before['group_constructions']+1
    assert after['fibre_assemblies'] == before['fibre_assemblies']
    capped = session.run(backend='galois', limits=DescentLimits(field_degree=1))
    assert capped.stop_reason == 'field_degree_limit' and capped.group is None
    assert session.cache_info() == after
    assert session.run(backend='galois').group is done.group
    shortcut = session.run(backend='auto', limits=DescentLimits(field_degree=0))
    assert shortcut.certified and shortcut.group.law == done.group.law
    assert shortcut.group is not done.group


def test_conditional_parent_is_not_upgraded_by_unconditional_relative_arithmetic():
    S = two_selmer((-1, 1), proof=False)
    session = DescentSession.from_two_selmer(S, proof=True)
    R = session.run()
    assert R.complete and R.finished and R.verify_evidence()
    assert not R.certified and R.assumptions and R.selmer is S
    assert R.completion_status == CompletionStatus.CONDITIONAL
    assert session.cache_info().get('arithmetic_contexts', 0) == 0
    with pytest.raises(ValueError): R.require_group(certified=True)
    with pytest.raises(ValueError): DescentSession.from_two_selmer(two_selmer((-1, 1), max_candidates=0))
    # A proof=False flag alone is not a fabricated GRH assumption on Q factors.
    split = four_descent((-1, 0), proof=False, until='two_selmer')
    assert split.selmer.certified and not split.assumptions and not split.complete


def test_model_and_basis_identity_cached_presentations_and_alternate_routes():
    change = WeierstrassChange((-1, 1), (1, 1, -1, 2))
    R = four_descent(change.target, minimal=False)
    assert R.certified and R.curve == change.target and R.selmer.curves.working == change.target
    assert all(v.covering.target == change.target for v in R.group.elements)
    S = two_selmer((-1, 0)); b0, b1 = S.basis_labels
    T = S.with_basis((b0^b1, b0))
    session = DescentSession.from_two_selmer(T)
    C = session.quartic_cover(1)
    assert C is session.quartic_cover(1) and C.selmer is T
    with pytest.raises(ValueError): fibre_representative(S, 1, S.torsion_data(), initial_cover=C)
    G = session.run().group
    H = DescentSession.from_two_selmer(S).run().group
    mapping = {g.key: next(h for h in H.elements
                if H.projection(h) == T.transport_to(S, G.projection(g))) for g in G.elements}
    for g in G.elements:
        for k in G.elements:
            assert mapping[(g+k).key] == mapping[g.key]+mapping[k.key]
    # Separate route keys: asking for the reference zero-fibre descent must not
    # simply return a cached transported zero fibre from the optimized route.
    session = DescentSession((-1, 1))
    optimized = session.run(until='fibres')
    reference = session.run(until='fibres', reuse_zero_fibre=False, prune=False)
    assert reference.fibres.complete and reference.verify_evidence()
    assert reference.fibres.fibre(0) is not optimized.fibres.fibre(0)
    assert reference.fibres.cardinality == optimized.fibres.cardinality
