"""API orchestration unit tests; stubs here are NOT mathematical verification.

The genuine Sage module exercises the same paths with actual descent data.
Existing exact-covering and portable-certificate suites remain enabled.
"""
from dataclasses import fields, replace
from types import SimpleNamespace

import pytest

import elliptic_descent.descent as api
from elliptic_descent import DescentLimits, DescentSession, DescentStage, FourDescentResult, FieldDegreeLimit
from elliptic_descent.status import CompletionStatus
from elliptic_descent.weierstrass import WeierstrassModel


def node(**kwargs):
    return SimpleNamespace(complete=True, assumptions=(), verify_evidence=lambda **kw: True, **kwargs)


def chain(conditional=False):
    S = node(dimension=1, curves=SimpleNamespace(original=WeierstrassModel((-1, 1))))
    B = node(selmer=S)
    F = node(selmer=S, pairing=B, certified=not conditional, known_classes=())
    G = node(fibres=F, cardinality=4, invariant_factors=(4,), galois_evidence=None)
    if conditional:
        for obj in (S, B, F, G):
            obj.assumptions = ('unproved arithmetic assumption for unit-test stub',)
    return S, B, F, G


def install_stages(monkeypatch, session, objects):
    S, B, F, G = objects
    calls = []
    def two(limits):
        calls.append('two_selmer'); session._selmer = S; return S
    def pairing(limits):
        calls.append('pairing'); session._pairing = B; return B
    def fibres(*args):
        calls.append('fibres'); return F
    def group(*args):
        calls.append('group'); return G
    for name, function in (('_two_selmer', two), ('_compute_pairing', pairing),
                            ('_compute_fibres', fibres), ('_compute_group', group)):
        monkeypatch.setattr(session, name, function)
    return calls


@pytest.mark.parametrize('name', [f.name for f in fields(DescentLimits)])
def test_all_budget_fields_are_exact_nonnegative_or_unbounded(name):
    assert getattr(DescentLimits(**{name: 0}), name) == 0
    assert getattr(DescentLimits(), name) is None
    for invalid in (-1, 0.0, True, '0'):
        with pytest.raises((ValueError, TypeError)):
            DescentLimits(**{name: invalid})
    L = DescentLimits(relative_candidates=0, relative_depth=2, relative_nodes=3)
    assert L.relative == dict(max_candidates=0, max_depth=2, max_nodes=3)


@pytest.mark.parametrize('until,index', [(s, i) for i, s in enumerate(DescentStage)])
def test_staged_completion_is_never_full_group_completion(monkeypatch, until, index):
    session = DescentSession((-1, 1))
    calls = install_stages(monkeypatch, session, chain())
    R = session.run(until=until)
    assert calls == [s.value for s in tuple(DescentStage)[:index+1]]
    assert R.stage == until and R.finished and R.verify_evidence()
    assert R.complete == R.certified == (until == DescentStage.GROUP)
    if index < 3:
        assert R.completion_status == CompletionStatus.INCOMPLETE
        with pytest.raises(ValueError): R.require_group()
        with pytest.raises(ValueError): _ = R.cardinality
    else:
        assert R.cardinality == 4 and R.invariant_factors == (4,)
        assert R.require_group(certified=True) is R.group
    assert session.last_result is R


@pytest.mark.parametrize('stage', range(3))
def test_incomplete_stage_stops_before_downstream_calls(monkeypatch, stage):
    session = DescentSession((-1, 1))
    objects = chain(); objects[stage].complete = False
    calls = install_stages(monkeypatch, session, objects)
    R = session.run()
    assert len(calls) == stage+1
    assert not R.complete and not R.finished and not R.certified and not R.verify_evidence()
    assert R.stop_reason == tuple(DescentStage)[stage].value+'_incomplete'
    assert R.group is None


def test_only_typed_resource_limit_is_returned_as_partial(monkeypatch):
    session = DescentSession((-1, 1))
    install_stages(monkeypatch, session, chain())
    def limited(*args): raise FieldDegreeLimit(16, 2)
    monkeypatch.setattr(session, '_compute_group', limited)
    R = session.run(limits=DescentLimits(field_degree=2))
    assert R.stop_reason == 'field_degree_limit'
    assert R.limit_detail == (('degree', 16), ('limit', 2))
    assert R.group is None and R.fibres.complete and R.verify_evidence()
    assert not R.complete and not R.certified
    for error in (RuntimeError('backend failure'), ArithmeticError('invalid mathematical evidence')):
        def fail(*args): raise error
        monkeypatch.setattr(session, '_compute_group', fail)
        with pytest.raises(type(error), match=str(error)): session.run()
        assert session.last_result.fibres is R.fibres
        assert session.last_result.group is None


def test_conditional_provenance_and_live_verification(monkeypatch):
    session = DescentSession((-1, 1))
    install_stages(monkeypatch, session, chain(conditional=True))
    R = session.run()
    assert R.complete and R.finished and R.verify_evidence() and not R.certified
    assert R.completion_status == CompletionStatus.CONDITIONAL
    assert R.assumptions
    with pytest.raises(ValueError): R.require_group(certified=True)
    S, B, F, G = chain()
    calls = []
    G.verify_evidence = lambda **kwargs: not calls.append(True)
    R = FourDescentResult(DescentStage.GROUP, DescentLimits(), S, B, F, G)
    assert R.certified and R.certified and len(calls) == 2
    G.verify_evidence = lambda **kwargs: False
    assert not R.certified and not R.verify_evidence()
    assert R.completion_status == CompletionStatus.INCOMPLETE
    assert not replace(R, pairing=node(selmer=S)).verify_evidence()
    assert not replace(R, fibres=None).complete


def test_on_demand_quartics_are_reused_only_within_the_bound_session():
    session = DescentSession((-1, 1))
    with pytest.raises(ValueError): session.quartic_cover(0)
    calls = []
    S = node()
    def coordinates(c):
        if c not in (0, 1): raise ValueError('outside basis')
        return c
    def create(c):
        calls.append(c); return object()
    S.ambient_label, S.quartic_cover = coordinates, create
    session._selmer = S
    C = session.quartic_cover(1)
    assert session.quartic_cover(1) is C and calls == [1]
    assert session.cache_info() == {'quartic_presentations': 1}
    other = DescentSession((-1, 1)); other._selmer = S
    assert other.quartic_cover(1) is not C
    with pytest.raises(ValueError): session.quartic_cover(2)
    with pytest.raises(TypeError): session.quartic_cover(0.0)
    info = session.cache_info(); info['quartic_presentations'] = 999
    assert session.cache_info()['quartic_presentations'] == 1


def test_pairing_cache_retries_unknowns_not_completed_entries():
    session = DescentSession((-1, 1)); S = node(dimension=3)
    session._selmer = S
    calls = []
    def entry(a, b, *, max_local_trials):
        calls.append((a, b, max_local_trials))
        return SimpleNamespace(selmer=S, left=a, right=b,
                               complete=(a, b) == (1, 2) or max_local_trials is None)
    session._builder = SimpleNamespace(entry=entry)
    limited = session._compute_pairing(DescentLimits(pairing_trials=0))
    assert not limited.complete and len(calls) == 3
    again = session._compute_pairing(DescentLimits(pairing_trials=0))
    assert not again.complete and len(calls) == 3
    full = session._compute_pairing(DescentLimits())
    assert full.complete and len(calls) == 5
    assert full.entries[0] is limited.entries[0]
    assert session._compute_pairing(DescentLimits(pairing_trials=0)) is full
    assert len(calls) == 5


def test_group_cache_respects_backend_identity_and_acceptance_degree(monkeypatch):
    session = DescentSession((-1, 1)); F = object(); calls = []
    def group(fibres, *, backend, max_field_degree):
        calls.append((fibres, backend, max_field_degree))
        return SimpleNamespace(galois_evidence=SimpleNamespace(torsion=SimpleNamespace(
            field=SimpleNamespace(degree=16))))
    monkeypatch.setattr(api, 'four_selmer_group', group)
    G = session._compute_group(F, 'galois', None)
    assert session._compute_group(F, 'galois', 16) is G and len(calls) == 1
    with pytest.raises(FieldDegreeLimit): session._compute_group(F, 'galois', 2)
    assert session._compute_group(F, 'auto', None) is not G and len(calls) == 2
    assert session._compute_group(object(), 'galois', None) is not G and len(calls) == 3


def test_input_errors_are_detected_before_backend_work(monkeypatch):
    for option in ('proof', 'minimal'):
        with pytest.raises(TypeError): DescentSession((-1, 1), **{option: 1})
    session = DescentSession((-1, 1))
    monkeypatch.setattr(session, '_two_selmer', lambda *args: pytest.fail('backend should not be called'))
    for options in ({'until': 'rank'}, {'limits': {}}, {'prune': 0},
                    {'reuse_zero_fibre': 1}, {'backend': 'guessed'}):
        with pytest.raises((ValueError, TypeError)): session.run(**options)
    assert session.last_result is None
    with pytest.raises(ValueError): DescentSession.from_two_selmer(object())
    with pytest.raises(TypeError): FieldDegreeLimit(16.0, 1)
    assert isinstance(FieldDegreeLimit(16, 1), RuntimeError)
    R = FourDescentResult(DescentStage.GROUP, DescentLimits())
    assert R.stage is None and not R.complete and not R.verify_evidence()
    assert R.summary()['cardinality'] is None
    with pytest.raises(ValueError): _ = R.curve
    with pytest.raises(ValueError): R.covering((0, 0))
    with pytest.raises(TypeError): R.require_group(certified=1)


def test_class_wrapper_caches_are_data_only_and_replacements_get_new_wrappers():
    from elliptic_descent.four_selmer_fibres import AssembledFibre
    from elliptic_descent import MarkedQuarticCover, construct_four_cover
    marking = MarkedQuarticCover.canonical((1, -1, 0, 0, 1))
    model = construct_four_cover(marking.source.f, (0, -1, 0, 0))
    lifts = SimpleNamespace(
        representative=SimpleNamespace(coordinates=0, marking=marking),
        marked_result=SimpleNamespace(relative_result=SimpleNamespace(coverings=(model,))))
    fibre = AssembledFibre(lifts, True, (), ((0,),))
    classes = fibre.classes
    assert fibre.classes is classes and classes[0].covering.verify()
    replacement = replace(fibre, orbits=((0,),))
    assert replacement.classes is not classes
    assert replacement.classes[0].covering == classes[0].covering
    empty = replace(fibre, orbits=())
    assert empty.classes == ()  # A previous cache cannot conceal changed accounting.
