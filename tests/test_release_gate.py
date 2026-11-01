"""The publication gate must fail closed; these are metadata/security unit tests.

No fixture below substitutes for a genuine Sage or independent mathematical run.
"""
import copy
import importlib.util
import io
from math import prod
from pathlib import Path
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import release_gate as gate
import publish_release as publisher
sys.path.pop(0)
SHA = 'a'*40


def runs():
    return [dict(id=i+1, run_attempt=1, path='.github/workflows/'+p, name=name,
                 head_sha=SHA, head_branch='main', event='push', status='completed', conclusion='success')
            for i, (p, name) in enumerate(gate.WORKFLOWS.items())]


def independent():
    data = []
    for name in gate.INDEPENDENT:
        n, labelled = prod(gate.GROUPS[name]), name in gate.LABELLED
        data.append(dict(version=gate.VERSION, revision=SHA, fixture=name,
                         verified=True, certified=True, complete=True, assumptions=[],
                         sha_two_and_four_trivial=True, classes=n,
                         expected_invariant_factors=gate.GROUPS[name], invariant_factors=gate.GROUPS[name],
                         labelled=labelled, ordered_sums_checked=n*n if labelled else 0,
                         pairwise_separations_checked=n*n if labelled else 0,
                         rational_point_witnesses=[{'key': [0, i], 'point': [1,0,0,0]} for i in range(n)] if labelled else [],
                         installed_module='/opt/sage/site-packages/elliptic_descent/__init__.py'))
    return data


def test_release_requires_every_exact_head_workflow():
    r = runs()
    assert set(gate.selected_runs(r, SHA)) == set(gate.WORKFLOWS)
    assert gate.selected_runs(r[:-1], SHA) is None
    r[0]['status'] = 'in_progress'
    assert gate.selected_runs(r, SHA) is None


@pytest.mark.parametrize('key,value', [('head_sha', 'b'*40), ('head_branch','other'),
    ('event','pull_request'), ('path','.github/workflows/fake.yml'), ('name','Fake tests')])
def test_unrelated_run_cannot_satisfy_release(key, value):
    r = runs()
    r[0][key] = value
    assert gate.selected_runs(r, SHA) is None


@pytest.mark.parametrize('outcome', ['failure','cancelled','skipped','timed_out',None])
def test_failed_latest_attempt_cannot_be_hidden_by_earlier_success(outcome):
    r = runs()
    newer = {**r[0], 'run_attempt': 2, 'conclusion': outcome}
    with pytest.raises(ValueError):
        gate.selected_runs([newer, *r], SHA)
    newer.update(conclusion='success')
    assert gate.selected_runs([newer, *r], SHA)['tests.yml']['run_attempt'] == 2


def test_independent_inventory_has_all_named_classes_and_sums():
    summary = gate.independent_records(independent(), gate.VERSION, SHA)
    assert summary == dict(curves=6, labelled_curves=4, labelled_classes=14, ordered_sums=52)


@pytest.mark.parametrize('key,value', [('revision','b'*40), ('version','0.4.0'),
    ('verified',False), ('certified',False), ('complete',False), ('assumptions',['GRH']),
    ('sha_two_and_four_trivial',False), ('classes',999), ('ordered_sums_checked',0),
    ('pairwise_separations_checked',0), ('labelled',False), ('rational_point_witnesses',[]),
    ('installed_module','/checkout/src/elliptic_descent/__init__.py')])
def test_bad_independent_record_blocks_publication(key,value):
    data = independent()
    row = next(r for r in data if r['fixture'] == 'rank_one')
    row[key] = value
    with pytest.raises(ValueError):
        gate.independent_records(data, gate.VERSION, SHA)


def test_missing_duplicate_and_changed_scope_independent_records_fail():
    data = independent()
    for bad in (data[:-1], data+[data[0]], [{**r, 'fixture':'unknown'} for r in data]):
        with pytest.raises(ValueError):
            gate.independent_records(bad, gate.VERSION, SHA)
    data = independent()
    row = next(r for r in data if r['fixture'] == 'rank_one')
    row['rational_point_witnesses'][1]['key'] = row['rational_point_witnesses'][0]['key']
    with pytest.raises(ValueError):
        gate.independent_records(data, gate.VERSION, SHA)


def junit_fixture(tmp_path):
    source, reports = tmp_path/'source', tmp_path/'reports'
    (source/'tests').mkdir(parents=True)
    reports.mkdir()
    for group in gate.SAGE_GROUPS:
        name = group.replace('-', '_')
        (source/'tests'/f'test_{name}.py').write_text(f'def test_{name}():\n    pass\n')
        (reports/f'sage-{group}.xml').write_text(
            f'<testsuites><testsuite><testcase classname="tests.test_{name}" name="test_{name}"/></testsuite></testsuites>')
    return source, reports


def test_junit_inventory_checks_source_functions_and_duplicates(tmp_path):
    source, reports = junit_fixture(tmp_path)
    assert gate.sage_reports(reports, source)['passed'] == len(gate.SAGE_GROUPS)
    (source/'tests'/'test_missing.py').write_text('def test_missing():\n    pass\n')
    with pytest.raises(ValueError):
        gate.sage_reports(reports, source)
    (source/'tests'/'test_missing.py').unlink()
    p = reports/'sage-core.xml'
    s = p.read_text()
    p.write_text(s.replace('/>', '><skipped/></testcase>'))
    with pytest.raises(ValueError):
        gate.sage_reports(reports, source)
    p.write_text(s.replace('</testsuite>', '<testcase classname="tests.test_core" name="test_core"/></testsuite>'))
    with pytest.raises(ValueError):
        gate.sage_reports(reports, source)


@pytest.mark.parametrize('name', ['../escape','/absolute','nested/../../escape','bad\\escape'])
def test_artifact_traversal_is_rejected(tmp_path, name):
    b = io.BytesIO()
    with zipfile.ZipFile(b, 'w') as z:
        z.writestr(name, b'bad')
    with pytest.raises(ValueError):
        publisher.extract_archive(b.getvalue(), tmp_path/'out')


def test_artifact_symlink_is_rejected_and_regular_bytes_extract(tmp_path):
    b = io.BytesIO()
    with zipfile.ZipFile(b, 'w') as z:
        info = zipfile.ZipInfo('link')
        info.external_attr = 0o120777 << 16
        z.writestr(info, '../outside')
    with pytest.raises(ValueError):
        publisher.extract_archive(b.getvalue(), tmp_path/'bad')
    b = io.BytesIO()
    with zipfile.ZipFile(b, 'w') as z:
        z.writestr('data/result.json', '{}')
    publisher.extract_archive(b.getvalue(), tmp_path/'good')
    assert (tmp_path/'good/data/result.json').read_text() == '{}'


def test_annotated_tags_resolve_without_rewriting_them(monkeypatch):
    def fake_api(path, optional=False):
        return {'object': {'type':'tag','sha':'b'*40}} if '/ref/' in path else {'object':{'type':'commit','sha':SHA}}
    monkeypatch.setattr(publisher, 'api', fake_api)
    assert publisher.tag_commit('owner/repo','v0.5.0') == SHA
    monkeypatch.setattr(publisher, 'api', lambda *a, **k: None)
    assert publisher.tag_commit('owner/repo','v0.5.0') is None


def test_release_workflow_does_not_publish_from_pull_request_code():
    text = (ROOT/'.github/workflows/release.yml').read_text()
    assert 'pull_request:' not in text and 'pull_request_target:' not in text
    assert "github.event.workflow_run.event == 'push'" in text
    assert 'branches: [main]' in text and 'persist-credentials: false' in text
    assert 'contents: write' in text and 'actions: read' in text
    assert 'cancel-in-progress: false' in text


def test_independent_positive_witness_helpers_are_exact():
    sys.path.insert(0, str(ROOT/'benchmarks'))
    try:
        spec = importlib.util.spec_from_file_location('independent_mw_audit', ROOT/'benchmarks/independent_mw.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    assert set(module.integer_roots(1,0,-4)) == {-2,2}
    assert module.integer_roots(2,0,-1) == ()
    assert module.integer_roots(0,2,-3) == ()
    assert module.integer_roots(0,2,-4) == (2,)
    assert module.integer_roots(0,0,0) is None
    import sympy as sp
    u = sp.symbols('u0:4')
    positive = sum(x*x for x in u)
    with pytest.raises(ArithmeticError, match='no insolubility claim'):
        module.small_point((positive, positive+u[0]**2), bound=1)


def test_same_named_python_reports_do_not_replace_sage_reports(tmp_path):
    python = tmp_path/'foundation-python-3.13'
    sage = tmp_path/'descent-sage-10.6'
    python.mkdir()
    sage.mkdir()
    (python/'etale-results.jsonl').write_text('not a Sage arithmetic report')
    with pytest.raises(ValueError, match='genuine-Sage'):
        gate.artifact_path(tmp_path, 'etale-results.jsonl')
    (sage/'etale-results.jsonl').write_text('genuine Sage report')
    assert gate.artifact_path(tmp_path, 'etale-results.jsonl') == sage/'etale-results.jsonl'
    (python/'independent-rank_one.jsonl').write_text('{}')
    (sage/'independent-rank_one.jsonl').write_text('{}')
    with pytest.raises(ValueError, match='exactly one'):
        gate.artifact_path(tmp_path, 'independent-rank_one.jsonl')
