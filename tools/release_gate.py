"""Fail-closed v0.5 artifact inventory; checks reporting, not mathematics anew.

Use the exact checkout and artifacts from one revision. Separate CI jobs must
have run the real evidence verifiers. No network access or publication here.
"""
import ast
import importlib.util
import json
from math import prod
from pathlib import Path
import re
import xml.etree.ElementTree as ET

VERSION = '0.5.0'
GROUPS = {
    'rank_one': (4,), 'split_two_torsion': (2, 2),
    'cyclic_two_torsion': (2,), 'cyclic_four_torsion': (4,),
    'full_two_and_four': (2, 4), 'zero_selmer': (),
    'nontrivial_sha_571a1': (2, 2), 'split_nontrivial_radical': (2, 2, 2, 2),
}
INDEPENDENT = set(GROUPS)-{'nontrivial_sha_571a1', 'split_nontrivial_radical'}
LABELLED = {'rank_one', 'split_two_torsion', 'cyclic_two_torsion', 'cyclic_four_torsion'}
API = {'rank_one', 'split_two_torsion', 'zero_selmer', 'nontrivial_sha_571a1'}
SEL2 = {'rank_one': 1, 'split_two_torsion': 2, 'cyclic_two_torsion': 1,
        'cyclic_four_torsion': 1, 'full_two_and_four': 2, 'zero_selmer': 0,
        'nontrivial_sha_571a1': 2}
PAIRING = {**{k: 0 for k in SEL2}, 'nontrivial_sha_571a1': 2,
           'rank_three': 0, 'split_nontrivial_radical': 2}
SAGE_GROUPS = {'core', 'api', 'pairing', 'radical', 'fibres-fixtures',
               'fibres-behaviour', 'fibres-reference', 'group-behaviour',
               *(f'group-curve-{i}' for i in range(8))}
WORKFLOWS = {
    'tests.yml': 'Tests', 'distribution.yml': 'Distribution audit',
    'group-benchmarks.yml': 'Four-Selmer group audit',
    'global-fibre-audit.yml': 'Global fibre audit', 'public-api.yml': 'Public API audit',
    'release-validation.yml': 'Independent release validation',
    'two-selmer-arithmetic.yml': 'Two-Selmer arithmetic audit',
    'two-quadratic-audit.yml': 'Two-quadratic arithmetic audit',
    'relative-arithmetic.yml': 'Relative arithmetic fixtures',
}


def require(value, message):
    if not value:
        raise ValueError(message)


def selected_runs(runs, sha):
    """Latest push/main run per exact workflow path; missing/pending means wait.

    A failed/cancelled latest run cannot be hidden by an older successful run.
    A same-name PR job or a run at another revision cannot satisfy this gate.
    """
    found = {}
    for r in runs:
        path = r.get('path', '').split('@')[0]
        name = path.removeprefix('.github/workflows/')
        if (path != '.github/workflows/'+name or name not in WORKFLOWS
                or r.get('name') != WORKFLOWS[name] or r.get('head_sha') != sha
                or r.get('head_branch') != 'main' or r.get('event') != 'push'):
            continue
        if name not in found or (r['id'], r.get('run_attempt', 1)) > (
                found[name]['id'], found[name].get('run_attempt', 1)):
            found[name] = r
    for name, r in found.items():
        if r.get('status') == 'completed':
            require(r.get('conclusion') == 'success', f'{name}: final run is not successful')
    if set(found) != set(WORKFLOWS) or any(r.get('status') != 'completed' for r in found.values()):
        return None
    return found


def records(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def clean(row, version, revision, *, complete='complete'):
    require(row.get('package_version', row.get('version')) == version, 'wrong benchmark version')
    require(row.get('revision') == revision, 'mixed benchmark revisions')
    require(row.get('verified') is True, 'unverified record')
    if complete:
        require(row.get(complete) is True and row.get('certified') is True, 'incomplete/uncertified record')
        require(row.get('assumptions') == [], 'conditional release record')


def unique(rows, expected):
    names = [r.get('fixture') for r in rows]
    require(len(names) == len(set(names)) and set(names) == set(expected), 'missing/duplicate/unexpected fixture')
    return {r['fixture']: r for r in rows}


def independent_records(rows, version, revision):
    data = unique(rows, INDEPENDENT)
    for name, r in data.items():
        clean(r, version, revision)
        n = prod(GROUPS[name])
        require(r.get('sha_two_and_four_trivial') is True, 'missing independent Sha bound')
        require(tuple(r.get('expected_invariant_factors', ())) == GROUPS[name]
                and tuple(r.get('invariant_factors', ())) == GROUPS[name]
                and type(r.get('classes')) is int and r['classes'] == n, 'independent group mismatch')
        labelled = name in LABELLED
        require(r.get('labelled') is labelled, 'changed independent audit scope')
        expected = n*n if labelled else 0
        require(r.get('ordered_sums_checked') == expected
                and r.get('pairwise_separations_checked') == expected, 'incomplete independent sums')
        witnesses = r.get('rational_point_witnesses', [])
        require(len(witnesses) == (n if labelled else 0), 'missing independent rational points')
        require(len({tuple(w['key']) for w in witnesses}) == len(witnesses), 'duplicate named class')
        require('site-packages' in r.get('installed_module', '') or
                'dist-packages' in r.get('installed_module', ''), 'independent check did not use installed wheel')
    return {'curves': len(data), 'labelled_curves': len(LABELLED),
            'labelled_classes': sum(prod(GROUPS[n]) for n in LABELLED),
            'ordered_sums': sum(prod(GROUPS[n])**2 for n in LABELLED)}


def sage_reports(directory, source):
    paths = list(Path(directory).rglob('sage-*.xml'))
    require({p.stem.removeprefix('sage-') for p in paths} == SAGE_GROUPS
            and len(paths) == len(SAGE_GROUPS), 'missing/duplicate Sage job report')
    identities, covered, total = set(), set(), 0
    for path in paths:
        root = ET.parse(path).getroot()
        cases = list(root.iter('testcase'))
        require(cases, 'empty Sage test report')
        for c in cases:
            require(not any(c.find(k) is not None for k in ('failure', 'error', 'skipped')), 'nonpassing Sage case')
            key = (c.get('classname'), c.get('name'))
            require(key not in identities, 'duplicate Sage test identity')
            identities.add(key)
            covered.add((key[0].replace('.', '/')+'.py', key[1].split('[')[0]))
            total += 1
    expected = set()
    for p in (Path(source)/'tests').glob('test_*.py'):
        for node in ast.parse(p.read_text()).body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith('test_'):
                expected.add((p.relative_to(source).as_posix(), node.name))
    require(expected <= covered, 'source test functions absent from genuine Sage reports')
    return {'passed': total, 'failures': 0, 'errors': 0, 'skips': 0,
            'source_functions': len(expected), 'modules': len({p for p, _ in expected})}


def artifact_path(directory, name):
    """Select Sage regression reports explicitly; Python jobs share filenames."""
    directory = Path(directory)
    regression = {'sage-results.jsonl', 'relative-results.jsonl',
                  'etale-results.jsonl', 'two-quadratic-results.jsonl',
                  'rational-root-results.jsonl', 'release-audit.json',
                  'marked-models-results.jsonl', 'two-selmer-results.jsonl',
                  'selmer-quartics-results.jsonl', 'cassels-tate-results.jsonl',
                  'pairing-replay.json'}
    if name in regression:
        path = directory/'descent-sage-10.6'/name
        require(path.is_file(), f'missing genuine-Sage regression report: {name}')
        return path
    paths = list(directory.rglob(name))
    require(len(paths) == 1, f'expected exactly one {name}')
    return paths[0]


def audit_directory(directory, source, revision, version=VERSION):
    directory, source = Path(directory), Path(source)
    require(re.fullmatch('[0-9a-f]{40}', revision) is not None, 'invalid revision')

    def one(name):
        return artifact_path(directory, name)

    def load(name):
        return records(one(name))

    spec = importlib.util.spec_from_file_location('legacy_release_audit', source/'benchmarks/release_audit.py')
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    corpora = {n: load(n) for n in legacy.CORPORA}
    arithmetic = load('etale-results.jsonl')
    old = legacy.audit_records(corpora, arithmetic, version=version)
    for rows in (*corpora.values(), arithmetic):
        require(all(r.get('revision') == revision for r in rows), 'mixed legacy benchmark revisions')

    sel2 = unique(load('two-selmer-results.jsonl'), SEL2)
    for name, r in sel2.items():
        clean(r, version, revision)
        require(r['selmer_dimension'] == r['eclib_selmer_dimension'] == SEL2[name], 'Selmer-2 comparison failed')
        require(r['unresolved'] == 0 and r['processed'] == r['total_candidates'], 'partial Selmer-2 ledger')
        require(all(r[k] is None for k in ('max_candidates', 'max_nodes', 'max_depth')), 'bounded release calculation')
    quartics = load('selmer-quartics-results.jsonl')
    seen = set()
    for r in quartics:
        clean(r, version, revision, complete='complete_parent')
        name, c = r['fixture'], r['coordinates']
        require(name in SEL2 and type(c) is int and 0 <= c < 2**SEL2[name], 'invalid quartic label')
        require((name, c) not in seen, 'duplicate quartic')
        seen.add((name, c))
    require(seen == {(n, c) for n, d in SEL2.items() for c in range(2**d)}, 'missing labelled quartic')
    pairing = unique(load('cassels-tate-results.jsonl'), PAIRING)
    for name, r in pairing.items():
        clean(r, version, revision)
        require(r['pairing_rank'] == r['pari_pairing_rank'] == PAIRING[name], 'pairing comparison failed')
        require(r['max_local_trials'] is None, 'bounded pairing')
    require(json.loads(one('pairing-replay.json').read_text()).get('passed') is True, 'pairing replay failed')
    marks = unique(load('marked-models-results.jsonl'), {'short','general','four_torsion','nonminimal'})
    for r in marks.values():
        clean(r, version, revision, complete=None)
        require(r['sage_normalization_verified'] is True, 'marked normalization failed')

    sums, witnesses, classes = 0, 0, 0
    for name, factors in GROUPS.items():
        n = prod(factors)
        g = unique(load(f'four-selmer-group-{name}.jsonl'), {name})[name]
        clean(g, version, revision)
        require(g['max_field_degree'] is None and tuple(g['invariant_factors']) == factors, 'wrong/bounded group')
        require(len(g['keys']) == n and len({tuple(k) for k in g['keys']}) == n, 'group labels incomplete')
        require(len(g['addition_table']) == n and all(len(row) == n for row in g['addition_table']), 'sum table incomplete')
        replay = unique(load(f'group-replay-{name}.jsonl'), {name})[name]
        require(replay.get('verified') is True and replay['classes'] == n
                and replay['sums_checked'] == n*n and tuple(replay['invariant_factors']) == factors, 'group replay failed')
        require(replay['parent_replay']['verified'] is True, 'group parent replay failed')
        f = unique(load(f'four-selmer-fibres-{name}.jsonl'), {name})[name]
        clean(f, version, revision)
        require(f['global_cardinality'] == f['predicted_cardinality'] == n and not f['unresolved'], 'fibre accounting failed')
        fr = json.loads(one(f'fibre-replay-{name}.json').read_text())
        require(fr['passed'] is True and fr['global_classes'] == n, 'fibre replay failed')
        sums += n*n
        classes += n
        witnesses += replay['parent_replay']['equivalence_witnesses']
    for name in API:
        r = unique(load(f'public-api-{name}.jsonl'), {name})[name]
        clean(r, version, revision)
        require(r['explicit_pipeline_equal'] is True and r['work']['arithmetic_contexts'] == 1
                and r['warm_work_delta'] and all(v == 0 for v in r['warm_work_delta'].values()), 'API reuse regression')
    independent = independent_records([r for n in INDEPENDENT for r in load(f'independent-{n}.jsonl')], version, revision)
    return {'version': version, 'revision': revision, 'passed': True,
            'sage': sage_reports(directory, source), 'legacy': old,
            'group_curves': len(GROUPS), 'group_classes': classes, 'ordered_sums': sums,
            'parent_equivalence_witnesses': witnesses, 'independent': independent,
            'selmer2_curves': len(sel2), 'quartics': len(quartics), 'pairings': len(pairing),
            'api_curves': len(API), 'scope': 'single-revision release inventory; not a new mathematical verifier'}
