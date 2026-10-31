"""The release-report validator rejects missing, duplicate and conditional records."""
import copy
import importlib.util
from pathlib import Path
import re

import pytest
from elliptic_descent import __version__

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('release_audit', ROOT/'benchmarks/release_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def corpus():
    records = {}
    for filename, cases in audit.CORPORA.items():
        records[filename] = []
        for name, (count, pattern) in cases.items():
            r = dict(fixture=name, repetition=0, package_version=__version__, complete=True,
                     certified=True, verified=True, proof='certified', assumptions=[],
                     unresolved=0, processed=4, total_model_classes=4,
                     max_nodes=None, max_depth=None, max_candidates=None,
                     factorization_type=pattern, solver_seconds=1.0, verification_seconds=0.5,
                     eclib_selmer_dimension=count.bit_length()-1 if count else None)
            if filename == 'rational-root-results.jsonl':
                r.update(coverings=count, orientation_multiplicity=1, marked_map_verified=True)
            else:
                r.update(oriented_coverings=count, unoriented_coverings=count//2)
            records[filename].append(r)
    arithmetic = [dict(fixture=n, repetition=0, package_version=__version__, verified=True,
                       proof='certified', assumptions=[], relative_descent_complete=False)
                  for n in sorted(audit.ARITHMETIC)]
    return records, arithmetic


def test_complete_reporting_matrix():
    records, arithmetic = corpus()
    report = audit.audit_records(records, arithmetic)
    assert report['passed'] and report['distinct_descent_fixtures'] == 22
    assert report['eclib_comparisons'] == 14 and len(report['factorization_types']) == 5
    assert report['arithmetic_records'] == 8


@pytest.mark.parametrize('key,value', [
    ('certified',False), ('verified',False), ('complete',False), ('unresolved',1),
    ('assumptions',['GRH']), ('package_version','0.3.0'), ('max_nodes',100),
    ('processed',3), ('oriented_coverings',999), ('solver_seconds',float('nan')),
])
def test_bad_records_fail_the_release_gate(key, value):
    records, arithmetic = corpus()
    records['relative-results.jsonl'][0][key] = value
    with pytest.raises(ValueError):
        audit.audit_records(records, arithmetic)


def test_missing_duplicate_or_mislabelled_corpora_are_not_passes():
    records, arithmetic = corpus()
    for bad in ('missing', 'duplicate', 'marking', 'eclib', 'arithmetic'):
        c, a = copy.deepcopy(records), copy.deepcopy(arithmetic)
        if bad == 'missing':
            c['relative-results.jsonl'].pop()
        elif bad == 'duplicate':
            c['relative-results.jsonl'].append(c['relative-results.jsonl'][0])
        elif bad == 'marking':
            c['rational-root-results.jsonl'][0]['marked_map_verified'] = False
        elif bad == 'eclib':
            c['rational-root-results.jsonl'][0]['eclib_selmer_dimension'] = None
        else:
            a[0]['relative_descent_complete'] = True
        with pytest.raises(ValueError):
            audit.audit_records(c,a)


def test_runtime_metadata_and_release_document_agree():
    version = re.search(r'^version = "([^"]+)"', (ROOT/'pyproject.toml').read_text(), re.M).group(1)
    assert version == __version__ == '0.4.0'
    assert '# v0.4.0' in (ROOT/'docs/v0.4.md').read_text()
