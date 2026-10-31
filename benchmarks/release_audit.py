"""Validate the full release benchmark corpus; never infer success from a count.

This checks reporting/accounting, not the underlying mathematics a second time.
Each benchmark must already run its exact evidence verifier. See docs/v0.4.md.
"""
import argparse
import json
import math
from pathlib import Path

from elliptic_descent import __version__

CORPORA = {
    'relative-results.jsonl': {
        'x4-x+1': (4, (4,)), 'handbook-positive': (2, (4,)),
        '571a1-first': (0, (4,)), '571a1-second': (0, (4,)),
        '571a1-third': (0, (4,)), 'x4+1': (4, (4,)),
    },
    'two-quadratic-results.jsonl': {
        'two_real_fields': (2, (2,2)), 'isomorphic_real_fields': (4, (2,2)),
        'nonmonic_positive': (4, (2,2)), 'two_imaginary_fields': (2, (2,2)),
        'mixed_signature': (2, (2,2)), 'affine_positive': (4, (2,2)),
        'real_obstruction': (0, (2,2)), 'two_adic_obstruction': (0, (2,2)),
    },
    'rational-root-results.jsonl': {
        'cubic_root_zero': (2, (3,1)), 'cubic_root_three': (2, (3,1)),
        'cubic_all_real': (2, (3,1)), 'quadratic_two_roots': (2, (2,1,1)),
        'imaginary_quadratic': (2, (2,1,1)), 'split_four_torsion': (4, (1,1,1,1)),
        'negative_cubic': (2, (3,1)), 'rational_split': (4, (1,1,1,1)),
    },
}
ARITHMETIC = {'irreducible', 'cubic_linear', 'two_quadratics', 'quadratic_two_linear',
              'split', 'isomorphic_quadratics', 'class_group_two_torsion', 'nonmonic'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit_records(corpora, arithmetic, *, version=__version__):
    require(set(corpora) == set(CORPORA), 'missing or unexpected benchmark corpus')
    rows, comparisons, patterns, times = 0, 0, set(), []
    for filename, expected in CORPORA.items():
        records = corpora[filename]
        seen = set()
        for r in records:
            name, repetition = r.get('fixture'), r.get('repetition')
            require(name in expected and type(repetition) is int and repetition >= 0,
                    f'{filename}: invalid fixture/repetition')
            require((name, repetition) not in seen, f'{filename}: duplicate result')
            seen.add((name, repetition))
            prefix = f'{filename}/{name}'
            require(r.get('package_version') == version, f'{prefix}: wrong package version')
            require(all(r.get(k) is True for k in ('complete','certified','verified')),
                    f'{prefix}: incomplete or unverified')
            require(r.get('proof') == 'certified' and r.get('assumptions') == [],
                    f'{prefix}: conditional arithmetic')
            require(r.get('unresolved') == 0 and r.get('processed') == r.get('total_model_classes')
                    and type(r.get('processed')) is int and r['processed'] >= 0,
                    f'{prefix}: unresolved or incomplete candidate accounting')
            require(all(k in r and r[k] is None for k in ('max_nodes','max_depth','max_candidates')),
                    f'{prefix}: release benchmarks must have no artificial limits')
            count, pattern = expected[name]
            actual = r.get('coverings', r.get('oriented_coverings'))
            require(type(actual) is int and actual == count, f'{prefix}: covering count mismatch')
            require(tuple(r.get('factorization_type', (4,))) == pattern,
                    f'{prefix}: factorization mismatch')
            if filename == 'rational-root-results.jsonl':
                require(r.get('orientation_multiplicity') == 1 and r.get('marked_map_verified') is True,
                        f'{prefix}: missing orientation or marked-map verification')
            else:
                require(type(r.get('unoriented_coverings')) is int and r['unoriented_coverings']*2 == actual,
                        f'{prefix}: orientation count mismatch')
            dimension = r.get('eclib_selmer_dimension')
            needs_eclib = filename == 'rational-root-results.jsonl' or (
                filename == 'two-quadratic-results.jsonl' and count > 0)
            if needs_eclib:
                require(type(dimension) is int and dimension >= 0 and 2**dimension == count,
                        f'{prefix}: missing or inconsistent eclib comparison')
                comparisons += 1
            for key in ('solver_seconds','verification_seconds'):
                value = r.get(key)
                require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
                        f'{prefix}: invalid timing')
            times.append(r['solver_seconds'])
            patterns.add(pattern)
            rows += 1
        require({name for name, _ in seen} == set(expected), f'{filename}: missing fixtures')
    seen = set()
    for r in arithmetic:
        key = (r.get('fixture'), r.get('repetition'))
        require(key[0] in ARITHMETIC and type(key[1]) is int and key[1] >= 0 and key not in seen,
                'invalid or duplicate arithmetic record')
        seen.add(key)
        require(r.get('package_version') == version and r.get('verified') is True,
                'unverified/wrong-version arithmetic')
        require(r.get('proof') == 'certified' and r.get('assumptions') == [],
                'conditional arithmetic fixture')
        require(r.get('relative_descent_complete') is False, 'arithmetic fixture overclaims descent')
    require({name for name, _ in seen} == ARITHMETIC, 'missing arithmetic fixtures')
    return {'package_version': version, 'passed': True, 'descent_records': rows,
            'distinct_descent_fixtures': sum(map(len, CORPORA.values())),
            'arithmetic_records': len(arithmetic), 'eclib_comparisons': comparisons,
            'factorization_types': sorted(patterns), 'total_solver_seconds': sum(times),
            'maximum_solver_seconds': max(times),
            'scope': 'verified corpus/accounting; not an independent class-group proof'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('.'))
    args = parser.parse_args()
    def load(name):
        return [json.loads(line) for line in (args.directory/name).read_text().splitlines() if line.strip()]
    summary = audit_records({name: load(name) for name in CORPORA}, load('etale-results.jsonl'))
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
