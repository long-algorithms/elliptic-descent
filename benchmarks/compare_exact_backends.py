"""Differential replay of actual covering data with Python and FLINT rationals.

Run after installing the optional fast extra. Each worker starts a fresh
process so SymPy's ground domain is selected before import. Both must verify
the same field, sixteen torsion points, labelled cocycles and every class sum.
No output from one backend is accepted on behalf of the other.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def worker():
    from dataclasses import replace
    import sympy
    from sympy.external.gmpy import GROUND_TYPES
    from elliptic_descent.galois_covering import evidence_from_record
    requested = os.environ['SYMPY_GROUND_TYPES']
    if GROUND_TYPES != requested:
        raise RuntimeError(f'expected {requested}, actually using {GROUND_TYPES}')
    record = json.loads((ROOT/'tests/data/group_split_full.json').read_text())
    evidence = evidence_from_record(record['galois_evidence'])
    if not evidence.verify():
        raise ArithmeticError('actual covering evidence did not verify')
    F, T = evidence.torsion.field, evidence.class_table
    direct = tuple(c.cocycle_direct(evidence.torsion) for c in evidence.coverings)
    if direct != T.cocycles:
        raise ArithmeticError('generator and direct geometric cocycles disagree')
    bad = replace(evidence.coverings[1], norm_witness=-evidence.coverings[1].norm_witness)
    if bad.verify(evidence.torsion):
        raise ArithmeticError('altered signed covering was accepted')
    result = {
        'field_degree': F.degree, 'keys': T.keys, 'cocycles': T.cocycles,
        'actions': evidence.torsion.module.actions, 'transitions': F.transitions,
        'origins': [[list(map(str, F.pack(v))) for v in c.origin(evidence.torsion)]
                    for c in evidence.coverings],
        'addition_table': [[T.add_indices(i,j) for j in range(len(T.keys))]
                           for i in range(len(T.keys))],
        'zero_index': T.zero_index, 'tamper_rejected': True,
    }
    print(json.dumps({'ground_types': GROUND_TYPES, 'sympy': sympy.__version__,
                      'verified': True, 'exact_data': result}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker()
        return
    reports, times = {}, {}
    for backend in ('python', 'flint'):
        env = dict(os.environ, SYMPY_GROUND_TYPES=backend)
        start = time.perf_counter()
        output = subprocess.check_output([sys.executable, __file__, '--worker'],
                                         env=env, text=True)
        reports[backend] = json.loads(output)
        times[backend] = time.perf_counter()-start
    if reports['python']['exact_data'] != reports['flint']['exact_data']:
        raise ArithmeticError('exact backends disagree on a field, cocycle or sum')
    data = reports['python']['exact_data']
    import flint
    print(json.dumps({'verified': True, 'field_degree': data['field_degree'],
                      'classes': len(data['keys']), 'sums_checked': len(data['keys'])**2,
                      'backends': list(reports), 'python_flint': flint.__version__,
                      'sympy': reports['python']['sympy'], 'seconds': times,
                      'scope': 'exact backend differential check on recorded geometric evidence'},
                     sort_keys=True))


if __name__ == '__main__':
    main()
