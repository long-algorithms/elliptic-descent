"""Labelled Cassels--Tate matrices with portable exact local-form evidence.

PARI supplies an independent implementation's pairing rank, not our matrix.
Exact forms, product square roots and local rational evaluation points can be
replayed with ordinary Python. Parent arithmetic remains a Sage dependency.
"""
import argparse
import json
import time

from foundation import metadata
from two_selmer import FIXTURES as TWO_SELMER_FIXTURES
from elliptic_descent import __version__, two_selmer, cassels_tate_pairing

FIXTURES = {**TWO_SELMER_FIXTURES,
            'rank_three': (0,0,1,-7,6), 'split_nontrivial_radical': (0,0,0,-289,0)}


def entry_record(e):
    record = {'left_coordinates': e.left, 'right_coordinates': e.right, 'bit': e.bit,
              'left_label': e.selmer.ambient_label(e.left),
              'right_label': e.selmer.ambient_label(e.right)}
    if e.evaluation is None:
        return {**record, 'reason': 'alternation_or_zero'}
    F = e.evaluation
    return {**record, 'reason': 'Fisher_Theorem_3_1',
            'normalized_quartics': [[str(c) for c in q.f] for q in F.triple.quartics],
            'product_square_root': [str(c) for c in F.triple.square_root],
            'gamma': [str(c) for c in F.triple.gamma],
            'support': F.support,
            'normalizations': [{
                'coordinates': N.labelled_cover.coordinates,
                'ambient_label': N.labelled_cover.ambient_label,
                'delta': [str(c) for c in N.labelled_cover.delta],
                'original_quartic': [str(c) for c in N.labelled_cover.quartic.f],
                'parameter_matrix': [[str(c) for c in row] for row in N.parameter_matrix],
                'theta_invariant': [str(c) for c in N.theta_invariant],
            } for N in e.normalized],
            'local': [{'place': r.place, 'kind': r.kind, 'bit': r.bit,
                       'trials': r.trials, 'point': None if r.point is None else [str(c) for c in r.point]}
                      for r in F.local]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=['all', *FIXTURES], default='all')
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    from sage.all import pari
    from sage.env import SAGE_VERSION
    names = tuple(FIXTURES) if args.fixture == 'all' else (args.fixture,)
    for name in names:
        for repetition in range(args.repeats):
            start = time.perf_counter()
            S = two_selmer(FIXTURES[name], proof=True)
            parent_seconds = time.perf_counter()-start
            start = time.perf_counter()
            P = cassels_tate_pairing(S)
            pairing_seconds = time.perf_counter()-start
            start = time.perf_counter()
            verified, pari_rank = False, None
            if args.verify:
                pari_rank = int(pari.ellinit(list(FIXTURES[name])).ellrank()[2])
                verified = P.certified and P.verify_evidence()
                if not verified or P.rank != pari_rank:
                    raise ArithmeticError(f'pairing verification failed: {name}')
                if not all(P.in_radical(e.coordinates) for e in S.torsion_data().entries):
                    raise ArithmeticError('torsion Kummer class outside radical')
            verification_seconds = time.perf_counter()-start
            print(json.dumps({
                **metadata(), 'package_version': __version__, 'sage': SAGE_VERSION,
                'fixture': name, 'repetition': repetition,
                'scope': 'labelled Cassels-Tate pairing on Selmer-2; not Selmer-4 assembly',
                'original_ainvs': [str(c) for c in S.curves.original.ainvs],
                'working_ainvs': [str(c) for c in S.curves.working.ainvs],
                'basis_labels': S.basis_labels, 'selmer_dimension': S.dimension,
                'matrix': P.matrix, 'pairing_rank': P.rank, 'radical_basis': P.radical_basis,
                'radical_ambient_basis': P.radical_ambient_basis, 'pari_pairing_rank': pari_rank,
                'complete': P.complete, 'certified': P.certified, 'verified': verified,
                'proof': S.arithmetic.arithmetic_proof.value, 'assumptions': P.assumptions,
                'max_local_trials': None, 'four_selmer_complete': False,
                'parent_selmer_seconds': parent_seconds, 'pairing_seconds': pairing_seconds,
                'verification_seconds': verification_seconds,
                'entries': [entry_record(e) for e in P.entries],
            }, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
