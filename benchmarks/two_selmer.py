"""Full labelled 2-Selmer benchmarks, with exact evidence and independent counts.

Run with sage -python. Enumeration is deliberately exhaustive, not a speed
comparison with eclib or PARI. No full 4-Selmer claim is made.
"""
import argparse
import json
import time

from foundation import metadata
from elliptic_descent import __version__
from elliptic_descent.two_selmer import two_selmer

FIXTURES = {
    'rank_one': (0,0,0,-1,1),
    'split_two_torsion': (0,0,0,-1,0),
    'cyclic_two_torsion': (0,0,0,1,0),
    'cyclic_four_torsion': (0,0,0,4,0),
    'full_two_and_four': (0,-5,0,-16,80),
    'zero_selmer': (0,-1,1,-10,-20),
    'nontrivial_sha_571a1': (0,-1,1,-929,-10595),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=['all',*FIXTURES], default='all')
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    from sage.all import EllipticCurve
    from sage.env import SAGE_VERSION
    names = tuple(FIXTURES) if args.fixture == 'all' else (args.fixture,)
    for name in names:
        for repetition in range(args.repeats):
            start = time.perf_counter()
            result = two_selmer(FIXTURES[name], proof=True)
            solver_seconds = time.perf_counter()-start
            start = time.perf_counter()
            data = result.torsion_data()
            verified, independent = False, None
            if args.verify:
                verified = result.certified and result.verify_evidence() and data.verify(result)
                independent = int(EllipticCurve(list(FIXTURES[name])).global_minimal_model()
                                  .mwrank_curve(verbose=False).selmer_rank())
                if not verified or result.dimension != independent:
                    raise ArithmeticError(f'labelled 2-Selmer audit failed: {name}')
                for c in range(result.cardinality):
                    if result.coordinates(result.representative(c)) != c:
                        raise ArithmeticError(f'coordinate audit failed: {name}')
            verification_seconds = time.perf_counter()-start
            print(json.dumps({
                **metadata(), 'package_version': __version__, 'sage': SAGE_VERSION,
                'fixture': name, 'repetition': repetition, 'scope': 'full labelled 2-Selmer over Q',
                'original_ainvs': [str(c) for c in result.curves.original.ainvs],
                'working_ainvs': [str(c) for c in result.curves.working.ainvs],
                'cubic': [str(c) for c in result.arithmetic.algebra.f],
                'factor_degrees': result.arithmetic.algebra.degrees,
                'support': result.arithmetic.support, 'ambient_dimension': result.arithmetic.dimension,
                'norm_kernel_dimension': len(result.candidate_space.basis),
                'total_candidates': result.total_candidates, 'processed': len(result.decisions),
                'unresolved': len(result.unresolved), 'selmer_dimension': result.dimension,
                'selmer_cardinality': result.cardinality, 'basis_labels': result.basis_labels,
                'basis_representatives': [[str(c) for c in result.representative(1<<i)]
                                          for i in range(result.dimension)],
                'torsion_two_order': len(data.entries), 'torsion_image_rank': data.image_rank,
                'torsion_kernel_order': len(data.kernel), 'rational_four_torsion_order': len(data.four_torsion),
                'complete': result.complete, 'certified': result.certified, 'verified': verified,
                'proof': result.arithmetic.arithmetic_proof.value, 'assumptions': result.assumptions,
                'eclib_selmer_dimension': independent, 'four_selmer_complete': False,
                'max_candidates': None, 'max_nodes': None, 'max_depth': None,
                'solver_seconds': solver_seconds, 'verification_seconds': verification_seconds,
                **dict(result.timings),
            }, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
