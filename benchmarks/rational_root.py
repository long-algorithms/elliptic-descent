"""Rational-root descents with independently verified eclib counts and marking.

Run with sage -python. No Magma execution or timing comparison is claimed.
One JSON line per fixture and repetition; all default searches are unbounded.
"""
import argparse
import json
import time

import sympy as sp

from foundation import metadata
from elliptic_descent import __version__, relative_four_descent
from elliptic_descent.rational_root import RationalRootIdentification

FIXTURES = {
    'cubic_root_zero': (0,1,-1,0,1),
    'cubic_root_three': (-3,4,-1,-3,1),
    'cubic_all_real': (0,1,-3,0,1),
    'quadratic_two_roots': (2,0,-3,0,1),
    'imaginary_quadratic': (-1,0,0,0,1),
    'split_four_torsion': (4,0,-5,0,1),
    'negative_cubic': (0,-1,1,0,-1),
    'rational_split': (sp.Rational(1,4),0,sp.Rational(-5,4),0,1),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=['all', *FIXTURES], default='all')
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    from sage.all import EllipticCurve, QQ
    from sage.env import SAGE_VERSION
    names = tuple(FIXTURES) if args.fixture == 'all' else (args.fixture,)
    for name in names:
        for repetition in range(args.repeats):
            start = time.perf_counter()
            result = relative_four_descent(FIXTURES[name], proof=True)
            solver_seconds = time.perf_counter()-start
            identification = RationalRootIdentification(result.source)
            verification_seconds, dimension, verified, marked = None, None, False, False
            if args.verify:
                start = time.perf_counter()
                marked = identification.check_marked_map()
                verified = result.certified and result.verify_evidence()
                E = EllipticCurve([QQ(str(c)) for c in identification.branch_cubic_ainvs])
                dimension = int(E.global_minimal_model().mwrank_curve(verbose=False).selmer_rank())
                if (not verified or not marked or not result.coverings
                        or len(result.coverings) != 2**dimension
                        or result.orientation.multiplicity != 1):
                    raise ArithmeticError(f'rational-root verification failed: {name}')
                verification_seconds = time.perf_counter()-start
            local = []
            for candidate in result.candidates:
                for report in candidate.local.finite:
                    certificate = report.certificate
                    local.append({'label': candidate.label, 'prime': report.place,
                                  'status': report.status.value, 'nodes': report.nodes,
                                  'depth': report.depth,
                                  'certificate': certificate.kind if certificate else None,
                                  'model_changes': certificate.change.steps
                                  if certificate and certificate.kind == 'model_change' else 0})
            print(json.dumps({
                **metadata(), 'package_version': __version__, 'sage': SAGE_VERSION,
                'fixture': name, 'repetition': repetition,
                'f_ascending': [str(c) for c in result.source.f],
                'root': str(identification.root),
                'branch_cubic_ainvs': [str(c) for c in identification.branch_cubic_ainvs],
                'scope': result.scope, 'factorization_type': result.factorization_type,
                'orientation_multiplicity': result.orientation.multiplicity,
                'reversing_multiplier': [str(c) for c in result.orientation.reversing_multiplier],
                'proof': result.arithmetic_proof.value, 'assumptions': result.assumptions,
                'support': result.support, 'norm_dimension': result.norm_dimension,
                'scalar_image_dimension': len(result.scalar_basis),
                'total_model_classes': result.total_model_classes, 'processed': len(result.candidates),
                'unresolved': len(result.unresolved), 'coverings': len(result.coverings),
                'complete': result.complete, 'certified': result.certified,
                'verified': verified, 'marked_map_verified': marked,
                'eclib_selmer_dimension': dimension, 'max_nodes': None,
                'max_depth': None, 'max_candidates': None,
                'solver_seconds': solver_seconds, 'verification_seconds': verification_seconds,
                **dict(result.timings), 'local_decisions': local,
            }, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
