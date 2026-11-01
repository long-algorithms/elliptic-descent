"""Staged and warm curve-level API audits against the explicit pipeline.

The two routes share mathematical kernels; this is an integration/reuse audit,
not an independent descent algorithm or a speed comparison with other software.
"""
import argparse
import json
from time import perf_counter

from foundation import metadata
from elliptic_descent import (
    __version__, DescentSession, DescentLimits, cassels_tate_pairing,
    assemble_four_selmer_fibres,
)

FIXTURES = {
    'rank_one': ((0, 0, 0, -1, 1), (4,)),
    'split_two_torsion': ((0, 0, 0, -1, 0), (2, 2)),
    'zero_selmer': ((0, -1, 1, -10, -20), ()),
    'nontrivial_sha_571a1': ((0, -1, 1, -929, -10595), (2, 2)),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=['all', *FIXTURES], default='all')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    names = tuple(FIXTURES) if args.fixture == 'all' else (args.fixture,)
    for name in names:
        a, expected = FIXTURES[name]
        session = DescentSession(a, proof=True)
        start = perf_counter()
        initial = session.run(limits=DescentLimits(two_selmer_candidates=0))
        if initial.complete or initial.pairing is not None:
            raise ArithmeticError('an empty candidate prefix became a complete result')
        intermediate = session.run(until='two_selmer')
        if not intermediate.finished or intermediate.complete:
            raise ArithmeticError('staged completion was confused with full group completion')
        result = session.run()
        cold_seconds = perf_counter()-start
        counts = session.cache_info()
        start = perf_counter()
        warm = session.run()
        warm_seconds = perf_counter()-start
        if (warm.group is not result.group or warm.selmer is not intermediate.selmer
                or counts != session.cache_info() or counts.get('arithmetic_contexts') != 1):
            raise ArithmeticError('repeated API access reran completed construction work')
        verified, reference_equal = False, False
        start = perf_counter()
        if args.verify:
            if not result.certified or not warm.verify_evidence() or result.invariant_factors != expected:
                raise ArithmeticError('integrated group evidence failed')
            B = cassels_tate_pairing(result.selmer)
            F = assemble_four_selmer_fibres(result.selmer, B)
            G = F.group()
            reference_equal = G.law == result.group.law and B.matrix == result.pairing.matrix
            if not reference_equal or not G.verify_evidence():
                raise ArithmeticError('integrated and explicit pipeline disagree')
            for v in result.group.elements:
                if result.covering(v.key) is not v.covering or v.covering.target != result.curve:
                    raise ArithmeticError('on-demand representative lost its original marking')
                if 2*v != result.group.include_two_selmer(result.group.projection(v)):
                    raise ArithmeticError('exact-sequence maps changed at the public interface')
            verified = True
        print(json.dumps({
            **metadata(), 'package_version': __version__, 'fixture': name,
            'scope': 'curve-level API integration and per-session reuse',
            'original_ainvs': list(a), 'cardinality': result.cardinality,
            'invariant_factors': result.invariant_factors,
            'complete': result.complete, 'certified': result.certified,
            'assumptions': result.assumptions, 'verified': verified,
            'explicit_pipeline_equal': reference_equal,
            'work': counts, 'warm_work_delta': {k: session.cache_info()[k]-v for k, v in counts.items()},
            'cold_staged_seconds': cold_seconds, 'warm_seconds': warm_seconds,
            'verification_and_reference_seconds': perf_counter()-start,
            'timings': dict(result.timings), 'group_method': result.group.method,
        }, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
