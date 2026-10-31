"""Model/marking benchmarks, not Selmer-group computations or speed claims."""
import argparse
import json
import time
from dataclasses import replace

from foundation import metadata
from elliptic_descent import __version__
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange
from elliptic_descent.cubic_conversion import IdentityTwoCover

FIXTURES = {
    'short': ((0,0,0,-1,1), (0,1)),
    'general': ((0,0,1,-7,6), (0,2)),
    'four_torsion': ((0,-5,0,-16,80), (8,12)),
    'nonminimal': ((0,0,0,-16,64), (0,8)),
}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats',type=int,default=1)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--sage',action='store_true')
    args=parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    for name,(a,P) in FIXTURES.items():
        for rep in range(args.repeats):
            start=time.perf_counter()
            E=WeierstrassModel(a)
            I=IdentityTwoCover(E)
            marked=replace(I.cover,translation=P)
            # Materialize generic maps so construction timing is not just wrappers.
            _=marked.map_xy
            construction=time.perf_counter()-start
            verification=None
            checked=False
            if args.verify:
                start=time.perf_counter()
                checked=(I.verify() and marked.verify()
                         and I.cover.evaluate(I.to_quartic(P))==E.double(P)
                         and marked.evaluate(I.to_quartic(P))==E.add(E.double(P),P))
                if not checked:
                    raise ArithmeticError(f'marking verification failed: {name}')
                verification=time.perf_counter()-start
            sage_check=None
            if args.sage:
                from elliptic_descent.backends.sage_models import rational_curve_context
                context=rational_curve_context(E)
                sage_check=context.to_working.verify() and context.to_original.verify()
                if not sage_check:
                    raise ArithmeticError('normalization check failed')
            print(json.dumps({**metadata(),'package_version':__version__,'fixture':name,
                'repetition':rep,'scope':'marked model interfaces; not Selmer computation',
                'verified':checked,'sage_normalization_verified':sage_check,
                'original_ainvs':[str(c) for c in E.ainvs],
                'quartic':[str(c) for c in I.cover.source.f],
                'construction_seconds':construction,'verification_seconds':verification,
                'curve_level_selmer_complete':False},sort_keys=True),flush=True)


if __name__=='__main__':
    main()
