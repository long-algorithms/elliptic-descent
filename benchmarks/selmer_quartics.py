"""All labelled classes -> marked quartics, with replayable exact conic evidence.

Construction includes its internal checks. A separate audit timing covers
label/local/marking checks. Selmer arithmetic is timed once per input curve.
No 4-Selmer group, live Magma comparison or speed superiority is claimed.
"""
import argparse
import json
import time

from foundation import metadata
from two_selmer import FIXTURES
from elliptic_descent import __version__, two_selmer
from elliptic_descent.etale import EtaleQuarticAlgebra


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture',choices=['all',*FIXTURES],default='all')
    parser.add_argument('--repeats',type=int,default=1)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.repeats<1:
        parser.error('--repeats must be positive')
    from sage.all import EllipticCurve
    from sage.env import SAGE_VERSION
    names=tuple(FIXTURES) if args.fixture=='all' else (args.fixture,)
    for name in names:
        for repetition in range(args.repeats):
            start=time.perf_counter()
            S=two_selmer(FIXTURES[name],proof=True)
            arithmetic_seconds=time.perf_counter()-start
            independent=None
            if args.verify:
                independent=int(EllipticCurve(list(FIXTURES[name])).global_minimal_model()
                                .mwrank_curve(verbose=False).selmer_rank())
                if not S.certified or not S.verify_evidence() or S.dimension!=independent:
                    raise ArithmeticError(f'Selmer arithmetic audit failed: {name}')
            labels=set()
            for coordinates in range(S.cardinality):
                start=time.perf_counter()
                C=S.quartic_cover(coordinates)
                construction_seconds=time.perf_counter()-start
                verified,verification_seconds=False,None
                if args.verify:
                    start=time.perf_counter()
                    verified=C.certified and C.verify()
                    if (not verified or C.ambient_label in labels
                            or S.coordinates(C.delta)!=coordinates or C.target!=S.curves.original):
                        raise ArithmeticError(f'labelled quartic audit failed: {name}/{coordinates}')
                    verification_seconds=time.perf_counter()-start
                labels.add(C.ambient_label)
                bits=max(max(abs(int(c.p)).bit_length(),int(c.q).bit_length()) for c in C.quartic.f)
                print(json.dumps({
                    **metadata(),'package_version':__version__,'sage':SAGE_VERSION,
                    'fixture':name,'repetition':repetition,'coordinates':coordinates,
                    'ambient_label':C.ambient_label,'selmer_basis_labels':S.basis_labels,
                    'selmer_cardinality':S.cardinality,'selmer_dimension':S.dimension,
                    'eclib_selmer_dimension':independent,
                    'original_ainvs':[str(a) for a in S.curves.original.ainvs],
                    'working_ainvs':[str(a) for a in S.curves.working.ainvs],
                    'working_to_original_urst':[str(a) for a in S.curves.to_original.urst],
                    'delta':[str(a) for a in C.delta],
                    'quartic_ascending':[str(a) for a in C.quartic.f],
                    'factorization_type':EtaleQuarticAlgebra(C.quartic.f).factorization_type,
                    'parametrization_matrix':[[str(a) for a in row] for row in C.presentation.parametrization.matrix],
                    'jacobian_to_original_urst':[str(a) for a in C.marking.change.urst],
                    'max_quartic_coefficient_bits':bits,'verified':verified,
                    'complete_parent':S.complete,'certified':C.certified,
                    'proof':S.arithmetic.arithmetic_proof.value,'assumptions':C.assumptions,
                    'scope':'marked quartic for one labelled 2-Selmer class',
                    'four_selmer_complete':False,'relative_descent_computed':False,
                    'selmer_arithmetic_seconds':arithmetic_seconds,
                    'construction_and_initial_checks_seconds':construction_seconds,
                    'verification_seconds':verification_seconds,
                },sort_keys=True),flush=True)
            if labels!=set(S.known_soluble_labels):
                raise ArithmeticError(f'missing or duplicate classes: {name}')


if __name__=='__main__':
    main()
