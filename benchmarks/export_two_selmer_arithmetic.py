"""Export exact certified cubic arithmetic for local-solver replay.

These files retain provenance but are not independently checkable class-group
certificates. No Selmer-completeness claim follows from a norm kernel alone.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter

from elliptic_descent.backends.sage_models import rational_curve_context
from elliptic_descent.backends.sage_two_selmer import cubic_square_classes
from elliptic_descent.cubic_arithmetic import kummer_algebra

FIXTURES = {
    'rank_one': (0,0,0,-1,1), 'split_two_torsion': (0,0,0,-1,0),
    'cyclic_two_torsion': (0,0,0,1,0), 'cyclic_four_torsion': (0,0,0,4,0),
    'full_two_and_four': (0,-5,0,-16,80), 'zero_selmer': (0,-1,1,-10,-20),
    'nontrivial_sha_571a1': (0,-1,1,-929,-10595),
}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    for name,curve in FIXTURES.items():
        start=perf_counter()
        curves=rational_curve_context(curve)
        A=cubic_square_classes(kummer_algebra(curves.working),proof=True)
        record={'fixture':name,'original_ainvs':curve,
                'working_ainvs':[str(c) for c in curves.working.ainvs],
                'cubic':[str(c) for c in A.algebra.f], 'support':A.support,
                'generators':[[str(c) for c in g] for g in A.generators],
                'norm_columns':A.norm_columns,'norm_kernel_basis':A.norm_kernel.basis,
                'norm_kernel_count':A.norm_kernel.cardinality,'backend':A.backend,
                'proof':A.arithmetic_proof.value,'selmer_complete':False,
                'seconds':perf_counter()-start}
        (args.output/(name+'.json')).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
        print(name,A.dimension,A.norm_kernel.cardinality,record['seconds'],flush=True)


if __name__=='__main__': main()
