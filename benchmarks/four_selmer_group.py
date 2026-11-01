"""Full group operations on actual global classes, with portable exact evidence.

The normal-field method is an expensive reference backend, not a performance
claim. Run with sage -python. Every fixture remains a required acceptance gate.
"""
import argparse
import json
import time

from foundation import metadata
from four_selmer_fibres import FIXTURES, geometry_record
from elliptic_descent import __version__, two_selmer, cassels_tate_pairing, assemble_four_selmer_fibres
from elliptic_descent.four_selmer_group import four_selmer_group
from elliptic_descent.galois_covering import evidence_record

EXPECTED={
    'rank_one':(4,), 'split_two_torsion':(2,2), 'cyclic_two_torsion':(2,),
    'cyclic_four_torsion':(4,), 'full_two_and_four':(2,4), 'zero_selmer':(),
    'nontrivial_sha_571a1':(2,2), 'split_nontrivial_radical':(2,2,2,2),
}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture',choices=['all',*FIXTURES],default='all')
    parser.add_argument('--backend',choices=['auto','galois'],default='auto')
    parser.add_argument('--repeats',type=int,default=1)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--max-field-degree',type=int)
    args=parser.parse_args()
    if args.repeats<1:
        parser.error('--repeats must be positive')
    from sage.env import SAGE_VERSION
    names=tuple(FIXTURES) if args.fixture=='all' else (args.fixture,)
    for name in names:
        for rep in range(args.repeats):
            start=time.perf_counter()
            S=two_selmer(FIXTURES[name],proof=True)
            B=cassels_tate_pairing(S)
            F=assemble_four_selmer_fibres(S,B,proof=True)
            parent_seconds=time.perf_counter()-start
            start=time.perf_counter()
            G=four_selmer_group(F,backend=args.backend,max_field_degree=args.max_field_degree)
            solver_seconds=time.perf_counter()-start
            start=time.perf_counter()
            verified=args.verify and G.certified and G.verify_evidence()
            if args.verify and (not verified or G.invariant_factors!=EXPECTED[name]):
                raise ArithmeticError(f'group verification failed: {name}')
            if args.verify:
                for a in G.elements:
                    if G.from_coordinates(a.coordinates)!=a or 2*a!=G.include_two_selmer(G.projection(a)):
                        raise ArithmeticError('group coordinate or doubling compatibility failed')
                    for b in G.elements:
                        if (a+b)-b!=a or G.projection(a+b)!=G.projection(a)^G.projection(b):
                            raise ArithmeticError('group inverse or projection failed')
            verification_seconds=time.perf_counter()-start
            fibre_record={**geometry_record(F),'fixture':name,
                    'original_ainvs':[str(c) for c in S.curves.original.ainvs],
                    'selmer_dimension':S.dimension,'pairing_matrix':B.matrix,'pairing_rank':B.rank,
                    'image_basis':F.image_basis,'image_classes':F.image_classes,
                    'kernel_dimension':F.kernel_dimension,'global_cardinality':F.cardinality,
                    'predicted_cardinality':F.predicted_cardinality,
                    'global_keys':[C.key for C in F.known_classes],'zero_class_key':F.zero_class.key}
            record={**metadata(),'package_version':__version__,'sage':SAGE_VERSION,
                    'fixture':name,'repetition':rep,'scope':'full group law on assembled marked 4-Selmer classes',
                    'method':G.method,'keys':G.law.keys,'addition_table':G.law.table,'zero_index':G.law.zero,
                    'invariant_factors':G.invariant_factors,'generator_keys':[a.key for a in G.generators],
                    'relations':G.relations,'complete':G.complete,'certified':G.certified,
                    'verified':verified,'assumptions':G.assumptions,'fibres':fibre_record,
                    'galois_evidence':None if G.galois_evidence is None else evidence_record(G.galois_evidence),
                    'parent_seconds':parent_seconds,'solver_seconds':solver_seconds,
                    'verification_seconds':verification_seconds,'max_field_degree':args.max_field_degree}
            print(json.dumps(record,sort_keys=True),flush=True)


if __name__=='__main__':
    main()
