"""Explicit global fibres, orbit witnesses and kernel maps; not Sel4 group addition.

Record complete exact equivalence geometry for replay without Sage. Parent
Selmer completeness and local enumeration remain separate certified dependencies.
"""
import argparse
import json
import time

import sympy as sp

from foundation import metadata
from two_selmer import FIXTURES as BASE
from elliptic_descent import __version__, two_selmer, cassels_tate_pairing, assemble_four_selmer_fibres
from elliptic_descent.cubic_arithmetic import _poly
from elliptic_descent.fibre_actions import zero_fibre_kummer
from elliptic_descent.backends.sage_pairing import cubic_product_square_root

FIXTURES = {**BASE, 'split_nontrivial_radical': (0,0,0,-289,0)}
EXPECTED = {'rank_one':4,'split_two_torsion':4,'cyclic_two_torsion':2,
            'cyclic_four_torsion':4,'full_two_and_four':8,'zero_selmer':1,
            'nontrivial_sha_571a1':4,'split_nontrivial_radical':16}


def vector(v):
    return [str(c) for c in v]


def point(P):
    return None if P is None else vector(P)


def matrix(M):
    return [vector(row) for row in M]


def norm_pair(D):
    return {'xi':vector(D.xi),'norm_witness':str(D.norm_witness)}


def class_square_root(arithmetic, value, representative):
    h = _poly(arithmetic.algebra.f)
    ratio = (_poly(value)*sp.invert(_poly(representative),h)).rem(h)
    return cubic_product_square_root(arithmetic,tuple(ratio.nth(i) for i in range(3)))


def geometry_record(F):
    S=F.selmer
    fibres=[]
    for f in F.fibres:
        C=f.representative
        actions=[]
        for a in f.actions:
            actions.append({'point':point(a.point),'working_point':point(a.action.point),
                'matrix':matrix(a.action.matrix),'y_scale':str(a.action.y_scale),
                'cubic_multiplier':vector(a.action.multiplier),'permutation':a.permutation,
                'witnesses':[{'raw':norm_pair(w.transported.raw),
                    'scalar':str(w.scalar_square.scalar),'multiplier':vector(w.scalar_square.multiplier),
                    'coordinate_matrix':matrix(w.coordinate_matrix.tolist())} for w in a.witnesses]})
        fibres.append({'coordinates':f.coordinates,'ambient_label':C.ambient_label,
            'delta':vector(C.delta),'conic_matrix':matrix(C.presentation.parametrization.matrix),
            'quartic':vector(C.quartic.f),'jacobian_to_original':vector(C.marking.change.urst),
            'total_relative_candidates':f.relative.total_model_classes,
            'relative_complete':f.lifts.complete,'relative_certified':f.lifts.certified,
            'relative_unresolved':len(f.lifts.unresolved),'expected_liftable':f.expected_liftable,
            'models':[norm_pair(D) for D in f.models],'actions':actions,'orbits':f.orbits})
    zero=F.fibre(0)
    kernel=[]
    for label,D in zip(F.kernel_map.labels,zero.models):
        delta=zero_fibre_kummer(zero.representative,D,S.arithmetic)
        representative=S.representative(label)
        root=class_square_root(S.arithmetic,delta,representative)
        kernel.append({'coordinates':label,'delta':vector(delta),
                       'representative':vector(representative),'square_root':vector(root)})
    return {'working_ainvs':vector(S.curves.working.ainvs),
            'working_to_original':vector(S.curves.to_original.urst),'fibres':fibres,
            'kernel_labels':F.kernel_map.labels,'kernel_classes':kernel,
            'kernel_quotient_basis':F.kernel_map.quotient.basis,
            'torsion_image_basis':F.torsion.image_basis,
            'torsion':[{'point':point(e.point),'coordinates':e.coordinates,
                        'halves':[point(P) for P in e.rational_halves]} for e in F.torsion.entries],
            'exclusions':[{'coordinates':e.coordinates,'test_vector':e.test_vector} for e in F.exclusions]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture',choices=['all',*FIXTURES],default='all')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    from sage.all import EllipticCurve, pari
    from sage.env import SAGE_VERSION
    names=tuple(FIXTURES) if args.fixture=='all' else (args.fixture,)
    for name in names:
        start=time.perf_counter()
        S=two_selmer(FIXTURES[name],proof=True)
        parent_seconds=time.perf_counter()-start
        start=time.perf_counter()
        B=cassels_tate_pairing(S)
        pairing_seconds=time.perf_counter()-start
        start=time.perf_counter()
        F=assemble_four_selmer_fibres(S,B)
        assembly_seconds=time.perf_counter()-start
        verified=False
        start=time.perf_counter()
        independent_dimension=independent_pairing_rank=None
        if args.verify:
            independent_dimension=int(EllipticCurve(list(FIXTURES[name])).global_minimal_model()
                                      .mwrank_curve(verbose=False).selmer_rank())
            independent_pairing_rank=int(pari.ellinit(list(FIXTURES[name])).ellrank()[2])
            verified=F.certified and F.verify_evidence()
            if (not verified or S.dimension!=independent_dimension or B.rank!=independent_pairing_rank
                    or F.cardinality!=EXPECTED[name]):
                raise ArithmeticError(f'global fibre verification failed: {name}')
        verification_seconds=time.perf_counter()-start
        record={**metadata(),'package_version':__version__,'sage':SAGE_VERSION,'fixture':name,
            'scope':'explicit global Selmer-4 fibres, projection and kernel; no global addition',
            'original_ainvs':vector(S.curves.original.ainvs),'selmer_dimension':S.dimension,
            'selmer_basis_labels':S.basis_labels,'pairing_matrix':B.matrix,'pairing_rank':B.rank,
            'image_basis':F.image_basis,'image_classes':F.image_classes,
            'kernel_dimension':F.kernel_dimension,'global_cardinality':F.cardinality,
            'predicted_cardinality':F.predicted_cardinality,'global_keys':[C.key for C in F.known_classes],
            'zero_class_key':F.zero_class.key,'complete':F.complete,'certified':F.certified,
            'verified':verified,'assumptions':F.assumptions,'unresolved':sum(len(f.lifts.unresolved) for f in F.fibres),
            'eclib_selmer_dimension':independent_dimension,'pari_pairing_rank':independent_pairing_rank,
            'max_fibres':None,'max_candidates':None,'max_nodes':None,'max_depth':None,
            'pruned':True,'group_operations_available':False,'parent_selmer_seconds':parent_seconds,
            'pairing_seconds':pairing_seconds,'assembly_seconds':assembly_seconds,
            'verification_seconds':verification_seconds,**geometry_record(F)}
        print(json.dumps(record,sort_keys=True),flush=True)


if __name__=='__main__':
    main()
