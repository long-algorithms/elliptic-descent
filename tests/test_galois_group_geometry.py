"""Portable exact certificates from actual PR5 covering models, not mocked Sage."""
from dataclasses import replace
from functools import lru_cache
import json
from pathlib import Path

import pytest
import sympy as sp

from elliptic_descent.galois_covering import (
    ExactGaloisField, evidence_from_record, evidence_record,
)
from elliptic_descent.four_selmer_group import FiniteClassLaw


@lru_cache(maxsize=1)
def record_and_evidence():
    r=json.loads((Path(__file__).parent/'data/group_zero_fibre.json').read_text())
    return r,evidence_from_record(r)


def test_real_covering_cocycles_reproduce_the_retained_torsion_orbits():
    r,E=record_and_evidence()
    assert E.verify()
    F,T=E.torsion.field,E.class_table
    assert F.degree==8 and T.verify() and len(T.keys)==2
    all_cocycles=tuple(tuple(c) for c in [item['cocycle'] for item in r['provenance']['relative_diagrams']])
    signatures=tuple(E.torsion.module.canonical(c) for c in all_cocycles)
    assert signatures[0]==signatures[2] and signatures[1]==signatures[3]
    assert signatures[0]!=signatures[1]
    for i,j in ((0,2),(1,3)):
        assert E.torsion.module.equivalence_witness(all_cocycles[i],all_cocycles[j]) is not None
    law=FiniteClassLaw(T.keys,tuple(tuple(T.add_indices(i,j) for j in range(2)) for i in range(2)),T.zero_index)
    assert law.verify() and law.decomposition[0]==(2,)
    roundtrip=evidence_from_record(json.loads(json.dumps(evidence_record(E))))
    assert roundtrip.verify() and roundtrip.class_table.cocycles==T.cocycles


def test_branch_origin_and_all_galois_conjugates_land_on_the_actual_four_torsion():
    _,E=record_and_evidence(); F=E.torsion.field
    assert E.torsion.verify()
    for C in E.coverings:
        Q=C.origin(E.torsion)
        assert C.verify(E.torsion) and C.evaluate(E.torsion,Q) is None
        for s in range(F.degree):
            point=C.evaluate(E.torsion,tuple(F.apply(s,v) for v in Q))
            assert point in E.torsion.point_indices
            assert E.torsion.curve.multiply(4,point) is None
    for v in E.torsion.root_values:
        assert F.unpack(F.pack(v))==v


def test_signed_trivialization_changes_are_checked_and_coboundaries_are_not_new_classes():
    _,E=record_and_evidence(); F=E.torsion.field; C=E.coverings[1]
    roots=C.roots_squared
    bad=replace(C,roots_squared=(F.pack(-F.unpack(roots[0])),*roots[1:]))
    assert not bad.verify(E.torsion)
    # Flipping two signs changes the origin by a 2-torsion translation.
    other=replace(C,roots_squared=(F.pack(-F.unpack(roots[0])),F.pack(-F.unpack(roots[1])),roots[2]))
    assert other.verify(E.torsion)
    assert E.torsion.module.equivalence_witness(C.cocycle(E.torsion),other.cocycle(E.torsion)) is not None
    assert not replace(C,norm_witness=-C.norm_witness).verify(E.torsion)
    assert not replace(C,xi=(0,0,0,0)).verify(E.torsion)
    assert not replace(C,jacobian_to_original=(1,0,0,0)).verify(E.torsion)


def test_field_and_torsion_data_cannot_omit_conjugates_or_fake_a_basis():
    _,E=record_and_evidence(); F=E.torsion.field
    assert not replace(F,automorphism_images=F.automorphism_images[:-1]).verify()
    assert not replace(F,generators=()).verify()
    assert not replace(F,automorphism_images=(F.automorphism_images[0],)*F.degree).verify()
    assert not ExactGaloisField((-1,0,1),((0,1),(0,-1)),(1,)).verify()
    with pytest.raises(ValueError): F.unpack((1,))
    assert not replace(E.torsion,basis=(E.torsion.basis[0],)*2).verify()
    assert not replace(E.torsion,roots=(E.torsion.roots[0],)*3).verify()
    assert not replace(E,coverings=(E.coverings[0],E.coverings[0])).verify()


def test_half_point_identity_symbolically_with_arbitrary_roots():
    x,w1,w2,w3=sp.symbols('x w1 w2 w3')
    roots=(x-w1*w1,x-w2*w2,x-w3*w3)
    a2=-sum(roots);a4=sum(roots[i]*roots[j] for i,j in ((0,1),(0,2),(1,2)));a6=-sp.prod(roots)
    pairs=w1*w2+w1*w3+w2*w3
    X=x+pairs; V=-w1*w2*w3+(w1+w2+w3)*pairs
    assert sp.expand(V*V-X**3-a2*X*X-a4*X-a6)==0
    slope=sp.cancel((3*X*X+2*a2*X+a4)/(2*V))
    double_x=sp.cancel(slope*slope-a2-2*X)
    double_y=sp.cancel(-V+slope*(X-double_x))
    assert sp.cancel(double_x-x)==0
    assert sp.cancel(double_y-w1*w2*w3)==0
