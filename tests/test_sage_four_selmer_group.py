"""Genuine Sage discovery/integration gates; unavailable Sage is not simulated."""
import os
from dataclasses import replace
from functools import lru_cache

import pytest

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE')=='1':
    import sage.all as sage
else:
    sage=pytest.importorskip('sage.all',reason='requires SageMath runtime')

from elliptic_descent import (
    two_selmer, cassels_tate_pairing, assemble_four_selmer_fibres, four_selmer_group,
    WeierstrassChange,
)

CASES=(
    ((0,0,0,-1,1),(4,)),
    ((0,0,0,-1,0),(2,2)),
    ((0,0,0,1,0),(2,)),
    ((0,0,0,4,0),(4,)),
    ((0,-5,0,-16,80),(2,4)),
    ((0,-1,1,-10,-20),()),
    ((0,-1,1,-929,-10595),(2,2)),
    ((0,0,0,-289,0),(2,2,2,2)),
)


@lru_cache(maxsize=16)
def fibres(a):
    S=two_selmer(a)
    return assemble_four_selmer_fibres(S,cassels_tate_pairing(S))


@pytest.mark.parametrize('a,expected',CASES,ids=[f'group_case_{i}' for i in range(len(CASES))])
def test_complete_group_on_the_actual_covering_classes(a,expected):
    F=fibres(a); G=F.group()
    assert G.complete and G.certified and G.verify_evidence()
    assert G.invariant_factors==expected and G.cardinality==F.cardinality
    assert G.zero.key==F.zero_class.key
    assert {v.key for v in G.elements}=={C.key for C in F.known_classes}
    for v in G.elements:
        assert v.covering.target==F.selmer.curves.original
        assert G.from_coordinates(v.coordinates)==v
        assert v-v==G.zero and 4*v==G.zero
        assert 2*v==G.include_two_selmer(G.projection(v))
        assert v.order() in (1,2,4)
        for w in G.elements:
            assert (v+w)-w==v and v+w==w+v
            assert G.projection(v+w)==G.projection(v)^G.projection(w)
            for z in G.elements:
                assert (v+w)+z==v+(w+z)
    for g,n,row in zip(G.generators,G.invariant_factors,G.relations):
        assert g.order()==n and n*g==G.zero
        assert len(row)==len(G.generators)
    for c in range(F.selmer.cardinality):
        for d in range(F.selmer.cardinality):
            assert G.include_two_selmer(c)+G.include_two_selmer(d)==G.include_two_selmer(c^d)
    if F.cardinality>4 and F.image_basis and F.kernel_dimension:
        assert G.method=='geometric_galois_cocycles'
        assert G.galois_evidence.verify()
    with pytest.raises(ValueError): G.element((999,0))
    with pytest.raises(TypeError): G.multiply(1.0,G.zero)


def test_forced_and_geometric_routes_agree_class_by_class():
    # Rational 2-torsion, not just the trivial class. The forced projection law
    # and the separately constructed splitting-field cocycles must agree.
    F=fibres((0,0,0,-1,0))
    direct=four_selmer_group(F,backend='galois')
    shortcut=F.group()
    assert direct.certified and direct.verify_evidence()
    assert direct.method=='geometric_galois_cocycles'
    assert direct.law==shortcut.law
    replay=F.group(evidence=direct.galois_evidence)
    assert replay.law==direct.law and replay.verify_evidence()
    with pytest.raises(RuntimeError): F.group(evidence=direct.galois_evidence,max_field_degree=1)
    with pytest.raises(ValueError): direct.zero+shortcut.zero
    bad=replace(direct.galois_evidence,coverings=tuple(reversed(direct.galois_evidence.coverings)))
    with pytest.raises(ValueError): F.group(evidence=bad)
    with pytest.raises(RuntimeError): four_selmer_group(F,backend='galois',max_field_degree=1)


def test_conditional_parent_partial_assembly_and_corrupt_operations():
    S=two_selmer((0,0,0,-1,1),proof=False)
    B=cassels_tate_pairing(S)
    F=assemble_four_selmer_fibres(S,B,proof=False)
    G=F.group()
    assert G.complete and G.assumptions and not G.certified
    assert G.verify_evidence()
    partial=assemble_four_selmer_fibres(S,B,max_fibres=0)
    with pytest.raises(ValueError): partial.group()
    law=replace(G.law,table=tuple((G.law.zero,)*G.cardinality for _ in G.elements))
    assert not replace(G,law=law).verify_evidence()
    assert not replace(G,method='unverified-count-fit').verify_evidence()
    with pytest.raises(ValueError): G.from_coordinates((0,0,0))


def test_nonminimal_target_and_changed_selmer_basis_keep_actual_maps():
    a=(0,0,0,-1,1)
    C=WeierstrassChange(a,(1,1,-1,2))
    S=two_selmer(C.target)
    G=assemble_four_selmer_fibres(S,cassels_tate_pairing(S)).group()
    assert G.certified and G.invariant_factors==(4,)
    assert all(g.covering.target==C.target for g in G.elements)
    S=two_selmer((0,0,0,-1,0))
    b0,b1=S.basis_labels
    T=S.with_basis((b0^b1,b0))
    G=assemble_four_selmer_fibres(S,cassels_tate_pairing(S)).group()
    H=assemble_four_selmer_fibres(T,cassels_tate_pairing(T)).group()
    # Here projection is an isomorphism, giving the natural class correspondence.
    mapping={g.key:next(h for h in H.elements if H.projection(h)==S.transport_to(T,G.projection(g)))
             for g in G.elements}
    for g in G.elements:
        for h in G.elements:
            assert mapping[(g+h).key]==mapping[g.key]+mapping[h.key]
