"""Actual global fibres, exact torsion quotients and kernel maps under Sage."""
import os
from dataclasses import replace
from functools import lru_cache
from itertools import product

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage=pytest.importorskip('sage.all',reason='requires SageMath runtime')

from elliptic_descent import two_selmer, WeierstrassChange
from elliptic_descent.cassels_tate import cassels_tate_pairing
from elliptic_descent.four_selmer_fibres import assemble_four_selmer_fibres

FIXTURES=(
    ('rank_one',(0,0,0,-1,1),4),
    ('split_two',(0,0,0,-1,0),4),
    ('cyclic_two',(0,0,0,1,0),2),
    ('cyclic_four',(0,0,0,4,0),4),
    ('full_two_and_four',(0,-5,0,-16,80),8),
    ('zero_selmer',(0,-1,1,-10,-20),1),
    ('sha_571a1',(0,-1,1,-929,-10595),4),
)


@lru_cache(maxsize=16)
def assembled(a):
    S=two_selmer(a,proof=True)
    B=cassels_tate_pairing(S)
    return assemble_four_selmer_fibres(S,B)


@pytest.mark.parametrize('name,a,size',FIXTURES)
def test_explicit_global_fibres_and_exact_sequence(name,a,size):
    F=assembled(a)
    assert F.complete and F.certified and F.verify_evidence(),name
    assert F.cardinality==size==F.predicted_cardinality
    S=F.selmer
    E=sage.EllipticCurve(list(a)).global_minimal_model()
    assert S.dimension==int(E.mwrank_curve(verbose=False).selmer_rank())
    assert F.pairing.rank==int(sage.pari.ellinit(list(a)).ellrank()[2])
    assert len(F.known_classes)==size
    assert {F.projection(C.key) for C in F.known_classes}==set(F.image_classes)
    assert len(set(C.key for C in F.known_classes))==size
    for C in F.known_classes:
        assert C.covering.target==S.curves.original and C.covering.verify()
    for f in F.fibres:
        assert len(f.models)==S.cardinality
        assert len(f.orbits)==1<<F.kernel_dimension
        for a,b in product(range(len(f.models)),repeat=2):
            evidence=f.equivalence(a,b)
            assert (evidence is not None)==any(a in o and b in o for o in f.orbits)
            if evidence is not None:
                assert evidence.verify() and evidence.coordinate_matrix.det()!=0
    I=F.kernel_map
    assert I.verify() and I.quotient.cardinality==1<<F.kernel_dimension
    for a,b in product(range(S.cardinality),repeat=2):
        same = I.apply(a).key==I.apply(b).key
        from elliptic_descent.quotient import normal_form
        assert same == (normal_form(a^b,F.torsion.image_basis)==0)
    assert all(F.projection(I.apply(c).key)==0 for c in range(S.cardinality))


def test_rational_four_torsion_changes_stabilizers_not_an_arbitrary_divisor():
    cyclic_two=assembled((0,0,0,1,0))
    cyclic_four=assembled((0,0,0,4,0))
    assert len(cyclic_two.torsion.entries)==len(cyclic_four.torsion.entries)==2
    assert cyclic_two.torsion.image_rank==1 and cyclic_four.torsion.image_rank==0
    assert all(len(o)==2 for f in cyclic_two.fibres for o in f.orbits)
    assert all(len(o)==1 for f in cyclic_four.fibres for o in f.orbits)
    mixed=assembled((0,-5,0,-16,80))
    assert len(mixed.torsion.entries)==4 and mixed.torsion.image_rank==1
    assert len(mixed.torsion.kernel)==2
    for f in mixed.fibres:
        for a,e in zip(f.actions,mixed.torsion.entries):
            assert (a.permutation==tuple(range(len(f.models)))) == (e.coordinates==0)


def test_pruned_and_unpruned_routes_agree_on_all_small_classes():
    pruned=assembled((0,0,0,-1,0))
    unpruned=assemble_four_selmer_fibres(pruned.selmer,pruned.pairing,prune=False)
    assert unpruned.certified and unpruned.verify_evidence()
    assert unpruned.cardinality==pruned.cardinality
    assert tuple((f.coordinates,len(f.orbits)) for f in unpruned.fibres)==tuple(
        (f.coordinates,len(f.orbits)) for f in pruned.fibres)
    negative=assembled((0,-1,1,-929,-10595))
    assert negative.image_classes==(0,) and len(negative.exclusions)==3
    assert all(e.verify(negative.pairing) for e in negative.exclusions)
    assert negative.kernel_map.quotient.cardinality==4


def test_rebased_parent_and_original_model_maps_preserve_kernel_and_projection():
    F=assembled((0,0,0,-1,0))
    B=F.pairing.with_basis(tuple(reversed(F.selmer.basis_labels)))
    other=assemble_four_selmer_fibres(B.selmer,B)
    assert other.certified and other.verify_evidence() and other.cardinality==F.cardinality
    for c in range(F.selmer.cardinality):
        d=B.selmer.coordinates_of_label(F.selmer.ambient_label(c))
        assert F.pairing.in_radical(c)==other.pairing.in_radical(d)
    E=assembled((0,0,0,-1,1))
    change=WeierstrassChange(E.selmer.curves.original,(sp.Rational(1,2),1,-1,3))
    S=two_selmer(change.target)
    moved=assemble_four_selmer_fibres(S,cassels_tate_pairing(S))
    assert moved.certified and moved.verify_evidence()
    assert moved.cardinality==E.cardinality
    assert all(C.covering.target==change.target for C in moved.known_classes)
    for c in range(E.selmer.cardinality):
        d=E.selmer.transport_to(S,c,change)
        assert moved.projection(moved.kernel_map.apply(d).key)==0


def test_partial_and_conditional_calculations_do_not_promote_global_claims():
    F=assembled((0,0,0,-1,1))
    for limits in ({'max_fibres':0},{'max_candidates':0},{'max_nodes':0}):
        partial=assemble_four_selmer_fibres(F.selmer,F.pairing,**limits)
        assert not partial.complete and not partial.certified and not partial.verify_evidence()
        assert partial.predicted_cardinality==4
        with pytest.raises(ValueError):_=partial.cardinality
    conditional=two_selmer((0,0,0,-1,1),proof=False)
    result=assemble_four_selmer_fibres(conditional,cassels_tate_pairing(conditional),proof=True)
    assert result.complete and result.verify_evidence() and result.assumptions and not result.certified
    assert result.cardinality==4
    with pytest.raises(ValueError):
        assemble_four_selmer_fibres(two_selmer((0,0,0,-1,1),max_candidates=0))
    with pytest.raises(ValueError):assemble_four_selmer_fibres(F.selmer,assembled((0,0,0,1,0)).pairing)


def test_corrupted_action_orbit_and_kernel_labels_fail_verification():
    F=assembled((0,0,0,-1,0)); f=F.fibres[0]
    duplicate=replace(f,orbits=((0,),(0,)))
    assert not duplicate.verify_structure(F.selmer,F.torsion)
    a=f.actions[1]
    bad=replace(a,permutation=tuple(range(len(f.models))))
    altered=replace(f,actions=(f.actions[0],bad,*f.actions[2:]))
    assert not altered.verify_structure(F.selmer,F.torsion)
    w=a.witnesses[0]
    wrong=replace(w,scalar_square=replace(w.scalar_square,scalar=2*w.scalar_square.scalar))
    assert not wrong.verify()
    badI=replace(F.kernel_map,labels=(0,)*len(f.models))
    assert not badI.verify()
    with pytest.raises(ValueError):F.projection((999,0))
    with pytest.raises(ValueError):F.projection((0,999))
    assert not replace(F,fibres=()).complete
