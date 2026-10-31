"""Group operations must agree with actual coverings, not just expected orders."""
import os
from dataclasses import replace
from functools import lru_cache
from itertools import product

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import two_selmer, cassels_tate_pairing, assemble_four_selmer_fibres
from elliptic_descent.four_selmer_group import four_selmer_group, SignedInverseWitness
from elliptic_descent.weierstrass import WeierstrassChange

FIXTURES = (
    ('rank_one', (0,0,0,-1,1), (4,)),
    ('split_two_torsion', (0,0,0,-1,0), (2,2)),
    ('cyclic_two_torsion', (0,0,0,1,0), (2,)),
    ('cyclic_four_torsion', (0,0,0,4,0), (4,)),
    ('full_two_and_four', (0,-5,0,-16,80), (2,4)),
    ('zero_selmer', (0,-1,1,-10,-20), ()),
    ('nontrivial_sha_571a1', (0,-1,1,-929,-10595), (2,2)),
    ('split_nontrivial_radical', (0,0,0,-289,0), (2,2,2,2)),
)


@lru_cache(maxsize=12)
def assembled(a):
    S = two_selmer(a)
    return assemble_four_selmer_fibres(S, cassels_tate_pairing(S))


@lru_cache(maxsize=12)
def group(a):
    return four_selmer_group(assembled(a))


@pytest.mark.parametrize('name,a,invariants', FIXTURES, ids=[f[0] for f in FIXTURES])
def test_group_arithmetic_on_every_constructed_class(name, a, invariants):
    G = group(a)
    F, S = G.fibres, G.fibres.selmer
    assert G.complete and G.certified and G.verify_evidence(), name
    assert G.invariant_factors == invariants
    assert G.cardinality == F.cardinality == sp.prod(invariants)
    assert {P.key for P in G.elements} == {P.key for P in F.known_classes}
    assert G.zero == F.zero_class
    assert tuple(G.order(P) for P in G.generators) == invariants
    for P in G.elements:
        assert P.covering.target == S.curves.original and P.covering.verify()
        assert G.element(G.coordinates(P)) == P
        assert G.add(P,G.zero) == P and G.add(P,G.negate(P)) == G.zero
        assert G.multiply(-1,P) == G.negate(P)
        assert G.multiply(G.order(P),P) == G.zero
        assert G.multiply(4,P) == G.zero
        assert G.multiply(2,P) == G.include_two_selmer(G.projection(P))
    for P,Q in product(G.elements, repeat=2):
        R = G.add(P,Q)
        assert G.add(P,Q) == G.add(Q,P)
        assert G.projection(R) == G.projection(P)^G.projection(Q)
        assert G.coordinates(R) == tuple((x+y)%n for x,y,n in
                                         zip(G.coordinates(P),G.coordinates(Q),invariants))
        for T in G.elements:
            assert G.add(G.add(P,Q),T) == G.add(P,G.add(Q,T))
    for a,b in product(range(S.cardinality),repeat=2):
        assert G.add(G.include_two_selmer(a),G.include_two_selmer(b)) == G.include_two_selmer(a^b)
    # Every row is an actual relation on the selected kernel/lift generators.
    generators = tuple(F.kernel_map.from_quotient(1<<i) for i in range(F.kernel_dimension))
    lookup = {P.key:P for P in F.known_classes}
    generators += tuple(lookup[key] for key in G.link.lift_keys)
    for row in G.relations:
        value = G.zero
        for coefficient,P in zip(row,generators):
            value = G.add(value,G.multiply(coefficient,P))
        assert value == G.zero


def test_cyclic_inverse_is_bound_to_the_signed_marking():
    G = group((0,0,0,-1,1))
    assert G.link.method == 'cyclic_signed_inverse'
    witness = G.link.inverse_witness
    assert witness.verify(G.fibres)
    P = next(P for P in G.elements if P.key == witness.source_key)
    assert G.negate(P).key == witness.target_key
    assert not replace(witness,target_key=witness.source_key).verify(G.fibres)
    other = four_selmer_group(G.fibres,lift_keys=(witness.target_key,))
    for P,Q in product(G.elements,repeat=2):
        assert G.add(P,Q) == other.add(P,Q)
    assert G.generators != other.generators  # Coordinates differ, geometric law does not.


def test_kernel_and_projection_shortcuts_are_exact_isomorphisms():
    for a,method in (((0,-1,1,-929,-10595),'kernel_inclusion_isomorphism'),
                     ((0,0,0,-1,0),'projection_isomorphism')):
        G = group(a)
        assert G.link.method == method and G.verify_evidence()
        assert all(G.order(P)<=2 for P in G.elements)
        assert G.link.context is None


def test_partial_conditional_and_wrong_class_inputs_fail_closed():
    F = assembled((0,0,0,-1,1))
    partial = assemble_four_selmer_fibres(F.selmer,F.pairing,max_fibres=0)
    with pytest.raises(ValueError): four_selmer_group(partial)
    G = group((0,0,0,-1,1))
    with pytest.raises(ValueError): G.add((99,0),G.zero)
    with pytest.raises(ValueError): G.element((4,))
    with pytest.raises(TypeError): G.multiply(1.0,G.zero)
    with pytest.raises(ValueError): four_selmer_group(F,lift_keys=((0,0),))
    with pytest.raises(ValueError): four_selmer_group(F,lift_keys=())
    corrupted = replace(G.link,normal_keys=(G.link.normal_keys[0],)*G.cardinality)
    bad = replace(G,link=corrupted)
    assert not bad.complete and not bad.certified and not bad.verify_evidence()
    with pytest.raises(ValueError): _=bad.invariant_factors
    S = two_selmer((0,0,0,-1,1),proof=False)
    conditional = four_selmer_group(assemble_four_selmer_fibres(S,cassels_tate_pairing(S),proof=True))
    assert conditional.complete and conditional.assumptions and not conditional.certified
    assert conditional.verify_evidence()


def test_original_nonminimal_model_retains_all_group_maps():
    F = assembled((0,0,0,-1,1))
    change = WeierstrassChange(F.selmer.curves.original,(sp.Rational(1,2),1,-1,3))
    S = two_selmer(change.target)
    moved = four_selmer_group(assemble_four_selmer_fibres(S,cassels_tate_pairing(S)))
    assert moved.certified and moved.invariant_factors == (4,) and moved.verify_evidence()
    assert all(P.covering.target==change.target for P in moved.elements)
    for P in moved.elements:
        assert moved.multiply(2,P)==moved.include_two_selmer(moved.projection(P))


def test_theta_relations_and_direct_sums_are_not_abstract_relabellings():
    G = group((0,-5,0,-16,80))
    assert G.link.method == 'exact_theta_cochains'
    C = G.link.context
    lookup = {P.key:values for P,values in zip(G.fibres.known_classes,G.link.cochains)}
    # Direct cochain multiplication does not use the extension carry formula.
    for P,Q in product(G.elements,repeat=2):
        R = G.add(P,Q)
        quotient = tuple(a*b/c for a,b,c in zip(lookup[P.key],lookup[Q.key],lookup[R.key]))
        witness = C.coboundary(quotient)
        assert witness is not None and witness.verify(C,quotient)
    # A deliberately altered geometric link must not report invariant factors.
    w = G.link.normal_witnesses[0]
    bad = replace(G,link=replace(G.link,normal_witnesses=(replace(w,origin_value=2*w.origin_value),
                                                        *G.link.normal_witnesses[1:])))
    assert not bad.complete and not bad.certified
    assert not bad.verify_evidence()
