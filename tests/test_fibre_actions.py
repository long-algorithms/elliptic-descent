"""Exact global-equivalence geometry, independent split arithmetic and kernel labels."""
from dataclasses import replace
from itertools import product
from types import SimpleNamespace

import pytest
import sympy as sp

from elliptic_descent.cubic_arithmetic import kummer_algebra, rational_two_torsion, kummer_value
from elliptic_descent.cubic_local import CubicCover
from elliptic_descent.selmer_quartic import quartic_presentation
from elliptic_descent.etale import EtaleQuarticAlgebra
from elliptic_descent.equivalence import orientation_certificate
from elliptic_descent.geometry import construct_four_cover, U
from elliptic_descent.arithmetic import square_class, quartic_support
from elliptic_descent.local_complete import complete_local_screen
from elliptic_descent.fibre_actions import quartic_deck_action, zero_fibre_kummer
from elliptic_descent.four_selmer_fibres import action_orbits
from elliptic_descent.backends.sage_fibres import oriented_class_matcher
from elliptic_descent.rational_root import RationalRootIdentification
from elliptic_descent.weierstrass import WeierstrassChange
from elliptic_descent.status import LocalStatus
from elliptic_descent.two_cover import X,Y,Z


@pytest.mark.parametrize('a', [(-1,0),(1,0),(4,0),(0,-5,0,-16,80)])
def test_deck_action_is_the_specified_rational_two_torsion_translation(a):
    P = quartic_presentation(CubicCover(a,(1,0,0)),conic_point=(1,0,0))
    E = P.auxiliary.curve
    torsion = rational_two_torsion(E)
    actions = {T:quartic_deck_action(P,T) for T in torsion}
    I = RationalRootIdentification(P.quartic)
    # Use the explicit degree-one identification, not invariants or point counts.
    for T,g in actions.items():
        assert g.verify()
        for Q in torsion:
            J = P.to_working.inverse().apply(Q)
            point = I.from_jacobian(J)
            got = P.to_working.apply(I.to_jacobian(g.evaluate(point)))
            assert got == E.add(Q,T)
            assert P.marking.evaluate(g.evaluate(point)) == P.marking.evaluate(point)
        for R,h in actions.items():
            composite = replace(g,matrix=tuple(map(tuple,(sp.Matrix(g.matrix)*sp.Matrix(h.matrix)).tolist())),
                                y_scale=g.y_scale*h.y_scale)
            assert composite.same_transformation(actions[E.add(T,R)])
    A=EtaleQuarticAlgebra(P.quartic.f)
    xi=A.embed_factor_unit(A.degrees.index(1),(A.f[4],))
    model=construct_four_cover(A.f,xi,norm_witness=1)
    for g in actions.values():
        c=g.transport(model)
        assert c.verify()
        assert not replace(c,raw=replace(c.raw,norm_witness=-c.raw.norm_witness)).verify()
        assert not replace(g,y_scale=-g.y_scale).verify()
    with pytest.raises(ValueError): actions[None].evaluate((0,0,0))


def test_deck_action_on_a_nonzero_intermediate_class_and_nonminimal_model():
    D=CubicCover((-1,0),kummer_value((-1,0),(0,0)))
    A=D.algebra
    parts=[(1,),(0,),(1,)]  # root order +1,0,-1
    P=quartic_presentation(D,conic_point=A.from_components(parts))
    for T in rational_two_torsion(D.curve):
        g=quartic_deck_action(P,T)
        assert g.verify()
        point=P.from_auxiliary((*A.from_components(parts),1))
        assert P.marking.evaluate(g.evaluate(point)) == P.marking.evaluate(point)
    E=WeierstrassChange((-1,0),(sp.Rational(1,3),2,-1,4)).target
    P=quartic_presentation(CubicCover(E,(1,0,0)),conic_point=(1,0,0))
    assert all(quartic_deck_action(P,T).verify() for T in rational_two_torsion(E))


def test_independent_split_orbits_and_actual_kernel_kummer_classes():
    # A genuinely split algebra, with exact rational sign/valuation coordinates;
    # no Sage or imported class-group data and no arbitrary permutation labels.
    P=quartic_presentation(CubicCover((-1,0),(1,0,0)),conic_point=(1,0,0))
    f=P.quartic.f
    A=EtaleQuarticAlgebra(f)
    assert A.degrees==(1,1,1,1)
    support=quartic_support(f)
    # Enumerate Q(S,2)^4 norm classes up to the diagonal scalar by fixing the
    # first component to one, and determine the fourth from the norm condition.
    reps=tuple(sp.Rational(sign*sp.prod(p for j,p in enumerate(support) if mask>>j & 1))
               for sign in (-1,1) for mask in range(1<<len(support)))
    byclass={square_class(v,support):v for v in reps}
    models=[]
    for a,b in product(reps,repeat=2):
        d=byclass[square_class(f[4]*a*b,support)]
        xi=A.from_components(((1,),(a,),(b,),(d,)))
        D=construct_four_cover(f,xi)
        local=complete_local_screen(D)
        assert local.status!=LocalStatus.UNKNOWN and local.verify(D)
        if local.status==LocalStatus.SOLUBLE: models.append(D)
    assert len(models)==4
    torsion=rational_two_torsion(P.auxiliary.curve)
    actions=[quartic_deck_action(P,T) for T in torsion]
    transported=[[g.transport(D) for D in models] for g in actions]
    matcher=oriented_class_matcher(f,models,tuple(t.raw for row in transported for t in row))
    matches=[[matcher.match(t) for t in row] for row in transported]
    permutations=[tuple(i for i,_ in row) for row in matches]
    assert action_orbits(permutations,4)==((0,1,2,3),)
    for row in matches:
        for _,w in row: assert w.verify() and w.coordinate_matrix.det()!=0
    # Minimal exact interface data for the GEOMETRIC zero-fibre formula, not
    # a mocked Selmer result or a claim of certified number-field arithmetic.
    E=P.auxiliary.curve
    C=SimpleNamespace(ambient_label=0,quartic=P.quartic,marking=P.marking,
                      selmer=SimpleNamespace(curves=SimpleNamespace(
                          working=E,to_working=WeierstrassChange(E))))
    arithmetic=SimpleNamespace(algebra=kummer_algebra(E))
    key=lambda d:tuple(square_class(c[0],(2,)) for c in arithmetic.algebra.to_components(d))
    labels=[key(zero_fibre_kummer(C,D,arithmetic)) for D in models]
    torsion_labels=[key(kummer_value(E,T)) for T in torsion]
    assert set(labels)==set(torsion_labels)
    for permutation,delta in zip(permutations,torsion_labels):
        for i,j in enumerate(permutation):
            assert labels[j]==tuple(a^b for a,b in zip(labels[i],delta))


def test_orbit_validator_rejects_invalid_permutations():
    assert action_orbits(((0,1,2,3),(1,0,3,2)),4)==((0,1),(2,3))
    for p in ((0,0),(0,2),(0,True),(1,),(-1,0)):
        with pytest.raises(ValueError):action_orbits((p,),2)
