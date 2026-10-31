"""Exact model/marking tests independent of Sage and numerical approximations."""
from dataclasses import replace
from random import Random

import pytest
import sympy as sp

from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange, x, y
from elliptic_descent.cubic_conversion import CubicToQuartic, IdentityTwoCover
from elliptic_descent.marked_cover import MarkedQuarticCover, function_equal
from elliptic_descent.geometry import construct_four_cover, U
from elliptic_descent.two_cover import QuarticTwoCover


@pytest.mark.parametrize('ainvs', [(-1,1), (0,0,1,-7,6), (1,2,3,4,5), (0,-5,0,-16,80)])
@pytest.mark.parametrize('urst', [(1,0,0,0), (2,3,4,5), (-2,-1,2,-3),
                                  (sp.Rational(1,2),sp.Rational(1,3),-1,sp.Rational(2,5))])
def test_model_changes_inverse_discriminant_and_composition(ainvs,urst):
    E = WeierstrassModel(ainvs)
    change = WeierstrassChange(E,urst)
    assert change.verify() and change.inverse().verify()
    assert change.then(change.inverse()) == WeierstrassChange(E)
    second = WeierstrassChange(change.target,(-3,2,-1,4))
    composite = change.then(second)
    assert composite.target == second.target and composite.verify()
    for a,b in zip(composite.maps,second.maps):
        assert sp.cancel(a-b.subs(dict(zip((x,y),change.maps)),simultaneous=True)) == 0
    XX,YY,ZZ = sp.symbols('XX YY ZZ')
    weighted = change.apply_weighted((XX,YY,ZZ))
    for affine,coordinate,degree in zip(change.maps,weighted[:2],(2,3)):
        assert sp.cancel(affine.subs({x:XX/ZZ**2,y:YY/ZZ**3},simultaneous=True)
                         -coordinate/weighted[2]**degree) == 0


def test_group_law_and_special_points():
    E = WeierstrassModel((-1,1)); P=(0,1)
    Q=E.double(P)
    assert Q == (sp.Rational(1,4),sp.Rational(-7,8))
    assert E.add(P,E.negate(P)) is None
    assert E.add(None,P) == E.add(P,None) == P
    assert E.double(None) is None
    change = WeierstrassChange(E,(-2,1,2,-3))
    assert change.target.double(change.apply(P)) == change.apply(Q)
    assert change.inverse().apply(change.apply(P)) == P
    assert change.apply(None) is None
    T = WeierstrassModel((-1,0))
    assert T.double((0,0)) is None
    # Nonzero a1,a3: the vertical inverse is not (x,-y).
    M=change.target; A=change.apply(P)
    assert M.add(A,M.negate(A)) is None


@pytest.mark.parametrize('f', [(1,-1,0,1), (0,-1,0,1), (1,2,3,4),
                              (sp.Rational(-2,3),1,sp.Rational(3,5),-2)])
def test_cubic_conversion_identity_and_infinite_branch(f):
    conversion=CubicToQuartic(f)
    assert conversion.verify() and conversion.quartic.f[4]
    assert conversion.to_quartic((1,0,0)) == (0,0,1)
    assert conversion.from_quartic((0,0,1)) == (1,0,0)
    assert conversion.quartic.check_map()
    if f == (0,-1,0,1):
        assert conversion.shift == 2  # Both shifts 0 and 1 are roots.


def test_cubic_conversion_affine_and_infinity_roundtrips():
    conversion=CubicToQuartic((1,-1,0,1),0)
    # x=shift goes to a quartic point at infinity rather than dividing by zero.
    assert conversion.to_quartic((0,1,1)) == (1,1,0)
    for P in [(0,1,1),(0,-1,1),(1,1,1),(1,0,0)]:
        assert conversion.from_quartic(conversion.to_quartic(P)) == P
    phi=conversion.transport_map((x,y))
    assert phi == (1/x,y/x**2)
    assert function_equal(phi[1]**2-(phi[0]**3-phi[0]+1),0,conversion.quartic.f)


@pytest.mark.parametrize('ainvs', [(-1,1),(-1,0),(0,0,1,-7,6),(1,2,3,4,5)])
def test_zero_selmer_class_is_doubling_not_degree_one(ainvs):
    I=IdentityTwoCover(ainvs)
    assert I.verify() and I.cover.target == WeierstrassModel(ainvs)
    assert I.from_quartic(I.to_quartic(None)) is None
    assert I.cover.evaluate(I.to_quartic(None)) is None
    assert not I.cover.verify_supplied_map(I.identification_maps)
    # The supplied cubic formula for [2] is transported to exactly this marking.
    a1,_,a3,_,_=I.curve.ainvs
    cubic_map=tuple(q.subs({x:x,y:(y-a1*x-a3)/2},simultaneous=True)
                    for q in I.curve.duplication_maps)
    C=I.conversion.marked_map(cubic_map,I.cover.change)
    assert C == I.cover
    with pytest.raises(ValueError):
        I.conversion.marked_map((x,(y-a1*x-a3)/2),I.cover.change)


def test_nonminimal_zero_cover_and_infinite_points():
    E=WeierstrassModel((-1,1)); change=WeierstrassChange(E,(sp.Rational(1,2),3,1,-2))
    I=IdentityTwoCover(change.target)
    assert I.verify()
    for P in [None,(0,1),(0,-1),E.double((0,1))]:
        Q=change.apply(P)
        assert I.from_quartic(I.to_quartic(Q)) == Q
        assert I.cover.evaluate(I.to_quartic(Q)) == I.curve.double(Q)
    assert IdentityTwoCover((-1,0)).cover.evaluate(IdentityTwoCover((-1,0)).to_quartic((0,0))) is None
    # Two rational 4-torsion signs double to the same nonzero 2-torsion point.
    I=IdentityTwoCover((0,-5,0,-16,80))
    for P in [(8,12),(8,-12)]:
        assert I.cover.evaluate(I.to_quartic(P)) == (5,0)


def test_translations_retargeting_and_full_marking_comparison():
    I=IdentityTwoCover((-1,1)); E=I.curve; T=(0,1)
    marked=replace(I.cover,translation=T)
    assert marked.verify()
    assert marked.evaluate(I.to_quartic(None)) == T
    assert not marked.verify_supplied_map(I.cover.map_xy)
    change=WeierstrassChange(E,(2,-1,1,3))
    new=marked.retarget(change)
    assert new.verify() and new.target==change.target
    for P in [None,T,E.double(T)]:
        assert new.evaluate(I.to_quartic(P)) == change.apply(marked.evaluate(I.to_quartic(P)))
    expected=tuple(q.subs(dict(zip((x,y),marked.map_xy)),simultaneous=True) for q in change.maps)
    assert new.verify_supplied_map(expected)
    # Constant maps can satisfy the target equation but are not marked 2-covers.
    assert E.equation(*T)==0
    assert not marked.verify_supplied_map(tuple(map(sp.Rational,T)))


def test_translation_exceptional_points_extend_without_denominator_errors():
    # On y^2=(x-5)(x-4)(x+4), P has order 4, and 2P=T has order 2.
    I=IdentityTwoCover((0,-5,0,-16,80)); T=(5,0)
    marked=replace(I.cover,translation=T)
    assert marked.verify()
    assert marked.evaluate(I.to_quartic((8,12))) is None
    assert marked.evaluate(I.to_quartic(None))==T
    assert marked.evaluate(I.to_quartic(T))==T


def test_composite_four_cover_keeps_original_target_and_checks_points():
    source=QuarticTwoCover((1,-1,0,0,1))
    change=WeierstrassChange(source.ainvs,(2,1,-1,3))
    marking=MarkedQuarticCover(source,change)
    D=construct_four_cover(source.f,(0,-1,0,0))
    marked=marking.compose(D)
    assert marked.verify() and marked.target==change.target
    assert marked.evaluate((1,0,0,0)) == change.apply(source.evaluate(D.evaluate((1,0,0,0))[1]))
    # Do not expand the degree-16 composite into huge multivariate fractions.
    # Stage identities are exact; independently check evaluation of the exported expressions.
    P=(1,0,0,0)
    assert tuple(sp.cancel(q.subs(dict(zip(U,P)),simultaneous=True)) for q in marked.map_xy) == marked.evaluate(P)
    with pytest.raises(ValueError):
        marked.evaluate((1,1,1,1))
    with pytest.raises(ValueError):
        marking.compose(construct_four_cover((4,0,-5,0,1),(1,0,0,0)))


def test_model_and_function_validation():
    for a in [(0,0),(1,2,3)]:
        with pytest.raises(ValueError): WeierstrassModel(a)
    with pytest.raises(TypeError): WeierstrassModel((1.0,2))
    E=WeierstrassModel((-1,1))
    for urst in [(0,0,0,0),(1,2,3)]:
        with pytest.raises(ValueError): WeierstrassChange(E,urst)
    with pytest.raises(ValueError): WeierstrassChange(E).then(WeierstrassChange((-2,1)))
    for P in [(1,2),(1,2,3)]:
        with pytest.raises(ValueError): E.point(P)
    with pytest.raises(ValueError): CubicToQuartic((0,-1,0,1),1)
    for f in [(1,2,3),(0,0,0,1),(1,2,3,0)]:
        with pytest.raises(ValueError): CubicToQuartic(f)
    for q in [sp.sqrt(x),sp.Symbol('other'),y/(y*y-(x**3-x+1)),sp.Float(1.5)*x]:
        assert not function_equal(q,q,(1,-1,0,1))
        with pytest.raises(ValueError): CubicToQuartic((1,-1,0,1)).transport_map((q,y))
    I=IdentityTwoCover(E)
    with pytest.raises(ValueError): MarkedQuarticCover(I.cover.source,WeierstrassChange(E))
    with pytest.raises(ValueError): replace(I.cover,translation=(1,2))


def test_seeded_model_changes_and_cubic_presentations():
    rng=Random(20261101)
    checked=0
    for _ in range(12):
        a=tuple(rng.randrange(-3,4) for i in range(5))
        try: E=WeierstrassModel(a)
        except ValueError: continue
        u=sp.Rational(rng.choice([-2,-1,1,2]),rng.randrange(1,4))
        change=WeierstrassChange(E,(u,*[rng.randrange(-2,3) for i in range(3)]))
        assert change.verify() and change.inverse().verify()
        assert IdentityTwoCover(change.target).verify()
        checked+=1
    assert checked>=10


def test_imported_markings_require_exact_targets_and_maps():
    from elliptic_descent.backends.sage import PariTwoCover
    from elliptic_descent.backends.sage_models import mark_pari_two_cover
    I=IdentityTwoCover((-1,1))
    a1,_,a3,_,_=I.curve.ainvs
    cubic_maps=tuple(q.subs({x:x,y:(y-a1*x-a3)/2},simultaneous=True)
                     for q in I.curve.duplication_maps)
    record=PariTwoCover(I.conversion.f,I.curve.ainvs,cubic_maps)
    marked,conversion=mark_pari_two_cover(record,I.cover.change)
    assert marked==I.cover and conversion==I.conversion
    record=PariTwoCover(I.cover.source.f,I.curve.ainvs,I.cover.map_xy)
    marked,conversion=mark_pari_two_cover(record,I.cover.change)
    assert marked==I.cover and conversion is None
    with pytest.raises(ValueError):
        mark_pari_two_cover(replace(record,ainvs=(0,0,0,-2,1)),I.cover.change)
    with pytest.raises(ValueError):
        mark_pari_two_cover(replace(record,map_xy=(x,y)),I.cover.change)
    with pytest.raises(ValueError):
        mark_pari_two_cover(record,I.cover.change,shift=0)
