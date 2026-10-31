"""Genuine Sage checks for model normalization, markings and relative descent."""
import os
from dataclasses import replace

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all',reason='requires SageMath runtime')

from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange
from elliptic_descent.cubic_conversion import IdentityTwoCover
from elliptic_descent.backends.sage_models import rational_curve_context


def coords(P):
    return None if P.is_zero() else tuple(sp.Rational(str(c)) for c in P.xy())


@pytest.mark.parametrize('a', [(-1,1),(0,0,1,-7,6),(1,2,3,4,5)])
@pytest.mark.parametrize('urst', [(2,3,4,5),(-2,-1,2,-3)])
def test_change_conventions_against_sage(a,urst):
    model=WeierstrassModel(a)
    E=sage.EllipticCurve([sage.QQ(str(c)) for c in model.ainvs])
    iso=WeierstrassIsomorphism(E,tuple(sage.QQ(c) for c in urst))
    change=WeierstrassChange(model,urst)
    assert change.target.ainvs == tuple(sp.Rational(str(c)) for c in iso.codomain().ainvs())
    assert change.source.discriminant == sp.Rational(str(E.discriminant()))
    assert change.target.discriminant == sp.Rational(str(iso.codomain().discriminant()))
    assert change.verify()
    inverse=~iso
    assert change.inverse().urst == tuple(sp.Rational(str(c)) for c in inverse.tuple())


def test_normalization_retains_nonminimal_original_model():
    change=WeierstrassChange((-1,1),(sp.Rational(1,6),2,-3,4))
    original=change.target
    context=rational_curve_context(original)
    E=sage.EllipticCurve([sage.QQ(str(c)) for c in original.ainvs])
    F=E.global_minimal_model()
    assert context.original==original
    assert context.working.ainvs == tuple(sp.Rational(str(c)) for c in F.ainvs())
    assert context.to_original.target==original and context.to_original.verify()
    marked=context.identity_cover()
    assert marked.target==original and marked.verify()
    fixed=rational_curve_context(E,minimal=False)
    assert fixed.original==fixed.working and fixed.to_original.urst==(1,0,0,0)
    with pytest.raises(ValueError):
        rational_curve_context(sage.EllipticCurve(sage.GF(101),[-1,1]))
    with pytest.raises(TypeError):
        rational_curve_context(original,minimal='true')


@pytest.mark.parametrize('a,P', [((-1,1),(0,1)),((0,0,1,-7,6),(0,2)),
                                ((0,-5,0,-16,80),(8,12))])
def test_zero_class_and_translated_markings_on_sage_multiples(a,P):
    I=IdentityTwoCover(a)
    E=sage.EllipticCurve([sage.QQ(str(c)) for c in I.curve.ainvs])
    point=E(*P)
    marking=replace(I.cover,translation=I.curve.point(P))
    assert I.verify() and marking.verify()
    for k in range(-3,6):
        Q=k*point
        source_point=I.to_quartic(coords(Q))
        assert I.from_quartic(source_point)==coords(Q)
        assert I.cover.evaluate(source_point)==coords(2*Q)
        assert marking.evaluate(source_point)==coords(2*Q+point)


def test_marked_relative_descent_does_not_change_the_original_curve():
    I=IdentityTwoCover((-1,1))
    original_change=WeierstrassChange(I.curve,(sp.Rational(1,2),1,-1,2))
    marking=I.cover.retarget(original_change)
    result=marking.descend(proof=True)
    assert result.complete and result.certified and result.verify_evidence()
    assert result.target==original_change.target
    E=sage.EllipticCurve([sage.QQ(str(c)) for c in I.curve.ainvs])
    assert len(result.coverings)==2**int(E.mwrank_curve(verbose=False).selmer_rank())
    assert all(D.target==result.target and D.verify() for D in result.coverings)
    partial=marking.descend(max_candidates=0)
    assert not partial.complete and not partial.certified and not partial.verify_evidence()
    assert partial.target==marking.target
    conditional=marking.descend(proof=False)
    assert conditional.complete and not conditional.certified
    assert conditional.verify_evidence()
