"""Split torsion parameters retain exact classes while avoiding artificial primes."""
import pytest
import sympy as sp

from elliptic_descent.cubic_local import CubicCover
from elliptic_descent.cubic_arithmetic import kummer_value
from elliptic_descent.fibre_presentations import split_torsion_parameters
from elliptic_descent.selmer_quartic import quartic_presentation
from elliptic_descent.arithmetic import quartic_support
from elliptic_descent.weierstrass import WeierstrassChange


@pytest.mark.parametrize('delta,T', [
    ((-1,-1,sp.Rational(1,289)),(0,0)),
    ((17,sp.Rational(-16,17),sp.Rational(1,289)),(17,0)),
    ((-17,sp.Rational(-18,17),sp.Rational(1,289)),(-17,0)),
])
def test_congruent_17_torsion_presentations_have_only_actual_support(delta,T):
    D=CubicCover((0,0,0,-289,0),delta)
    point,change=split_torsion_parameters(D,T)
    C=quartic_presentation(D,conic_point=point,parameter_change=change)
    assert C.verify() and D.evaluate((*point,1))==T
    assert C.quartic.f[1]==C.quartic.f[3]==0
    assert quartic_support(C.quartic.f)==(2,17)
    assert C.marking.target==D.curve


def test_split_torsion_parameters_on_changed_original_models():
    change=WeierstrassChange((-1,0),(sp.Rational(1,3),2,-1,4))
    E=change.target
    for T in ((0,0),(1,0),(-1,0)):
        P=change.apply(T)
        D=CubicCover(E,kummer_value(E,P))
        point,M=split_torsion_parameters(D,P)
        C=quartic_presentation(D,conic_point=point,parameter_change=M)
        assert C.verify() and C.to_working.target==E
        assert D.evaluate((*point,1))==P


def test_invalid_or_wrong_kummer_classes_are_not_used_as_shortcuts():
    D=CubicCover((-1,0),(1,0,0))
    with pytest.raises(ValueError):split_torsion_parameters(D,None)
    with pytest.raises(ValueError):split_torsion_parameters(D,(0,0))
    with pytest.raises(ValueError):split_torsion_parameters(CubicCover((1,0),(1,0,0)),(0,0))


def test_same_support_smaller_coefficients_are_not_preferred():
    # The biquadratic presentation is smaller but makes the local lifting
    # tree much worse. Retain the original conic presentation on support ties.
    from elliptic_descent.conic_parametrization import ConicParametrization
    from elliptic_descent.selmer_quartic import QuarticPresentation
    from elliptic_descent.fibre_presentations import _score, _prefer
    D = CubicCover((0,1,0,-24,36), (-3,sp.Rational(1,6),sp.Rational(1,6)))
    original = QuarticPresentation(D, ConicParametrization(
        ((6,-84,186),(7,-74,217),(1,-14,55))))
    point, M = split_torsion_parameters(D,(2,0))
    smaller = quartic_presentation(D,conic_point=point,parameter_change=M)
    assert original.verify() and smaller.verify()
    assert _score(original)[0] == _score(smaller)[0] == 2
    assert _score(smaller)[1] < _score(original)[1]
    assert not _prefer(_score(smaller), _score(original))
    assert _prefer((1,100,(2,)), _score(original))
