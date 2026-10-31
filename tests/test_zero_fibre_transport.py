"""Zero-fibre reuse is an exact isomorphism, including its signed marked map."""
from dataclasses import replace
from types import SimpleNamespace

import pytest
import sympy as sp

from elliptic_descent.cubic_arithmetic import kummer_value, kummer_algebra
from elliptic_descent.cubic_local import CubicCover, cubic_local_solubility
from elliptic_descent.selmer_quartic import quartic_presentation
from elliptic_descent.zero_fibre import ZeroBranchData, transport_zero_cover, TransferredLocalEvidence
from elliptic_descent.fibre_actions import zero_fibre_kummer
from elliptic_descent.weierstrass import WeierstrassChange
from elliptic_descent.status import LocalStatus


@pytest.mark.parametrize('a,P', [((-1,1),None),((-1,1),(0,1)),((-1,0),(0,0)),
                                ((1,0),(0,0)),((4,0),(2,4)),
                                ((0,-5,0,-16,80),(8,12)),((0,0,1,-7,6),(0,2))])
def test_linear_transport_and_doubling_diagram(a,P):
    auxiliary=CubicCover(a,kummer_value(a,P))
    E=auxiliary.curve
    presentation=quartic_presentation(CubicCover(E,(1,0,0)),conic_point=(1,0,0))
    transport=transport_zero_cover(presentation,auxiliary)
    assert transport.verify() and transport.coordinate_matrix.det()!=0
    # Actual map-to-E comparison, not just equality of defining equations.
    if P is None:
        point=(1,0,0,0)
    else:
        parts=[(sp.S.One,)+(sp.S.Zero,)*(d-1) for d in auxiliary.algebra.degrees]
        for i,g in enumerate(auxiliary.algebra.factors):
            if len(g)==2 and P[0]==-g[0]: parts[i]=(sp.S.Zero,)
        point=(*auxiliary.algebra.from_components(parts),sp.S.One)
    relative_point=transport.from_auxiliary(point)
    assert transport.to_auxiliary(relative_point)==point
    equations,Cpoint=transport.model.evaluate(relative_point)
    assert not any(equations)
    assert presentation.marking.evaluate(Cpoint)==E.double(auxiliary.evaluate(point))
    # The derived kernel label recovers exactly this delta, not merely a square class.
    C=SimpleNamespace(ambient_label=0,quartic=presentation.quartic,marking=presentation.marking,
        selmer=SimpleNamespace(curves=SimpleNamespace(working=E,to_working=WeierstrassChange(E))))
    assert zero_fibre_kummer(C,transport.model,SimpleNamespace(algebra=kummer_algebra(E)))==auxiliary.delta
    local=cubic_local_solubility(auxiliary)
    assert local.status==LocalStatus.SOLUBLE
    evidence=TransferredLocalEvidence(transport,local)
    assert evidence.verify(transport.model)
    assert not replace(transport,model=replace(transport.model,norm_witness=-transport.model.norm_witness)).verify()
    assert not replace(transport,equation_matrix=((0,0),(0,0))).verify()
    assert not evidence.verify(replace(transport.model,norm_witness=-transport.model.norm_witness))


@pytest.mark.parametrize('delta', [(1,0,0),(-9827,-292,16),(-36815,-1096,60)])
def test_sha_zero_fibre_reuses_real_models_and_exact_local_evidence(delta):
    E=(0,-1,1,-929,-10595)
    auxiliary=CubicCover(E,delta)
    presentation=quartic_presentation(CubicCover(E,(1,0,0)),conic_point=(1,0,0))
    transport=transport_zero_cover(presentation,auxiliary)
    assert transport.verify()
    local=cubic_local_solubility(auxiliary)
    assert local.status==LocalStatus.SOLUBLE and local.verify(auxiliary)
    assert TransferredLocalEvidence(transport,local).verify(transport.model)


def test_nonminimal_models_and_different_zero_quartic_parameters():
    change=WeierstrassChange((-1,1),(sp.Rational(1,2),1,-1,3))
    E=change.target
    auxiliary=CubicCover(E,kummer_value(E,change.apply((0,1))))
    for parameters in (None,((1,1),(1,2))):
        P=quartic_presentation(CubicCover(E,(1,0,0)),conic_point=(1,0,0),parameter_change=parameters)
        assert transport_zero_cover(P,auxiliary).verify()
    with pytest.raises(ValueError):
        transport_zero_cover(P,CubicCover((-1,1),(1,0,0)))
