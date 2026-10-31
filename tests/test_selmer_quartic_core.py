"""Exact conic, inverse-chart, arithmetic-class and quartic-marking checks."""
from dataclasses import replace
from itertools import product
from random import Random

import pytest
import sympy as sp

from elliptic_descent.conic_parametrization import ConicParametrization, parametrize_conic
from elliptic_descent.cubic_arithmetic import kummer_value, cubic_multiply
from elliptic_descent.cubic_local import CubicCover
from elliptic_descent.geometry import U
from elliptic_descent.rational_root import RationalRootIdentification
from elliptic_descent.selmer_quartic import QuarticPresentation, quartic_presentation
from elliptic_descent.two_cover import X, Y, Z
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange, x, y
from elliptic_descent.marked_cover import function_equal


def same_projective(a, b):
    return any(a) and any(b) and all(a[i]*b[j] == a[j]*b[i]
                                    for i in range(len(a)) for j in range(len(a)))


def same_weighted(a, b):
    X1, Y1, Z1 = a
    X2, Y2, Z2 = b
    if X1*Z2 != X2*Z1:
        return False
    return Y1*Z2**2 == Y2*Z1**2 if Z1 and Z2 else Y1*X2**2 == Y2*X1**2


@pytest.mark.parametrize('q,P', [
    (U[0]*U[2]-U[1]**2, (1,0,0)),
    (U[0]**2+U[1]**2-U[2]**2, (0,1,1)),
    (2*U[0]**2+3*U[1]**2-5*U[2]**2, (1,1,1)),
])
def test_conic_inverse_has_no_missing_tangent_or_infinite_point(q,P):
    C = parametrize_conic(q,P)
    assert C.verify(q) and sp.Matrix(C.matrix).det() != 0
    for a,b in ((1,0),(0,1),(1,1),(-2,3),(sp.Rational(2,3),sp.Rational(-4,5))):
        v = tuple(p.subs({X:a,Z:b},simultaneous=True) for p in C.forms)
        for scalar in (1,-3,sp.Rational(2,7)):
            scaled = tuple(scalar*c for c in v)
            s,t,scale = C.inverse_coordinates(scaled)
            assert a*t == b*s
            assert tuple(p.subs({X:s,Z:t},simultaneous=True) for p in C.forms) == tuple(scale*c for c in scaled)
    # Also check the original point, including the tangent-line parameter.
    a,b,scale = C.inverse_coordinates(P)
    assert tuple(p.subs({X:a,Z:b},simultaneous=True) for p in C.forms) == tuple(scale*c for c in P)


@pytest.mark.parametrize('a', [(-1,1),(-1,0),(4,0),(0,0,1,-7,6),(1,2,3,4,5)])
def test_zero_class_quartics_and_signed_markings(a):
    E = WeierstrassModel(a)
    D = CubicCover(E,(1,0,0))
    C = quartic_presentation(D,conic_point=(1,0,0))
    assert C.verify() and C.quartic.f[4] != 0
    assert C.to_working.target == E and C.to_working.verify()
    assert C.marking.verify_supplied_map(C.map_working_xy)
    branch = C.from_auxiliary((1,0,0,0))
    assert C.marking.evaluate(branch) is None
    # This is [2] composed with an identification, not the degree-one identity.
    root = branch[0]/branch[2]  # Infinity was moved away from the branch divisor.
    identification = RationalRootIdentification(C.quartic, root)
    assert identification.check_marked_map()
    assert not C.marking.verify_supplied_map(tuple(q.subs(
        dict(zip((x,y),tuple(p.subs({X:x,Y:y,Z:1},simultaneous=True)
                           for p in identification.map_xy))),simultaneous=True)
        for q in C.to_working.maps))


@pytest.mark.parametrize('a,P', [((-1,1),(0,1)),((-1,1),(1,1)),
                                ((-1,0),(0,0)),((-1,0),(1,0)),
                                ((4,0),(2,4)),((0,-5,0,-16,80),(8,12))])
def test_kummer_class_point_roundtrips_and_original_map(a,P):
    E = WeierstrassModel(a)
    D = CubicCover(E,kummer_value(E,P))
    parts = [(sp.S.One,)+(sp.S.Zero,)*(d-1) for d in D.algebra.degrees]
    for i,g in enumerate(D.algebra.factors):
        if len(g)==2 and P[0]==-g[0]:
            parts[i]=(sp.S.Zero,)
    u = D.algebra.from_components(parts)
    C = quartic_presentation(D,conic_point=u)
    assert C.verify()
    A = (*u,sp.S.One)
    point = C.from_auxiliary(A)
    assert same_projective(C.to_auxiliary(point),A)
    assert C.marking.evaluate(point) == D.evaluate(A)
    for k in (-2,sp.Rational(2,3)):
        scaled = tuple(k*v for v in A)
        assert same_weighted(C.from_auxiliary(scaled),point)
    # Direct x-theta identity in the function field: no label-by-invariants shortcut.
    forms = tuple(p.subs({X:x,Z:1},simultaneous=True)/y for p in C.parametrization.forms)
    square = cubic_multiply(forms,forms,D.algebra.f)
    value = cubic_multiply(D.delta,square,D.algebra.f)
    assert function_equal(value[0],C.map_working_xy[0],C.quartic.f)
    assert function_equal(value[1],-1,C.quartic.f)
    assert function_equal(value[2],0,C.quartic.f)


def test_parameter_changes_and_square_rescaled_arithmetic_values():
    E=WeierstrassModel((-1,1))
    D=CubicCover(E,(-1,0,1))
    C=quartic_presentation(D,conic_point=(0,1,0))
    for change in (((0,1),(1,0)),((1,1),(1,2)),((sp.Rational(1,2),-1),(2,3))):
        other=quartic_presentation(D,conic_point=(0,1,0),parameter_change=change)
        assert other.verify() and other.quartic.f[4]
        assert other.auxiliary.delta == C.auxiliary.delta
    # A conic point from the same curve gives a different but verified parametrization.
    P=tuple(p.subs({X:1,Z:2},simultaneous=True) for p in C.parametrization.forms)
    other=quartic_presentation(D,conic_point=P)
    assert other.verify()
    v=(1,1,0)
    delta=D.algebra.multiply(D.delta,D.algebra.multiply(v,v))
    changed=CubicCover(E,delta)
    # v*u_new=u_old determines a rational conic point on the changed model.
    columns=[D.algebra.multiply(v,tuple(int(i==j) for j in range(3))) for i in range(3)]
    matrix=sp.Matrix.hstack(*(sp.Matrix(c) for c in columns))
    point=tuple(matrix.inv()*sp.Matrix((0,1,0)))
    assert quartic_presentation(changed,conic_point=point).verify()


def test_all_exceptional_point_evaluators_extend_the_maps():
    E=WeierstrassModel((0,-5,0,-16,80))
    C=quartic_presentation(CubicCover(E,(1,0,0)),conic_point=(1,0,0))
    branch=C.from_auxiliary((1,0,0,0))
    identification=RationalRootIdentification(C.quartic,branch[0]/branch[2])
    for P in (None,(5,0),(4,0),(-4,0),(8,12),(8,-12)):
        J=C.to_working.inverse().apply(P)
        Q=identification.from_jacobian(J)
        A=C.to_auxiliary(Q)
        assert same_weighted(C.from_auxiliary(A),Q)
        assert C.marking.evaluate(Q)==Dvalue(C,A)==E.double(P)


def Dvalue(C,P):
    return C.auxiliary.evaluate(P)


def test_invalid_conic_and_point_data_fail_instead_of_omitting_classes():
    u,v,w=U[:3]
    for q,P in ((u*u,(0,1,0)),(u*w-v*v,(0,0,0)),(u*w-v*v,(1,1,0))):
        with pytest.raises(ValueError):parametrize_conic(q,P)
    with pytest.raises(TypeError):parametrize_conic(u*w-v*v,(1.0,0,0))
    with pytest.raises(ValueError):ConicParametrization(((0,0,0),)*3)
    C=quartic_presentation(CubicCover((-1,1),(1,0,0)))
    with pytest.raises(ValueError):C.parametrization.reparameterize(((1,0),(0,0)))
    with pytest.raises(ValueError):C.to_auxiliary((0,0,0))
    with pytest.raises(ValueError):C.to_auxiliary((1,0,0))
    with pytest.raises(ValueError):C.from_auxiliary((0,0,0,0))
    with pytest.raises(ValueError):C.from_auxiliary((1,1,1,1))
    with pytest.raises(ValueError):C.parametrization.inverse_coordinates((0,0,0))
    with pytest.raises(ValueError):quartic_presentation(C.auxiliary,conic_point=(0,1,0))
    with pytest.raises(TypeError):quartic_presentation((1,2,3))
    bad=ConicParametrization(((1,0,0),(0,1,0),(0,0,1)))
    assert not bad.verify(C.auxiliary.q[2])
    with pytest.raises(ValueError):QuarticPresentation(C.auxiliary,bad)


def test_seeded_nonminimal_general_weierstrass_markings():
    rng=Random(50303)
    E=WeierstrassModel((-1,1))
    for _ in range(6):
        u=sp.Rational(rng.choice((-3,-2,1,2)),rng.randrange(1,4))
        change=WeierstrassChange(E,(u,*[rng.randrange(-2,3) for i in range(3)]))
        target=change.target
        D=CubicCover(target,(1,0,0))
        C=quartic_presentation(D,conic_point=(1,0,0))
        assert C.verify() and C.to_working.target==target
        assert C.to_working.target.ainvs[0:3] == target.ainvs[0:3]
