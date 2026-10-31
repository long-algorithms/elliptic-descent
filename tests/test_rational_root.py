"""Exact rational-root identifications, stabilizers and real-isolation regressions."""
from dataclasses import replace

import pytest
import sympy as sp

from elliptic_descent.etale import EtaleQuarticAlgebra
from elliptic_descent.equivalence import orientation_certificate, scalar_square_change
from elliptic_descent.geometry import construct_four_cover
from elliptic_descent.rational_root import RationalRootIdentification
from elliptic_descent.real import real_squareclass_solubility
from elliptic_descent.status import LocalStatus
from elliptic_descent.two_cover import X, Y, Z

FIXTURES = (
    (0, 1, -1, 0, 1),          # rational root zero, cubic field
    (-3, 4, -1, -3, 1),        # translated rational root
    (0, 1, -3, 0, 1),          # three real conjugates in the cubic factor
    (-1, 0, 0, 0, 1),          # imaginary quadratic component
    (2, 0, -3, 0, 1),          # two rational roots and a quadratic field
    (4, 0, -5, 0, 1),          # split, Jacobian has rational 4-torsion
    (0, -1, 1, 0, -1),         # negative leading coefficient
    (sp.Rational(1, 4), 0, sp.Rational(-5, 4), 0, 1),
)


def same_weighted_point(p, q):
    x, y, z = p
    a, b, c = q
    if x*c != a*z:
        return False
    if z and c:
        return y*c*c == b*z*z
    return x and a and y*a*a == b*x*x


@pytest.mark.parametrize('f', FIXTURES)
def test_each_rational_root_gives_the_marked_jacobian(f):
    algebra = EtaleQuarticAlgebra(f)
    for g in algebra.factors:
        if len(g) != 2:
            continue
        identification = RationalRootIdentification(f, -g[0])
        assert identification.check_marked_map()
        origin = (-g[0], 0, 1)
        assert identification.to_jacobian(origin) is None
        assert same_weighted_point(identification.from_jacobian(None), origin)
        # Every rational branch point is covered, not just the chosen origin.
        for h in algebra.factors:
            if len(h) == 2:
                P = (-h[0], 0, 1)
                assert same_weighted_point(identification.from_jacobian(
                    identification.to_jacobian(P)), P)


@pytest.mark.parametrize('f', FIXTURES)
def test_orientation_reversal_commutes_with_the_fixed_target(f):
    algebra = EtaleQuarticAlgebra(f)
    certificate = orientation_certificate(algebra)
    assert certificate.multiplicity == 1 and certificate.verify(f)
    v = certificate.reversing_multiplier
    assert algebra.multiply(v, v) == algebra.one and algebra.norm(v) == -1
    # If a is negative/nonsquare, put a in a rational component: Norm(xi)=a.
    index = algebra.degrees.index(1)
    xi = algebra.embed_factor_unit(index, (algebra.f[4],))
    D = construct_four_cover(f, xi, norm_witness=1)
    opposite = replace(D, norm_witness=-D.norm_witness)
    witness = certificate.identification(D)
    assert witness.verify(D, opposite)
    assert witness.coordinate_matrix(f).det() == -1
    changed, actual = scalar_square_change(D, 1, v)
    assert changed == opposite and actual == witness
    assert witness.verify(opposite, D)
    assert not replace(certificate, reversing_multiplier=algebra.one).verify(f)
    assert not replace(certificate, reversing_multiplier=algebra.idempotents[index]).verify(f)
    assert not replace(certificate, reversing_multiplier=None).verify(f)
    assert not replace(certificate, reversing_multiplier=(1.0, 0, 0, 0)).verify(f)


def test_all_linear_factor_choices_identify_the_same_two_orientations():
    f = (4, 0, -5, 0, 1)
    algebra = EtaleQuarticAlgebra(f)
    D = construct_four_cover(f, algebra.one)
    certificate = orientation_certificate(f)
    for i in range(4):
        flip = algebra.embed_factor_unit(i, (-1,))
        alternative = replace(certificate, reversing_multiplier=flip)
        assert alternative.verify(f)
        assert alternative.identification(D).verify(D, replace(D, norm_witness=-1))
    # Flipping two factors is a deck transformation, not an orientation reversal.
    v = algebra.multiply(algebra.embed_factor_unit(0, (-1,)), algebra.embed_factor_unit(1, (-1,)))
    assert not replace(certificate, reversing_multiplier=v).verify(f)


@pytest.mark.parametrize('f', [(1,-1,0,0,1), (6,0,-5,0,1)])
def test_no_accidental_orientation_collapse_without_a_rational_root(f):
    certificate = orientation_certificate(f)
    assert certificate.multiplicity == 2 and certificate.verify(f)
    D = construct_four_cover(f, (1,0,0,0))
    assert certificate.identification(D) is None
    assert not replace(certificate, reversing_multiplier=(1,0,0,0)).verify(f)
    with pytest.raises(ValueError):
        RationalRootIdentification(f)


@pytest.mark.parametrize('f', [(2,0,-3,0,1), (-12,10,4,-5,1)])
def test_touching_rational_root_isolation_is_closed_and_disjoint(f):
    # Second polynomial is (t-2)(t-3)(t^2-2).
    D = construct_four_cover(f, (1,0,0,0))
    result = real_squareclass_solubility(D)
    assert result.status == LocalStatus.SOLUBLE and result.verify(D)
    assert all(a[1] < b[0] for a,b in zip(result.intervals, result.intervals[1:]))
    polynomial = sp.Poly.from_list(list(reversed(f)), X, domain=sp.QQ)
    for a,b in result.intervals:
        assert a == b or (polynomial.eval(a) != 0 and polynomial.eval(b) != 0)


def test_real_decisions_are_unchanged_for_algebra_squares():
    f = (2,0,-3,0,1)
    algebra = EtaleQuarticAlgebra(f)
    v = (1,2,3,4)
    square = algebra.multiply(v,v)
    D = construct_four_cover(f, square, norm_witness=algebra.norm(v))
    result = real_squareclass_solubility(D)
    assert result.status == LocalStatus.SOLUBLE and result.verify(D)
    forged = replace(result, intervals=((-2,-1),)+result.intervals[1:])
    assert not forged.verify(D)


def test_point_roundtrips_at_infinity_and_rational_four_torsion():
    identification = RationalRootIdentification((4,0,-5,0,1), root=2)
    for P in ((1,1,0), (1,-1,0), (0,2,1), (0,-2,1), (1,0,1), (2,0,1)):
        assert same_weighted_point(identification.from_jacobian(identification.to_jacobian(P)), P)
        x,y,z = P
        scaled = (-3*x, 9*y, -3*z)
        assert identification.to_jacobian(P) == identification.to_jacobian(scaled)
    # E: y^2=(x-5)(x-4)(x+4); (8,12) doubles to (5,0).
    for P in (None, (8,12), (8,-12), (5,0), (4,0), (-4,0)):
        assert identification.to_jacobian(identification.from_jacobian(P)) == P
    Cpoint = identification.from_jacobian((8,12))
    assert identification.source.evaluate(Cpoint) == (5,0)


def test_invalid_identification_inputs_and_points():
    f = (4,0,-5,0,1)
    with pytest.raises(ValueError):
        RationalRootIdentification(f, 0)
    with pytest.raises(TypeError):
        RationalRootIdentification(f, 1.0)
    identification = RationalRootIdentification(f)
    for P in ((0,0,0), (0,1,0), (0,0,1), (1,2)):
        with pytest.raises(ValueError):
            identification.to_jacobian(P)
    with pytest.raises(TypeError):
        identification.to_jacobian((1.0,1,0))
    for P in ((0,0), (1,2,3)):
        with pytest.raises(ValueError):
            identification.from_jacobian(P)
    with pytest.raises(TypeError):
        identification.from_jacobian((8.0,12))


def test_seeded_rational_roots_preserve_marking_under_nonmonic_presentations():
    from random import Random
    rng = Random(20261031)
    t = sp.Symbol('t')
    checked = 0
    for _ in range(12):
        root = sp.Rational(rng.randrange(-3, 4), rng.randrange(1, 4))
        cubic = t**3+sum(rng.randrange(-3, 4)*t**i for i in range(3))
        scale = sp.Rational(rng.choice((-2, -1, 1, 2)), rng.randrange(1, 4))
        polynomial = sp.Poly(scale*(t-root)*cubic, t, domain=sp.QQ)
        if polynomial.discriminant() == 0:
            continue
        f = tuple(polynomial.nth(i) for i in range(5))
        identification = RationalRootIdentification(f, root)
        assert identification.check_marked_map()
        assert identification.to_jacobian((root, 0, 1)) is None
        certificate = orientation_certificate(f)
        assert certificate.verify(f) and certificate.multiplicity == 1
        checked += 1
    assert checked >= 8
