"""Exact map-equivalence tests, including the special 2+2 stabilizer."""
from dataclasses import replace
from itertools import product

import pytest
import sympy as sp

from elliptic_descent.etale import EtaleQuarticAlgebra
from elliptic_descent.equivalence import (
    OrientationCertificate, orientation_certificate,
    ScalarSquareEquivalence, scalar_square_change,
)
from elliptic_descent.geometry import construct_four_cover
from elliptic_descent.relative import relative_four_descent


@pytest.mark.parametrize("f", [
    (1, -1, 0, 0, 1), (6, 0, -5, 0, 1), (16, 0, -10, 0, 1),
    (2, 0, 3, 0, 1), (-6, 0, 1, 0, 1), (4, 0, 5, 0, 1),
    (12, 0, -10, 0, 2),
])
def test_supported_orientation_hypotheses(f):
    proof = orientation_certificate(f)
    assert proof.verify(f) and proof.multiplicity == 2
    assert proof.factorization_type in ((4,), (2, 2))
    assert not replace(proof, factors=((0, 1),)*4).verify(f)


@pytest.mark.parametrize("f", [(-3, 4, -1, -3, 1), (2, 0, -3, 0, 1), (4, 0, -5, 0, 1)])
def test_rational_root_orientation_rule_has_a_negative_norm_witness(f):
    algebra = EtaleQuarticAlgebra(f)
    linear_index = algebra.degrees.index(1)
    flip = algebra.embed_factor_unit(linear_index, (-1,))
    assert algebra.multiply(flip, flip) == algebra.one
    assert algebra.norm(flip) == -1
    certificate = orientation_certificate(algebra)
    assert certificate.verify(f) and certificate.multiplicity == 1
    assert algebra.norm(certificate.reversing_multiplier) == -1


@pytest.mark.parametrize("f", [(6, 0, -5, 0, 1), (16, 0, -10, 0, 1), (2, 0, 3, 0, 1)])
def test_component_sign_stabilizers_preserve_orientation(f):
    algebra = EtaleQuarticAlgebra(f)
    source = construct_four_cover(f, (1, 0, 0, 0))
    for first, second in product((-1, 1), repeat=2):
        v = algebra.from_components(((first, 0), (second, 0)))
        target, witness = scalar_square_change(source, 1, v)
        assert algebra.norm(v) == 1
        assert target.xi == source.xi and target.norm_witness == source.norm_witness
        assert witness.verify(source, target)
        assert not witness.verify(source, replace(target, norm_witness=-target.norm_witness))


@pytest.mark.parametrize("f,components,c", [
    ((16, 0, -10, 0, 1), ((0, 1), (0, 2)), sp.Rational(1, 8)),
    ((4, 0, 5, 0, 1), ((0, 2), (0, 1)), sp.Rational(-1, 4)),
])
def test_nonsquare_scalar_stabilizer_with_isomorphic_factors(f, components, c):
    algebra = EtaleQuarticAlgebra(f)
    source = construct_four_cover(f, (1, 0, 0, 0))
    v = algebra.from_components(components)
    assert algebra.multiply(v, v) == algebra.scalar(1/c)
    assert c*c*algebra.norm(v) == 1
    target, witness = scalar_square_change(source, c, v)
    assert target.xi == source.xi and target.norm_witness == source.norm_witness
    assert witness.verify(source, target)


def test_signed_norm_in_a_nonstabilizing_equivalence_and_inverse():
    f = (6, 0, -5, 0, 1)
    algebra = EtaleQuarticAlgebra(f)
    source = construct_four_cover(f, (1, 0, 0, 0))
    v = algebra.embed_factor_unit(0, (0, 1))
    assert algebra.norm(v) < 0
    target, witness = scalar_square_change(source, sp.Rational(-2, 3), v)
    assert target.norm_witness < 0
    assert witness.verify(source, target)
    assert not witness.verify(source, replace(target, norm_witness=-target.norm_witness))
    recovered, inverse = scalar_square_change(target, sp.Rational(-3, 2), algebra.inverse(v))
    assert recovered == source
    assert inverse.verify(target, recovered)
    assert not replace(witness, scalar=1).verify(source, target)
    assert not replace(witness, multiplier=(0, 0, 0, 0)).verify(source, target)
    assert witness.coordinate_matrix(f).det() == algebra.norm(v)


def test_zero_divisors_and_inexact_equivalence_data_are_rejected():
    source = construct_four_cover((6, 0, -5, 0, 1), (1, 0, 0, 0))
    algebra = EtaleQuarticAlgebra(source.f)
    for c, v in [(0, algebra.one), (1, algebra.idempotents[0])]:
        with pytest.raises(ValueError):
            scalar_square_change(source, c, v)
        assert not ScalarSquareEquivalence(c, v).verify(source, source)
    with pytest.raises(TypeError):
        scalar_square_change(source, 1.0, algebra.one)
    assert not ScalarSquareEquivalence(1.0, algebra.one).verify(source, source)
    assert not OrientationCertificate(((1, 0, 1), (3, 0, 1))).verify(source.f)
