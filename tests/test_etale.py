"""Independent exact checks of CRT arithmetic in every quartic splitting type."""
from dataclasses import FrozenInstanceError
from random import Random

import pytest
import sympy as sp

from elliptic_descent import EtaleQuarticAlgebra
from elliptic_descent.geometry import _mul, _norm, construct_four_cover

T = sp.Symbol("t")
EXPRESSIONS = (
    T**4-T+1,
    (T-3)*(T**3-T+1),
    (T*T-2)*(T*T-3),
    (T*T-2)*(T*T-1),
    (T*T-1)*(T*T-4),
    (T*T-2)*(T*T-8),  # Isomorphic quadratic fields, distinct labelled factors.
    T*(T**3-T+1),      # A rational zero branch point.
    (2*T-1)*(3*T+2)*(T*T-2)/7,
)


def coefficients(expression):
    poly = sp.Poly(expression, T, domain=sp.QQ)
    return tuple(poly.nth(i) for i in range(5))


@pytest.mark.parametrize("expression", EXPRESSIONS)
def test_factorization_idempotents_and_basis_roundtrips(expression):
    A = EtaleQuarticAlgebra(coefficients(expression))
    assert A.verify_decomposition()
    assert sum(A.degrees) == 4
    assert sum(sp.Poly.from_list(list(reversed(g)), T).degree() for g in A.factors) == 4
    assert all(sp.Poly.from_list(list(reversed(g)), T).is_irreducible for g in A.factors)
    for i, e in enumerate(A.idempotents):
        parts = A.to_components(e)
        for j, part in enumerate(parts):
            assert part == (sp.Integer(i == j),)+(sp.S.Zero,)*(A.degrees[j]-1)
    with pytest.raises(FrozenInstanceError):
        A.f = (1, 0, 0, 0, 1)


@pytest.mark.parametrize("expression", EXPRESSIONS)
def test_arithmetic_against_independent_global_operations(expression):
    A = EtaleQuarticAlgebra(coefficients(expression))
    rng = Random(1357911)
    for _ in range(12):
        a, b = (tuple(sp.Rational(rng.randrange(-9, 10), rng.randrange(1, 8)) for _ in range(4))
                for _ in range(2))
        parts = A.to_components(a)
        assert A.from_components(parts) == a
        # Multiply using the previous fixed-degree reduction implementation.
        assert A.multiply(a, b) == _mul(a, b, A.f)
        assert A.norm(a) == _norm(a, A.f)
        assert A.norm(A.multiply(a, b)) == A.norm(a)*A.norm(b)
        assert A.norm(a) == sp.prod(A.component_norms(a))
        if A.is_unit(a):
            assert A.multiply(a, A.inverse(a)) == A.one
            assert A.norm(A.inverse(a))*A.norm(a) == 1
    assert A.norm(A.zero) == 0
    assert A.norm(A.one) == 1
    assert not A.is_unit(A.zero)
    with pytest.raises(ZeroDivisionError):
        A.inverse(A.zero)


@pytest.mark.parametrize("expression", EXPRESSIONS[1:])
def test_factorwise_unit_embedding_and_zero_divisors(expression):
    A = EtaleQuarticAlgebra(coefficients(expression))
    for i, degree in enumerate(A.degrees):
        g = (sp.Rational(3, 5),)+(sp.S.Zero,)*(degree-1)
        embedded = A.embed_factor_unit(i, g)
        assert A.is_unit(embedded)
        assert A.norm(embedded) == sp.Rational(3, 5)**degree
        for j, component in enumerate(A.to_components(embedded)):
            assert component == (g if i == j else (1,)+(0,)*(A.degrees[j]-1))
        assert not A.is_unit(A.idempotents[i])
        assert A.norm(A.idempotents[i]) == 0
        with pytest.raises(ZeroDivisionError):
            A.inverse(A.idempotents[i])
        # This is arithmetic groundwork for the later orientation PR, not an
        # assertion that the complete covering-equivalence problem is solved.
        flip = A.embed_factor_unit(i, (-1,)+(0,)*(degree-1))
        assert A.multiply(flip, flip) == A.one
        assert A.norm(flip) == (-1)**degree


def test_component_multiplication_and_norms_with_isomorphic_factors():
    A = EtaleQuarticAlgebra(coefficients((T*T-2)*(T*T-8)))
    assert len(A.factors) == 2 and A.degrees == (2, 2)
    for i in range(2):
        x = A.embed_factor_unit(i, (0, 1))
        assert A.to_components(x)[1-i] == (1, 0)
        assert A.component_norms(x)[1-i] == 1
        assert A.component_norms(x)[i] == A.factors[i][0]
    assert A.embed_factor_unit(0, (0, 1)) != A.embed_factor_unit(1, (0, 1))


@pytest.mark.parametrize("scale", [sp.Rational(2), sp.Rational(-3, 5), sp.Rational(1, 49)])
def test_nonmonic_geometry_is_preserved(scale):
    f = coefficients((T*T-2)*(T*T-3))
    A = EtaleQuarticAlgebra(f)
    B = EtaleQuarticAlgebra(tuple(scale*c for c in f))
    assert B.f[4] == scale
    assert A.factors == B.factors and A.idempotents == B.idempotents
    a = (1, 2, 3, 4)
    assert A.to_components(a) == B.to_components(a)
    assert A.norm(a) == B.norm(a)
    assert B.verify_decomposition()


def test_norm_witness_and_geometry_interoperability():
    for f in (coefficients((T*T-2)*(T*T-3)), coefficients((T*T-1)*(T*T-4))):
        A = EtaleQuarticAlgebra(f)
        # A square in the algebra has square norm; use it as a construction fixture.
        v = tuple(map(sp.Rational, (1, 2, 3, 4)))
        xi = A.multiply(v, v)
        D = construct_four_cover(f, xi, norm_witness=A.norm(v))
        assert D.check_construction()


def test_exact_input_and_shape_validation():
    for f in ((1, 0, 0, 1), (1, 0, 0, 0, 0), (1, 0, -2, 0, 1)):
        with pytest.raises(ValueError):
            EtaleQuarticAlgebra(f)
    with pytest.raises(TypeError):
        EtaleQuarticAlgebra((1.0, -1, 0, 0, 1))
    A = EtaleQuarticAlgebra((6, 0, -5, 0, 1))
    for value in ((1, 2, 3), (1, 2, 3, 4, 5)):
        with pytest.raises(ValueError):
            A.to_components(value)
    with pytest.raises(TypeError):
        A.to_components((1.2, 0, 0, 0))
    for parts in (((1,),), ((1, 0, 0), (2, 0)), ((1, 0), (2,))):
        with pytest.raises(ValueError):
            A.from_components(parts)
    with pytest.raises(TypeError):
        A.from_components(((1, 0), (1.0, 0)))
    with pytest.raises(ValueError):
        A.embed_factor_unit(0, (0, 0))
    for i in (-1, 2):
        with pytest.raises(IndexError):
            A.embed_factor_unit(i, (1, 0))
    for i in (True, 0.0):
        with pytest.raises(TypeError):
            A.embed_factor_unit(i, (1, 0))


def test_linear_factors_are_evaluation_not_formal_generators():
    A = EtaleQuarticAlgebra(coefficients((2*T-1)*(3*T+2)*(T*T-2)/7))
    element = (1, 2, 3, 4)
    for g, part in zip(A.factors, A.to_components(element)):
        if len(g) == 2:
            root = -g[0]
            assert part == (sum(c*root**i for i, c in enumerate(element)),)
