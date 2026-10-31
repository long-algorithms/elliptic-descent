from dataclasses import replace
from itertools import product
from random import Random

import pytest
import sympy as sp

from elliptic_descent.geometry import U, construct_four_cover
from elliptic_descent.local import (
    QuadraticSystem, bad_primes, everywhere_local_solubility,
    padic_solubility, real_solubility,
)
from elliptic_descent.status import LocalStatus

u, v, w, z = U
POSITIVE = (u*u+v*v+w*w+z*z, u*u+2*v*v+3*w*w+4*z*z)
SINGULAR_REDUCTION = (u*u-17*v*v, w*w-17*z*z)


def test_two_adic_obstruction_with_solutions_mod_two():
    first = padic_solubility(POSITIVE, 2, max_depth=1)
    assert first.status == LocalStatus.UNKNOWN
    result = padic_solubility(POSITIVE, 2)
    assert result.status == LocalStatus.INSOLUBLE and result.depth == 2
    assert result.verify(POSITIVE)
    assert not result.verify(POSITIVE, max_nodes=1)
    forged = replace(result, certificate=replace(result.certificate, depth=1))
    assert not forged.verify(POSITIVE)


def test_two_adic_singular_hensel_lift():
    assert padic_solubility(SINGULAR_REDUCTION, 2, max_depth=1).status == LocalStatus.UNKNOWN
    result = padic_solubility(SINGULAR_REDUCTION, 2)
    assert result.status == LocalStatus.SOLUBLE
    assert result.certificate.kind == "hensel" and result.depth > 1
    assert result.verify(SINGULAR_REDUCTION)


@pytest.mark.parametrize("p", [3, 5, 7, 11])
def test_odd_prime_hensel_witnesses(p):
    result = padic_solubility(POSITIVE, p)
    assert result.status == LocalStatus.SOLUBLE
    assert result.verify(POSITIVE)
    assert not replace(result, status=LocalStatus.INSOLUBLE).verify(POSITIVE)
    forged = replace(result, certificate=replace(result.certificate, point=(1, 0, 0, 0)))
    assert not forged.verify(POSITIVE)


def test_budget_and_large_prime_never_mean_insoluble():
    for kwargs in ({"max_nodes": 0}, {"max_depth": 0}, {"max_nodes": 1}):
        result = padic_solubility(POSITIVE, 2, **kwargs)
        assert result.status == LocalStatus.UNKNOWN and result.certificate is None
    result = padic_solubility(POSITIVE, 1_000_000_007, max_nodes=1)
    assert result.status == LocalStatus.UNKNOWN and result.nodes == 1
    for p in (0, 1, 4, 9):
        with pytest.raises(ValueError):
            padic_solubility(POSITIVE, p)


def test_all_projective_charts_and_rational_scaling():
    equations = (u*u+v*v+w*w, u*u+2*v*v+3*w*w)
    result = padic_solubility(equations, 2)
    assert result.status == LocalStatus.SOLUBLE
    assert result.certificate.point == (0, 0, 0, 1)
    assert result.verify(equations)
    scaled = (sp.Rational(7, 12)*equations[0], -sp.Rational(5, 18)*equations[1])
    assert result.verify(scaled)


def test_real_definite_obstruction_and_unknown():
    result = real_solubility(POSITIVE)
    assert result.status == LocalStatus.INSOLUBLE and result.verify(POSITIVE)
    assert real_solubility(SINGULAR_REDUCTION).status == LocalStatus.UNKNOWN
    forged = replace(result, certificate=replace(result.certificate, pencil=(0, 0)))
    assert not forged.verify(POSITIVE)


def test_real_rational_witness_and_global_conclusion():
    D = construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0))
    result = everywhere_local_solubility(D, point=(1, 0, 0, 0))
    assert result.status == LocalStatus.SOLUBLE
    assert all(r.verify(D) for r in result.reports)
    assert everywhere_local_solubility(POSITIVE).status == LocalStatus.INSOLUBLE
    with pytest.raises(ValueError):
        real_solubility(POSITIVE, point=(1, 0, 0, 0))


def test_bad_prime_bound_and_degenerate_pencil():
    assert bad_primes(POSITIVE).primes == (2, 3)
    D = construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0))
    assert bad_primes(D).primes == (2, 229)
    with pytest.raises(NotImplementedError):
        bad_primes(SINGULAR_REDUCTION)


def test_quadratic_input_validation():
    for equations in ((u+1, v*v), (0, v*v), (u**3, v*v), (u*u,)):
        with pytest.raises(ValueError):
            QuadraticSystem.from_quadrics(equations)
    with pytest.raises(TypeError):
        QuadraticSystem.from_quadrics((1.2*u*u, v*v))


def test_depth_one_against_independent_primitive_enumeration():
    rng = Random(112233)
    for p in (2, 3, 5):
        for _ in range(12):
            q = tuple(sum(rng.randrange(-3, 4)*x*x for x in U) for _ in range(2))
            if any(expr == 0 for expr in q):
                continue
            # Independent enumeration: all primitive vectors, no chart normalization.
            system = QuadraticSystem.from_quadrics(q)
            any_point = any(any(a) and all(value % p == 0 for value in system.values(a))
                            for a in product(range(p), repeat=4))
            result = padic_solubility(q, p, max_depth=1)
            assert (result.status == LocalStatus.INSOLUBLE) == (not any_point)
            if result.status != LocalStatus.UNKNOWN:
                assert result.verify(q)
