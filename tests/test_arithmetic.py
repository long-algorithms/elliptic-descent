from itertools import product
from random import Random

import pytest
import sympy as sp

from elliptic_descent.arithmetic import (
    FiniteNormSpace, affine_f2, prime, prime_support, quartic_support, square_class,
)
from elliptic_descent.geometry import _norm
from elliptic_descent.status import ProofStatus

F = (1, -1, 0, 0, 1)


def xor_image(columns, x):
    result = 0
    for i, c in enumerate(columns):
        if x >> i & 1:
            result ^= c
    return result


@pytest.mark.parametrize("n", range(7))
def test_affine_solve_against_exhaustive_search(n):
    rng = Random(700+n)
    for _ in range(40):
        columns = tuple(rng.randrange(32) for _ in range(n))
        target = rng.randrange(32)
        solution = affine_f2(columns, target)
        expected = {x for x in range(1 << n) if xor_image(columns, x) == target}
        assert set(solution.elements()) == expected
        assert solution.cardinality == len(expected)


def test_empty_and_inconsistent_affine_systems():
    assert list(affine_f2([], 0).elements()) == [0]
    assert list(affine_f2([], 1).elements()) == []
    assert set(affine_f2([0, 0, 0], 0).elements()) == set(range(8))
    with pytest.raises(ValueError):
        affine_f2([-1], 0)


@pytest.mark.parametrize("x", [1, -1, 2, -6, sp.Rational(-3, 8), sp.Rational(49, 25)])
def test_rational_square_classes(x):
    assert square_class(x, (2, 3)) == square_class(x*sp.Rational(121, 169), (2, 3))
    assert square_class(x*x, (2, 3)) == 0


def test_square_class_validation(monkeypatch):
    assert square_class(-6, (2, 3)) == 7
    assert square_class(sp.Rational(2, 3), (2, 3)) == 6
    for x in (0, 5):
        with pytest.raises(ValueError):
            square_class(x, (2, 3))
    with pytest.raises(TypeError):
        square_class(1.0, (2,))
    with pytest.raises(ValueError):
        square_class(1, (2, 2))
    # Exercise the explicit missing-backend path even inside the Sage job.
    import sys
    from elliptic_descent.arithmetic import _proven_prime
    _proven_prime.cache_clear()
    with monkeypatch.context() as m:
        m.setitem(sys.modules, "sage.all", None)
        with pytest.raises(NotImplementedError):
            prime(2**64+13)
    with pytest.raises(ValueError):
        prime(9)


def test_support_bounds_and_nonmonic_input():
    assert quartic_support(F) == (2, 229)
    assert prime_support((sp.Rational(-15, 14), 25)) == (2, 3, 5, 7)
    assert 7 in quartic_support((sp.Rational(1, 7), -1, 0, 0, 4))
    with pytest.raises(ValueError):
        prime_support((0,))


def test_norm_candidates_and_truncation():
    # Supplied generators do not claim to be a complete independent K(S,2) basis.
    space = FiniteNormSpace(F, ((-1, 0, 0, 0), (0, -1, 0, 0)), (2, 229))
    batch = space.enumerate()
    assert batch.norm_fibre_exhausted and batch.total_count == 4
    assert batch.arithmetic_proof == ProofStatus.UNVERIFIED
    assert not space.enumerate(limit=2).norm_fibre_exhausted
    assert not space.enumerate(limit=0).norm_fibre_exhausted
    for candidate in batch.candidates:
        assert _norm(candidate.xi, F) == candidate.norm_witness**2
    with pytest.raises(ValueError):
        space.enumerate(limit=-1)


def test_nontrivial_affine_norm_target_and_empty_norm_fibre():
    # Norm(theta)=constant/lc=2; lc=2, so theta (but not 1) satisfies the target.
    f = (4, 0, 0, 0, 2)
    space = FiniteNormSpace(f, ((0, 1, 0, 0),), (2,))
    assert list(space.norm_fibre().elements()) == [1]
    assert space.enumerate().candidates[0].norm_witness == 1
    empty = FiniteNormSpace(f, (), (2,)).enumerate()
    assert empty.norm_fibre_exhausted and empty.total_count == 0


def test_assumption_bookkeeping():
    conditional = FiniteNormSpace(F, (), (), ProofStatus.CONDITIONAL, ("GRH",))
    assert conditional.enumerate().assumptions == ("GRH",)
    with pytest.raises(ValueError):
        FiniteNormSpace(F, (), (), ProofStatus.CONDITIONAL)
    with pytest.raises(ValueError):
        FiniteNormSpace(F, (), (), ProofStatus.CERTIFIED, ("GRH",))
