"""Group arithmetic is tested separately from geometric Selmer identification."""
from itertools import product
from random import Random

import pytest

from elliptic_descent.exponent_four import ExponentFourPresentation


@pytest.mark.parametrize('k,doubles,invariants', [
    (0, (), ()), (0, (0, 0), (2, 2)), (1, (), (2,)),
    (1, (1,), (4,)), (2, (1,), (2, 4)), (2, (1, 2), (4, 4)),
    (1, (1, 1), (2, 4)), (2, (3, 3, 0), (2, 2, 2, 4)),
])
def test_orders_invariant_coordinates_and_all_group_laws(k, doubles, invariants):
    P = ExponentFourPresentation(k, doubles)
    assert P.invariant_factors == invariants
    values = tuple(product(range(1 << k), range(1 << len(doubles))))
    assert P.cardinality == len(values)
    recovered = set()
    for a in values:
        assert P.add(a, (0, 0)) == a
        assert P.add(a, P.negate(a)) == (0, 0)
        assert P.multiply(4, a) == (0, 0)
        assert P.multiply(-1, a) == P.negate(a)
        assert P.multiply(P.order(a), a) == (0, 0)
        coordinates = P.to_invariants(a)
        assert P.from_invariants(coordinates) == a
        recovered.add(coordinates)
        for b in values:
            assert P.add(a, b) == P.add(b, a)
            expected = tuple((x+y) % n for x, y, n in zip(coordinates, P.to_invariants(b), invariants))
            assert P.to_invariants(P.add(a, b)) == expected
            for c in values:
                assert P.add(P.add(a, b), c) == P.add(a, P.add(b, c))
    assert len(recovered) == P.cardinality


def test_seeded_doubling_matrices_and_carry_relations():
    rng = Random(50601)
    for _ in range(100):
        k, r = rng.randrange(5), rng.randrange(5)
        P = ExponentFourPresentation(k, tuple(rng.randrange(1 << k) for j in range(r)))
        for j, d in enumerate(P.doubles):
            assert P.add((0, 1 << j), (0, 1 << j)) == (d, 0)
        for _ in range(10):
            a = (rng.randrange(1 << k), rng.randrange(1 << r))
            assert P.from_invariants(P.to_invariants(a)) == a


def test_invalid_coordinates_are_not_silently_reduced():
    with pytest.raises(ValueError): ExponentFourPresentation(1, (2,))
    P = ExponentFourPresentation(1, (1,))
    for bad in ((2, 0), (0, 2), (-1, 0), (0, 0, 0)):
        with pytest.raises(ValueError): P.point(bad)
    with pytest.raises(TypeError): P.add((0.0, 0), (0, 0))
    for bad in ((4,), (-1,), (), (0, 0)):
        with pytest.raises(ValueError): P.from_invariants(bad)
    with pytest.raises(TypeError): P.multiply(1.0, (0, 1))
