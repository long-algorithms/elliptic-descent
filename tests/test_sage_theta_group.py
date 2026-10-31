"""Exact algebraic theta data; no numerical recognition or assumed class law."""
import os
from dataclasses import replace
from functools import lru_cache

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import construct_four_cover, MarkedQuarticCover
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange
from elliptic_descent.backends.sage_theta_group import ThetaContext, flex_trivialization


@lru_cache(maxsize=4)
def context(ainvs):
    return ThetaContext(WeierstrassModel(ainvs))


@pytest.mark.parametrize('ainvs', [(0, 0, 0, -1, 0), (0, 0, 0, 0, 1), (0, 0, 0, -1, 1)])
def test_exact_theta_representation_characters_and_rational_descent(ainvs):
    C = context(ainvs)
    assert len(C.points) == len(C.matrices) == len(C.characters) == 16
    assert len(set(C.characters)) == 16
    values = tuple(sage.QQbar(1) for _ in range(16))
    witness = C.coboundary(values)
    assert witness is not None and witness.verify(C, values)
    for i in (1, 7, 15):
        # A character is a coboundary in theta-cochain coordinates.
        witness = C.coboundary(C.characters[i])
        assert witness is not None and witness.verify(C, C.characters[i])
        assert not replace(witness, origin_value=2).verify(C, C.characters[i])
    with pytest.raises(ValueError): C.coboundary((0,)*16)


@pytest.mark.parametrize('f,xi', [
    ((4, 0, -5, 0, 1), (1, 0, 0, 0)),
    ((6, 0, -5, 0, 1), (1, 0, 0, 0)),
    ((1, -1, 0, 0, 1), (0, -1, 0, 0)),
])
def test_flex_trivializations_preserve_signed_markings(f, xi):
    D = construct_four_cover(f, xi)
    marking = MarkedQuarticCover.canonical(f)
    T = flex_trivialization(marking.compose(D))
    assert T.verify()
    assert not replace(T, point=(sage.QQbar(0),)*4).verify()
    changed = WeierstrassChange(marking.target, (sp.Rational(1, 2), 1, -1, 3))
    moved = flex_trivialization(marking.retarget(changed).compose(D))
    assert moved.verify() and moved.covering.target == changed.target


def test_opposite_norm_witnesses_give_inverse_global_classes():
    f = (6, 0, -5, 0, 1)
    marking = MarkedQuarticCover.canonical(f)
    D = construct_four_cover(f, (1, 0, 0, 0))
    C = context(marking.target.ainvs)
    plus = flex_trivialization(marking.compose(D))
    minus = flex_trivialization(marking.compose(replace(D, norm_witness=-D.norm_witness)))
    a, b = C.cochain(plus), C.cochain(minus)
    product = tuple(x*y for x, y in zip(a, b))
    witness = C.coboundary(product)
    assert witness is not None and witness.verify(C, product)
