"""Origin torsors, contracted sums and equivariant (not just algebra) isomorphisms."""
import os
from dataclasses import replace
from functools import lru_cache
import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage=pytest.importorskip('sage.all',reason='requires SageMath runtime')

from elliptic_descent import construct_four_cover, MarkedQuarticCover, QuarticTwoCover, WeierstrassChange
from elliptic_descent.backends.sage_finite_torsors import OriginContext, contracted_sum, torsor_isomorphism

@lru_cache(maxsize=1)
def fixture():
    # Explicit zero-fibre transport of the unit cubic Kummer class for y^2=x^3-x.
    # A unit relative norm pair alone would NOT guarantee a trivial 4-cover.
    f=(0,4,24,44,24)
    D=construct_four_cover(f,(1,sp.Rational(23,4),11,6),norm_witness=sp.Rational(-1,96))
    Q=QuarticTwoCover(f)
    M=MarkedQuarticCover(Q,WeierstrassChange(Q.ainvs,(-2,-8,0,0)))
    assert M.target.ainvs==(0,0,0,-1,0)
    C=OriginContext(M.compose(D))
    A=C.torsor(M.compose(D))
    return C,A

def test_rational_origin_action_and_sum_with_zero():
    C,A=fixture()
    assert C.verify() and A.verify()
    identity=torsor_isomorphism(A,A)
    assert identity is not None and identity.verify(A,A)
    summed=contracted_sum(A,A)
    witness=torsor_isomorphism(summed,A)
    assert summed.verify() and witness is not None and witness.verify(summed,A)
    assert not replace(witness,coefficients=(0,)*16).verify(summed,A)
    assert not replace(A,denominator_inverse=A.tensor.ring.zero()).verify()
