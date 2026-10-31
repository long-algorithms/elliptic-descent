"""Integral quadric equations can have a rational target-equation residual."""
import pytest
import sympy as sp

from elliptic_descent.cubic_local import CubicCover
from elliptic_descent.geometry import U


@pytest.mark.parametrize('delta', [(4, 0, 0), (16, 0, 0), (0, 0, 4)])
def test_integral_quadrics_reduce_rational_residual_over_qq(delta):
    # These square-equivalent presentations reproduce the Sage failure. The
    # unscaled unit class does not: its quadrics already force QQ inference.
    D = CubicCover((0, 0, 1, -7, 6), delta)
    assert all(c.is_Integer for q in D.quadrics for c in sp.Poly(q, *U).coeffs())
    X, Y, Z = D.map_weighted
    a1, a2, a3, a4, a6 = D.curve.ainvs
    residual = sp.expand(Y**2 + a1*X*Y*Z + a3*Y*Z**3 - X**3
                         - a2*X**2*Z**2 - a4*X*Z**4 - a6*Z**6)
    assert any(c.q != 1 for c in sp.Poly(residual, *U, domain=sp.QQ).coeffs())
    assert D.verify()
