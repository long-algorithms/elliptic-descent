"""Exact square-class to degree-4 genus-one model construction.

Reference: Womack, *Explicit Descent on Elliptic Curves* (2003),
and Fisher, *Some improvements to 4-descent on an elliptic curve* (2008).

This module constructs *candidate* intersections of two quadrics.  It does not
compute Selmer groups, impose local conditions, minimize, or prove a covering
is everywhere locally soluble.  All coordinates use ascending powers of theta.
"""

from dataclasses import dataclass
from itertools import combinations
from math import isqrt
from typing import Sequence

import sympy as sp

U = sp.symbols("u0:4")


def _rat(value):
    """Reject non-rational or inexact scalar inputs."""
    if isinstance(value, float):
        raise TypeError("floating-point input is not allowed")
    v = sp.sympify(value)
    if v.is_Rational is not True:
        raise TypeError(f"expected an exact rational, got {value!r}")
    return sp.Rational(v)


def _four(values, name):
    if len(values) != 4:
        raise ValueError(f"{name} must have four ascending coefficients")
    return tuple(_rat(v) for v in values)


def _five(values):
    if len(values) != 5:
        raise ValueError("quartic must have five ascending coefficients")
    values = tuple(_rat(v) for v in values)
    if not values[4]:
        raise ValueError("quartic leading coefficient must be nonzero")
    t = sp.Symbol("t")
    f = sum(values[i] * t**i for i in range(5))
    if sp.degree(sp.gcd(f, sp.diff(f, t)), t) > 0:
        raise ValueError("quartic must be squarefree")
    return values


def _mul(a, b, f):
    """Multiply in Q[t]/(f), allowing polynomial coefficients in U."""
    out = [sp.S.Zero] * 7
    for i in range(4):
        for j in range(4):
            out[i + j] += a[i] * b[j]
    leading = f[4]
    for k in range(6, 3, -1):
        factor = out[k] / leading
        for i in range(4):
            out[k - 4 + i] -= factor * f[i]
        out[k] = sp.S.Zero
    return tuple(sp.expand(x) for x in out[:4])


def _norm(v, f):
    basis = [(sp.Integer(1) if j == i else sp.Integer(0)) for i in range(4) for j in range(4)]
    columns = [_mul(v, basis[4 * i:4 * i + 4], f) for i in range(4)]
    return sp.expand(sp.Matrix.hstack(*(sp.Matrix(c) for c in columns)).det(method="berkowitz"))


def _sqrt_rational(v):
    v = _rat(v)
    if v < 0:
        raise ValueError("norm divided by the leading coefficient is not a square")
    p, q = map(int, sp.fraction(v))
    a, b = isqrt(p), isqrt(q)
    if a * a != p or b * b != q:
        raise ValueError("norm divided by the leading coefficient is not a square")
    return sp.Rational(a, b)


@dataclass(frozen=True)
class FourCoverCandidate:
    """A pair of quadrics and a weighted-projective map to y^2=F(x,z)."""

    f: tuple
    xi: tuple
    norm_witness: sp.Rational
    q: tuple
    norm_u: sp.Expr

    @property
    def quadrics(self):
        return self.q[2], self.q[3]

    @property
    def map_xyz(self):
        """(X,Y,Z), weighted degrees (2,4,2), respectively."""
        return -self.q[0], self.f[4] * self.norm_witness * self.norm_u, self.q[1]

    def evaluate(self, point):
        """Evaluate equations and (X,Y,Z) at rational projective coordinates.

        Does not imply local solubility or verify smoothness of the model.
        """
        if len(point) != 4:
            raise ValueError("point must have four coordinates")
        coords = tuple(map(_rat, point))
        if all(c == 0 for c in coords):
            raise ValueError("projective point cannot be zero")
        s = dict(zip(U, coords))
        equations = tuple(sp.expand(Q.subs(s)) for Q in self.quadrics)
        image = tuple(sp.expand(v.subs(s)) for v in self.map_xyz)
        return equations, image

    def check_algebra(self):
        """Independently check the universal norm identity in four formal variables.

        This checks the norm witnesses and the quartic identity in the quotient
        algebra.  It is NOT a certification of local solubility or completeness.
        """
        if _norm(self.xi, self.f) != self.f[4] * self.norm_witness**2:
            return False
        symbols = sp.symbols("v0:4")
        v = tuple(symbols)
        Q = _mul(self.xi, _mul(v, v, self.f), self.f)
        if not all(sp.expand(a - b) == 0 for a, b in zip(Q, (q.subs(dict(zip(U, v))) for q in self.q))):
            return False
        a, b = sp.symbols("a b")
        linear = (a, b, sp.S.Zero, sp.S.Zero)
        direct = sp.expand(_norm(linear, self.f))
        formula = sp.expand(sum(self.f[i] * (-a)**i * b**(4 - i) for i in range(5)) / self.f[4])
        return sp.expand(direct - formula) == 0

    def check_smoothness(self):
        """Certify smoothness over Qbar via the projective Jacobian criterion.

        The singular locus of an intersection of two quadrics is cut out
        by the quadrics and all 2x2 minors of their 2x4 Jacobian.
        Exhaustively test the four standard projective affine charts using
        exact Groebner bases over Q.  This is deliberately an audit routine,
        not the fastest method for a large-coefficient pencil.
        """
        jac = sp.Matrix(2, 4, lambda i, j: sp.diff(self.quadrics[i], U[j]))
        minors = tuple(sp.expand(jac[0, i]*jac[1, j] - jac[0, j]*jac[1, i])
                       for i, j in combinations(range(4), 2))
        for k in range(4):
            coordinates = U[:k] + U[k+1:]
            specialization = {U[k]: sp.Integer(1)}
            equations = tuple(sp.expand(q.subs(specialization))
                              for q in (*self.quadrics, *minors))
            ideal = sp.groebner(equations, *coordinates, order="grevlex")
            if sp.expand(ideal.reduce(sp.Integer(1))[1]) != 0:
                return False
        return True

    def check_map_basepoint_free(self):
        """Certify the coordinate map has no projective common zero on D.

        It suffices to test q0=q1=q2=q3=0, because the map has
        X=-q0 and Z=q1.  This checks all four projective charts over
        the algebraic closure, not merely over Q-rational points.
        """
        for k in range(4):
            coordinates = U[:k] + U[k+1:]
            specialization = {U[k]: sp.Integer(1)}
            equations = tuple(sp.expand(q.subs(specialization)) for q in self.q)
            ideal = sp.groebner(equations, *coordinates, order="grevlex")
            if sp.expand(ideal.reduce(sp.Integer(1))[1]) != 0:
                return False
        return True

    def check_map_identity(self):
        """Verify Y^2-F(X,Z) vanishes modulo the ideal (q2,q3).

        A direct exact polynomial check, independent of the norm-law argument.
        This does not certify the map is everywhere defined or its degree.
        """
        x, y, z = self.map_xyz
        residual = sp.expand(y**2 - sum(self.f[i] * x**i * z**(4-i) for i in range(5)))
        ideal = sp.groebner(self.quadrics, *U, order="grevlex")
        return sp.expand(ideal.reduce(residual)[1]) == 0


def construct_four_cover(f: Sequence, xi: Sequence, *, norm_witness=None):
    """Build quadrics from xi with Norm(xi)=lc(f)*norm_witness**2.

    The candidate is {q2=q3=0} in P^3. On it, the map
    (X,Y,Z)=(-q0,lc(f)*r*Norm(U),q1) satisfies Y^2=F(X,Z).
    Neither global nor local Selmer membership is tested here.
    """
    f, xi = _five(f), _four(xi, "xi")
    nx = _norm(xi, f)
    if not nx:
        raise ValueError("xi must be invertible in the quartic algebra")
    r = _sqrt_rational(nx / f[4]) if norm_witness is None else _rat(norm_witness)
    if sp.expand(nx - f[4] * r**2) != 0:
        raise ValueError("incorrect norm witness")
    u = tuple(U)
    q = _mul(xi, _mul(u, u, f), f)
    return FourCoverCandidate(f=f, xi=xi, norm_witness=r, q=q, norm_u=_norm(u, f))
