"""Canonical quartic 2-cover and its Jacobian map.

Fisher--Sills, Lemma 3.1, specialized to y^2=F(x,z). Coordinates on the
quartic are (X,Y,Z) of weights (1,2,1); those on E have weights (2,3,1).
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .geometry import _five, _rat

X, Y, Z = sp.symbols("X Y Z")


@dataclass(frozen=True)
class QuarticTwoCover:
    f: tuple

    def __post_init__(self):
        object.__setattr__(self, "f", _five(self.f))

    @cached_property
    def binary_form(self):
        return sum(c*X**i*Z**(4-i) for i, c in enumerate(self.f))

    @property
    def ainvs(self):
        e, d, c, b, a = self.f
        return (sp.S.Zero, c, sp.S.Zero, b*d-4*a*e,
                a*d*d+b*b*e-4*a*c*e)

    @cached_property
    def map_weighted(self):
        g = self.binary_form
        h = sp.expand((sp.diff(g, X, Z)**2-sp.diff(g, X, 2)*sp.diff(g, Z, 2))/3)
        k = sp.expand((sp.diff(g, Z)*sp.diff(h, X)-sp.diff(g, X)*sp.diff(h, Z))/12)
        return sp.expand((h-4*self.f[2]*g)/3), k, 2*Y

    def check_map(self):
        x, y, z = self.map_weighted
        _, a2, _, a4, a6 = self.ainvs
        residual = sp.expand(y*y-x**3-a2*x*x*z*z-a4*x*z**4-a6*z**6)
        remainder = sp.Poly(residual, Y).rem(sp.Poly(Y*Y-self.binary_form, Y))
        return remainder.is_zero

    def evaluate(self, point):
        """Return affine E-coordinates, or None for the point at infinity."""
        if len(point) != 3:
            raise ValueError("expected weighted quartic coordinates (X,Y,Z)")
        point = tuple(map(_rat, point))
        if point[0] == point[2] == 0:
            raise ValueError("invalid weighted projective point")
        substitutions = dict(zip((X, Y, Z), point))
        if point[1]**2 != self.binary_form.subs(substitutions):
            raise ValueError("point is not on the quartic")
        x, y, z = (coordinate.subs(substitutions) for coordinate in self.map_weighted)
        return None if not z else (sp.cancel(x/z**2), sp.cancel(y/z**3))

    def compose(self, covering):
        """Weighted equations for the degree-16 map D -> C -> E."""
        if tuple(covering.f) != self.f:
            raise ValueError("covering has a different target quartic")
        substitutions = dict(zip((X, Y, Z), covering.map_xyz))
        return tuple(v.subs(substitutions, simultaneous=True) for v in self.map_weighted)
