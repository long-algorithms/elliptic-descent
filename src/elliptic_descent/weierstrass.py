"""Exact marked Weierstrass models over Q; no Selmer or rank computation.

The (u,r,s,t) convention agrees with Sage's WeierstrassIsomorphism:
source x=u^2*x'+r, y=u^3*y'+s*u^2*x'+t. Maps run source -> target.
See docs/marked-coverings.md for directions, weights and trust boundaries.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .geometry import _rat

x, y = sp.symbols('x y')


@dataclass(frozen=True)
class WeierstrassModel:
    """A nonsingular rational model; None denotes its origin in point APIs."""
    ainvs: tuple

    def __post_init__(self):
        values = tuple(_rat(a) for a in self.ainvs)
        if len(values) == 2:
            values = (sp.S.Zero,)*3 + values
        if len(values) != 5:
            raise ValueError('expected two or five exact Weierstrass coefficients')
        object.__setattr__(self, 'ainvs', values)
        if not self.discriminant:
            raise ValueError('singular Weierstrass model')

    @property
    def b_invariants(self):
        a1,a2,a3,a4,a6 = self.ainvs
        return (a1*a1+4*a2, a1*a3+2*a4, a3*a3+4*a6,
                a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4)

    @property
    def discriminant(self):
        b2,b4,b6,b8 = self.b_invariants
        return -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6

    def equation(self, X=x, Y=y):
        a1,a2,a3,a4,a6 = self.ainvs
        return Y*Y+a1*X*Y+a3*Y-X**3-a2*X*X-a4*X-a6

    def point(self, P):
        """Validate an affine rational point or the origin, without coercion loss."""
        if P is None:
            return None
        if len(P) != 2:
            raise ValueError('expected an affine point (x,y), or None for the origin')
        P = tuple(_rat(a) for a in P)
        if self.equation(*P) != 0:
            raise ValueError('point is not on this Weierstrass model')
        return P

    def negate(self, P):
        P = self.point(P)
        if P is None:
            return None
        X,Y = P
        return X,-Y-self.ainvs[0]*X-self.ainvs[2]

    def add(self, P, Q):
        """Exact chord/tangent addition, including vertical tangents and infinity."""
        P,Q = self.point(P),self.point(Q)
        if P is None:
            return Q
        if Q is None:
            return P
        X,Y = P
        U,V = Q
        a1,a2,a3,a4,_ = self.ainvs
        if X == U:
            if Q == self.negate(P):
                return None
            slope = (3*X*X+2*a2*X+a4-a1*Y)/(2*Y+a1*X+a3)
        else:
            slope = (V-Y)/(U-X)
        new_x = slope*slope+a1*slope-a2-X-U
        new_y = -Y-slope*(new_x-X)-a1*new_x-a3
        return self.point((sp.cancel(new_x),sp.cancel(new_y)))

    def double(self, P):
        return self.add(P,P)

    def translation_maps(self, P):
        """Generic rational maps for +P. Use add() to evaluate exceptional points."""
        P = self.point(P)
        if P is None:
            return x,y
        U,V = P
        a1,a2,a3,_,_ = self.ainvs
        slope = (y-V)/(x-U)
        X = slope*slope+a1*slope-a2-x-U
        Y = -y-slope*(X-x)-a1*X-a3
        return sp.cancel(X),sp.cancel(Y)

    @property
    def duplication_maps(self):
        a1,a2,a3,a4,_ = self.ainvs
        slope = (3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
        X = slope*slope+a1*slope-a2-2*x
        Y = -y-slope*(X-x)-a1*X-a3
        return sp.cancel(X),sp.cancel(Y)


@dataclass(frozen=True)
class WeierstrassChange:
    """Origin-preserving source -> target isomorphism, with exact inverse.

    target is derived, never inferred from j-invariants. To connect to an
    externally specified target compare the complete target.ainvs tuple.
    """
    source: WeierstrassModel
    urst: tuple = (1,0,0,0)

    def __post_init__(self):
        source = self.source if isinstance(self.source,WeierstrassModel) else WeierstrassModel(self.source)
        values = tuple(_rat(a) for a in self.urst)
        if len(values) != 4 or not values[0]:
            raise ValueError('expected (u,r,s,t) with u nonzero')
        object.__setattr__(self,'source',source)
        object.__setattr__(self,'urst',values)

    @cached_property
    def target(self):
        u,r,s,t = self.urst
        a1,a2,a3,a4,a6 = self.source.ainvs
        return WeierstrassModel(((a1+2*s)/u,
            (a2-a1*s+3*r-s*s)/u**2,
            (a3+a1*r+2*t)/u**3,
            (a4+2*a2*r-a1*(r*s+t)-a3*s+3*r*r-2*s*t)/u**4,
            (a6-a1*r*t+a2*r*r-a3*t+a4*r+r**3-t*t)/u**6))

    @property
    def maps(self):
        u,r,s,t = self.urst
        return (x-r)/u**2,(y-s*(x-r)-t)/u**3

    def apply(self, P):
        P = self.source.point(P)
        if P is None:
            return None
        values = dict(zip((x,y),P))
        return self.target.point(tuple(sp.cancel(q.subs(values,simultaneous=True)) for q in self.maps))

    def inverse(self):
        u,r,s,t = self.urst
        return WeierstrassChange(self.target,(1/u,-r/u**2,-s/u,(r*s-t)/u**3))

    def then(self, other):
        """Return other o self, with endpoints checked (not guessed)."""
        if not isinstance(other,WeierstrassChange) or self.target != other.source:
            raise ValueError('model-change endpoints do not match')
        u,r,s,t = self.urst
        U,R,S,T = other.urst
        return WeierstrassChange(self.source,(u*U,r+u*u*R,s+u*S,t+u*u*s*R+u**3*T))

    def apply_weighted(self, coordinates):
        """Transport (X,Y,Z) of weights (2,3,1), with no point coercion."""
        if len(coordinates) != 3:
            raise ValueError('expected three weighted coordinates')
        X,Y,Z = coordinates
        u,r,s,t = self.urst
        return (X-r*Z**2)/u**2,(Y-s*(X-r*Z**2)*Z-t*Z**3)/u**3,Z

    def verify(self):
        u,_,_,_ = self.urst
        transformed = self.target.equation(*self.maps)
        return (sp.cancel(u**6*transformed-self.source.equation()) == 0
                and self.source.discriminant == u**12*self.target.discriminant
                and self.inverse().target == self.source)
