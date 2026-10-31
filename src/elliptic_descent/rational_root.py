"""Exact identification of a rational-root quartic with its marked Jacobian.

For a root alpha put A=f'(alpha), B=f''(alpha)/2, C=f'''(alpha)/6.
The reciprocal change produces v^2=u^3+B*u^2+A*C*u+A^2*lc(f).
After x=u+(B-f_2)/3 and y=-v, the original canonical quartic map is [2].
These identities are verified exactly; no point search or rank computation
is used. The descent engine retains its unified norm-pair construction.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .etale import EtaleQuarticAlgebra
from .geometry import _rat
from .two_cover import QuarticTwoCover, X, Y, Z


@dataclass(frozen=True)
class RationalRootIdentification:
    """C -> Jac(C), sending the chosen rational branch point to the origin.

    source is a QuarticTwoCover or five ascending rational coefficients.
    root=None chooses the first labelled rational factor. Both finite and
    infinite points of C are supported in weighted coordinates (X,Y,Z).
    Elliptic points use affine (x,y), with None denoting the origin.
    """
    source: QuarticTwoCover
    root: object = None

    def __post_init__(self):
        source = self.source if isinstance(self.source, QuarticTwoCover) else QuarticTwoCover(tuple(self.source))
        if self.root is None:
            algebra = EtaleQuarticAlgebra(source.f)
            roots = tuple(-g[0] for g in algebra.factors if len(g) == 2)
            if not roots:
                raise ValueError("the quartic has no rational root")
            root = roots[0]
        else:
            root = _rat(self.root)
        if sum(c*root**i for i, c in enumerate(source.f)) != 0:
            raise ValueError("the specified rational number is not a root")
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "root", root)

    @cached_property
    def taylor_coefficients(self):
        """A,B,C,D in f(alpha+h)=A*h+B*h^2+C*h^3+D*h^4."""
        h = sp.Symbol("h")
        p = sp.Poly(sum(c*(self.root+h)**i for i, c in enumerate(self.source.f)), h)
        return tuple(p.nth(i) for i in range(1, 5))

    @property
    def shift(self):
        return (self.taylor_coefficients[1]-self.source.f[2])/3

    @property
    def branch_cubic_ainvs(self):
        """Derivative-based elliptic model, independent of the quartic invariants."""
        a, b, c, d = self.taylor_coefficients
        return (sp.S.Zero, b, sp.S.Zero, a*c, a*a*d)

    @property
    def map_xy(self):
        """Rational affine coordinates on the canonical Jacobian (away from root)."""
        a = self.taylor_coefficients[0]
        w = X-self.root*Z
        return self.shift+a*Z/w, -a*Y/w**2

    def to_jacobian(self, point):
        if len(point) != 3:
            raise ValueError("expected weighted quartic coordinates (X,Y,Z)")
        x, y, z = map(_rat, point)
        if not (x or z):
            raise ValueError("invalid weighted projective point")
        if y*y != self.source.binary_form.subs({X: x, Z: z}):
            raise ValueError("point is not on the source quartic")
        w = x-self.root*z
        if not w:
            return None  # Extension of the birational map at the chosen root.
        a = self.taylor_coefficients[0]
        return (self.shift+a*z/w, -a*y/w**2)

    def from_jacobian(self, point):
        """Inverse map, including the origin and quartic points at infinity."""
        if point is None:
            return (self.root, sp.S.Zero, sp.S.One)
        if len(point) != 2:
            raise ValueError("expected an affine Jacobian point (x,y) or None")
        x, y = map(_rat, point)
        _, a2, _, a4, a6 = self.source.ainvs
        if y*y != x**3+a2*x*x+a4*x+a6:
            raise ValueError("point is not on the specified canonical Jacobian")
        a = self.taylor_coefficients[0]
        u = x-self.shift
        return (self.root*u+a, -a*y, u)

    def check_marked_map(self):
        """Prove pi_C = [2] o identification by exact function-field identities.

        Checks the derivative model, its isomorphism, the inverse map, and the
        duplication identity. Both maps extend uniquely on the smooth proper
        curves; point evaluators handle the excluded branch point explicitly.
        """
        if not self.source.check_map():
            return False
        a, b, c, d = self.taylor_coefficients
        if not a:  # Squarefreeness makes the root simple.
            return False
        _, a2, _, a4, a6 = self.source.ainvs
        u, v = sp.symbols("u v")
        cubic = u**3+b*u*u+a*c*u+a*a*d
        translated = (u+self.shift)**3+a2*(u+self.shift)**2+a4*(u+self.shift)+a6
        if sp.expand(cubic-translated) != 0:
            return False
        F = self.source.binary_form
        ideal = sp.groebner((Y*Y-F.subs(Z, 1),), Y, X, domain=sp.QQ)

        def vanishes(value):
            numerator, denominator = sp.fraction(sp.cancel(value))
            return (ideal.reduce(denominator)[1] != 0
                    and ideal.reduce(numerator)[1] == 0)

        x, y = (sp.cancel(value.subs(Z, 1)) for value in self.map_xy)
        if not vanishes(y*y-x**3-a2*x*x-a4*x-a6):
            return False
        slope = (3*x*x+2*a2*x+a4)/(2*y)
        doubled_x = slope*slope-a2-2*x
        doubled_y = -y+slope*(x-doubled_x)
        px, py, pz = (q.subs(Z, 1) for q in self.source.map_weighted)
        if not (vanishes(doubled_x-px/pz**2) and vanishes(doubled_y-py/pz**3)):
            return False
        # The inverse polynomial map has no basepoint: Z=0 makes X=A!=0.
        inverse = {X: self.root*(u-self.shift)+a, Y: -a*v, Z: u-self.shift}
        residual = sp.expand((Y*Y-F).subs(inverse, simultaneous=True))
        target_ideal = sp.groebner((v*v-u**3-a2*u*u-a4*u-a6,), v, u, domain=sp.QQ)
        if target_ideal.reduce(residual)[1] != 0:
            return False
        inverse_xy = tuple(sp.cancel(q.subs(inverse, simultaneous=True)) for q in self.map_xy)
        return inverse_xy == (u, v)
