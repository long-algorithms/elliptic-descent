"""Exact cubic-to-quartic presentation conversion, with all projective points.

For y^2=f(x) of degree three and f(s)!=0, set x=s+1/t, y=w/t^2.
The new curve is w^2=t^4*f(s+1/t), a genuine squarefree quartic. No degree-three
polynomial is passed to the quartic-algebra constructor by padding with zero.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .geometry import _rat
from .two_cover import QuarticTwoCover, X, Y, Z
from .weierstrass import WeierstrassModel, WeierstrassChange, x, y
from .marked_cover import MarkedQuarticCover, function_equal


@dataclass(frozen=True)
class CubicToQuartic:
    """A presentation isomorphism, not a Selmer computation or a covering map.

    Points on both double covers use weights (1,2,1). The cubic branch point
    at infinity (1,0,0) is included and maps to (0,0,1) on the quartic.
    """
    f: tuple
    shift: object = None

    def __post_init__(self):
        f = tuple(_rat(c) for c in self.f)
        if len(f) != 4 or not f[3]:
            raise ValueError('expected four ascending coefficients of a genuine cubic')
        p = sp.Poly(sum(c*x**i for i,c in enumerate(f)),x,domain=sp.QQ)
        if not p.discriminant():
            raise ValueError('the cubic must be squarefree')
        s = self.shift
        if s is None:
            s = next(sp.Integer(k) for k in range(4) if p.eval(k) != 0)
        s = _rat(s)
        if p.eval(s) == 0:
            raise ValueError('the reciprocal shift must not be a root')
        object.__setattr__(self,'f',f)
        object.__setattr__(self,'shift',s)

    @cached_property
    def binary_form(self):
        return sum(c*X**i*Z**(4-i) for i,c in enumerate(self.f))

    @cached_property
    def quartic(self):
        s = self.shift
        c0,c1,c2,c3 = self.f
        return QuarticTwoCover((0,c3,c2+3*c3*s,c1+2*c2*s+3*c3*s*s,
                                 c0+c1*s+c2*s*s+c3*s**3))

    def _point(self, point, form):
        if len(point) != 3:
            raise ValueError('expected weighted coordinates (X,Y,Z)')
        a,b,c = tuple(_rat(v) for v in point)
        if not (a or c) or b*b != form.subs({X:a,Z:c},simultaneous=True):
            raise ValueError('invalid point on the specified double-cover model')
        return a,b,c

    def to_quartic(self, point):
        a,b,c = self._point(point,self.binary_form)
        return c,b,a-self.shift*c

    def from_quartic(self, point):
        a,b,c = self._point(point,self.quartic.binary_form)
        return self.shift*a+c,b,a

    def verify(self):
        # Inverse matrix [[0,1],[1,-s]] cancels [[s,1],[1,0]] exactly.
        forward = sp.Matrix([[0,1],[1,-self.shift]])
        backward = sp.Matrix([[self.shift,1],[1,0]])
        converted = self.binary_form.subs({X:self.shift*X+Z,Z:X},simultaneous=True)
        return (forward*backward == sp.eye(2)
                and sp.expand(converted-self.quartic.binary_form) == 0)

    def transport_map(self, coordinates):
        """Precompose supplied cubic affine functions with the inverse conversion.

        Checks rational-function domains only. It does NOT establish that an
        arbitrary input map is a 2-covering; marked_map() requires that evidence.
        """
        if len(coordinates) != 2:
            raise ValueError('expected two affine coordinate functions')
        if not all(function_equal(q,q,self.f) for q in coordinates):
            raise ValueError('undefined or nonrational function on the cubic')
        return tuple(sp.cancel(q.subs({x:self.shift+1/x,y:y/x**2},simultaneous=True))
                     for q in map(sp.sympify,coordinates))

    def marked_map(self, coordinates, change, translation=None):
        """Accept a supplied cubic map only after checking its full marking.

        change starts at the converted quartic's canonical Jacobian. It and
        translation must be supplied explicitly; matching j-invariants is not
        sufficient and no marking/translation is guessed.
        """
        cover = MarkedQuarticCover(self.quartic,change,translation)
        if not cover.verify_supplied_map(self.transport_map(coordinates)):
            raise ValueError('supplied map does not equal the declared marked 2-covering')
        return cover


@dataclass(frozen=True)
class IdentityTwoCover:
    """The zero 2-Selmer class: [2]:E->E, presented by a quartic.

    This is multiplication by TWO, not the identity degree-one morphism.
    An arbitrary rational (possibly nonminimal) input Weierstrass model is kept.
    """
    curve: WeierstrassModel
    shift: object = None

    def __post_init__(self):
        E = self.curve if isinstance(self.curve,WeierstrassModel) else WeierstrassModel(self.curve)
        object.__setattr__(self,'curve',E)
        # Validate the requested shift eagerly.
        object.__setattr__(self,'shift',self.conversion.shift)

    @cached_property
    def conversion(self):
        b2,b4,b6,_ = self.curve.b_invariants
        return CubicToQuartic((b6,2*b4,b2,4),self.shift)

    @cached_property
    def cover(self):
        a1,_,a3,_,_ = self.curve.ainvs
        C = self.conversion.quartic
        change = WeierstrassChange(C.ainvs,(-2,-4*self.shift,-a1,-4*a3))
        if change.target != self.curve:
            raise ArithmeticError('completed-square marking has the wrong target')
        return MarkedQuarticCover(C,change)

    def to_quartic(self, point):
        P = self.curve.point(point)
        if P is None:
            return self.conversion.to_quartic((1,0,0))
        a1,_,a3,_,_ = self.curve.ainvs
        a,b = P
        return self.conversion.to_quartic((a,2*b+a1*a+a3,1))

    def from_quartic(self, point):
        a,b,c = self.conversion.from_quartic(point)
        if not c:
            return None
        a1,_,a3,_,_ = self.curve.ainvs
        return self.curve.point((a/c,(b/c**2-a1*a/c-a3)/2))

    @property
    def identification_maps(self):
        """Degree-one map from the quartic to the original E, for auditing."""
        a1,_,a3,_,_ = self.curve.ainvs
        return self.shift+1/x,(y/x**2-a1*(self.shift+1/x)-a3)/2

    def verify(self):
        if not self.conversion.verify() or not self.cover.verify():
            return False
        a,b = self.identification_maps
        if not function_equal(self.curve.equation(a,b),0,self.cover.source.f):
            return False
        double = tuple(q.subs({x:a,y:b},simultaneous=True) for q in self.curve.duplication_maps)
        return self.cover.verify_supplied_map(double)
