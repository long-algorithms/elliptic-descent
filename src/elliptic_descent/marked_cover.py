"""Explicit markings of quartic 2-coverings on the user's Weierstrass model.

A marking is canonical quartic map, then an origin-preserving isomorphism,
then a rational translation. This construction certifies covering type, not
merely a map identity or a degree. Arbitrary supplied maps must match it.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .two_cover import QuarticTwoCover, X, Y, Z
from .weierstrass import WeierstrassModel, WeierstrassChange, x, y


def reduce_function(expression, coefficients):
    """Normalize to A(x)+B(x)*y in the quadratic function field, exactly.

    Reducing before composing translations avoids large spurious powers of y.
    Denominators are tested in the function field, not just as polynomials.
    """
    expression = sp.sympify(expression)
    if expression.has(sp.Float):
        raise ValueError('inexact function coefficients')
    f = sum(c*x**i for i,c in enumerate(coefficients))
    n,d = sp.fraction(sp.cancel(expression))
    field = sp.QQ.frac_field(x)
    modulus = sp.Poly(y*y-f,y,domain=field)
    n = sp.Poly(n,y,domain=field).rem(modulus)
    d = sp.Poly(d,y,domain=field).rem(modulus)
    d0,d1 = d.nth(0),d.nth(1)
    norm = sp.cancel(d0*d0-f*d1*d1)
    if norm == 0:
        raise ValueError('identically zero function-field denominator')
    A = sp.cancel((n.nth(0)*d0-f*n.nth(1)*d1)/norm)
    B = sp.cancel((n.nth(1)*d0-n.nth(0)*d1)/norm)
    return A+B*y


def function_equal(left, right, coefficients):
    """Equality in Q(x,y)/(y^2-f); reject nonrational/undefined expressions."""
    try:
        a = reduce_function(left,coefficients)
        b = reduce_function(right,coefficients)
        return sp.cancel(a-b) == 0
    except (sp.PolynomialError, ValueError, TypeError, sp.CoercionFailed):
        return False


@dataclass(frozen=True)
class MarkedQuarticCover:
    """A degree-four 2-cover of target, with its marking retained exactly.

    change runs Jac(source) -> target. translation belongs to target(Q).
    Translations change the covering class in general, so are never discarded.
    """
    source: QuarticTwoCover
    change: WeierstrassChange
    translation: tuple | None = None

    def __post_init__(self):
        source = self.source if isinstance(self.source,QuarticTwoCover) else QuarticTwoCover(self.source)
        if not isinstance(self.change,WeierstrassChange):
            raise TypeError('an explicit WeierstrassChange is required')
        if self.change.source != WeierstrassModel(source.ainvs):
            raise ValueError('marking must start at the canonical Jacobian of this quartic')
        object.__setattr__(self,'source',source)
        object.__setattr__(self,'translation',self.target.point(self.translation))

    @property
    def target(self):
        return self.change.target

    @classmethod
    def canonical(cls, f):
        source = f if isinstance(f,QuarticTwoCover) else QuarticTwoCover(f)
        return cls(source,WeierstrassChange(source.ainvs))

    @cached_property
    def map_xy(self):
        a,b,c = (q.subs({X:x,Y:y,Z:1},simultaneous=True) for q in self.source.map_weighted)
        reduce = lambda q: reduce_function(q,self.source.f)
        values = {x:reduce(a/c**2),y:reduce(b/c**3)}
        changed = tuple(reduce(q.subs(values,simultaneous=True)) for q in self.change.maps)
        if self.translation is None:
            return changed
        A,B = changed
        u,v = self.translation
        a1,a2,a3,_,_ = self.target.ainvs
        slope = reduce((B-v)/(A-u))
        new_x = reduce(slope*slope+a1*slope-a2-A-u)
        new_y = reduce(-B-slope*(new_x-A)-a1*new_x-a3)
        return new_x,new_y

    def verify(self):
        return (self.source.check_map() and self.change.verify()
                and function_equal(self.target.equation(*self.map_xy),0,self.source.f))

    def verify_supplied_map(self, coordinates):
        """Compare a claimed marking, not merely its target equation or j-invariant."""
        return (len(coordinates) == 2 and self.verify()
                and all(function_equal(a,b,self.source.f)
                        for a,b in zip(coordinates,self.map_xy)))

    def evaluate(self, point):
        """Evaluate through proper maps; handles rational-map exceptional points."""
        return self.target.add(self.change.apply(self.source.evaluate(point)),self.translation)

    def retarget(self, change):
        """Postcompose with a specified target-model isomorphism."""
        if not isinstance(change,WeierstrassChange) or change.source != self.target:
            raise ValueError('retargeting must start at the current target')
        return MarkedQuarticCover(self.source,self.change.then(change),change.apply(self.translation))

    def compose(self, covering):
        return MarkedFourCover(covering,self)

    def descend(self, **kwargs):
        """Relative descent above this source, with maps to the original target.

        This is not a curve-level Selmer computation. The same source/marking
        is retained even for partial or conditional results.
        """
        from .relative import relative_four_descent
        return MarkedRelativeResult(relative_four_descent(self.source,**kwargs),self)


@dataclass(frozen=True)
class MarkedFourCover:
    model: object
    marking: MarkedQuarticCover

    def __post_init__(self):
        if tuple(self.model.f) != self.marking.source.f:
            raise ValueError('4-covering has a different intermediate quartic')

    @property
    def target(self):
        return self.marking.target

    @cached_property
    def map_xy(self):
        a,b,c = self.model.map_xyz
        return tuple(q.subs({x:a/c,y:b/c**2},simultaneous=True)
                     for q in self.marking.map_xy)

    def verify(self):
        """Stagewise exact proof of the composite, not a new local/Selmer test."""
        return self.model.check_construction() and self.marking.verify()

    def evaluate(self, point):
        equations,image = self.model.evaluate(point)
        if any(equations):
            raise ValueError('point is not on the quadric intersection')
        return self.marking.evaluate(image)


@dataclass(frozen=True)
class MarkedRelativeResult:
    relative_result: object
    marking: MarkedQuarticCover

    def __post_init__(self):
        if self.relative_result.source != self.marking.source:
            raise ValueError('relative result and marking have different sources')

    @property
    def target(self):
        return self.marking.target

    @property
    def coverings(self):
        return tuple(self.marking.compose(D) for D in self.relative_result.coverings)

    @property
    def complete(self):
        return self.relative_result.complete

    @property
    def certified(self):
        return self.relative_result.certified and self.marking.verify()

    @property
    def unresolved(self):
        return self.relative_result.unresolved

    def verify_evidence(self, **kwargs):
        return self.marking.verify() and self.relative_result.verify_evidence(**kwargs)
