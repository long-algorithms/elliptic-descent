"""Cubic etale algebra and the classical x-theta Kummer map for 2-descent.

Use the monic completed-square equation v^2=h(x), v=y+(a1*x+a3)/2.
Unlike the quartic relative-descent presentation, H^1(Q,E[2]) is a norm
kernel in L*/L*2: there is NO rational-scalar quotient. See docs/two-selmer.md.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .arithmetic import prime_support
from .geometry import _rat, _sqrt_rational
from .weierstrass import WeierstrassModel

T = sp.Symbol('theta')


def _vector(values, length):
    if len(values) != length:
        raise ValueError(f'expected {length} ascending exact coefficients')
    return tuple(_rat(c) for c in values)


def _poly(values):
    return sp.Poly.from_list(list(reversed(values)), T, domain=sp.QQ)


def cubic_multiply(left, right, f):
    """Internal exact reduction; also allows formal polynomial coefficients."""
    out = [sp.S.Zero]*5
    for i in range(3):
        for j in range(3):
            out[i+j] += left[i]*right[j]
    for k in (4, 3):
        for i in range(3):
            out[k-3+i] -= out[k]*f[i]
    return tuple(sp.expand(c) for c in out[:3])


def cubic_norm(element, f):
    columns = [cubic_multiply(element, tuple(sp.Integer(i == j) for j in range(3)), f)
               for i in range(3)]
    return sp.expand(sp.Matrix.hstack(*(sp.Matrix(c) for c in columns)).det(method='berkowitz'))


@dataclass(frozen=True)
class CubicEtaleAlgebra:
    """Labelled factor algebra of a monic squarefree rational cubic."""
    f: tuple

    def __post_init__(self):
        f = _vector(self.f, 4)
        if f[3] != 1 or not _poly(f).discriminant():
            raise ValueError('expected a monic squarefree cubic')
        object.__setattr__(self, 'f', f)

    @cached_property
    def factors(self):
        factors = tuple(sorted((tuple(g.monic().nth(i) for i in range(g.degree()+1))
                                for g, _ in _poly(self.f).factor_list()[1]),
                               key=lambda g: (len(g), g)))
        product = sp.Poly(1, T, domain=sp.QQ)
        for g in factors:
            product *= _poly(g)
        if product != _poly(self.f):
            raise ArithmeticError('cubic factorization reconstruction failed')
        return factors

    @property
    def degrees(self):
        return tuple(len(g)-1 for g in self.factors)

    @property
    def one(self):
        return (sp.S.One, sp.S.Zero, sp.S.Zero)

    @property
    def support(self):
        return prime_support((2, _poly(self.f).discriminant(), *(c.q for c in self.f)))

    @cached_property
    def idempotents(self):
        modulus = _poly(self.f)
        result = []
        for g in map(_poly, self.factors):
            h = modulus.exquo(g)
            e = (h*sp.invert(h, g)).rem(modulus)
            result.append(tuple(e.nth(i) for i in range(3)))
        return tuple(result)

    def to_components(self, element):
        p = _poly(_vector(element, 3))
        return tuple(tuple(p.rem(_poly(g)).nth(i) for i in range(len(g)-1)) for g in self.factors)

    def from_components(self, components):
        if len(components) != len(self.factors):
            raise ValueError('incorrect component count')
        p = sp.Poly(0, T, domain=sp.QQ)
        for values, degree, e in zip(components, self.degrees, self.idempotents):
            p += _poly(_vector(values, degree))*_poly(e)
        p = p.rem(_poly(self.f))
        return tuple(p.nth(i) for i in range(3))

    def embed_unit(self, index, values):
        if type(index) is not int or not 0 <= index < len(self.factors):
            raise ValueError('invalid factor index')
        values = _vector(values, self.degrees[index])
        if not any(values):
            raise ValueError('zero is not a factor unit')
        components = [(sp.S.One,)+(sp.S.Zero,)*(d-1) for d in self.degrees]
        components[index] = values
        return self.from_components(components)

    def multiply(self, left, right):
        return cubic_multiply(_vector(left, 3), _vector(right, 3), self.f)

    def norm(self, element):
        """Component resultants, independent of the multiplication-matrix norm."""
        return sp.Rational(sp.prod(_poly(g).resultant(_poly(c))
                                  for g, c in zip(self.factors, self.to_components(element))))

    def verify(self):
        for i in range(3):
            e = tuple(sp.Integer(i == j) for j in range(3))
            if self.from_components(self.to_components(e)) != e:
                return False
        for i, e in enumerate(self.idempotents):
            for j, v in enumerate(self.idempotents):
                if self.multiply(e, v) != (e if i == j else (0, 0, 0)):
                    return False
        return tuple(map(sum, zip(*self.idempotents))) == self.one


def kummer_algebra(curve):
    curve = curve if isinstance(curve, WeierstrassModel) else WeierstrassModel(curve)
    b2, b4, b6, _ = curve.b_invariants
    return CubicEtaleAlgebra((b6/4, b4/2, b2/4, sp.S.One))


def kummer_value(curve, point):
    """L*/L*2 representative; O -> 1 and a rational root uses h'(root).

    This is an arithmetic class, not a claim that a point has been divided by
    two. The derivative replaces the zero component at a nonzero 2-torsion point.
    """
    curve = curve if isinstance(curve, WeierstrassModel) else WeierstrassModel(curve)
    P = curve.point(point)
    algebra = kummer_algebra(curve)
    if P is None:
        return algebra.one
    x, _ = P
    parts = list(algebra.to_components((x, -1, 0)))
    for i, g in enumerate(algebra.factors):
        if len(g) == 2 and x == -g[0]:
            parts[i] = (_poly(algebra.f).diff().eval(x),)
    result = algebra.from_components(parts)
    if not algebra.norm(result):
        raise ArithmeticError('Kummer value is not a unit')
    _sqrt_rational(algebra.norm(result))
    return result


def rational_two_torsion(curve):
    curve = curve if isinstance(curve, WeierstrassModel) else WeierstrassModel(curve)
    a1, _, a3, _, _ = curve.ainvs
    return (None,)+tuple(curve.point((-g[0], (a1*g[0]-a3)/2))
                        for g in kummer_algebra(curve).factors if len(g) == 2)


def rational_halves_of_two_torsion(curve, point):
    """Exhaust all rational halves by the exact duplication polynomial.

    No height search. A rational half has rational x among this quartic's roots;
    rational square roots supply its two y choices, checked by the group law.
    """
    curve = curve if isinstance(curve, WeierstrassModel) else WeierstrassModel(curve)
    P = curve.point(point)
    if P is None:
        return rational_two_torsion(curve)
    if curve.double(P) is not None:
        raise ValueError('point is not 2-torsion')
    h = _poly(kummer_algebra(curve).f)
    A = h.nth(2)
    q = h.diff()**2-4*h*sp.Poly(2*T+A+P[0], T, domain=sp.QQ)
    a1, _, a3, _, _ = curve.ainvs
    result = []
    for x in sorted(q.ground_roots()):
        if not h.eval(x):
            continue
        try:
            v = _sqrt_rational(h.eval(x))
        except ValueError:
            continue
        for sign in (-1, 1):
            Q = curve.point((x, sign*v-(a1*x+a3)/2))
            if curve.double(Q) != P:
                raise ArithmeticError('duplication root failed point verification')
            if Q not in result:
                result.append(Q)
    return tuple(result)
