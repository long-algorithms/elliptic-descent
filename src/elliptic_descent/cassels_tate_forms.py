"""Fisher's binary-quartic formula, Theorem 3.1 of arXiv:2208.14977.

Only exact polynomial identities live here. These objects do not establish
local solubility, identify Selmer labels, or assert a complete pairing matrix.
All coefficient tuples are ascending; the cubic generator is phi, not theta.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .cubic_arithmetic import CubicEtaleAlgebra, cubic_multiply, T
from .geometry import _five, _rat
from .two_cover import X, Z


def invariants(f):
    e, d, c, b, a = _five(f)
    return 12*a*e-3*b*d+c*c, 72*a*c*e-27*a*d*d-27*b*b*e+9*b*c*d-2*c**3


def binary_value(coefficients, x, z):
    n = len(coefficients)-1
    return sum(c*x**i*z**(n-i) for i, c in enumerate(coefficients))


def change_parameters(f, matrix):
    M = sp.Matrix(matrix)
    if M.shape != (2, 2) or not M.det():
        raise ValueError('expected an invertible rational parameter matrix')
    M = M.applyfunc(_rat)
    u, v = M*sp.Matrix((X, Z))
    form = sp.Poly(sp.expand(binary_value(_five(f), u, v)), X, Z, domain=sp.QQ)
    return tuple(form.coeff_monomial(X**i*Z**(4-i)) for i in range(5))


def substitute_generator(value, scale, shift):
    p = sp.Poly(sum(c*(scale*T+shift)**i for i, c in enumerate(value)), T, domain=sp.QQ)
    return tuple(p.nth(i) for i in range(3))


@dataclass(frozen=True)
class FisherQuartic:
    f: tuple

    def __post_init__(self):
        object.__setattr__(self, 'f', _five(self.f))

    @property
    def invariants(self):
        return invariants(self.f)

    @cached_property
    def algebra(self):
        I, J = self.invariants
        return CubicEtaleAlgebra((J, -3*I, 0, 1))

    @property
    def z(self):
        e, d, c, b, a = self.f
        return ((3*b*b-8*a*c)/3, 4*a/3, sp.S.Zero)

    @cached_property
    def hessian(self):
        e, d, c, b, a = self.f
        return (3*d*d-8*c*e, 4*(c*d-6*b*e), 2*(2*c*c-24*a*e-3*b*d),
                4*(b*c-6*a*d), 3*b*b-8*a*c)

    @cached_property
    def H(self):
        """Three binary quadratic forms, coefficients of 1,phi,phi^2."""
        I, _ = self.invariants
        g = binary_value(self.f, X, Z)
        h = binary_value(self.hessian, X, Z)
        return (sp.expand(sp.diff(h, X, 2)/36+2*I*Z*Z/9),
                sp.expand(sp.diff(g, X, 2)/9), -sp.Rational(2, 9)*Z*Z)

    def verify_identity(self):
        g = binary_value(self.f, X, Z)
        h = binary_value(self.hessian, X, Z)
        G = (h/3, 4*g/3, sp.S.Zero)
        left = cubic_multiply(self.z, G, self.algebra.f)
        right = cubic_multiply(self.H, self.H, self.algebra.f)
        return all(sp.expand(a-b) == 0 for a, b in zip(left, right))


def inverse_element(algebra, value):
    p = sp.Poly(sum(c*T**i for i, c in enumerate(value)), T, domain=sp.QQ)
    modulus = sp.Poly(sum(c*T**i for i, c in enumerate(algebra.f)), T, domain=sp.QQ)
    if not algebra.norm(value):
        raise ValueError('zero divisor has no multiplicative inverse')
    inv = sp.invert(p, modulus)
    return tuple(inv.nth(i) for i in range(3))


@dataclass(frozen=True)
class FisherTriple:
    """Three same-invariant forms and an exact product square root.

    No reduction of z(g) by a square is permitted: that would change gamma.
    Different legitimate square roots in an etale product are allowed.
    """
    quartics: tuple
    square_root: tuple

    def __post_init__(self):
        Q = tuple(q if isinstance(q, FisherQuartic) else FisherQuartic(q) for q in self.quartics)
        if len(Q) != 3 or len({q.invariants for q in Q}) != 1:
            raise ValueError('expected three quartics with the same I and J')
        m = tuple(_rat(c) for c in self.square_root)
        if len(m) != 3:
            raise ValueError('expected three cubic-algebra coefficients')
        object.__setattr__(self, 'quartics', Q)
        object.__setattr__(self, 'square_root', m)
        A = Q[0].algebra
        if any(not A.norm(q.z) for q in Q) or not A.norm(m):
            raise ValueError('Fisher cubic invariants and square root must be units')
        product = A.multiply(A.multiply(Q[0].z, Q[1].z), Q[2].z)
        if A.multiply(m, m) != product:
            raise ValueError('the supplied product square root is incorrect')
        if not any(self.gamma):
            raise ArithmeticError('the Fisher quadratic is identically zero')

    @cached_property
    def gamma(self):
        q1, q2, q3 = self.quartics
        A = q1.algebra
        multiplier = A.multiply(A.multiply(q2.z, q3.z), inverse_element(A, self.square_root))
        form = cubic_multiply(multiplier, q1.H, A.f)[2]
        p = sp.Poly(form, X, Z, domain=sp.QQ)
        return tuple(p.coeff_monomial(X**i*Z**(2-i)) for i in range(3))

    @property
    def hilbert_parameter(self):
        return self.quartics[1].f[4]

    def verify(self):
        Q = self.quartics
        A = Q[0].algebra
        return (len({q.invariants for q in Q}) == 1
                and all(q.verify_identity() and A.norm(q.z) for q in Q)
                and A.multiply(self.square_root, self.square_root) ==
                    A.multiply(A.multiply(Q[0].z, Q[1].z), Q[2].z)
                and any(self.gamma))
