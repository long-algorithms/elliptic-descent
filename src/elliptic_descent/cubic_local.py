"""Auxiliary local-test models for labelled cubic-algebra 2-Selmer classes.

For N(delta)=r^2, write delta*(u0+u1*theta+u2*theta^2)^2=Q0+Q1*theta+Q2*theta^2.
The smooth intersection Q2=0, Q1+u3^2=0 is the 2-cover of the supplied E.
These models decide Selmer membership; binary-quartic presentation is separate.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .cubic_arithmetic import T, _vector, kummer_algebra, cubic_multiply, cubic_norm
from .geometry import U, _sqrt_rational
from .local import bad_primes
from .local_complete import padic_solubility_complete
from .real import _interval_value
from .status import LocalStatus
from .weierstrass import WeierstrassModel


@dataclass(frozen=True)
class CubicCover:
    curve: WeierstrassModel
    delta: tuple

    def __post_init__(self):
        E = self.curve if isinstance(self.curve, WeierstrassModel) else WeierstrassModel(self.curve)
        object.__setattr__(self, 'curve', E)
        object.__setattr__(self, 'delta', _vector(self.delta, 3))
        if not self.norm_witness:
            raise ValueError('a descent class must be a unit')

    @cached_property
    def algebra(self):
        return kummer_algebra(self.curve)

    @cached_property
    def norm_witness(self):
        return _sqrt_rational(self.algebra.norm(self.delta))

    @cached_property
    def q(self):
        return cubic_multiply(self.delta, cubic_multiply(U[:3], U[:3], self.algebra.f), self.algebra.f)

    @property
    def quadrics(self):
        return self.q[2], self.q[1]+U[3]**2

    @cached_property
    def norm_u(self):
        return cubic_norm(U[:3], self.algebra.f)

    @property
    def map_weighted(self):
        a1, _, a3, _, _ = self.curve.ainvs
        z = U[3]
        return self.q[0], self.norm_witness*self.norm_u-(a1*self.q[0]*z+a3*z**3)/2, z

    def verify(self):
        """Exact construction and polynomial map identity, no local assertion."""
        f = sum(c*T**i for i, c in enumerate(self.algebra.f))
        u = sum(c*T**i for i, c in enumerate(U[:3]))
        if sp.expand(sp.resultant(f, u, T)-self.norm_u) != 0:
            return False
        if cubic_norm(self.delta, self.algebra.f) != self.norm_witness**2:
            return False
        X, Y, Z = self.map_weighted
        a1, a2, a3, a4, a6 = self.curve.ainvs
        residual = sp.expand(Y**2+a1*X*Y*Z+a3*Y*Z**3-X**3-a2*X**2*Z**2-a4*X*Z**4-a6*Z**6)
        ideal = sp.groebner(self.quadrics, *U, order='grevlex', domain=sp.QQ)
        return ideal.reduce(residual)[1] == 0

    def evaluate(self, point):
        values = _vector(point, 4)
        if not any(values):
            raise ValueError('zero is not a projective point')
        sub = dict(zip(U, values))
        if any(q.subs(sub) != 0 for q in self.quadrics):
            raise ValueError('point is not on the auxiliary covering')
        X, Y, Z = (q.subs(sub) for q in self.map_weighted)
        return None if not Z else self.curve.point((X/Z**2, Y/Z**3))


def _allowed(signs):
    # No global +/- ambiguity: coefficient of theta is -z^2, not an
    # arbitrary scalar. x-theta is positive before x and negative after x.
    return any(signs == (1,)*k+(-1,)*(len(signs)-k) for k in range(len(signs)+1))


@dataclass(frozen=True)
class CubicRealEvidence:
    f: tuple
    delta: tuple
    intervals: tuple
    signs: tuple
    status: LocalStatus

    def verify(self, model):
        if self.f != model.algebra.f or self.delta != model.delta:
            return False
        p = sp.Poly(sum(c*T**i for i, c in enumerate(self.f)), T, domain=sp.QQ)
        if len(self.signs) != len(self.intervals) or len(self.signs) != p.count_roots(-sp.oo, sp.oo):
            return False
        previous = -sp.oo
        for interval, sign in zip(self.intervals, self.signs):
            try:
                lo, hi = _vector(interval, 2)
            except (TypeError, ValueError):
                return False
            if not previous < lo <= hi or p.count_roots(lo, hi) != 1:
                return False
            previous = hi
            a, b = _interval_value(self.delta, lo, hi)
            if not ((sign == 1 and a > 0) or (sign == -1 and b < 0)):
                return False
        expected = LocalStatus.SOLUBLE if _allowed(self.signs) else LocalStatus.INSOLUBLE
        return self.status == expected


def cubic_real_solubility(model):
    p = sp.Poly(sum(c*T**i for i, c in enumerate(model.algebra.f)), T, domain=sp.QQ)
    intervals, signs, previous = [], [], -sp.oo
    for (lo, hi), _ in p.intervals():
        while True:
            a, b = _interval_value(model.delta, lo, hi)
            closed = lo == hi or (p.eval(lo) != 0 and p.eval(hi) != 0)
            if previous < lo and closed and (a > 0 or b < 0):
                break
            if lo == hi:
                raise ArithmeticError('nonunit class or inconsistent isolation')
            lo, hi = p.refine_root(lo, hi, steps=1)
        intervals.append((lo, hi))
        signs.append(1 if a > 0 else -1)
        previous = hi
    signs = tuple(signs)
    return CubicRealEvidence(model.algebra.f, model.delta, tuple(intervals), signs,
                             LocalStatus.SOLUBLE if _allowed(signs) else LocalStatus.INSOLUBLE)


@dataclass(frozen=True)
class CubicLocalEvidence:
    real: CubicRealEvidence
    finite: tuple
    good_reduction: object

    @property
    def status(self):
        if self.real.status == LocalStatus.INSOLUBLE or any(r.status == LocalStatus.INSOLUBLE for r in self.finite):
            return LocalStatus.INSOLUBLE
        if any(r.certificate and r.certificate.kind == 'rational_point' for r in self.finite):
            return LocalStatus.SOLUBLE
        if self.good_reduction is not None and self.real.status == LocalStatus.SOLUBLE:
            soluble = {r.place for r in self.finite if r.status == LocalStatus.SOLUBLE}
            if soluble >= set(self.good_reduction.primes):
                return LocalStatus.SOLUBLE
        return LocalStatus.UNKNOWN

    def verify(self, model, *, max_nodes=1_000_000):
        if self.status == LocalStatus.UNKNOWN or not self.real.verify(model):
            return False
        if self.real.status == LocalStatus.INSOLUBLE:
            return not self.finite and self.good_reduction is None
        return (self.good_reduction == bad_primes(model)
                and all(r.verify(model, max_nodes=max_nodes) for r in self.finite
                        if r.status != LocalStatus.UNKNOWN))


def cubic_local_solubility(model, *, max_depth=None, max_nodes=None):
    real = cubic_real_solubility(model)
    if real.status == LocalStatus.INSOLUBLE:
        return CubicLocalEvidence(real, (), None)
    bound = bad_primes(model)
    finite = []
    for p in bound.primes:
        r = padic_solubility_complete(model, p, max_depth=max_depth, max_nodes=max_nodes,
                                      _smooth_checked=True)
        finite.append(r)
        if r.status == LocalStatus.INSOLUBLE or (r.certificate and r.certificate.kind == 'rational_point'):
            break
    return CubicLocalEvidence(real, tuple(finite), bound)
