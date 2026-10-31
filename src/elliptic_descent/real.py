"""Complete real-place test for square-class quadric constructions.

In R[t]/(f), complex factors impose no square-root condition. At ordered
real roots alpha_i the signs of xi(alpha_i) must be a threshold pattern,
up to a common sign. Root isolation and interval signs are exact.
"""
from dataclasses import dataclass

import sympy as sp

from .geometry import _five, _four, _norm, _rat
from .status import LocalStatus


def _interval_value(coefficients, lo, hi):
    a = b = sp.S.Zero
    for coefficient in reversed(coefficients):
        products = (a*lo, a*hi, b*lo, b*hi)
        a, b = min(products)+coefficient, max(products)+coefficient
    return a, b


def _patterns(signs):
    """Possible xi signs: common sign times sign(x-alpha_i)."""
    return tuple((k, s) for k in range(len(signs)+1) for s in (1, -1)
                 if all(sign == s*(1 if i < k else -1)
                        for i, sign in enumerate(signs)))


@dataclass(frozen=True)
class RealSquareClassResult:
    f: tuple
    xi: tuple
    intervals: tuple
    signs: tuple
    status: LocalStatus
    place: str = "real"

    def verify(self, model):
        if tuple(model.f) != self.f or tuple(model.xi) != self.xi:
            return False
        t = sp.Symbol("t")
        f = sp.Poly.from_list(list(reversed(self.f)), t, domain=sp.QQ)
        if len(self.intervals) != len(self.signs) or len(self.intervals) != f.count_roots(-sp.oo, sp.oo):
            return False
        previous = -sp.oo
        for interval, sign in zip(self.intervals, self.signs):
            try:
                lo, hi = (_rat(value) for value in interval)
            except (TypeError, ValueError):
                return False
            if not previous < lo <= hi:
                return False
            previous = hi
            if f.count_roots(lo, hi) != 1:
                return False
            low, high = _interval_value(self.xi, lo, hi)
            if not ((sign == 1 and low > 0) or (sign == -1 and high < 0)):
                return False
        expected = LocalStatus.SOLUBLE if _patterns(self.signs) else LocalStatus.INSOLUBLE
        return self.status == expected


def real_squareclass_solubility(model):
    """Decide R-solubility for D_xi, not for an arbitrary quadric intersection."""
    f, xi = _five(model.f), _four(model.xi, "xi")
    if not _norm(xi, f):
        raise ValueError("xi must be invertible")
    t = sp.Symbol("t")
    polynomial = sp.Poly.from_list(list(reversed(f)), t, domain=sp.QQ)
    intervals, signs = [], []
    previous = -sp.oo
    for (lo, hi), multiplicity in polynomial.intervals():
        if multiplicity != 1:
            raise ValueError("quartic must be squarefree")
        while True:
            low, high = _interval_value(xi, lo, hi)
            # SymPy isolates roots in open intervals: an endpoint can itself
            # be a different (rational) root. Close and separate each interval
            # before storing it, or a following singleton cannot be separated.
            closed_isolation = (lo == hi or
                                (polynomial.eval(lo) != 0 and polynomial.eval(hi) != 0))
            if previous < lo and closed_isolation and (low > 0 or high < 0):
                break
            if lo == hi:
                raise ArithmeticError("inconsistent root isolation or nonunit xi")
            lo, hi = polynomial.refine_root(lo, hi, steps=1)
        intervals.append((lo, hi))
        signs.append(1 if low > 0 else -1)
        previous = hi
    status = LocalStatus.SOLUBLE if _patterns(tuple(signs)) else LocalStatus.INSOLUBLE
    return RealSquareClassResult(f, xi, tuple(intervals), tuple(signs), status)
