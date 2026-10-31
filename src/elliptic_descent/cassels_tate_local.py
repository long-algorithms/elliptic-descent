"""Exact local evaluations and finite support for Fisher's pairing formula.

A rational parameter where g is a nonzero Q_p square is already an exact
existence witness for a Q_p point. Thus no approximate p-adic root is exported.
Unknown searches never contribute a guessed zero to the pairing.
"""
from dataclasses import dataclass
from itertools import count
from math import gcd, lcm

import sympy as sp

from .arithmetic import nonnegative_integer, prime, prime_support
from .cassels_tate_forms import FisherTriple, binary_value
from .geometry import _rat
from .hilbert import local_square, hilbert_bit


def pairing_support(triple):
    """Fisher Remark 3.3: small primes, bad reduction, denominators and content.

    Away from these primes, g1 and gamma are integral, disc(g1) and content
    of gamma are units, and the first Hilbert argument is a unit. For p>=11
    the Hasse bound supplies a smooth residue point outside gamma=0.
    """
    g = triple.quartics[0].f
    gamma = triple.gamma
    t = sp.Symbol('t')
    disc = sp.discriminant(sum(c*t**i for i, c in enumerate(g)), t)
    denominator = lcm(*(int(c.q) for c in gamma))
    content = gcd(*(abs(int(c*denominator)) for c in gamma))
    return prime_support((2, 3, 5, 7, disc, triple.hilbert_parameter,
                          denominator, content, *(c.q for c in g)))


def _real_samples(g):
    """Rational samples in every sign interval, using exact root isolation."""
    t = sp.Symbol('t')
    p = sp.Poly(sum(c*t**i for i, c in enumerate(g)), t, domain=sp.QQ)
    intervals = [list(I) for I, _ in p.intervals()]
    for i in range(len(intervals)-1):
        while intervals[i][1] >= intervals[i+1][0]:
            for j in (i, i+1):
                a, b = intervals[j]
                if a != b:
                    intervals[j] = list(p.refine_root(a, b, steps=1))
    if not intervals:
        gaps = [(-1, 1)]
    else:
        gaps = [(intervals[0][0]-2, intervals[0][0]-1)]
        gaps += [(a[1], b[0]) for a, b in zip(intervals, intervals[1:])]
        gaps += [(intervals[-1][1]+1, intervals[-1][1]+2)]
    # Three distinct samples avoid the at most two zeros of a nonzero gamma.
    for a, b in gaps:
        for n in (1, 2, 3):
            yield sp.Rational((4-n)*a+n*b, 4), sp.S.One


def _finite_samples(p):
    """A dense rational sequence on P^1(Q_p), favoring small heights.

    A short low-height rational sweep reaches p-adic neighborhoods which may
    require a very large integer in either of the two standard affine charts.
    The subsequent unbounded chart traversal is retained, so this is purely
    an optimization and not a probabilistic or bounded local test.
    """
    yield sp.S.One, sp.S.Zero
    for n in range(-8, 9):
        yield sp.Integer(n), sp.S.One
        if n:
            yield sp.S.One, sp.Integer(p*n)
    # A small Farey-style sweep: every nonzero primitive rational denominator
    # up to 64 occurs; no unsupported local-solubility claim is made by this.
    for height in range(1, 65):
        for b in range(1, height+1):
            if b == height:
                numerators = range(-height, height+1)
            else:
                numerators = (-height, height)
            for a in numerators:
                if gcd(a, b) == 1:
                    yield sp.Integer(a), sp.Integer(b)
    # Unbounded original coverage of the affine and infinity p-adic charts.
    for n in count(9):
        yield sp.Integer(n), sp.S.One
        yield sp.Integer(-n), sp.S.One
        yield sp.S.One, sp.Integer(p*n)
        yield sp.S.One, sp.Integer(-p*n)

@dataclass(frozen=True)
class PairingLocalEvidence:
    place: object
    point: tuple | None
    bit: int | None
    trials: int
    kind: str

    def verify(self, triple):
        try:
            if type(self.trials) is not int or self.trials < 0 or type(self.bit) is not int or self.bit not in (0, 1):
                return False
            place = self.place if self.place == 'real' else prime(self.place)
            parameter = triple.hilbert_parameter
            if self.kind == 'square_parameter':
                return self.point is None and self.bit == 0 and local_square(parameter, place)
            if self.kind != 'local_point' or self.point is None or len(self.point) != 2:
                return False
            a, b = map(_rat, self.point)
            if not (a or b):
                return False
            value = binary_value(triple.quartics[0].f, a, b)
            gamma = binary_value(triple.gamma, a, b)
            return (bool(value) and bool(gamma) and local_square(value, place)
                    and hilbert_bit(parameter, gamma, place) == self.bit)
        except (TypeError, ValueError):
            return False


def local_pairing(triple, place, *, max_trials=None):
    if not isinstance(triple, FisherTriple):
        raise TypeError('expected a verified FisherTriple construction')
    place = place if place == 'real' else prime(place)
    limit = None if max_trials is None else nonnegative_integer(max_trials, 'max_trials')
    if limit == 0:
        return PairingLocalEvidence(place, None, None, 0, 'unknown')
    if local_square(triple.hilbert_parameter, place):
        return PairingLocalEvidence(place, None, 0, 0, 'square_parameter')
    samples = _real_samples(triple.quartics[0].f) if place == 'real' else _finite_samples(place)
    trials = 0
    for a, b in samples:
        if limit is not None and trials >= limit:
            return PairingLocalEvidence(place, None, None, trials, 'unknown')
        trials += 1
        value = binary_value(triple.quartics[0].f, a, b)
        gamma = binary_value(triple.gamma, a, b)
        if value and gamma and local_square(value, place):
            bit = hilbert_bit(triple.hilbert_parameter, gamma, place)
            return PairingLocalEvidence(place, (a, b), bit, trials, 'local_point')
    # Only the real iterator is finite. For a locally soluble smooth quartic
    # at least one open real interval exists; failure is not a pairing value.
    raise ArithmeticError('no real point off the evaluation divisor; check the input Selmer claim')


@dataclass(frozen=True)
class FisherEvaluation:
    triple: FisherTriple
    support: tuple
    local: tuple

    @property
    def complete(self):
        places = ('real', *self.support)
        return (tuple(r.place for r in self.local) == places
                and all(type(r.bit) is int and r.bit in (0, 1) for r in self.local))

    @property
    def bit(self):
        if not self.complete:
            raise ValueError('unresolved local evaluations have no pairing value')
        return sum(r.bit for r in self.local) % 2

    def verify(self):
        return (self.complete and self.triple.verify()
                and self.support == pairing_support(self.triple)
                and all(r.verify(self.triple) for r in self.local))


def evaluate_fisher(triple, *, max_local_trials=None):
    if not triple.verify():
        raise ValueError('incorrect invariant identities')
    support = pairing_support(triple)
    reports = tuple(local_pairing(triple, p, max_trials=max_local_trials)
                    for p in ('real', *support))
    return FisherEvaluation(triple, support, reports)
