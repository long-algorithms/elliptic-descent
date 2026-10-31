"""Exact finite-support square classes and affine norm equations.

A norm fibre in K(S,2) is NOT a relative 4-Selmer set. Rational-scalar
quotients, orientations and covering equivalence remain separate obligations.
See docs/arithmetic.md and Fisher (2008), Section 2.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import islice
from math import isqrt

import sympy as sp

from .geometry import _five, _four, _mul, _norm, _rat, _sqrt_rational
from .status import ProofStatus


def nonnegative_integer(value, name):
    value = _rat(value)
    if value.q != 1 or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def prime(value):
    """Prove primality; large primes use Sage/PARI with proof=True.

    Pure-Python use is supported below 2**64. The complete descent already
    requires Sage, so it has no prime-size restriction. Failed certification
    raises instead of treating a probable prime as a proven prime.
    """
    return _proven_prime(nonnegative_integer(value, "prime"))


@lru_cache(maxsize=2048)
def _proven_prime(p):
    if p < 2**64:
        valid = sp.isprime(p)
    else:
        try:
            from sage.all import ZZ
        except ImportError as error:
            raise NotImplementedError("primality proofs above 64 bits require SageMath") from error
        valid = ZZ(p).is_prime(proof=True)
    if not valid:
        raise ValueError("expected a prime")
    return p


def prime_support(values):
    """Factor exact rationals, checking reconstruction and every prime factor."""
    support = set()
    for value in values:
        value = _rat(value)
        if not value:
            raise ValueError("zero has no finite prime support")
        for n in (abs(int(value.p)), int(value.q)):
            factors = sp.factorint(n)
            if sp.prod(p**e for p, e in factors.items()) != n:
                raise ArithmeticError("factorization reconstruction failed")
            support.update(prime(p) for p in factors)
    return tuple(sorted(support))


def quartic_support(f):
    """Conservative bad-presentation support, NOT a Selmer completeness flag.

    Outside S: 2, lc(f), all coefficient denominators and disc(f) are units.
    See the necessary finite-support lemma in docs/arithmetic.md.
    """
    f = _five(f)
    t = sp.Symbol("t")
    disc = sp.discriminant(sum(c*t**i for i, c in enumerate(f)), t)
    return prime_support((2, f[4], disc, *(c.q for c in f)))


def square_class(value, support):
    """Bits for sign and S-valuations; reject odd valuation outside S.

    No factorization of value is required: after stripping S, the remaining
    numerator and denominator must both be integer squares.
    """
    value = _rat(value)
    if not value:
        raise ValueError("zero is not a square class")
    support = tuple(prime(p) for p in support)
    if len(set(support)) != len(support):
        raise ValueError("support contains duplicate primes")
    n, d = abs(int(value.p)), int(value.q)
    bits = 1 if value < 0 else 0
    for i, p in enumerate(support, 1):
        parity = 0
        while n % p == 0:
            n //= p
            parity ^= 1
        while d % p == 0:
            d //= p
            parity ^= 1
        bits |= parity << i
    if isqrt(n)**2 != n or isqrt(d)**2 != d:
        raise ValueError("square class has odd valuation outside the declared support")
    return bits


@dataclass(frozen=True)
class AffineF2:
    """Solutions are particular + span(basis), encoded as nonnegative bitsets."""
    nvariables: int
    particular: int | None
    basis: tuple[int, ...]

    @property
    def cardinality(self):
        return 0 if self.particular is None else 1 << len(self.basis)

    def elements(self):
        if self.particular is None:
            return
        for mask in range(self.cardinality):
            v = self.particular
            for i, b in enumerate(self.basis):
                if mask >> i & 1:
                    v ^= b
            yield v


def affine_f2(columns, target):
    """Solve XOR_j(columns[j]*x_j)=target by exact elimination over F_2."""
    columns = tuple(nonnegative_integer(c, "column") for c in columns)
    target = nonnegative_integer(target, "target")
    n = len(columns)
    width = max((target.bit_length(), *(c.bit_length() for c in columns)))
    rows = [sum(((c >> k) & 1) << j for j, c in enumerate(columns))
            | (((target >> k) & 1) << n) for k in range(width)]
    pivots = []
    for j in range(n):
        k = next((i for i in range(len(pivots), len(rows)) if rows[i] >> j & 1), None)
        if k is None:
            continue
        r = len(pivots)
        rows[r], rows[k] = rows[k], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> j & 1:
                rows[i] ^= rows[r]
        pivots.append(j)
    if any(row == 1 << n for row in rows):
        return AffineF2(n, None, ())
    particular = sum(((rows[i] >> n) & 1) << j for i, j in enumerate(pivots))
    basis = tuple((1 << j) | sum(((rows[i] >> j) & 1) << k
                                for i, k in enumerate(pivots))
                  for j in range(n) if j not in pivots)
    return AffineF2(n, particular, basis)


@dataclass(frozen=True)
class NormCandidate:
    label: int
    xi: tuple
    norm_witness: sp.Rational


@dataclass(frozen=True)
class CandidateBatch:
    candidates: tuple[NormCandidate, ...]
    total_count: int
    arithmetic_proof: ProofStatus
    assumptions: tuple[str, ...]

    @property
    def norm_fibre_exhausted(self):
        """Only the specified finite norm fibre, never the relative Selmer set."""
        return len(self.candidates) == self.total_count


@dataclass(frozen=True)
class FiniteNormSpace:
    """A declared basis of K(S,2); completeness is inherited from the backend.

    Directly constructing this object does not certify basis independence.
    Use sage_square_classes() for a checked number-field backend.
    """
    f: tuple
    generators: tuple
    support: tuple[int, ...]
    arithmetic_proof: ProofStatus = ProofStatus.UNVERIFIED
    assumptions: tuple[str, ...] = ()
    backend: str = "supplied generators (basis not certified)"

    def __post_init__(self):
        object.__setattr__(self, "f", _five(self.f))
        object.__setattr__(self, "generators", tuple(_four(g, "generator") for g in self.generators))
        object.__setattr__(self, "support", tuple(sorted(prime(p) for p in self.support)))
        if len(set(self.support)) != len(self.support):
            raise ValueError("duplicate support primes")
        if self.arithmetic_proof == ProofStatus.CONDITIONAL and not self.assumptions:
            raise ValueError("conditional arithmetic requires explicit assumptions")
        if self.arithmetic_proof == ProofStatus.CERTIFIED and self.assumptions:
            raise ValueError("certified arithmetic cannot carry unproved assumptions")

    def norm_fibre(self):
        columns = [square_class(_norm(g, self.f), self.support) for g in self.generators]
        return affine_f2(columns, square_class(self.f[4], self.support))

    def enumerate(self, *, limit=4096):
        """Enumerate at most limit candidates; truncation is explicit in the batch."""
        limit = nonnegative_integer(limit, "limit")
        fibre = self.norm_fibre()
        result = []
        for label in islice(fibre.elements(), limit):
            xi = (sp.S.One, sp.S.Zero, sp.S.Zero, sp.S.Zero)
            for j, g in enumerate(self.generators):
                if label >> j & 1:
                    xi = _mul(xi, g, self.f)
            r = _sqrt_rational(_norm(xi, self.f) / self.f[4])
            result.append(NormCandidate(label, xi, r))
        return CandidateBatch(tuple(result), fibre.cardinality,
                              self.arithmetic_proof, self.assumptions)
