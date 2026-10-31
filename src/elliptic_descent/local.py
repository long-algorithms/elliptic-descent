"""Exact, bounded local tests for intersections of two quadrics in P^3.

Reference path: projective residue enumeration and strong Hensel lifting.
This is not the optimized Fisher--Sills algorithm. See docs/local.md.
Every resource exhaustion is UNKNOWN, never an insolubility certificate.
"""
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm

import sympy as sp

from .arithmetic import nonnegative_integer, prime, prime_support
from .geometry import U, _rat
from .status import LocalStatus

MONOMIALS = tuple(combinations_with_replacement(range(4), 2))


@dataclass(frozen=True)
class QuadraticSystem:
    coefficients: tuple

    @classmethod
    def from_quadrics(cls, quadrics):
        if isinstance(quadrics, cls):
            return quadrics
        if hasattr(quadrics, "quadrics"):
            quadrics = quadrics.quadrics
        if len(quadrics) != 2:
            raise ValueError("expected two homogeneous quadrics in U")
        rows = []
        for q in quadrics:
            p = sp.Poly(q, *U)
            if p.is_zero or any(sum(m) != 2 for m, _ in p.terms()):
                raise ValueError("each equation must be a nonzero homogeneous quadric")
            coeffs = [_rat(p.coeff_monomial(U[i]*U[j])) for i, j in MONOMIALS]
            scale = lcm(*(int(c.q) for c in coeffs))
            row = tuple(int(c*scale) for c in coeffs)
            content = gcd(*row)
            sign = 1 if next(c for c in row if c) > 0 else -1
            rows.append(tuple(sign*c//content for c in row))
        return cls(tuple(rows))

    @property
    def key(self):
        return sha256(repr(self.coefficients).encode("ascii")).hexdigest()

    def values(self, point):
        return tuple(sum(c*point[i]*point[j] for c, (i, j) in zip(row, MONOMIALS))
                     for row in self.coefficients)

    def jacobian(self, point):
        rows = []
        for row in self.coefficients:
            gradient = [0]*4
            for c, (i, j) in zip(row, MONOMIALS):
                gradient[i] += c*point[j]
                gradient[j] += c*point[i]
            rows.append(gradient)
        return rows

    def matrices(self):
        """Integral Hessians: 2Q(u)=u^T H u over Z; p=2 uses coefficients."""
        result = []
        for row in self.coefficients:
            H = sp.zeros(4)
            for c, (i, j) in zip(row, MONOMIALS):
                H[i, j] += c
                H[j, i] += c
            result.append(H)
        return tuple(result)


def _point(point):
    if len(point) != 4:
        raise ValueError("point must have four rational coordinates")
    point = tuple(_rat(c) for c in point)
    if not any(point):
        raise ValueError("zero is not a projective point")
    scale = lcm(*(int(c.q) for c in point))
    point = tuple(int(c*scale) for c in point)
    divisor = gcd(*point)
    return tuple(c//divisor for c in point)


def _valuation(n, p):
    if not n:
        return float("inf")  # Used only for comparisons, never arithmetic proof.
    n, v = abs(int(n)), 0
    while n % p == 0:
        n //= p
        v += 1
    return v


@dataclass(frozen=True)
class LocalCertificate:
    system_key: str
    place: int | str
    kind: str
    point: tuple = ()
    pivot: int = 0
    minor: tuple = ()
    depth: int = 0
    pencil: tuple = ()


@dataclass(frozen=True)
class LocalResult:
    place: int | str
    status: LocalStatus
    certificate: LocalCertificate | None = None
    nodes: int = 0
    depth: int = 0
    reason: str = ""

    def verify(self, quadrics, *, max_nodes=1_000_000):
        return verify_local_result(quadrics, self, max_nodes=max_nodes)


def _hensel(system, point, p, pivot):
    values = system.values(point)
    if values == (0, 0):
        return LocalCertificate(system.key, p, "rational_point", tuple(point))
    v = min(_valuation(a, p) for a in values)
    J = system.jacobian(point)
    for i, j in combinations((k for k in range(4) if k != pivot), 2):
        determinant = J[0][i]*J[1][j] - J[0][j]*J[1][i]
        if determinant and v > 2*_valuation(determinant, p):
            return LocalCertificate(system.key, p, "hensel", tuple(point), pivot, (i, j))
    return None


def _digits(base, length):
    # Unlike itertools.product(range(base), ...), this does not materialize
    # a potentially enormous pool before the node budget can be checked.
    for number in range(base**length):
        digits = [0]*length
        for i in range(length-1, -1, -1):
            number, digits[i] = divmod(number, base)
        yield tuple(digits)


def _initial_points(p):
    # Disjoint charts: the first unit coordinate is normalized to 1.
    for pivot in range(4):
        for suffix in _digits(p, 3-pivot):
            yield pivot, (0,)*pivot + (1,) + suffix


def _lifts(frontier, modulus, p):
    for pivot, point in frontier:
        free = tuple(i for i in range(4) if i != pivot)
        for digits in _digits(p, 3):
            lifted = list(point)
            for i, digit in zip(free, digits):
                lifted[i] += modulus*digit
            yield pivot, tuple(lifted)


def padic_solubility(quadrics, p, *, max_depth=6, max_nodes=100_000):
    """Decide soluble/insoluble when a certificate is found; otherwise UNKNOWN.

    max_nodes counts all tested residue points across all depths. Normalizing
    each equation to primitive integral coefficients is safe over Q_p.
    """
    system, p = QuadraticSystem.from_quadrics(quadrics), prime(p)
    max_depth = nonnegative_integer(max_depth, "max_depth")
    max_nodes = nonnegative_integer(max_nodes, "max_nodes")
    frontier, nodes, previous_modulus = [], 0, 1
    for depth in range(1, max_depth+1):
        modulus = previous_modulus*p
        points = _initial_points(p) if depth == 1 else _lifts(frontier, previous_modulus, p)
        next_frontier = []
        for pivot, point in points:
            if nodes >= max_nodes:
                return LocalResult(p, LocalStatus.UNKNOWN, nodes=nodes, depth=depth,
                                   reason="node budget exhausted")
            nodes += 1
            if any(v % modulus for v in system.values(point)):
                continue
            certificate = _hensel(system, point, p, pivot)
            if certificate:
                return LocalResult(p, LocalStatus.SOLUBLE, certificate, nodes, depth)
            next_frontier.append((pivot, point))
        if not next_frontier:
            certificate = LocalCertificate(system.key, p, "empty_modulus", depth=depth)
            return LocalResult(p, LocalStatus.INSOLUBLE, certificate, nodes, depth)
        frontier, previous_modulus = next_frontier, modulus
    return LocalResult(p, LocalStatus.UNKNOWN, nodes=nodes, depth=max_depth,
                       reason="depth budget exhausted")


def _definite(H):
    minors = [H[:i, :i].det() for i in range(1, 5)]
    return all(d > 0 for d in minors) or all((-1)**i*d > 0 for i, d in enumerate(minors, 1))


def real_solubility(quadrics, *, point=None, pencil_bound=3):
    """Exact rational witnesses or definite-pencil obstructions; otherwise unknown.

    The finite pencil search is not a complete real-solubility algorithm.
    """
    system = QuadraticSystem.from_quadrics(quadrics)
    bound = nonnegative_integer(pencil_bound, "pencil_bound")
    points = [_point(point)] if point is not None else product(range(-1, 2), repeat=4)
    for candidate in points:
        if any(candidate) and system.values(candidate) == (0, 0):
            certificate = LocalCertificate(system.key, "real", "rational_point", tuple(candidate))
            return LocalResult("real", LocalStatus.SOLUBLE, certificate)
    if point is not None:
        raise ValueError("supplied point does not lie on the intersection")
    A, B = system.matrices()
    for a, b in product(range(-bound, bound+1), repeat=2):
        if (a or b) and gcd(a, b) == 1 and _definite(a*A+b*B):
            certificate = LocalCertificate(system.key, "real", "definite_pencil", pencil=(a, b))
            return LocalResult("real", LocalStatus.INSOLUBLE, certificate)
    return LocalResult("real", LocalStatus.UNKNOWN, reason="no real certificate found")


def verify_local_result(quadrics, result, *, max_nodes=1_000_000):
    """Check witnesses; negative finite-place evidence uses independent enumeration.

    For empty_modulus this directly enumerates P^3(Z/p^k), rather than calling
    the lifting routine. Insufficient verifier budget returns False, not proof.
    """
    system = QuadraticSystem.from_quadrics(quadrics)
    c = result.certificate
    if c is None or c.system_key != system.key or c.place != result.place:
        return False
    if c.kind == "rational_point":
        try:
            point = _point(c.point)
        except (TypeError, ValueError):
            return False
        return result.status == LocalStatus.SOLUBLE and system.values(point) == (0, 0)
    if c.kind == "definite_pencil":
        if c.place != "real" or result.status != LocalStatus.INSOLUBLE or len(c.pencil) != 2:
            return False
        a, b = map(_rat, c.pencil)
        A, B = system.matrices()
        return _definite(a*A+b*B)
    if c.place == "real":
        return False
    p = prime(c.place)
    if c.kind == "hensel":
        if result.status != LocalStatus.SOLUBLE or len(c.point) != 4 or len(c.minor) != 2:
            return False
        if not (0 <= c.pivot < 4) or c.point[c.pivot] != 1:
            return False
        if any(_rat(x).q != 1 for x in c.point):
            return False
        i, j = c.minor
        if not (0 <= i < j < 4) or c.pivot in (i, j):
            return False
        J = system.jacobian(c.point)
        determinant = J[0][i]*J[1][j] - J[0][j]*J[1][i]
        return bool(determinant and min(_valuation(v, p) for v in system.values(c.point))
                    > 2*_valuation(determinant, p))
    if c.kind != "empty_modulus" or result.status != LocalStatus.INSOLUBLE or c.depth < 1:
        return False
    modulus = p**c.depth
    # The independent checker is exhaustive; reject over-budget certificates
    # before materializing any residue ranges.
    if p**(3*(c.depth-1))*(p**3+p**2+p+1) > max_nodes:
        return False
    nodes = 0
    for pivot in range(4):
        ranges = [range(0, modulus, p) if i < pivot else range(modulus)
                  for i in range(4) if i != pivot]
        free = tuple(i for i in range(4) if i != pivot)
        for coordinates in product(*ranges):
            if nodes >= max_nodes:
                return False
            nodes += 1
            point = [0]*4
            point[pivot] = 1
            for i, value in zip(free, coordinates):
                point[i] = value
            if not any(v % modulus for v in system.values(point)):
                return False
    return True


@dataclass(frozen=True)
class GoodReductionBound:
    primes: tuple[int, ...]
    shift: int
    determinant_polynomial: tuple
    discriminant: int


def bad_primes(quadrics):
    """Conservative finite local set from a regular squarefree determinant pencil.

    At primes outside this set the intersection is a smooth genus-one curve;
    the Hasse bound and Hensel lifting imply Q_p-solubility. Always includes 2.
    A singular pencil is unsupported, not an empty descent.
    """
    A, B = QuadraticSystem.from_quadrics(quadrics).matrices()
    t = sp.Symbol("t")
    for shift in range(5):
        M = A + shift*B
        if M.det():
            g = sp.Poly((t*M+B).det(), t)
            d = int(g.discriminant())
            if not d:
                break
            # Including lc(g) makes the finite-field diagonalization argument
            # valid without having to find a new regular member modulo p.
            primes = prime_support((2, d, g.LC()))
            return GoodReductionBound(primes, shift, tuple(reversed(g.all_coeffs())), d)
    raise NotImplementedError("singular or nonregular determinant pencil")


@dataclass(frozen=True)
class LocalScreen:
    reports: tuple[LocalResult, ...]
    good_reduction: GoodReductionBound | None

    @property
    def status(self):
        if any(r.status == LocalStatus.INSOLUBLE for r in self.reports):
            return LocalStatus.INSOLUBLE
        if any(r.certificate and r.certificate.kind == "rational_point" for r in self.reports):
            return LocalStatus.SOLUBLE
        if self.good_reduction is not None:
            required = {"real", *self.good_reduction.primes}
            if {r.place for r in self.reports if r.status == LocalStatus.SOLUBLE} >= required:
                return LocalStatus.SOLUBLE
        return LocalStatus.UNKNOWN


def everywhere_local_solubility(quadrics, *, point=None, max_depth=6, max_nodes=100_000):
    """Certified local screening of one model, NOT exhaustive global descent."""
    real = real_solubility(quadrics, point=point)
    if real.status == LocalStatus.INSOLUBLE or (real.certificate and real.certificate.kind == "rational_point"):
        return LocalScreen((real,), None)
    bound = bad_primes(quadrics)
    reports = [real]
    for p in bound.primes:
        report = padic_solubility(quadrics, p, max_depth=max_depth, max_nodes=max_nodes)
        reports.append(report)
        if report.status == LocalStatus.INSOLUBLE or (report.certificate and report.certificate.kind == "rational_point"):
            break
    return LocalScreen(tuple(reports), bound)
