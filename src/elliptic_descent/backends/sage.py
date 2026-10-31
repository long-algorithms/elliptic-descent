"""Lazy Sage/PARI adapters using documented public APIs.

K(S,2) includes the class-group contribution; it is not just an S-unit group.
proof=True is passed explicitly, never inherited from Sage's global setting.
"""
from dataclasses import dataclass

import sympy as sp

from ..arithmetic import FiniteNormSpace, prime, quartic_support
from ..geometry import _five, _rat
from ..status import ProofStatus


def _sage():
    try:
        import sage.all as sage
    except ImportError as error:
        raise ImportError("SageMath is required; run with 'sage -python'") from error
    return sage


def number_field_square_classes(K, support, *, proof=True):
    """K(S,2), with S containing ALL primes above each supplied rational prime."""
    if type(proof) is not bool:
        raise TypeError("proof must be an explicit boolean")
    support = tuple(sorted(set(prime(p) for p in support)))
    S = [P for p in support for P in K.primes_above(p)]
    V, generators, from_V, to_V = K.selmer_space(S, 2, proof=proof)
    if len(generators) != V.dimension():
        raise ArithmeticError("backend returned an inconsistent square-class basis")
    for g, e in zip(generators, V.basis()):
        if to_V(g) != e:
            raise ArithmeticError("square-class coordinate round trip failed")
    return V, generators, from_V, to_V


@dataclass(frozen=True)
class SageSquareClasses:
    space: FiniteNormSpace
    field: object
    vector_space: object
    from_vector: object
    to_vector: object

    def coordinates(self, xi):
        """Bit coordinates of an element in this K(S,2); rejects outsiders."""
        K = self.field
        value = sum(K(str(_rat(c)))*K.gen()**i for i, c in enumerate(xi))
        vector = self.to_vector(value)
        return sum(int(v) << i for i, v in enumerate(vector))


def sage_square_classes(f, support=None, *, proof=True):
    """Build a finite norm space for an irreducible quartic over Q.

    Normalizing the field polynomial to monic does not alter the geometric f.
    Reducible quartics require an etale-product backend (not yet implemented).
    """
    sage = _sage()
    f = _five(f)
    support = quartic_support(f) if support is None else tuple(support)
    R = sage.PolynomialRing(sage.QQ, "t")
    g = R([sage.QQ(str(c)) for c in f])
    if not g.is_irreducible():
        raise NotImplementedError("Sage arithmetic adapter currently requires an irreducible quartic")
    K = sage.NumberField(g.monic(), "theta")
    V, gens, from_V, to_V = number_field_square_classes(K, support, proof=proof)
    coefficients = tuple(tuple(sp.Rational(str(c)) for c in (a.list() + [0]*4)[:4]) for a in gens)
    from sage.env import SAGE_VERSION
    space = FiniteNormSpace(f, coefficients, tuple(support),
                            ProofStatus.CERTIFIED if proof else ProofStatus.CONDITIONAL,
                            () if proof else ("GRH for number-field class-group computations",),
                            f"Sage {SAGE_VERSION}: selmer_space(proof={proof})")
    return SageSquareClasses(space, K, V, from_V, to_V)


@dataclass(frozen=True)
class PariTwoCover:
    """PARI basis entry, with map identity checked independently in SymPy.

    No certification of the entire PARI basis is inferred from this identity.
    Degree below four is preserved, not silently converted to a quartic.
    """
    f: tuple
    ainvs: tuple
    map_xy: tuple
    source: str = "PARI ell2cover (basis metadata not independently certified)"

    def check_map(self):
        x, y = sp.symbols("x y")
        f = sum(c*x**i for i, c in enumerate(self.f))
        a1, a2, a3, a4, a6 = self.ainvs
        X, Y = self.map_xy
        ideal = sp.groebner((y*y-f,), y, x, domain=sp.QQ)
        for coordinate in (X, Y):
            denominator = sp.fraction(sp.cancel(coordinate))[1]
            if ideal.reduce(denominator)[1] == 0:
                return False
        residual = sp.cancel(Y*Y+a1*X*Y+a3*Y-X**3-a2*X**2-a4*X-a6)
        numerator = sp.fraction(residual)[0]
        return ideal.reduce(numerator)[1] == 0


def pari_two_covers(ainvs):
    """Import PARI's basis quartics and maps; this does not enumerate all classes."""
    sage = _sage()
    ainvs = tuple(_rat(a) for a in ainvs)
    if len(ainvs) == 2:
        ainvs = (sp.S.Zero,)*3 + ainvs
    if len(ainvs) != 5:
        raise ValueError("expected two or five Weierstrass coefficients")
    E = sage.EllipticCurve([sage.QQ(str(a)) for a in ainvs])
    raw = sage.pari(list(E.ainvs())).ellinit().ell2cover()
    x, y = sp.symbols("x y")
    result = []
    for R, P in raw:
        f = tuple(sp.Rational(str(R.polcoef(i))) for i in range(int(R.poldegree())+1))
        coordinates = tuple(sp.sympify(str(P[i]).replace("^", "**"), locals={"x": x, "y": y})
                            for i in range(2))
        cover = PariTwoCover(f, ainvs, coordinates)
        if not cover.check_map():
            raise ArithmeticError("PARI covering map failed exact verification")
        result.append(cover)
    return tuple(result)
