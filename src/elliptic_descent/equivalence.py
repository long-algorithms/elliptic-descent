"""Fixed-target scalar-square equivalences and supported orientation policies.

For degrees (4) and (2,2), a scalar-square stabilizer cannot reverse the norm
witness. Rational-root patterns need a different rule and are rejected here.
See docs/two-quadratic-descent.md for the proof and the precise quotient.
"""
from dataclasses import dataclass
from typing import Sequence

import sympy as sp

from .etale import EtaleQuarticAlgebra
from .geometry import FourCoverCandidate, U, _four, _rat, construct_four_cover


@dataclass(frozen=True)
class OrientationCertificate:
    """A checkable factorization, with a documented stabilizer theorem.

    This verifies the hypothesis of the orientation lemma, not local solubility
    or the completeness of any candidate enumeration.
    """
    factors: tuple

    @property
    def factorization_type(self):
        return tuple(sorted((len(g)-1 for g in self.factors), reverse=True))

    @property
    def multiplicity(self):
        if self.factorization_type not in ((4,), (2, 2)):
            raise NotImplementedError("orientation rule requires type 4 or 2+2")
        return 2

    def verify(self, f):
        try:
            algebra = EtaleQuarticAlgebra(f)
            return (algebra.factors == self.factors
                    and self.factorization_type in ((4,), (2, 2)))
        except (TypeError, ValueError):
            return False


def orientation_certificate(f) -> OrientationCertificate:
    algebra = f if isinstance(f, EtaleQuarticAlgebra) else EtaleQuarticAlgebra(f)
    if algebra.factorization_type not in ((4,), (2, 2)):
        raise NotImplementedError(
            "complete relative descent currently supports irreducible quartics "
            "and products of two irreducible quadratics; rational-root cases "
            "require different covering-equivalence bookkeeping")
    return OrientationCertificate(algebra.factors)


@dataclass(frozen=True)
class ScalarSquareEquivalence:
    """Witness xi' = c*v^2*xi and r' = c^2*N(v)*r over a fixed quartic.

    The coordinate change is U_old = v*U_new. Its determinant is N(v), and
    the weighted map coordinates scale by (c,c^2,c), not by a common scalar.
    """
    scalar: sp.Rational
    multiplier: tuple

    def coordinate_matrix(self, f):
        algebra = EtaleQuarticAlgebra(f)
        v = _four(self.multiplier, "multiplier")
        if not algebra.is_unit(v):
            raise ValueError("equivalence multiplier must be a unit")
        columns = [algebra.multiply(v, tuple(sp.Integer(i == j) for j in range(4)))
                   for i in range(4)]
        return sp.ImmutableMatrix.hstack(*(sp.ImmutableMatrix(c) for c in columns))

    def verify(self, source: FourCoverCandidate, target: FourCoverCandidate) -> bool:
        try:
            c, v = _rat(self.scalar), _four(self.multiplier, "multiplier")
            if not c or source.f != target.f:
                return False
            algebra = EtaleQuarticAlgebra(source.f)
            if not algebra.is_unit(v):
                return False
            expected = tuple(c*a for a in algebra.multiply(algebra.multiply(v, v), source.xi))
            if target.xi != expected:
                return False
            norm = algebra.norm(v)
            if target.norm_witness != c*c*norm*source.norm_witness:
                return False
            if not source.check_construction() or not target.check_construction():
                return False
            matrix = self.coordinate_matrix(source.f)
            if matrix.det() != norm:
                return False
            substitutions = dict(zip(U, matrix*sp.Matrix(U)))
            if any(sp.expand(new-c*old.subs(substitutions, simultaneous=True)) != 0
                   for new, old in zip(target.q, source.q)):
                return False
            return all(sp.expand(new-c**weight*old.subs(substitutions, simultaneous=True)) == 0
                       for new, old, weight in zip(target.map_xyz, source.map_xyz, (1, 2, 1)))
        except (TypeError, ValueError, ZeroDivisionError):
            return False


def scalar_square_change(source: FourCoverCandidate, scalar, multiplier: Sequence):
    """Construct an equivalent oriented model and its exact witness.

    Preserve the signed norm witness. Replacing it by the positive square root
    could change the map to the fixed quartic when N(multiplier) is negative.
    """
    c, v = _rat(scalar), _four(multiplier, "multiplier")
    algebra = EtaleQuarticAlgebra(source.f)
    if not c or not algebra.is_unit(v):
        raise ValueError("nonzero scalar and invertible multiplier are required")
    xi = tuple(c*a for a in algebra.multiply(algebra.multiply(v, v), source.xi))
    witness = c*c*algebra.norm(v)*source.norm_witness
    target = construct_four_cover(source.f, xi, norm_witness=witness)
    return target, ScalarSquareEquivalence(c, v)
