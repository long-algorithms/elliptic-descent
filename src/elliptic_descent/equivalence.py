"""Fixed-target scalar-square equivalences and supported orientation policies.

For degrees (4) and (2,2), a scalar-square stabilizer cannot reverse the norm
witness. With a rational root, a linear-factor sign change identifies the two
orientations. Its CRT multiplier is retained as exact evidence. See
docs/rational-root-descent.md for the proof and the precise quotient.
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
    reversing_multiplier: tuple | None = None

    @property
    def factorization_type(self):
        return tuple(sorted((len(g)-1 for g in self.factors), reverse=True))

    @property
    def multiplicity(self):
        if self.factorization_type in ((4,), (2, 2)):
            return 2
        if self.factorization_type in ((3, 1), (2, 1, 1), (1, 1, 1, 1)):
            return 1
        raise NotImplementedError("expected a squarefree quartic factorization")

    def verify(self, f):
        try:
            algebra = EtaleQuarticAlgebra(f)
            if algebra.factors != self.factors:
                return False
            if self.multiplicity == 2:
                return self.reversing_multiplier is None
            v = _four(self.reversing_multiplier, "orientation multiplier")
            return (algebra.multiply(v, v) == algebra.one
                    and algebra.norm(v) == -1)
        except (TypeError, ValueError, NotImplementedError):
            return False

    def identification(self, model):
        """Witness identifying (xi,r) with (xi,-r), or None if inequivalent.

        Call witness.verify(model, replace(model, norm_witness=-r)) to audit
        the actual polynomial maps. This is not an identification of the
        target's points: it is an isomorphism between two diagrams over C.
        """
        if not self.verify(model.f):
            raise ValueError("orientation evidence does not match the source quartic")
        if self.multiplicity == 2:
            return None
        return ScalarSquareEquivalence(sp.S.One, self.reversing_multiplier)


def orientation_certificate(f) -> OrientationCertificate:
    algebra = f if isinstance(f, EtaleQuarticAlgebra) else EtaleQuarticAlgebra(f)
    if 1 not in algebra.degrees:
        return OrientationCertificate(algebra.factors)
    # Flip only one rational factor. The multiplier has square 1 and norm -1.
    # Keep factors labelled, including separately embedded isomorphic fields.
    index = algebra.degrees.index(1)
    v = algebra.embed_factor_unit(index, (-1,))
    result = OrientationCertificate(algebra.factors, v)
    if not result.verify(algebra.f):
        raise ArithmeticError("orientation-reversing CRT witness failed")
    return result


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
