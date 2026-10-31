"""Certified factorwise L(S,2) for the monic cubic used in full 2-descent.

This is a finite ambient space, NOT yet the Selmer group. Its entire norm
kernel is screened locally by two_selmer.py. No PARI ell2cover basis is assumed.
"""
from dataclasses import dataclass

import sympy as sp

from .sage import _sage, number_field_square_classes
from .sage_etale import SageFactorSquareClasses
from ..arithmetic import affine_f2, nonnegative_integer, square_class
from ..cubic_arithmetic import CubicEtaleAlgebra, cubic_norm
from ..status import ProofStatus


@dataclass(frozen=True)
class CubicSquareClasses:
    algebra: CubicEtaleAlgebra
    support: tuple
    factors: tuple
    offsets: tuple
    generators: tuple
    norm_columns: tuple
    arithmetic_proof: ProofStatus
    assumptions: tuple
    backend: str

    @property
    def dimension(self):
        return len(self.generators)

    @property
    def norm_kernel(self):
        return affine_f2(self.norm_columns, 0)

    def coordinates(self, element):
        parts = self.algebra.to_components(element)
        if not all(any(c) for c in parts):
            raise ValueError('square classes require an invertible cubic-algebra element')
        label = 0
        for f, c, offset in zip(self.factors, parts, self.offsets):
            v = f.to_vector(f.as_field_element(c))
            label |= sum(int(bit) << (offset+i) for i, bit in enumerate(v))
        return label

    def representative(self, label):
        label = nonnegative_integer(label, 'ambient label')
        if label >= 1 << self.dimension:
            raise ValueError('label exceeds ambient dimension')
        result = self.algebra.one
        for i, g in enumerate(self.generators):
            if label >> i & 1:
                result = self.algebra.multiply(result, g)
        return result

    def norm_label(self, label):
        label = nonnegative_integer(label, 'ambient label')
        if label >= 1 << self.dimension:
            raise ValueError('label exceeds ambient dimension')
        result = 0
        for i, column in enumerate(self.norm_columns):
            if label >> i & 1:
                result ^= column
        return result

    def verify_basis(self):
        """Exact basis-coordinate/norm checks, not a second class-group proof."""
        if not self.algebra.verify() or self.offsets[0] != 0 or self.offsets[-1] != self.dimension:
            return False
        if sum(f.dimension for f in self.factors) != self.dimension:
            return False
        for i, g in enumerate(self.generators):
            if self.coordinates(g) != 1 << i:
                return False
            norm = self.algebra.norm(g)
            if cubic_norm(g, self.algebra.f) != norm or square_class(norm, self.support) != self.norm_columns[i]:
                return False
        return True


def cubic_square_classes(algebra, *, proof=True):
    if type(proof) is not bool:
        raise TypeError('proof must be a boolean')
    if not isinstance(algebra, CubicEtaleAlgebra):
        algebra = CubicEtaleAlgebra(algebra)
    sage = _sage()
    R = sage.PolynomialRing(sage.QQ, 'theta')
    support = algebra.support
    factors, generators, norms, offsets = [], [], [], [0]
    for i, polynomial in enumerate(algebra.factors):
        degree = len(polynomial)-1
        if degree == 1:
            K = sage.QQ
            V, gens, from_V, to_V = K.selmer_space([sage.ZZ(p) for p in support], sage.ZZ(2), proof=proof)
        else:
            K = sage.NumberField(R([sage.QQ(str(c)) for c in polynomial]), f'root_{i}')
            V, gens, from_V, to_V = number_field_square_classes(K, support, proof=proof)
        assumptions = () if proof or degree == 1 else ('GRH for number-field class-group computations',)
        factor = SageFactorSquareClasses(i, polynomial, K, V, tuple(gens), from_V, to_V,
            ProofStatus.CONDITIONAL if assumptions else ProofStatus.CERTIFIED, assumptions)
        if len(gens) != factor.dimension:
            raise ArithmeticError('inconsistent cubic factor basis')
        for g, basis in zip(gens, V.basis()):
            if to_V(g) != basis or to_V(from_V(basis)) != basis:
                raise ArithmeticError('cubic factor coordinate round trip failed')
            embedded = algebra.embed_unit(i, factor.as_coefficients(g))
            norm = sp.Rational(str(g if degree == 1 else g.norm()))
            if algebra.norm(embedded) != norm or cubic_norm(embedded, algebra.f) != norm:
                raise ArithmeticError('factor/CRT/determinant norms disagree')
            generators.append(embedded)
            norms.append(square_class(norm, support))
        factors.append(factor)
        offsets.append(offsets[-1]+factor.dimension)
    assumptions = tuple(dict.fromkeys(a for f in factors for a in f.assumptions))
    from sage.env import SAGE_VERSION
    context = CubicSquareClasses(algebra, support, tuple(factors), tuple(offsets), tuple(generators),
                                  tuple(norms), ProofStatus.CONDITIONAL if assumptions else ProofStatus.CERTIFIED,
                                  assumptions, f'Sage {SAGE_VERSION}: cubic factorwise selmer_space(proof={proof})')
    if not context.verify_basis():
        raise ArithmeticError('cubic square-class basis verification failed')
    return context
