"""Factorwise Sage square classes in a quartic etale algebra.

Computes A(S,2), not the relative Selmer set or an orientation quotient.
The relative-descent layer adds equivalences and local decisions separately.
The legacy sage_square_classes API remains irreducible-only. See
docs/etale-arithmetic.md for the arithmetic-only contract.
"""
from dataclasses import dataclass
from typing import Callable

import sympy as sp

from ..arithmetic import (FiniteNormSpace, affine_f2, nonnegative_integer,
                          prime, quartic_support, square_class)
from ..etale import EtaleQuarticAlgebra
from ..geometry import _norm
from ..status import ProofStatus
from .sage import _sage, number_field_square_classes


@dataclass(frozen=True)
class SageFactorSquareClasses:
    """One labelled factor; degree-one factors use QQ rather than a fake field."""
    index: int
    polynomial: tuple
    field: object
    vector_space: object
    generators: tuple
    from_vector: Callable
    to_vector: Callable
    arithmetic_proof: ProofStatus
    assumptions: tuple[str, ...]

    @property
    def degree(self):
        return len(self.polynomial)-1

    @property
    def dimension(self):
        return int(self.vector_space.dimension())

    def as_field_element(self, coefficients):
        if len(coefficients) != self.degree:
            raise ValueError("incorrect factor coordinate count")
        # A linear factor coordinate already is the value at its rational root.
        if self.degree == 1:
            return self.field(str(coefficients[0]))
        return sum((self.field(str(c))*self.field.gen()**i
                    for i, c in enumerate(coefficients)), self.field.zero())

    def as_coefficients(self, element):
        element = self.field(element)
        if self.degree == 1:
            return (sp.Rational(str(element)),)
        values = element.list()
        return tuple(sp.Rational(str(c)) for c in values+[0]*(self.degree-len(values)))


@dataclass(frozen=True)
class SageEtaleSquareClasses:
    """Direct sum of factor square-class spaces with global power-basis outputs.

    Bit coordinates concatenate factor bases in algebra.factors order.
    representatives -> coordinates is exact on classes; the other round trip
    is only modulo A* squared, not equality of field elements.
    """
    algebra: EtaleQuarticAlgebra
    space: FiniteNormSpace
    factors: tuple[SageFactorSquareClasses, ...]
    offsets: tuple[int, ...]
    norm_columns: tuple[int, ...]

    @property
    def dimension(self):
        return self.offsets[-1]

    @property
    def factor_dimensions(self):
        return tuple(factor.dimension for factor in self.factors)

    def coordinates(self, element):
        """Class coordinates; reject zero divisors and odd valuations outside S."""
        components = self.algebra.to_components(element)
        if not all(any(c) for c in components):
            raise ValueError("square classes are defined only for units of the etale algebra")
        label = 0
        for factor, coefficients, offset in zip(self.factors, components, self.offsets):
            value = factor.as_field_element(coefficients)
            vector = factor.to_vector(value)
            label |= sum(int(bit) << (offset+j) for j, bit in enumerate(vector))
        return label

    def _label(self, label):
        label = nonnegative_integer(label, "square-class label")
        if label >= 1 << self.dimension:
            raise ValueError("square-class label exceeds the product-space dimension")
        return label

    def representative(self, label):
        """A unit in the original quartic power basis representing a bit vector."""
        label = self._label(label)
        components = []
        for factor, offset in zip(self.factors, self.offsets):
            vector = factor.vector_space([(label >> (offset+j)) & 1
                                          for j in range(factor.dimension)])
            components.append(factor.as_coefficients(factor.from_vector(vector)))
        return self.algebra.from_components(components)

    def norm_coordinates(self, label):
        """Square-class norm in QQ(S,2), encoded as sign then support-prime bits."""
        label = self._label(label)
        result = 0
        for j, column in enumerate(self.norm_columns):
            if label >> j & 1:
                result ^= column
        return result

    def norm_fibre(self):
        """The affine norm constraint only; no local or covering-equivalence test."""
        return affine_f2(self.norm_columns, square_class(self.algebra.f[4], self.space.support))


def sage_etale_square_classes(f, support=None, *, proof=True):
    """Compute A(S,2) for every squarefree quartic factorization type over QQ.

    f is five ascending rational coefficients or an EtaleQuarticAlgebra.
    support is a finite set of rational primes; every prime above it is used
    in each number-field factor. None uses quartic_support(f), while () really
    means empty support. Duplicate or nonprime entries are rejected.

    proof=True requests certified number-field arithmetic. proof=False retains
    GRH dependencies for nonrational factors; wholly split rational arithmetic
    is unconditional even with proof=False. No global Sage proof flag is used.

    No rational-scalar, target-automorphism or orientation quotient is applied.
    Use the separate relative_four_descent interface for local conditions and
    factorization-dependent covering equivalences.
    """
    if type(proof) is not bool:
        raise TypeError("proof must be an explicit boolean")
    algebra = f if isinstance(f, EtaleQuarticAlgebra) else EtaleQuarticAlgebra(tuple(f))
    support = quartic_support(algebra.f) if support is None else tuple(support)
    support = tuple(sorted(prime(p) for p in support))
    if len(set(support)) != len(support):
        raise ValueError("support contains duplicate primes")
    if not algebra.verify_decomposition():
        raise ArithmeticError("etale algebra decomposition verification failed")
    sage = _sage()
    R = sage.PolynomialRing(sage.QQ, "t")
    factors, generators, norm_columns, offsets = [], [], [], [0]
    for index, polynomial in enumerate(algebra.factors):
        degree = len(polynomial)-1
        if degree == 1:
            K = sage.QQ
            V, gens, from_V, to_V = K.selmer_space(
                [sage.ZZ(p) for p in support], sage.ZZ(2), proof=proof)
        else:
            K = sage.NumberField(R([sage.QQ(str(c)) for c in polynomial]), f"theta_{index}")
            V, gens, from_V, to_V = number_field_square_classes(K, support, proof=proof)
        assumptions = (() if proof or degree == 1 else
                       ("GRH for number-field class-group computations",))
        factor = SageFactorSquareClasses(
            index, polynomial, K, V, tuple(gens), from_V, to_V,
            ProofStatus.CONDITIONAL if assumptions else ProofStatus.CERTIFIED,
            assumptions,
        )
        if len(gens) != factor.dimension:
            raise ArithmeticError("factor square-class basis has inconsistent dimension")
        for generator, basis in zip(gens, V.basis()):
            if to_V(generator) != basis or to_V(from_V(basis)) != basis:
                raise ArithmeticError("factor square-class coordinate round trip failed")
            embedded = algebra.embed_factor_unit(index, factor.as_coefficients(generator))
            # Cross-check three paths: native factor norm, CRT/resultants, and
            # the existing quartic multiplication-matrix determinant.
            native_norm = generator if degree == 1 else generator.norm()
            norm = sp.Rational(str(native_norm))
            if algebra.norm(embedded) != norm or _norm(embedded, algebra.f) != norm:
                raise ArithmeticError("embedded square-class generator norm mismatch")
            generators.append(embedded)
            norm_columns.append(square_class(norm, support))
        factors.append(factor)
        offsets.append(offsets[-1]+factor.dimension)
    assumptions = tuple(dict.fromkeys(a for factor in factors for a in factor.assumptions))
    from sage.env import SAGE_VERSION
    space = FiniteNormSpace(
        algebra.f, tuple(generators), support,
        ProofStatus.CONDITIONAL if assumptions else ProofStatus.CERTIFIED,
        assumptions, f"Sage {SAGE_VERSION}: factorwise selmer_space(proof={proof})",
    )
    context = SageEtaleSquareClasses(algebra, space, tuple(factors), tuple(offsets), tuple(norm_columns))
    # This also guards against accidentally embedding a generator additively,
    # or conflating factors which happen to define isomorphic fields.
    for j, generator in enumerate(generators):
        if context.coordinates(generator) != 1 << j:
            raise ArithmeticError("product square-class coordinate round trip failed")
    return context
