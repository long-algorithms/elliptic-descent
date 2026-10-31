"""Exact quartic etale algebras, with explicit Chinese-remainder maps.

The geometric quartic is retained unchanged. Its monic irreducible factors
label the components, including distinct factors defining isomorphic fields.
This module requires only SymPy; no Selmer or descent completeness is claimed.
"""
from dataclasses import dataclass
from functools import cached_property
from typing import Sequence

import sympy as sp

from .geometry import _five, _four, _rat

RationalVector = tuple[sp.Rational, ...]
_T = sp.Symbol("t")


def _polynomial(coefficients):
    return sp.Poly.from_list(list(reversed(coefficients)), _T, domain=sp.QQ)


def _coefficients(polynomial, length):
    return tuple(sp.Rational(polynomial.nth(i)) for i in range(length))


@dataclass(frozen=True)
class EtaleQuarticAlgebra:
    """Q[t]/(f) and its labelled product of fields, for any squarefree quartic.

    All coefficient vectors use ascending powers. Global vectors have length
    four; a component vector has length equal to that factor's degree. Linear
    factors have one rational coordinate (the value at the root).

    Components are ordered by degree, then lexicographically by the ascending
    coefficients of their monic factors. Scalar rescaling of f preserves this
    ordering, but the original geometric f is never replaced by a monic curve.
    """

    f: tuple

    def __post_init__(self):
        object.__setattr__(self, "f", _five(self.f))

    @cached_property
    def _modulus(self):
        return _polynomial(self.f).monic()

    @cached_property
    def factors(self) -> tuple[RationalVector, ...]:
        """Monic irreducible factors; equal-degree isomorphic fields stay separate."""
        _, factorization = self._modulus.factor_list()
        if any(exponent != 1 for _, exponent in factorization):
            raise ArithmeticError("quartic factorization was not squarefree")
        factors = tuple(sorted(
            (_coefficients(g.monic(), int(g.degree())+1) for g, _ in factorization),
            key=lambda coefficients: (len(coefficients), coefficients),
        ))
        product = sp.Poly(1, _T, domain=sp.QQ)
        for g in factors:
            product *= _polynomial(g)
        if product != self._modulus:
            raise ArithmeticError("quartic factorization reconstruction failed")
        return factors

    @cached_property
    def _factor_polynomials(self):
        return tuple(_polynomial(g) for g in self.factors)

    @property
    def degrees(self) -> tuple[int, ...]:
        return tuple(len(g)-1 for g in self.factors)

    @property
    def factorization_type(self) -> tuple[int, ...]:
        """Descending degree partition, independent of component ordering."""
        return tuple(sorted(self.degrees, reverse=True))

    @property
    def zero(self) -> RationalVector:
        return (sp.S.Zero,)*4

    @property
    def one(self) -> RationalVector:
        return (sp.S.One, sp.S.Zero, sp.S.Zero, sp.S.Zero)

    def scalar(self, value) -> RationalVector:
        return (_rat(value), sp.S.Zero, sp.S.Zero, sp.S.Zero)

    @cached_property
    def idempotents(self) -> tuple[RationalVector, ...]:
        """CRT idempotent e_i is 1 in component i and 0 in every other one."""
        values = []
        for g in self._factor_polynomials:
            h = self._modulus.exquo(g)
            e = (h*sp.invert(h, g)).rem(self._modulus)
            values.append(_coefficients(e, 4))
        return tuple(values)

    def to_components(self, element: Sequence) -> tuple[RationalVector, ...]:
        polynomial = _polynomial(_four(element, "element"))
        return tuple(_coefficients(polynomial.rem(g), int(g.degree()))
                     for g in self._factor_polynomials)

    def from_components(self, components: Sequence[Sequence]) -> RationalVector:
        if len(components) != len(self.factors):
            raise ValueError("expected one coordinate vector for each labelled factor")
        polynomial = sp.Poly(0, _T, domain=sp.QQ)
        for values, degree, e in zip(components, self.degrees, self.idempotents):
            if len(values) != degree:
                raise ValueError(f"component must have exactly {degree} ascending coefficients")
            values = tuple(_rat(c) for c in values)
            polynomial += _polynomial(values)*_polynomial(e)
        return _coefficients(polynomial.rem(self._modulus), 4)

    def multiply(self, left: Sequence, right: Sequence) -> RationalVector:
        """Polynomial-remainder multiplication, independent of geometry._mul."""
        a = _polynomial(_four(left, "left operand"))
        b = _polynomial(_four(right, "right operand"))
        return _coefficients((a*b).rem(self._modulus), 4)

    def is_unit(self, element: Sequence) -> bool:
        return all(any(component) for component in self.to_components(element))

    def inverse(self, element: Sequence) -> RationalVector:
        components = self.to_components(element)
        if not all(any(component) for component in components):
            raise ZeroDivisionError("element is not a unit in the product algebra")
        inverses = tuple(_coefficients(sp.invert(_polynomial(c), g), int(g.degree()))
                         for c, g in zip(components, self._factor_polynomials))
        return self.from_components(inverses)

    def component_norms(self, element: Sequence) -> RationalVector:
        """Exact norms via resultants in each monic factor, including rational ones."""
        return tuple(sp.Rational(g.resultant(_polynomial(c)))
                     for g, c in zip(self._factor_polynomials, self.to_components(element)))

    def norm(self, element: Sequence) -> sp.Rational:
        """Product of component norms, independent of the global determinant norm."""
        return sp.Rational(sp.prod(self.component_norms(element)))

    def embed_factor_unit(self, index: int, element: Sequence) -> RationalVector:
        """Embed a factor unit multiplicatively: (1,...,element,...,1).

        The additive lift e_i*element is a zero divisor in a nontrivial product
        and must NOT be used as a square-class generator.
        """
        if isinstance(index, bool) or not isinstance(index, (int, sp.Integer)):
            raise TypeError("factor index must be an integer")
        if not 0 <= index < len(self.factors):
            raise IndexError("factor index out of range")
        components = [(sp.S.One,)+(sp.S.Zero,)*(degree-1) for degree in self.degrees]
        if len(element) != self.degrees[index]:
            raise ValueError("incorrect number of factor coordinates")
        element = tuple(_rat(c) for c in element)
        if not any(element):
            raise ValueError("a square-class generator must be a unit")
        components[index] = element
        return self.from_components(components)

    def verify_decomposition(self) -> bool:
        """Audit factorization, idempotents and the two exact rational linear maps."""
        if sum(self.degrees) != 4 or any(not g.is_irreducible for g in self._factor_polynomials):
            return False
        if tuple(map(sum, zip(*self.idempotents))) != self.one:
            return False
        for i, e in enumerate(self.idempotents):
            for j, other in enumerate(self.idempotents):
                if self.multiply(e, other) != (e if i == j else self.zero):
                    return False
        # Four basis vectors test each composition as a Q-linear map.
        for i in range(4):
            basis = tuple(sp.Integer(i == j) for j in range(4))
            if self.from_components(self.to_components(basis)) != basis:
                return False
        for i, degree in enumerate(self.degrees):
            for j in range(degree):
                basis = [tuple(sp.S.Zero for _ in range(d)) for d in self.degrees]
                basis[i] = tuple(sp.Integer(k == j) for k in range(degree))
                if self.to_components(self.from_components(basis)) != tuple(basis):
                    return False
        return True
