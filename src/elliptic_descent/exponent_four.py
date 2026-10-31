"""Exact arithmetic in a witnessed abelian extension of exponent at most four.

This module alone does NOT identify a group with a Selmer group. Its caller
must link the normal forms and generator relations to actual covering classes.
All bit vectors are local to the recorded kernel and quotient bases.
"""
from dataclasses import dataclass
from functools import cached_property

from .arithmetic import affine_f2, nonnegative_integer
from .quotient import row_basis
from .two_selmer import F2Labels


@dataclass(frozen=True)
class ExponentFourPresentation:
    """Kernel generators e_i, lifts f_j, relations 2e_i=0, 2f_j=d_j.

    A normal form is (a,b), where a is a kernel bitset and b a lift bitset.
    Carries, not XOR alone, implement addition. No geometry is inferred here.
    """
    kernel_dimension: int
    doubles: tuple

    def __post_init__(self):
        k = nonnegative_integer(self.kernel_dimension, 'kernel dimension')
        doubles = tuple(nonnegative_integer(v, 'doubling column') for v in self.doubles)
        if any(v >= 1 << k for v in doubles):
            raise ValueError('doubling column exceeds the kernel dimension')
        object.__setattr__(self, 'kernel_dimension', k)
        object.__setattr__(self, 'doubles', doubles)

    @property
    def image_dimension(self):
        return len(self.doubles)

    @property
    def cardinality(self):
        return 1 << (self.kernel_dimension + self.image_dimension)

    def point(self, value):
        if len(value) != 2:
            raise ValueError('expected a (kernel bitset, image bitset) normal form')
        a, b = (nonnegative_integer(v, 'normal-form coordinate') for v in value)
        if a >= 1 << self.kernel_dimension or b >= 1 << self.image_dimension:
            raise ValueError('normal-form coordinate is outside this presentation')
        return a, b

    def doubling_column(self, b):
        _, b = self.point((0, b))
        result = 0
        for j, d in enumerate(self.doubles):
            if b >> j & 1:
                result ^= d
        return result

    def add(self, left, right):
        a, b = self.point(left)
        c, d = self.point(right)
        return a ^ c ^ self.doubling_column(b & d), b ^ d

    def negate(self, value):
        a, b = self.point(value)
        return a ^ self.doubling_column(b), b

    def multiply(self, scalar, value):
        if type(scalar) is not int:
            raise TypeError('scalar must be an exact integer')
        value = self.point(value)
        result = (0, 0)
        for _ in range(scalar % 4):
            result = self.add(result, value)
        return result

    def order(self, value):
        value = self.point(value)
        if value == (0, 0):
            return 1
        return 4 if self.doubling_column(value[1]) else 2

    @cached_property
    def invariant_generators(self):
        # Choose lifts with independent doubles, then a basis of ker(doubling).
        pivots, images = [], []
        for j, d in enumerate(self.doubles):
            if len(row_basis((*images, d))) > len(images):
                pivots.append(j)
                images.append(d)
        image = list(images)
        complement = []
        for j in range(self.kernel_dimension):
            v = 1 << j
            if len(row_basis((*image, v))) > len(image):
                image.append(v)
                complement.append(v)
        nullspace = affine_f2(self.doubles, 0).basis
        order_two = tuple((v, 0) for v in complement) + tuple((0, b) for b in nullspace)
        order_four = tuple((0, 1 << j) for j in pivots)
        result = order_two + order_four
        if any(self.order(g) != n for g, n in zip(result, self.invariant_factors)):
            raise ArithmeticError('invariant generators have inconsistent orders')
        return result

    @property
    def invariant_factors(self):
        rank = len(row_basis(self.doubles))
        return (2,) * (self.kernel_dimension + self.image_dimension - 2*rank) + (4,) * rank

    def from_invariants(self, coordinates):
        values = tuple(nonnegative_integer(v, 'invariant coordinate') for v in coordinates)
        if len(values) != len(self.invariant_factors) or any(
                v >= n for v, n in zip(values, self.invariant_factors)):
            raise ValueError('invariant-factor coordinates are not canonical')
        result = (0, 0)
        for c, g in zip(values, self.invariant_generators):
            result = self.add(result, self.multiply(c, g))
        return result

    def to_invariants(self, value):
        a, b = self.point(value)
        generators = self.invariant_generators
        quotient_indices = tuple(i for i, g in enumerate(generators) if g[1])
        quotient_basis = F2Labels(self.image_dimension, tuple(generators[i][1] for i in quotient_indices))
        low_bits = quotient_basis.coordinates(b)
        coordinates = [0]*len(generators)
        used = (0, 0)
        for j, i in enumerate(quotient_indices):
            coordinates[i] = (low_bits >> j) & 1
            if coordinates[i]:
                used = self.add(used, generators[i])
        residue = self.add((a, b), self.negate(used))
        if residue[1]:
            raise ArithmeticError('quotient elimination failed')
        kernel_indices = tuple(i for i, (g, n) in enumerate(zip(generators, self.invariant_factors))
                               if not g[1] or n == 4)
        kernel_vectors = tuple(generators[i][0] if self.invariant_factors[i] == 2
                               else self.doubling_column(generators[i][1]) for i in kernel_indices)
        high_bits = F2Labels(self.kernel_dimension, kernel_vectors).coordinates(residue[0])
        for j, i in enumerate(kernel_indices):
            coordinates[i] += ((high_bits >> j) & 1)*(2 if self.invariant_factors[i] == 4 else 1)
        result = tuple(coordinates)
        if self.from_invariants(result) != (a, b):
            raise ArithmeticError('invariant coordinate round trip failed')
        return result
