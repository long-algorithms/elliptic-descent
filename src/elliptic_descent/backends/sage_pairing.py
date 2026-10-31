"""Square-root discovery for Fisher's exact cubic product invariant.

Reuse the parent's labelled factor fields. Never extend the fields to force a
square root; the returned root is independently squared in the rational algebra.
The product root is construction evidence, not a new class-group computation.
"""
from ..geometry import _sqrt_rational


def cubic_product_square_root(arithmetic, value):
    A = arithmetic.algebra
    if not A.norm(value):
        raise ValueError('the product invariant must be a unit')
    components = []
    for factor, part in zip(arithmetic.factors, A.to_components(value)):
        if factor.degree == 1:
            coordinates = (_sqrt_rational(part[0]),)
        else:
            element = factor.as_field_element(part)
            root = element.sqrt(extend=False)
            if root*root != element:
                raise ArithmeticError('number-field square root failed exact verification')
            coordinates = factor.as_coefficients(root)
        # Fix the otherwise irrelevant independent signs reproducibly.
        coordinates = min(coordinates, tuple(-c for c in coordinates))
        components.append(coordinates)
    result = A.from_components(components)
    if A.multiply(result, result) != tuple(value):
        raise ArithmeticError('CRT product square root failed rational-algebra verification')
    return result
