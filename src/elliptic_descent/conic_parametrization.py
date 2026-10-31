"""Exact conic/Veronese parametrizations with algebraically checked inverses.

A rational point on a nonsingular conic gives three binary quadratic forms.
Their coefficient matrix is invertible, so the parametrization is an
isomorphism of the conic with P^1, not an unverified dominant map.
"""
from dataclasses import dataclass
from math import gcd, lcm

import sympy as sp

from .geometry import U, _rat
from .two_cover import X, Z


def _matrix(rows, size):
    M = sp.Matrix(rows)
    if M.shape != (size, size):
        raise ValueError(f'expected a {size} by {size} rational matrix')
    return sp.ImmutableMatrix(size, size, [_rat(c) for c in M])


def _primitive_matrix(M):
    """Common rational scaling changes a conic presentation, not its image."""
    denominator = lcm(*(int(c.q) for c in M))
    integers = [int(c*denominator) for c in M]
    content = gcd(*integers)
    if not content:
        raise ValueError('zero parametrization')
    if next(c for c in integers if c) < 0:
        content = -content
    return tuple(tuple(sp.Rational(c*denominator, content) for c in row)
                 for row in M.tolist())


def _from_forms(forms):
    M = sp.Matrix([[sp.Poly(p, X, Z, domain=sp.QQ).coeff_monomial(m)
                    for m in (X*X, X*Z, Z*Z)] for p in forms])
    return ConicParametrization(_primitive_matrix(M))


@dataclass(frozen=True)
class ConicParametrization:
    """u=M*(X^2,X*Z,Z^2), with M invertible over Q.

    The inverse is obtained from M^-1*u=(a,b,c), ac=b^2. The two charts
    (a:b) and (b:c) cover the conic, including the chosen tangent point.
    """
    matrix: tuple

    def __post_init__(self):
        M = _matrix(self.matrix, 3)
        if not M.det():
            raise ValueError('the Veronese coefficient matrix must be invertible')
        object.__setattr__(self, 'matrix', tuple(map(tuple, M.tolist())))

    @property
    def forms(self):
        return tuple(sp.expand(p) for p in sp.Matrix(self.matrix)*sp.Matrix((X*X, X*Z, Z*Z)))

    def verify(self, quadratic):
        q = sp.Poly(quadratic, *U[:3], domain=sp.QQ)
        if any(sum(m) != 2 for m, _ in q.terms()):
            return False
        H = sp.hessian(q.as_expr(), U[:3])
        if not H.det():
            return False
        substitution = dict(zip(U[:3], self.forms))
        return sp.expand(q.as_expr().subs(substitution, simultaneous=True)) == 0

    def reparameterize(self, change):
        """Precompose an exact GL2 change; do not silently allow a degree drop."""
        N = _matrix(change, 2)
        if not N.det():
            raise ValueError('the parameter change must be invertible')
        a, b = N*sp.Matrix((X, Z))
        return _from_forms(tuple(p.subs({X:a, Z:b}, simultaneous=True) for p in self.forms))

    def inverse_coordinates(self, values):
        """Return (X,Z,scale) with p(X,Z)=scale*values; no square roots."""
        if len(values) != 3:
            raise ValueError('expected three conic coordinates')
        values = sp.Matrix([_rat(c) for c in values])
        if not any(values):
            raise ValueError('zero is not a projective conic point')
        a, b, c = sp.Matrix(self.matrix).inv()*values
        if a*c != b*b:
            raise ValueError('point is not on the parametrized conic')
        return (a, b, a) if a else (b, c, c)


def parametrize_conic(quadratic, point):
    """Construct and verify the line-through-point parametrization exactly.

    With q(P)=0, d=X*V+Z*W, use q(d)*P-(P^T*H*d)*d. The complements V,W
    are coordinate vectors completing P to a basis. A determinant check
    excludes a basepoint or a multiple parametrization.
    """
    if len(point) != 3:
        raise ValueError('expected three rational conic-point coordinates')
    P = sp.Matrix([_rat(c) for c in point])
    q = sp.Poly(quadratic, *U[:3], domain=sp.QQ)
    if not any(P) or any(sum(m) != 2 for m, _ in q.terms()):
        raise ValueError('expected a nonzero point and a homogeneous conic')
    H = sp.hessian(q.as_expr(), U[:3])
    if not H.det() or (P.T*H*P)[0] != 0:
        raise ValueError('point is not on a nonsingular conic')
    pivot = next(i for i, c in enumerate(P) if c)
    i, j = [i for i in range(3) if i != pivot]
    d = X*sp.eye(3)[:, i]+Z*sp.eye(3)[:, j]
    forms = (d.T*H*d)[0]*P/2-(P.T*H*d)[0]*d
    result = _from_forms(forms)
    if not result.verify(quadratic):
        raise ArithmeticError('conic parametrization failed its exact identity')
    return result
