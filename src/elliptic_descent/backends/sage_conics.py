"""A rational conic point from PARI qfsolve, checked in exact rational arithmetic.

The conic underlying an everywhere locally soluble 2-cover is locally soluble
and hence rationally soluble by Hasse--Minkowski. A failed backend response
raises; it is never interpreted as an empty Selmer set or an omitted class.
"""
import sympy as sp

from .sage import _sage
from ..geometry import U, _rat


def rational_conic_point(quadratic):
    H = sp.hessian(quadratic, U[:3])
    if not H.det():
        raise ValueError('a nonsingular ternary quadratic form is required')
    # Cheap exact witnesses are sufficient, but failure here is not an obstruction.
    for i in range(3):
        P = tuple(sp.Integer(i == j) for j in range(3))
        if (sp.Matrix(P).T*H*sp.Matrix(P))[0] == 0:
            return P
    sage = _sage()
    from sage.quadratic_forms.qfsolve import qfsolve
    matrix = sage.matrix(sage.QQ, [[sage.QQ(str(c)) for c in row] for row in H.tolist()])
    solution = qfsolve(matrix)
    try:
        P = tuple(_rat(c) for c in solution)
    except TypeError as error:
        raise ArithmeticError('qfsolve did not return a rational conic point') from error
    if len(P) != 3 or not any(P) or (sp.Matrix(P).T*H*sp.Matrix(P))[0] != 0:
        raise ArithmeticError('qfsolve conic point failed exact verification')
    return P
