"""Exact rational origin torsors of marked 4-coverings.

The fibre above O is a degree-sixteen etale E[4]-torsor. Work in its rational
coordinate algebra, rather than a common algebraic-closure splitting field.
All discovered points and maps are verified as exact polynomial identities.
This module is a reference backend, not a practical complexity guarantee.
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations, count

import sympy as sp

from .sage import _sage
from ..geometry import U


def _q(c):
    return _sage().QQ(str(c))


def _eval(expr, values):
    """Evaluate a rational polynomial without string parsing or lost exactness."""
    result = 0
    for powers, coefficient in sp.Poly(expr, *U, domain=sp.QQ).terms():
        term = _q(coefficient)
        for v, n in zip(values, powers):
            term *= v**n
        result += term
    return result


def _coefficients(poly, n):
    return tuple(poly[i] for i in range(n))


class EtaleAlgebra:
    """A monogenic finite etale Q-algebra, with exact power-basis coordinates."""
    def __init__(self, coefficients):
        s = _sage()
        self.polynomial_ring = s.PolynomialRing(s.QQ, 't')
        self.polynomial = self.polynomial_ring([_q(c) for c in coefficients]).monic()
        self.degree = int(self.polynomial.degree())
        if self.degree < 1 or not self.polynomial.is_squarefree():
            raise ValueError('the algebra polynomial must be nonconstant and separable')
        self.ring = self.polynomial_ring.quotient(self.polynomial, 'z')
        self.generator = self.ring.gen()

    def vector(self, value):
        return _sage().vector(_sage().QQ, _coefficients(self.ring(value).lift(), self.degree))

    def element(self, coordinates):
        return self.ring(self.polynomial_ring(list(coordinates)))

    @cached_property
    def factors(self):
        return tuple(g.monic() for g, _ in self.polynomial.factor())

    def component(self, factor):
        s = _sage()
        if factor.degree() == 1:
            return s.QQ, -factor[0]
        K = s.NumberField(factor, 'a')
        return K, K.gen()

    def evaluate(self, value, root):
        return self.ring(value).lift()(root)

    def reconstruct(self, values):
        """CRT of factor polynomials; no component ordering or field merging."""
        if len(values) != len(self.factors):
            raise ValueError('component count does not match factorization')
        result = self.polynomial_ring.zero()
        for g, value in zip(self.factors, values):
            h = self.polynomial // g
            e = h*h.inverse_mod(g)
            result += self.polynomial_ring(value)*e
        return self.ring(result)

    def multiplication_matrix(self, value):
        s = _sage()
        return s.matrix(s.QQ, [self.vector(value*self.generator**j)
                              for j in range(self.degree)]).transpose()

    @cached_property
    def trace_pairing(self):
        s = _sage()
        traces = [self.multiplication_matrix(self.generator**i).trace()
                  for i in range(2*self.degree-1)]
        return s.matrix(s.QQ, self.degree, self.degree, lambda i, j: traces[i+j])


def _origin_algebra(covering):
    """Find an affine chart and primitive element of the full geometric O-fibre."""
    s = _sage()
    if covering.marking.translation is not None or not covering.verify():
        raise ValueError('an exactly verified zero-translation marking is required')
    D = covering.model
    R = s.PolynomialRing(s.QQ, names=('v1', 'v2', 'v3'), order='degrevlex')
    v = R.gens()
    for k in count():
        ell = tuple(s.QQ(a) for a in (1, k, k*k, k*k*k))
        values = (1-sum(a*b for a, b in zip(ell[1:], v)), *v)
        equations = tuple(R(_eval(q, values)) for q in (*D.quadrics, D.norm_u))
        I = R.ideal(equations)
        if I.dimension() != 0 or I.vector_space_dimension() != 16:
            continue
        monomials = tuple(I.normal_basis())
        def vector(p):
            p = I.reduce(p)
            return s.vector(s.QQ, [p.monomial_coefficient(m) for m in monomials])
        for j in count(1):
            primitive = v[0]+j*v[1]+j*j*v[2]
            powers, term = [], R.one()
            for _ in range(17):
                powers.append(vector(term))
                term = I.reduce(term*primitive)
            B = s.matrix(s.QQ, powers[:16]).transpose()
            if not B.det():
                continue
            relation = B.solve_right(powers[16])
            algebra = EtaleAlgebra(tuple(-a for a in relation)+(s.QQ(1),))
            coordinates = tuple(algebra.element(B.solve_right(vector(p))) for p in values)
            return OriginGeometry(covering, algebra, ell, (0, 1, j, j*j), coordinates)


def _jets(matrices, point, K):
    s = _sage()
    pivot = next(i for i, c in enumerate(point) if c)
    Q = point/point[pivot]
    J = s.matrix(K, [list(H*Q) for H in matrices])
    others = tuple(i for i in range(4) if i != pivot)
    dependent = next(pair for pair in combinations(others, 2)
                     if J.matrix_from_columns(pair).det())
    free = next(i for i in others if i not in dependent)
    inverse = J.matrix_from_columns(dependent).inverse()
    linear = s.vector(K, 4); linear[free] = 1
    for i, a in zip(dependent, -inverse*J.column(free)):
        linear[i] = a
    jets = [Q, linear]
    for n in range(2, 5):
        rhs = s.vector(K, [sum((jets[i]*H*jets[n-i] for i in range(1, n)), K(0))/2
                           for H in matrices])
        v = s.vector(K, 4)
        for i, a in zip(dependent, -inverse*rhs):
            v[i] = a
        jets.append(v)
    return tuple(jets)


def _solve(columns, rhs, R, degree, K):
    s = _sage()
    mons = [R.monomial(i, j, k, degree-i-j-k) for i in range(degree+1)
            for j in range(degree-i+1) for k in range(degree-i-j+1)]
    M = s.matrix(K, [[p.monomial_coefficient(m) for p in columns] for m in mons])
    if M.rank() != len(columns):
        raise ArithmeticError('ambiguous flex polynomial identity')
    answer = M.solve_right(s.vector(K, [rhs.monomial_coefficient(m) for m in mons]))
    if sum((a*p for a, p in zip(answer, columns)), R.zero()) != rhs:
        raise ArithmeticError('flex identity is inconsistent')
    return answer


def _embedded(ainvs, z):
    a1, a2, a3, a4, a6 = ainvs
    z0, z1, z2, z3 = z
    return (z0*z3-z1*z1, z2*z2+a1*z1*z2+a3*z0*z2-z1*z3
            -a2*z0*z3-a4*z0*z1-a6*z0*z0)


def _flex_data(covering, point, K):
    """Riemann--Roch at one finite-algebra component; no splitting field needed."""
    s = _sage(); D = covering.model
    H = tuple(s.matrix(K, [[_q(c) for c in row] for row in sp.hessian(q, U).tolist()])
              for q in D.quadrics)
    jets = _jets(H, s.vector(K, point), K)
    Q, v1, v2 = jets[:3]
    q0, q1 = (_eval(q, Q) for q in D.q[:2])
    if not q1:
        raise ArithmeticError('a degree-four branch point cannot be at infinity')
    alpha = -q0/q1
    if sum(K(_q(c))*alpha**i for i, c in enumerate(D.f)):
        raise ArithmeticError('origin point did not map to a branch root')
    H01 = [s.matrix(K, [[_q(c) for c in row] for row in sp.hessian(q, U).tolist()])
           for q in D.q[:2]]
    Hbranch = -H01[0]-alpha*H01[1]
    w2 = Q*Hbranch*v2 + (v1*Hbranch*v1)/2
    # Differentiate the rational norm form directly, avoiding its four roots.
    y1 = K(_q(D.f[4]*D.norm_witness))*sum(
        (_eval(sp.diff(D.norm_u, U[i]), Q)*v1[i] for i in range(4)), K(0))
    if not w2 or not y1:
        raise ArithmeticError('origin point does not have the expected simple differential')
    desired = K(_q(covering.marking.change.urst[0]))*q1*w2/(2*y1)
    return jets, desired


def field_flex_matrix(covering, point, K):
    """Return B with [4]B=pi and B(point)=O, using exact field arithmetic."""
    s = _sage()
    jets, desired = _flex_data(covering, point, K)
    kernel = lambda vs: s.matrix(K, [list(v) for v in vs]).right_kernel().basis()
    base = kernel(jets[:4])
    if len(base) != 1 or not base[0]*jets[4]:
        raise ArithmeticError('origin is not a degree-four flex')
    L0 = s.vector(K, list(base[0]))
    Lx = s.vector(K, list(next(v for v in kernel(jets[:2]) if v*jets[2])))
    Ly = s.vector(K, list(next(v for v in kernel(jets[:1]) if v*jets[1])))
    R = s.PolynomialRing(K, names=('u0','u1','u2','u3'))
    z = R.gens()
    linear = lambda v: sum((R(v[i])*z[i] for i in range(4)), R.zero())
    f0, fx, fy = map(linear, (L0, Lx, Ly))
    quadrics = tuple(R(_eval(q, z)) for q in covering.model.quadrics)
    terms = (fx**3, fx*fy*f0, fx*fx*f0, fy*f0*f0, fx*f0*f0, f0**3)
    ideal = tuple(q*v for q in quadrics for v in z)
    b, B, C, D, E, F = _solve((*terms, *ideal), fy*fy*f0, R, 3, K)[:6]
    if not b:
        raise ArithmeticError('singular flex cubic')
    ainvs = (-B, C, -b*D, b*E, b*b*F)
    last = _solve(tuple(f0*v for v in z)+quadrics, b*b*fx*fx, R, 2, K)
    B0 = s.matrix(K, [list(L0), list(b*Lx), list(b*Ly), list(last[:4])])
    actual = -(Lx*jets[2])/(Ly*jets[1])
    u = desired/actual
    target = tuple(K(_q(a)) for a in covering.target.ainvs)
    a1, a2, a3, _, _ = ainvs
    ss = (u*target[0]-a1)/2
    r = (u*u*target[1]-a2+a1*ss+ss*ss)/3
    t = (u**3*target[2]-a3-a1*r)/2
    change = s.matrix(K, [[1,0,0,0],[-r/u**2,1/u**2,0,0],
                          [(ss*r-t)/u**3,-ss/u**3,1/u**3,0],
                          [r*r/u**4,-2*r/u**4,0,1/u**4]])
    B = change*B0
    B /= next(v for v in B.list() if v)
    if not verify_field_flex(covering, point, K, B):
        raise ArithmeticError('finite-field-component flex map failed verification')
    return B


def verify_field_flex(covering, point, K, B):
    s = _sage()
    if B.nrows()!=4 or B.ncols()!=4 or not B.det():
        return False
    R = s.PolynomialRing(K, names=('u0','u1','u2','u3'))
    z = s.vector(R, R.gens())
    quadrics = tuple(R(_eval(q, z)) for q in covering.model.quadrics)
    Q = s.vector(K, point)
    if not any(Q) or any(_eval(q, Q) for q in covering.model.quadrics):
        return False
    image = B*Q
    if any(image[i] for i in range(3)) or not image[3]:
        return False
    try:
        for q in _embedded(tuple(K(_q(a)) for a in covering.target.ainvs), B*z):
            _solve(quadrics, q, R, 2, K)
        jets, desired = _flex_data(covering, Q, K)
        if any(B.row(0)*v for v in jets[:4]) or not B.row(0)*jets[4]:
            return False
        if any(B.row(1)*v for v in jets[:2]) or not B.row(1)*jets[2]:
            return False
        if B.row(2)*jets[0] or not B.row(2)*jets[1]:
            return False
        return -(B.row(1)*jets[2])/(B.row(2)*jets[1]) == desired
    except (ValueError, ArithmeticError, ZeroDivisionError):
        return False


@dataclass
class OriginGeometry:
    covering: object
    algebra: EtaleAlgebra
    affine_form: tuple
    separating_form: tuple
    point: tuple

    def verify(self):
        A, P = self.algebra, self.point
        if A.degree != 16 or len(P)!=4:
            return False
        if sum(a*p for a,p in zip(self.affine_form,P)) != 1:
            return False
        if sum(a*p for a,p in zip(self.separating_form,P)) != A.generator:
            return False
        return not any(_eval(q,P) for q in (*self.covering.model.quadrics,self.covering.model.norm_u))

    @cached_property
    def flexes(self):
        result = []
        for factor in self.algebra.factors:
            K, root = self.algebra.component(factor)
            point = tuple(self.algebra.evaluate(p, root) for p in self.point)
            B = field_flex_matrix(self.covering, point, K)
            result.append((factor, K, root, point, B))
        return tuple(result)

    @cached_property
    def inverse_flex_matrix(self):
        s = _sage(); A = self.algebra
        inverses = [B.inverse() for _,_,_,_,B in self.flexes]
        def coefficients(value, factor):
            return [value] if factor.degree()==1 else value.list()
        return s.matrix(A.ring, 4, 4, lambda i,j: A.reconstruct([
            coefficients(B[i,j], factor) for B, factor in zip(inverses,A.factors)]))


def origin_geometry(covering):
    result = _origin_algebra(covering)
    if not result.verify():
        raise ArithmeticError('rational origin fibre failed its exact equations')
    return result
