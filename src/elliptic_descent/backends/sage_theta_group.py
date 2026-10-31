"""Exact theta-cochain comparisons of marked 4-coverings, over Sage QQbar.

A geometric trivialization B:C->E with [4]B=pi determines a theta cochain.
Products encode cohomological addition. Rational interpolation on E[4] tests
coboundaries up to ALL sixteen Weil characters. No group is guessed from its
order or pairing rank; no numerical tolerance is used for a decision.

This reference backend matches an existing exhaustive set of coverings rather
than splitting a new central simple algebra to construct another model.
See docs/four-selmer-group.md and Fisher (2008), Section 8.
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations, count

import sympy as sp

from .sage import _sage
from ..geometry import U


def _field():
    return _sage().QQbar


def _rational(value):
    return _sage().QQ(str(value))


def _coefficients(expression, degree):
    p = sp.Poly(expression, *U, domain=sp.QQ)
    return {tuple(m): _rational(c) for m, c in p.terms()}


def _polynomial(expression, ring):
    return sum((ring.base_ring()(_rational(c))*ring.monomial(*m)
                for m, c in sp.Poly(expression, *U, domain=sp.QQ).terms()), ring.zero())


def _monomials(ring, degree):
    return tuple(ring.monomial(i, j, k, degree-i-j-k)
                 for i in range(degree+1) for j in range(degree-i+1)
                 for k in range(degree-i-j+1))


def _solve_polynomials(columns, rhs, degree):
    sage = _sage()
    monomials = _monomials(rhs.parent(), degree)
    M = sage.matrix(_field(), [[q.monomial_coefficient(m) for q in columns] for m in monomials])
    v = sage.vector(_field(), [rhs.monomial_coefficient(m) for m in monomials])
    if M.rank() != len(columns):
        raise ArithmeticError('polynomial identity has an ambiguous coefficient solution')
    answer = M.solve_right(v)
    if sum((c*q for c, q in zip(answer, columns)), rhs.parent().zero()) != rhs:
        raise ArithmeticError('polynomial coefficient identity failed')
    return answer


def _normalized_matrix(M):
    a = next((v for v in M.list() if v != 0), None)
    if a is None:
        raise ValueError('zero projective matrix')
    return M/a


def _projective_equal(A, B):
    return _normalized_matrix(A) == _normalized_matrix(B)


def _sort_algebraics(values):
    return tuple(sorted(values, key=lambda a: (a.real(), a.imag())))


def _hessians(model):
    sage = _sage()
    K = _field()
    return tuple(sage.matrix(K, [[_rational(c) for c in row] for row in sp.hessian(q, U).tolist()])
                 for q in model.quadrics)


def _jets(matrices, point):
    """Local affine Hensel expansion at a smooth point, through order four."""
    sage, K = _sage(), _field()
    pivot = next(i for i, c in enumerate(point) if c)
    Q = point/point[pivot]
    J = sage.matrix(K, [list(H*Q) for H in matrices])
    others = tuple(i for i in range(4) if i != pivot)
    dependent = next((pair for pair in combinations(others, 2)
                      if J.matrix_from_columns(pair).det()), None)
    if dependent is None:
        raise ArithmeticError('the chosen point is singular')
    free = next(i for i in others if i not in dependent)
    inverse = J.matrix_from_columns(dependent).inverse()
    linear = sage.vector(K, 4)
    linear[free] = 1
    solved = -inverse*J.column(free)
    for i, a in zip(dependent, solved):
        linear[i] = a
    jets = [Q, linear]
    for n in range(2, 5):
        rhs = sage.vector(K, [sum((jets[i]*H*jets[n-i] for i in range(1, n)), K.zero())/2
                               for H in matrices])
        v = sage.vector(K, 4)
        for i, a in zip(dependent, -inverse*rhs):
            v[i] = a
        jets.append(v)
    return tuple(jets)


def _root_point(model):
    """An exact geometric point over a branch point, hence a 4-flex."""
    sage, K = _sage(), _field()
    R = sage.PolynomialRing(K, 'a')
    f = R([_rational(c) for c in model.f])
    roots = _sort_algebraics(f.roots(multiplicities=False))
    if len(roots) != 4:
        raise ArithmeticError('the quartic is not separable')
    # Rational branch points keep the common cases in much smaller fields.
    from ..etale import EtaleQuarticAlgebra
    linear = [g for g in EtaleQuarticAlgebra(model.f).factors if len(g) == 2]
    alpha = K(_rational(-linear[0][0])) if linear else roots[0]
    xi = R([_rational(c) for c in model.xi])
    values = sage.vector(K, [((a-alpha)/xi(a)).sqrt() for a in roots])
    evaluation = sage.matrix(K, [[a**i for i in range(4)] for a in roots])
    Q = evaluation.solve_right(values)
    if not any(Q):
        raise ArithmeticError('zero cannot be a flex point')
    H = _hessians(model)
    if any(Q*M*Q for M in H):
        raise ArithmeticError('constructed flex point failed the source equations')
    return alpha, roots, Q


def _intermediate_differential(covering, alpha, roots, jets):
    """Linear term of -x/y for I_alpha o (D->C) on the original E.

    I_alpha sends the branch point to O and pi_C=[2]I_alpha. Therefore the
    marked trivialization D->E must have HALF this invariant differential.
    """
    sage, K = _sage(), _field()
    D = covering.model
    if covering.marking.translation is not None:
        raise ValueError('theta trivialization requires the retained zero-translation marking')
    Q, v1, v2 = jets[:3]
    H0, H1 = tuple(sage.matrix(K, [[_rational(c) for c in row]
                       for row in sp.hessian(q, U).tolist()]) for q in D.q[:2])
    H = -H0-alpha*H1
    w2 = Q*H*v2 + (v1*H*v1)/2
    q1 = (Q*H1*Q)/2
    evaluation = lambda v, a: sum((c*a**i for i, c in enumerate(v)), K.zero())
    norm_linear = evaluation(v1, alpha)
    for a in roots:
        if a != alpha:
            norm_linear *= evaluation(Q, a)
    y1 = _rational(D.f[4])*_rational(D.norm_witness)*norm_linear
    if not q1 or not w2 or not y1:
        raise ArithmeticError('exceptional flex differential is not simple')
    return _rational(covering.marking.change.urst[0])*q1*w2/y1


def _changed_ainvs(ainvs, urst):
    a1, a2, a3, a4, a6 = ainvs
    u, r, s, t = urst
    return ((a1+2*s)/u, (a2-a1*s+3*r-s*s)/u**2, (a3+a1*r+2*t)/u**3,
            (a4+2*a2*r-a1*(r*s+t)-a3*s+3*r*r-2*s*t)/u**4,
            (a6-a1*r*t+a2*r*r-a3*t+a4*r+r**3-t*t)/u**6)


def _embedded_equations(ainvs, z):
    a1, a2, a3, a4, a6 = ainvs
    z0, z1, z2, z3 = z
    return (z0*z3-z1*z1,
            z2*z2+a1*z1*z2+a3*z0*z2-z1*z3-a2*z0*z3-a4*z0*z1-a6*z0*z0)


@dataclass(frozen=True)
class FlexTrivialization:
    """A checked projective C->E isomorphism whose marked map is [4]."""
    covering: object
    alpha: object
    point: tuple
    matrix: object

    def verify(self):
        sage, K = _sage(), _field()
        D, B = self.covering.model, self.matrix
        if B.nrows() != 4 or B.ncols() != 4 or not B.det():
            return False
        R = sage.PolynomialRing(K, names=('u0', 'u1', 'u2', 'u3'))
        quadrics = tuple(_polynomial(q, R) for q in D.quadrics)
        Q = sage.vector(K, self.point)
        if not any(Q) or any(q(*Q) for q in quadrics):
            return False
        image = B*Q
        if any(image[i] for i in range(3)) or not image[3]:
            return False
        try:
            transformed = B*sage.vector(R, R.gens())
            E = tuple(K(_rational(c)) for c in self.covering.target.ainvs)
            for q in _embedded_equations(E, transformed):
                _solve_polynomials(quadrics, q, 2)
            jets = _jets(_hessians(D), Q)
            # Orders (4,2,1,0) at the flex fix the differential and signed marking.
            if any(B.row(0)*v for v in jets[:4]) or not B.row(0)*jets[4]:
                return False
            if any(B.row(1)*v for v in jets[:2]) or not B.row(1)*jets[2]:
                return False
            if B.row(2)*jets[0] or not B.row(2)*jets[1]:
                return False
            P = sage.PolynomialRing(K, 't')([_rational(c) for c in D.f])
            roots = _sort_algebraics(P.roots(multiplicities=False))
            if P(self.alpha) or sum((Q[i]*self.alpha**i for i in range(4)), K.zero()):
                return False
            expected = _intermediate_differential(self.covering, self.alpha, roots, jets)
            actual = -(B.row(1)*jets[2])/(B.row(2)*jets[1])
            return 2*actual == expected
        except (ValueError, ArithmeticError, ZeroDivisionError):
            return False


def flex_trivialization(covering):
    """Build B by Riemann--Roch pole orders at one exact 4-flex.

    This uses linear algebra and polynomial identities, not numerical curve
    recognition, and works equally for special j-invariants and rational torsion.
    """
    sage, K = _sage(), _field()
    alpha, roots, point = _root_point(covering.model)
    jets = _jets(_hessians(covering.model), point)
    Q = jets[0]
    kernel = lambda vs: sage.matrix(K, [list(v) for v in vs]).right_kernel().basis()
    base = kernel(jets[:4])
    if len(base) != 1 or not base[0]*jets[4]:
        raise ArithmeticError('point is not an ordinary degree-four flex')
    L0 = base[0]
    Lx = next(v for v in kernel(jets[:2]) if v*jets[2])
    Ly = next(v for v in kernel(jets[:1]) if v*jets[1])
    R = sage.PolynomialRing(K, names=('u0', 'u1', 'u2', 'u3'))
    z = sage.vector(R, R.gens())
    s, x, y = (v*z for v in (L0, Lx, Ly))
    quadrics = tuple(_polynomial(q, R) for q in covering.model.quadrics)
    terms = (x**3, x*y*s, x*x*s, y*s*s, x*s*s, s**3)
    ideal = tuple(q*v for q in quadrics for v in R.gens())
    solution = _solve_polynomials((*terms, *ideal), y*y*s, 3)
    b, B, C, D, E, F = solution[:6]
    if not b:
        raise ArithmeticError('degenerate cubic leading coefficient')
    ainvs = (-B, C, -b*D, b*E, b*b*F)
    last = _solve_polynomials(tuple(s*v for v in R.gens())+quadrics, b*b*x*x, 2)
    B0 = sage.matrix(K, [list(L0), list(b*Lx), list(b*Ly), list(last[:4])])
    actual = -(Lx*jets[2])/(Ly*jets[1])
    desired = _intermediate_differential(covering, alpha, roots, jets)/2
    u = desired/actual
    target = tuple(K(_rational(a)) for a in covering.target.ainvs)
    a1, a2, a3, _, _ = ainvs
    s0 = (u*target[0]-a1)/2
    r = (u*u*target[1]-a2+a1*s0+s0*s0)/3
    t = (u**3*target[2]-a3-a1*r)/2
    if _changed_ainvs(ainvs, (u, r, s0, t)) != target:
        raise ArithmeticError('signed differential does not give the original elliptic model')
    change = sage.matrix(K, [[1, 0, 0, 0], [-r/u**2, 1/u**2, 0, 0],
                [(s0*r-t)/u**3, -s0/u**3, 1/u**3, 0],
                [r*r/u**4, -2*r/u**4, 0, 1/u**4]])
    B = _normalized_matrix(change*B0)
    result = FlexTrivialization(covering, alpha, tuple(Q), B)
    if not result.verify():
        raise ArithmeticError('projective flex trivialization failed exact verification')
    return result


def _torsion_points(E):
    """Enumerate E[4] by exact halves of E[2], including special j-invariants."""
    sage, K = _sage(), _field()
    R = sage.PolynomialRing(K, 'x')
    x = R.gen()
    a1, a2, a3, a4, a6 = E.a_invariants()
    h = x**3+(a1*a1/4+a2)*x*x+(a1*a3/2+a4)*x+(a3*a3/4+a6)
    points = [E(0)]
    for root in _sort_algebraics(h.roots(multiplicities=False)):
        T = E(root, -(a1*root+a3)/2)
        points.append(T)
        polynomial = h.derivative()**2-4*h*(2*x+h[2]+root)
        for z in _sort_algebraics(polynomial.roots(multiplicities=False)):
            v = h(z).sqrt()
            for sign in (-1, 1):
                P = E(z, sign*v-(a1*z+a3)/2)
                if 2*P != T:
                    raise ArithmeticError('algebraic halving failed the group-law check')
                if P not in points:
                    points.append(P)
    if len(points) != 16 or any(4*P != E(0) for P in points):
        raise ArithmeticError('4-torsion enumeration is not complete')
    key = lambda P: (P[0].real(), P[0].imag(), P[1].real(), P[1].imag())
    return (E(0),)+tuple(sorted(points[1:], key=key))


def _translation_matrix(E, point):
    """Exact linear action on |4O|, proved from the chord formula in Qbar(E)."""
    sage, K = _sage(), _field()
    if point.is_zero():
        return sage.identity_matrix(K, 4)
    R = sage.PolynomialRing(K, 'x')
    x = R.gen()
    S = sage.PolynomialRing(R, 'y')
    y = S.gen()
    a1, a2, a3, a4, a6 = E.a_invariants()
    modulus = y*y+(a1*x+a3)*y-x**3-a2*x*x-a4*x-a6
    reduce = lambda f: S(f).quo_rem(modulus)[1]
    tx, ty = point.xy()
    d, m = x-tx, y-ty
    nx = reduce(m*m+a1*m*d-(a2+x+tx)*d*d)
    ny = reduce(-(y+a3)*d**3-m*(nx-x*d*d)-a1*nx*d)
    base = (S(1), S(x), y, S(x*x))
    # Three function identities simultaneously determine the 16 matrix entries.
    relations = []
    for row, numerator, denominator in ((1, nx, d*d), (2, ny, d**3), (3, reduce(nx*nx), d**4)):
        columns = [S(0)]*16
        for j, z in enumerate(base):
            columns[j] = -reduce(z*numerator)
            columns[4*row+j] = reduce(z*denominator)
        for power_y in (0, 1):
            deg = max([int(c[power_y].degree()) for c in columns]+[0])
            for power_x in range(deg+1):
                relations.append([c[power_y][power_x] for c in columns])
    kernel = sage.matrix(K, relations).right_kernel().basis()
    if len(kernel) != 1:
        raise ArithmeticError('translation does not have a unique projective linear action')
    M = _normalized_matrix(sage.matrix(K, 4, 4, list(kernel[0])))
    if not M.det():
        raise ArithmeticError('singular torsion translation matrix')
    forms = M*sage.vector(S, base)
    for row, numerator, denominator in ((1, nx, d*d), (2, ny, d**3), (3, reduce(nx*nx), d**4)):
        if reduce(forms[row]*denominator-forms[0]*numerator):
            raise ArithmeticError('torsion action failed the exact chord identity')
    return M


@dataclass(frozen=True)
class ThetaCoboundary:
    """q(T)*character_S(T)=r(T), with r in the rational torsion etale algebra."""
    character_index: int
    origin_value: object
    polynomial: tuple

    def verify(self, context, values):
        sage = _sage()
        if (type(self.character_index) is not int or not 0 <= self.character_index < 16
                or len(values) != 16 or len(self.polynomial) != 15):
            return False
        try:
            origin = sage.QQ(self.origin_value)
            polynomial = tuple(sage.QQ(c) for c in self.polynomial)
        except (TypeError, ValueError):
            return False
        if not origin or values[0]*context.characters[self.character_index][0] != origin:
            return False
        for z, value, char in zip(context.nodes, values[1:], context.characters[self.character_index][1:]):
            evaluation = sum((c*z**i for i, c in enumerate(polynomial)), _field().zero())
            if not evaluation or value*char != evaluation:
                return False
        return True


class ThetaContext:
    """The same rational, Galois-equivariant E[4] representation for every cover.

    Two exact chord-formula matrices generate all sixteen translations. Their
    fourth powers and commutator certify every projective multiplication law.
    A rational resultant constructs the interpolation polynomial directly,
    avoiding gratuitous arithmetic in the full torsion splitting field.
    """
    def __init__(self, curve):
        sage, K = _sage(), _field()
        self.curve = curve
        self.elliptic = sage.EllipticCurve(K, [_rational(a) for a in curve.ainvs])
        self.points = _torsion_points(self.elliptic)
        P = next(P for P in self.points if not (2*P).is_zero())
        for Q in self.points:
            words = tuple(a*P+b*Q for a in range(4) for b in range(4))
            if len(set(words)) == 16:
                break
        else:
            raise ArithmeticError('no basis of the exact 4-torsion module')
        self.generator_indices = self.points.index(P), self.points.index(Q)
        self.coordinates = tuple(divmod(words.index(T), 4) for T in self.points)
        self.addition = tuple(tuple(self.coordinates.index(((a+c)%4,(b+d)%4))
                              for c,d in self.coordinates) for a,b in self.coordinates)
        A, B = _translation_matrix(self.elliptic,P), _translation_matrix(self.elliptic,Q)
        identity = sage.identity_matrix(K,4)
        if not _projective_equal(A**4,identity) or not _projective_equal(B**4,identity):
            raise ArithmeticError('theta generators do not have projective order four')
        AB, BA = A*B, B*A
        index = next(i for i,v in enumerate(BA.list()) if v)
        scalar = AB.list()[index]/BA.list()[index]
        if AB != scalar*BA or scalar*scalar != -1:
            raise ArithmeticError('theta-generator commutator is not primitive of order four')
        matrices = tuple(_normalized_matrix(A**a*B**b) for a,b in self.coordinates)
        for M in matrices:
            M.set_immutable()
        self.matrices = matrices
        self.characters = tuple(tuple(scalar**(a*d-b*c) for c,d in self.coordinates)
                                for a,b in self.coordinates)
        if len(set(self.characters)) != 16:
            raise ArithmeticError('theta commutator pairing is degenerate')
        a1, _, a3, _, _ = self.elliptic.a_invariants()
        for shift in count():
            nodes = tuple(P[0]+shift*(P[1]+(a1*P[0]+a3)/2) for P in self.points[1:])
            if len(set(nodes)) == 15:
                self.separating_shift, self.nodes = shift, nodes
                break
        # Derive the degree-six x-polynomial of exact-order-four points from
        # h(x(2P))=0. The degree-twelve numerator is its square. All roots and
        # simplicity are checked against the actual torsion points below.
        x,z = sp.symbols('x z')
        b2,b4,b6,_ = curve.b_invariants
        h = x**3+(b2/4)*x*x+(b4/2)*x+b6/4
        n = sp.expand(sp.diff(h,x)**2-4*h*(2*x+b2/4))
        numerator = sp.Poly(n**3+(b2/4)*n*n*(4*h)+(b4/2)*n*(4*h)**2+(b6/4)*(4*h)**3,x,domain=sp.QQ)
        q = numerator.sqf_part().monic()
        if q.degree()!=6 or numerator.monic()!=q*q:
            raise ArithmeticError('duplication did not give the six exact-order-four x coordinates')
        polynomial = sp.Poly(h.subs(x,z)*sp.resultant(q.as_expr(),(z-x)**2-shift*shift*h,x),z,domain=sp.QQ).monic()
        if polynomial.degree()!=15 or polynomial.discriminant()==0:
            raise ArithmeticError('torsion separating polynomial has the wrong degree or repeated roots')
        self.node_polynomial = tuple(sage.QQ(str(c)) for c in reversed(polynomial.all_coeffs()))
        R = sage.PolynomialRing(K,'z')
        polynomial = R(self.node_polynomial)
        if any(polynomial(a)!=0 for a in self.nodes):
            raise ArithmeticError('rational torsion polynomial does not vanish at the actual nodes')
        z = R.gen()
        derivative = polynomial.derivative()
        self.interpolants = tuple((polynomial.quo_rem(z-a)[0]/derivative(a)).list() for a in self.nodes)

    def verify(self):
        """Recheck the actual chord actions, generators, characters and nodes."""
        sage,K = _sage(),_field()
        if len(self.points)!=16 or len(set(self.points))!=16 or not self.points[0].is_zero():
            return False
        if any(4*P!=self.elliptic(0) for P in self.points):
            return False
        i,j = self.generator_indices
        A = _translation_matrix(self.elliptic,self.points[i])
        B = _translation_matrix(self.elliptic,self.points[j])
        if self.matrices[i]!=A or self.matrices[j]!=B:
            return False
        identity=sage.identity_matrix(K,4)
        if not _projective_equal(A**4,identity) or not _projective_equal(B**4,identity):
            return False
        AB,BA=A*B,B*A
        index=next(k for k,v in enumerate(BA.list()) if v)
        scalar=AB.list()[index]/BA.list()[index]
        if AB!=scalar*BA or scalar*scalar!=-1:
            return False
        if len(self.coordinates)!=16 or set(self.coordinates)!=set((a,b) for a in range(4) for b in range(4)):
            return False
        for T,(a,b),M,row in zip(self.points,self.coordinates,self.matrices,self.characters):
            if T!=a*self.points[i]+b*self.points[j] or M!=_normalized_matrix(A**a*B**b):
                return False
            if row!=tuple(scalar**(a*d-b*c) for c,d in self.coordinates):
                return False
        a1,_,a3,_,_=self.elliptic.a_invariants()
        nodes=tuple(P[0]+self.separating_shift*(P[1]+(a1*P[0]+a3)/2) for P in self.points[1:])
        if self.nodes!=nodes or len(set(nodes))!=15:
            return False
        R=sage.PolynomialRing(K,'z'); polynomial=R(self.node_polynomial)
        if polynomial.degree()!=15 or not polynomial.is_monic() or any(polynomial(a)!=0 for a in nodes):
            return False
        return all(R(row)(b)==int(i==j) for i,row in enumerate(self.interpolants) for j,b in enumerate(nodes))

    def cochain(self, trivialization):
        if trivialization.covering.target != self.curve or not trivialization.verify():
            raise ValueError('trivialization does not have this exact elliptic target')
        B = trivialization.matrix
        inverse = B.inverse()
        values = []
        for M in self.matrices:
            action = inverse*M*B
            scale = next(v for v in action.list() if v)
            values.append(1/scale)
        return tuple(values)

    def coboundary(self, values):
        """Decide equality of cohomology classes by all sixteen exact characters."""
        sage, K = _sage(), _field()
        values = tuple(K(c) for c in values)
        if len(values) != 16 or not all(values):
            raise ValueError('a theta cochain must be a unit at every torsion point')
        # A character times a rational constant is frequent (including every
        # identity test). Check all these exact cases before any interpolation.
        for index,character in enumerate(self.characters):
            twisted=tuple(a*b for a,b in zip(values,character))
            if all(v==twisted[0] for v in twisted):
                try:
                    rational=sage.QQ(twisted[0])
                except (TypeError,ValueError):
                    continue
                result=ThetaCoboundary(index,rational,(rational,)+(sage.QQ(0),)*14)
                if not result.verify(self,values):
                    raise ArithmeticError('constant theta-cochain witness failed')
                return result
        for index, character in enumerate(self.characters):
            twisted = tuple(a*b for a, b in zip(values, character))
            try:
                origin = sage.QQ(twisted[0])
                coefficients = []
                for degree in range(15):
                    c = sum((v*poly[degree] for v, poly in zip(twisted[1:], self.interpolants)), K.zero())
                    coefficients.append(sage.QQ(c))
            except (TypeError, ValueError):
                continue  # Exact nonrationality, never a numerical rejection.
            result = ThetaCoboundary(index, origin, tuple(coefficients))
            if not result.verify(self, values):
                raise ArithmeticError('rational theta-cochain interpolation failed')
            return result
        return None
