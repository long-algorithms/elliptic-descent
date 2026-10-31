"""Exact rational 2-torsion actions on relative covering diagrams.

Global 4-cover equivalence allows translations of the intermediate 2-cover.
Retain the translations, their norm-pair action, and actual equivalence matrices.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .cubic_arithmetic import T, _poly, kummer_algebra
from .etale import EtaleQuarticAlgebra
from .geometry import U, _rat, construct_four_cover, _sqrt_rational
from .rational_root import RationalRootIdentification
from .two_cover import X, Y, Z


def _matrix_tuple(M):
    return tuple(tuple(_rat(v) for v in row) for row in M.tolist())


def _veronese_pair(values):
    a, b, c = values
    if a*c != b*b or not any(values):
        raise ArithmeticError('invalid projective conic image')
    return sp.Matrix((a, b) if a else (b, c))


def _mobius_from_veronese(B):
    """Recover a PGL2 map from three rational conic points, then verify it."""
    p = _veronese_pair(B*sp.Matrix((1, 0, 0)))
    q = _veronese_pair(B*sp.Matrix((0, 0, 1)))
    r = _veronese_pair(B*sp.Matrix((1, 1, 1)))
    coefficients = sp.Matrix.hstack(p, q).inv()*r
    M = sp.Matrix.hstack(coefficients[0]*p, coefficients[1]*q)
    a, b = M*sp.Matrix((X, Z))
    transformed = B*sp.Matrix((X*X, X*Z, Z*Z))
    proposed = sp.Matrix((a*a, a*b, b*b))
    old, new = (sp.Poly(v, X, Z) for v in (transformed[0], proposed[0]))
    monomial, coefficient = new.terms()[0]
    scale = _rat(old.coeff_monomial(monomial)/coefficient)
    if not M.det() or not scale or any(sp.expand(v-scale*w) for v,w in zip(transformed, proposed)):
        raise ArithmeticError('Veronese automorphism does not match the recovered Mobius map')
    return M, scale


def torsion_multiplier(presentation, point):
    """Even signs for T in E[2]: + at its root, - at the other roots."""
    D = presentation.auxiliary
    P = D.curve.point(point)
    if P is None:
        return D.algebra.one
    if D.curve.double(P) is not None:
        raise ValueError('the intermediate translation must be rational 2-torsion')
    index = next((i for i,g in enumerate(D.algebra.factors) if len(g)==2 and -g[0]==P[0]), None)
    if index is None:
        raise ArithmeticError('rational 2-torsion root is absent from the labelled cubic')
    return tuple(2*c-int(j==0) for j,c in enumerate(D.algebra.idempotents[index]))


@dataclass(frozen=True)
class QuarticDeckAction:
    """Translation on the intermediate torsor, over its marked E.

    matrix acts on (X,Z); y_scale acts on Y of weight two. point belongs to
    the WORKING E. The cubic sign action retains its actual E[2] label.
    """
    presentation: object
    point: object
    matrix: tuple
    y_scale: sp.Rational
    multiplier: tuple

    @property
    def source(self):
        return self.presentation.quartic

    def verify(self):
        try:
            P = self.presentation
            A = P.auxiliary.algebra
            w = torsion_multiplier(P, self.point)
            M, k = sp.Matrix(self.matrix), _rat(self.y_scale)
            if (tuple(self.multiplier) != w or M.shape != (2,2) or not M.det() or not k
                    or any(v.is_Rational is not True for v in M)
                    or A.multiply(w,w) != A.one or A.norm(w) != 1):
                return False
            W = sp.Matrix.hstack(*(sp.Matrix(A.multiply(w,tuple(int(i==j) for j in range(3))))
                                  for i in range(3)))
            a,b = M*sp.Matrix((X,Z))
            changed = tuple(sp.expand(f.subs({X:a,Z:b},simultaneous=True)) for f in P.parametrization.forms)
            original = W*sp.Matrix(P.parametrization.forms)
            if any(sp.expand(v-k*z) for v,z in zip(changed,original)):
                return False
            # z on the auxiliary curve is unchanged: p(M)=k*W*p and Y'=kY.
            return (P.parametrization.verify(P.auxiliary.q[2])
                    and sp.expand(self.source.binary_form.subs({X:a,Z:b},simultaneous=True)
                                  -k*k*self.source.binary_form) == 0)
        except (TypeError, ValueError, ArithmeticError):
            return False

    def same_transformation(self, other):
        """Weighted projective equality, allowing scalar-equivalent GL2 matrices."""
        if self.source != other.source:
            return False
        A,B = sp.Matrix(self.matrix),sp.Matrix(other.matrix)
        i = next(i for i,c in enumerate(B) if c)
        ratio = A[i]/B[i]
        return ratio != 0 and A == ratio*B and self.y_scale == ratio**2*other.y_scale

    def evaluate(self, point):
        if len(point)!=3:
            raise ValueError('expected weighted quartic coordinates')
        a,b,c = map(_rat,point)
        if not (a or c) or b*b != self.source.binary_form.subs({X:a,Z:c},simultaneous=True):
            raise ValueError('point is not on the intermediate quartic')
        u,v = sp.Matrix(self.matrix)*sp.Matrix((a,c))
        return u,self.y_scale*b,v

    @cached_property
    def algebra_data(self):
        A = EtaleQuarticAlgebra(self.source.f)
        a,b,c,d = sp.Matrix(self.matrix)
        t = (a,-c,sp.S.Zero,sp.S.Zero)
        theta = A.multiply((-b,d,0,0),A.inverse(t))
        powers = [A.one]
        for _ in range(3):
            powers.append(A.multiply(powers[-1],theta))
        S = sp.Matrix.hstack(*(sp.Matrix(v) for v in powers))
        if not S.det() or A.norm(t) != self.y_scale**2:
            raise ArithmeticError('branch algebra action or its norm is invalid')
        residual, power = [sp.S.Zero]*4, A.one
        for coefficient in self.source.f:
            residual = [v+coefficient*w for v,w in zip(residual,power)]
            power = A.multiply(power,theta)
        if any(residual):
            raise ArithmeticError('branch automorphism does not preserve the quartic polynomial')
        L = sp.Matrix.hstack(*(sp.Matrix(A.multiply(t,v)) for v in powers))
        return A,t,S,L

    def transport(self, model):
        if model.f != self.source.f or not self.verify():
            raise ValueError('norm pair and deck action have different or invalid sources')
        A,t,S,_ = self.algebra_data
        xi = A.multiply(t,tuple(S*sp.Matrix(model.xi)))
        raw = construct_four_cover(model.f,xi,norm_witness=self.y_scale*model.norm_witness)
        return NormPairAction(self,model,raw)


def quartic_deck_action(presentation, point):
    D = presentation.auxiliary
    point = D.curve.point(point)
    w = torsion_multiplier(presentation,point)
    W = sp.Matrix.hstack(*(sp.Matrix(D.algebra.multiply(w,tuple(int(i==j) for j in range(3))))
                          for i in range(3)))
    P = sp.Matrix(presentation.parametrization.matrix)
    M,scale = _mobius_from_veronese(P.inv()*W*P)
    result = QuarticDeckAction(presentation,point,_matrix_tuple(M),1/scale,w)
    if not result.verify():
        raise ArithmeticError('rational 2-torsion translation failed exact verification')
    return result


@dataclass(frozen=True)
class NormPairAction:
    """raw(V) = deck(source(U)) with V=sigma(U), including the signed norm."""
    action: QuarticDeckAction
    source: object
    raw: object

    def verify(self):
        try:
            g,D,R = self.action,self.source,self.raw
            if not g.verify() or D.f != g.source.f or R.f != D.f:
                return False
            A,t,S,L = g.algebra_data
            if (R.xi != A.multiply(t,tuple(S*sp.Matrix(D.xi)))
                    or R.norm_witness != g.y_scale*D.norm_witness
                    or not D.check_construction() or not R.check_construction()):
                return False
            substitutions = dict(zip(U,S*sp.Matrix(U)))
            if any(sp.expand(q.subs(substitutions,simultaneous=True)-r)
                   for q,r in zip(R.q,L*sp.Matrix(D.q))):
                return False
            a,b,c,d = sp.Matrix(g.matrix)
            if L[:,:2] != sp.Matrix([[a,-b],[-c,d],[0,0],[0,0]]) or not L[2:,2:].det():
                return False
            return sp.expand(R.map_xyz[1].subs(substitutions,simultaneous=True)
                             -g.y_scale*D.map_xyz[1]) == 0
        except (TypeError, ValueError, ArithmeticError):
            return False


@dataclass(frozen=True)
class FibreEquivalence:
    """An explicit target -> source isomorphism over the original elliptic E.

    The intermediate map is a recorded rational 2-torsion translation.
    scalar_square identifies its raw norm pair with the selected target pair.
    """
    transported: NormPairAction
    target: object
    scalar_square: object

    @property
    def coordinate_matrix(self):
        S = self.transported.action.algebra_data[2]
        return sp.ImmutableMatrix(S.inv()*self.scalar_square.coordinate_matrix(self.target.f))

    def verify(self):
        return (self.transported.verify()
                and self.scalar_square.verify(self.transported.raw,self.target)
                and self.coordinate_matrix.det() != 0)


def zero_fibre_kummer(labelled_cover, relative_model, arithmetic):
    """Actual kernel-map class in the cubic algebra, not a cardinality label.

    For the zero intermediate class choose a rational branch root alpha. Put
    A=f'(alpha), s=(B-f2)/3, j=u^2*theta+r and beta=alpha+A/(j-s). Then
    delta=(s-j)*xi(beta)/xi(alpha), with x_work-theta equal to delta times
    (U(beta)/(u*U(alpha)))^2. Both marked maps and the branch substitution agree.
    """
    C = labelled_cover
    if C.ambient_label != 0 or C.marking.translation is not None or relative_model.f != C.quartic.f:
        raise ValueError('kernel coordinates require the marked zero intermediate class')
    I = RationalRootIdentification(C.quartic)
    change = C.marking.change.then(C.selmer.curves.to_working)
    if change.source.ainvs != I.source.ainvs or change.target != C.selmer.curves.working:
        raise ValueError('zero-cover identification has inconsistent endpoints')
    if arithmetic.algebra != kummer_algebra(change.target):
        raise ValueError('wrong cubic arithmetic for the zero fibre')
    A = arithmetic.algebra
    u,r,_,_ = change.urst
    j = sp.Poly(u*u*T+r,T,domain=sp.QQ)
    h = _poly(A.f)
    derivative = I.taylor_coefficients[0]
    beta = (sp.Poly(I.root,T,domain=sp.QQ)+derivative*sp.invert(j-I.shift,h)).rem(h)
    def evaluate(coefficients, argument):
        result = sp.Poly(0,T,domain=sp.QQ)
        for coefficient in reversed(coefficients):
            result = (result*argument+coefficient).rem(h)
        return result
    if not evaluate(C.quartic.f,beta).is_zero:
        raise ArithmeticError('zero-fibre branch substitution does not satisfy f')
    xi = relative_model.xi
    scalar = sum(c*I.root**i for i,c in enumerate(xi))
    if not scalar:
        raise ValueError('relative class is not a branch-algebra unit')
    delta = ((-j+I.shift)*evaluate(xi,beta)).mul_ground(1/scalar).rem(h)
    result = tuple(delta.nth(i) for i in range(3))
    _sqrt_rational(A.norm(result))
    return result
