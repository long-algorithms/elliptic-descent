"""Contracted products and exact equivariant isomorphisms of origin torsors.

Adding E[4]-torsors means quotienting (P,Q) by (P+T,Q-T), not multiplying
unlabelled algebra isomorphism types. Every discovered isomorphism is checked
against the retained E[4] action. Discovery uses exact number-field roots;
verification uses rational polynomial and matrix identities.
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import count

from .sage import _sage
from .sage_origin_torsors import (EtaleAlgebra, origin_geometry, field_flex_matrix,
                                 verify_field_flex, _q)


class TensorAlgebra:
    """Q[a,b]/(f(a),g(b)), in lexicographic tensor power-basis coordinates."""
    def __init__(self, left, right):
        s = _sage()
        self.left_algebra, self.right_algebra = left, right
        self.ring = s.PolynomialRing(s.QQ, names=('a','b'), order='degrevlex')
        self.a, self.b = self.ring.gens()
        self.ideal = self.ring.ideal(left.polynomial(self.a), right.polynomial(self.b))
        self.monomials = tuple(self.a**i*self.b**j for i in range(left.degree) for j in range(right.degree))

    def reduce(self, p):
        return self.ideal.reduce(self.ring(p))

    def left(self, p):
        return self.ring(self.left_algebra.ring(p).lift()(self.a))

    def right(self, p):
        return self.ring(self.right_algebra.ring(p).lift()(self.b))

    def multiply(self, p, q):
        return self.reduce(p*q)

    def powers(self, p, length):
        answer, current = [], self.ring.one()
        for _ in range(length):
            answer.append(current)
            current = self.multiply(current,p)
        return tuple(answer)

    def vector(self, p):
        p = self.reduce(p)
        return _sage().vector(_sage().QQ, [p.monomial_coefficient(m) for m in self.monomials])

    def matrix(self, p):
        s = _sage()
        return s.matrix(s.QQ, self.left_algebra.degree, self.right_algebra.degree, self.vector(p))

    def element(self, vector):
        if len(vector)!=len(self.monomials):
            raise ValueError('incorrect tensor coordinate count')
        return sum((c*m for c,m in zip(vector,self.monomials)), self.ring.zero())

    def substitute(self, p, first, second):
        """Evaluate in the quotient; reduce throughout rather than expand huge powers."""
        p = self.ring(p)
        pa = self.powers(first,int(p.degree(self.a))+1)
        pb = self.powers(second,int(p.degree(self.b))+1)
        return self.reduce(sum((c*self.multiply(pa[i],pb[j]) for (i,j),c in p.dict().items()),self.ring.zero()))

    def inverse(self, p):
        """Factorwise Euclidean inverses; do not invert a 256-by-256 rational matrix."""
        s = _sage(); A,B=self.left_algebra,self.right_algebra
        M = self.matrix(p)
        components=[]
        for factor in A.factors:
            K,root = A.component(factor)
            R=s.PolynomialRing(K,'x')
            coefficients=[sum((K(M[i,j])*root**i for i in range(A.degree)),K(0)) for j in range(B.degree)]
            inverse=R(coefficients).inverse_mod(R(B.polynomial.list()))
            components.append([[inverse[j]] if factor.degree()==1 else inverse[j].list() for j in range(B.degree)])
        result=sum((self.left(A.reconstruct([part[j] for part in components]))*self.b**j
                    for j in range(B.degree)),self.ring.zero())
        if self.multiply(p,result)!=1:
            raise ArithmeticError('tensor inverse failed exact multiplication')
        return result


class OriginContext:
    """A rational E[4] model obtained from the verified zero covering's O-fibre."""
    def __init__(self, zero_covering):
        s=_sage()
        self.geometry=origin_geometry(zero_covering)
        self.algebra=self.geometry.algebra
        factor=next((g for g in self.algebra.factors if g.degree()==1),None)
        if factor is None:
            raise ValueError('the zero covering must have a rational origin lift')
        self.identity_root=-factor[0]
        self.identity_point=tuple(self.algebra.evaluate(p,self.identity_root) for p in self.geometry.point)
        self.trivialization=field_flex_matrix(zero_covering,self.identity_point,s.QQ)
        self.elliptic_coordinates=tuple(sum((self.trivialization[i,j]*self.geometry.point[j] for j in range(4)),self.algebra.ring.zero()) for i in range(4))
        z0,z1,z2,z3=self.elliptic_coordinates
        a1,_,a3,_,_=zero_covering.target.ainvs
        negative=(z0,z1,-z2-_q(a1)*z1-_q(a3)*z0,z3)
        inverse=self.trivialization.inverse()
        Q=tuple(sum((inverse[i,j]*negative[j] for j in range(4)),self.algebra.ring.zero()) for i in range(4))
        denominator=sum(a*q for a,q in zip(self.geometry.affine_form,Q))
        self.antipode=sum(a*q for a,q in zip(self.geometry.separating_form,Q))/denominator
        self.antipode_matrix=s.matrix(s.QQ,[self.algebra.vector(self.antipode**i) for i in range(16)]).transpose()
        if not self.verify():
            raise ArithmeticError('rational E[4] context failed exact verification')

    def verify(self):
        A=self.algebra
        expected=tuple(sum((self.trivialization[i,j]*self.geometry.point[j]
                            for j in range(4)),A.ring.zero()) for i in range(4))
        if expected!=self.elliptic_coordinates:
            return False
        return (self.geometry.verify() and A.degree==16
                and verify_field_flex(self.geometry.covering,self.identity_point,_sage().QQ,self.trivialization)
                and A.polynomial(self.antipode)==0
                and A.evaluate(self.antipode,self.antipode)==A.generator
                and A.evaluate(self.antipode,self.identity_root)==self.identity_root)

    def torsor(self, covering):
        if covering.target!=self.geometry.covering.target:
            raise ValueError('origin torsors must have the same marked elliptic target')
        geometry=self.geometry if covering==self.geometry.covering else origin_geometry(covering)
        T=TensorAlgebra(geometry.algebra,self.algebra)
        inverse=geometry.inverse_flex_matrix
        Q=tuple(sum((T.left(inverse[i,j])*T.right(self.elliptic_coordinates[j]) for j in range(4)),T.ring.zero()) for i in range(4))
        denominator=sum(a*q for a,q in zip(geometry.affine_form,Q))
        numerator=sum(a*q for a,q in zip(geometry.separating_form,Q))
        action=T.multiply(numerator,T.inverse(denominator))
        result=FiniteTorsor(self,geometry.algebra,action,geometry=geometry,
                            denominator_inverse=T.inverse(denominator))
        if not result.verify():
            raise ArithmeticError('origin torsor action failed exact verification')
        return result


@dataclass
class FiniteTorsor:
    context: OriginContext
    algebra: EtaleAlgebra
    action: object
    geometry: object=None
    contraction: object=None
    denominator_inverse: object=None

    @cached_property
    def tensor(self):
        return TensorAlgebra(self.algebra,self.context.algebra)

    @cached_property
    def action_powers(self):
        return self.tensor.powers(self.tensor.ring(self.action),16)

    @cached_property
    def action_matrices(self):
        return tuple(self.tensor.matrix(p) for p in self.action_powers)

    def verify(self):
        T=self.tensor; A=self.algebra
        if A.degree!=16 or self.context.algebra.degree!=16:
            return False
        action=T.ring(self.action)
        value=sum((c*p for c,p in zip(A.polynomial.list(),T.powers(action,17))),T.ring.zero())
        if T.reduce(value):
            return False
        # Identity must be exactly the retained rational origin in E[4].
        identity_vector=T.matrix(action)*_sage().vector(_sage().QQ,
                           [self.context.identity_root**i for i in range(16)])
        if A.element(identity_vector)!=A.generator:
            return False
        # One geometric component identifies this torsor with E[4] over its
        # residue field. This determinant also checks freeness/transitivity.
        factor=min(A.factors,key=lambda g:g.degree())
        K,root=A.component(factor)
        orbit=_sage().matrix(K,16,16,[sum((K(M[i,j])*root**i for i in range(16)),K(0))
                   for M in self.action_matrices for j in range(16)])
        if not orbit.det():
            return False
        if self.geometry is not None:
            if not self.geometry.verify() or self.geometry.algebra is not A:
                return False
            inverse=self.geometry.inverse_flex_matrix
            if tuple(g for g,_,_,_,_ in self.geometry.flexes)!=A.factors:
                return False
            for g,K,root,P,B in self.geometry.flexes:
                if g(root)!=0 or P!=tuple(A.evaluate(p,root) for p in self.geometry.point):
                    return False
                actual=_sage().matrix(K,4,4,lambda i,j:A.evaluate(inverse[i,j],root))
                if actual*B!=_sage().identity_matrix(K,4) or not verify_field_flex(self.geometry.covering,P,K,B):
                    return False
            Q=tuple(sum((T.left(inverse[i,j])*T.right(self.context.elliptic_coordinates[j])
                         for j in range(4)),T.ring.zero()) for i in range(4))
            denominator=sum(a*q for a,q in zip(self.geometry.affine_form,Q))
            numerator=sum(a*q for a,q in zip(self.geometry.separating_form,Q))
            return (self.denominator_inverse is not None
                    and T.multiply(denominator,T.ring(self.denominator_inverse))==1
                    and T.multiply(numerator,T.ring(self.denominator_inverse))==action)
        return self.contraction is not None and self.contraction.verify(self)


def _projector(left,right):
    s=_sage(); H=left.context.algebra
    if left.context is not right.context:
        raise ValueError('contracted torsors must share the same E[4] context')
    gram=H.trace_pairing/16
    inverse=left.context.antipode_matrix
    right_matrices=tuple(M*inverse.transpose() for M in right.action_matrices)
    columns=[(A*gram*B.transpose()).list() for A in left.action_matrices for B in right_matrices]
    return s.matrix(s.QQ,columns).transpose()


@dataclass
class ContractedProduct:
    left: FiniteTorsor
    right: FiniteTorsor
    primitive: object
    power_matrix: object

    @cached_property
    def tensor(self):
        return TensorAlgebra(self.left.algebra,self.right.algebra)

    def coaction(self):
        """Translate the first factor; quotient by the anti-diagonal action."""
        s=_sage(); G=self.tensor.matrix(self.primitive)
        M=self.left.action_matrices
        columns=[]
        for h in range(16):
            transformation=s.matrix(s.QQ,16,16,lambda i,j:M[j][i,h])
            image=s.vector(s.QQ,(transformation*G).list())
            columns.append(self.power_matrix.solve_right(image))
        return s.matrix(s.QQ,columns).transpose()

    def verify(self,result):
        T=self.tensor; M=self.power_matrix
        if (result.context is not self.left.context or self.right.context is not result.context
                or not self.left.verify() or not self.right.verify()):
            return False
        if M.nrows()!=256 or M.ncols()!=16 or M.rank()!=16:
            return False
        powers=T.powers(T.ring(self.primitive),17)
        actual=_sage().matrix(_sage().QQ,[T.vector(p) for p in powers[:16]]).transpose()
        if actual!=M:
            return False
        projector=_projector(self.left,self.right)
        if projector.rank()!=16 or projector*M!=M:
            return False
        if T.reduce(sum((c*p for c,p in zip(result.algebra.polynomial.list(),powers)),T.ring.zero())):
            return False
        return result.tensor.matrix(result.action)==self.coaction()


def contracted_sum(left,right):
    """The actual H^1 sum: (P x Q)/(T,-T), with its residual E[4] action."""
    s=_sage(); T=TensorAlgebra(left.algebra,right.algebra)
    projector=_projector(left,right)
    if projector.rank()!=16:
        raise ArithmeticError('anti-diagonal invariants do not have dimension sixteen')
    basis=projector.column_space().basis()
    for n in count(1):
        v=sum((s.QQ(n)**i*b for i,b in enumerate(basis)),s.vector(s.QQ,256))
        g=T.element(v)
        powers=T.powers(g,17)
        M=s.matrix(s.QQ,[T.vector(p) for p in powers[:16]]).transpose()
        if M.rank()!=16:
            continue
        relation=M.solve_right(T.vector(powers[16]))
        algebra=EtaleAlgebra(tuple(-c for c in relation)+(s.QQ(1),))
        evidence=ContractedProduct(left,right,g,M)
        action_matrix=evidence.coaction()
        tensor=TensorAlgebra(algebra,left.context.algebra)
        result=FiniteTorsor(left.context,algebra,tensor.element(action_matrix.list()),contraction=evidence)
        if not result.verify():
            raise ArithmeticError('contracted product failed exact algebra/action identities')
        return result


@dataclass(frozen=True)
class TorsorIsomorphism:
    """Pullback of the source primitive element as a polynomial on the target."""
    coefficients: tuple

    def verify(self,source,target):
        s=_sage()
        if source.context is not target.context or len(self.coefficients)!=16:
            return False
        try:
            coefficients=tuple(s.QQ(c) for c in self.coefficients)
        except (TypeError,ValueError):
            return False
        A=target.algebra; value=A.element(coefficients)
        if source.algebra.polynomial(value):
            return False
        M=s.matrix(s.QQ,[A.vector(value**i) for i in range(16)]).transpose()
        if not M.det():
            return False
        left=sum((c*p for c,p in zip(coefficients,target.action_powers)),target.tensor.ring.zero())
        source_matrix=source.tensor.matrix(source.action)
        right=target.tensor.element((M*source_matrix).list())
        return target.tensor.reduce(left-right)==0


def torsor_isomorphism(source,target):
    """Complete exact discovery; algebra isomorphism alone is NEVER sufficient.

    Over one residue field K of target, its chosen point identifies target_K
    with E[4]_K. Each root of the source primitive polynomial gives a candidate
    equivariant map. It descends to Q iff all reconstructed coefficients are
    rational. This avoids exponentially combining choices in different factors.
    """
    s=_sage()
    if source.context is not target.context:
        raise ValueError('torsors use different E[4] contexts')
    A=target.algebra
    factor=min(A.factors,key=lambda g:g.degree())
    K,root=A.component(factor)
    orbit=s.matrix(K,16,16,lambda j,h:sum((K(target.action_matrices[j][i,h])*root**i
                                                 for i in range(16)),K(0))).transpose()
    if not orbit.det():
        raise ArithmeticError('target action does not identify its geometric orbit')
    R=s.PolynomialRing(K,'x')
    polynomial=R(source.algebra.polynomial.list())
    action=source.tensor.matrix(source.action)
    for value in polynomial.roots(multiplicities=False):
        rhs=s.vector(K,[sum((K(action[i,h])*value**i for i in range(16)),K(0)) for h in range(16)])
        coordinates=orbit.solve_right(rhs)
        try:
            rational=tuple(s.QQ(c) for c in coordinates)
        except (TypeError,ValueError):
            continue
        witness=TorsorIsomorphism(rational)
        if not witness.verify(source,target):
            raise ArithmeticError('rational equivariant map failed full polynomial verification')
        return witness
    return None
