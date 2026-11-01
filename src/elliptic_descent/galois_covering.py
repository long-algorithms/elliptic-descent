"""Portable exact Galois certificates for the actual marked 4-coverings.

A splitting-field discovery backend is separate. All data here use rational
power-basis coefficients, checked without Sage. The curve is trivialized using
a branch root and the classical half-point formula, retaining the signed norm.
Its flex-origin cocycle is then computed on every field automorphism.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .geometry import _rat, construct_four_cover
from .two_cover import QuarticTwoCover
from .weierstrass import WeierstrassModel, WeierstrassChange
from .finite_cohomology import FiniteFourModule, CohomologyClassTable

_t = sp.Symbol('a')


def _rational_tuple(values):
    return tuple(_rat(c) for c in values)


def polynomial_value(coefficients, x, zero):
    result = zero
    for c in reversed(coefficients):
        result = result*x+c
    return result


@dataclass(frozen=True)
class ExactGaloisField:
    """Q[a], a monic irreducible polynomial, all conjugates and generators.

    Having degree-many distinct roots in this field proves normality in
    characteristic zero. Automorphisms are verified as substitution maps.
    The first image is a itself. No numerical root matching is used in replay.
    """
    polynomial: tuple
    automorphism_images: tuple
    generators: tuple

    def __post_init__(self):
        f = _rational_tuple(self.polynomial)
        if len(f)<2 or f[-1] != 1:
            raise ValueError('expected a positive-degree monic field polynomial')
        object.__setattr__(self,'polynomial',f)
        object.__setattr__(self,'automorphism_images',tuple(_rational_tuple(v) for v in self.automorphism_images))
        object.__setattr__(self,'generators',tuple(self.generators))

    @property
    def degree(self):
        return len(self.polynomial)-1

    @cached_property
    def domain(self):
        p = sp.Poly.from_list(list(reversed(self.polynomial)),_t,domain=sp.QQ)
        if not p.is_irreducible:
            raise ValueError('field polynomial is reducible')
        # The formal root fixes a power basis; no floating-point embedding is used.
        return sp.QQ.algebraic_field((p,_t))

    def rational(self, value):
        return self.domain.convert(_rat(value))

    @property
    def zero(self):
        return self.domain.zero

    @property
    def one(self):
        return self.domain.one

    @property
    def generator(self):
        return self.domain.unit

    def unpack(self, coefficients):
        coefficients = _rational_tuple(coefficients)
        if len(coefficients) != self.degree:
            raise ValueError('field element has the wrong power-basis length')
        return polynomial_value(tuple(self.rational(c) for c in coefficients),self.generator,self.zero)

    def pack(self, value):
        values = tuple(sp.Rational(c.numerator,c.denominator) for c in reversed(value.to_list()))
        return values+(sp.S.Zero,)*(self.degree-len(values))

    @cached_property
    def images(self):
        return tuple(self.unpack(v) for v in self.automorphism_images)

    def apply(self, index, value):
        if type(index) is not int or not 0 <= index < len(self.images):
            raise ValueError('invalid field automorphism index')
        return polynomial_value(tuple(self.rational(c) for c in self.pack(value)),self.images[index],self.zero)

    @cached_property
    def transitions(self):
        lookup = {a:i for i,a in enumerate(self.images)}
        return tuple(tuple(lookup[self.apply(s,self.images[g])] for g in self.generators)
                     for s in range(len(self.images)))

    def verify(self):
        try:
            n = self.degree
            if (len(self.images)!=n or len(set(self.images))!=n or self.images[0]!=self.generator
                    or len(set(self.generators))!=len(self.generators)
                    or any(type(g) is not int or not 0<g<n for g in self.generators)):
                return False
            f = tuple(self.rational(c) for c in self.polynomial)
            if any(polynomial_value(f,a,self.zero) for a in self.images):
                return False
            seen,todo = {0},[0]
            while todo:
                for j in self.transitions[todo.pop()]:
                    if j not in seen:
                        seen.add(j);todo.append(j)
            return len(seen)==n
        except (TypeError,ValueError,KeyError):
            return False


class FieldCurve:
    """Exact general Weierstrass point operations in the certificate's field."""
    def __init__(self, field, ainvs):
        self.field = field
        self.model = WeierstrassModel(ainvs)
        self.a = tuple(field.rational(a) for a in self.model.ainvs)

    def point(self, P):
        if P is None:
            return None
        x,y = P
        a1,a2,a3,a4,a6 = self.a
        if y*y+a1*x*y+a3*y-x**3-a2*x*x-a4*x-a6:
            raise ValueError('algebraic point is not on the retained elliptic curve')
        return x,y

    def negate(self, P):
        if P is None:
            return None
        x,y = P
        return x,-y-self.a[0]*x-self.a[2]

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x,y = P;v,w = Q
        a1,a2,a3,a4,_ = self.a
        if x == v:
            if Q == self.negate(P):
                return None
            slope = (3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
        else:
            slope = (w-y)/(v-x)
        X = slope*slope+a1*slope-a2-x-v
        return self.point((X,-y-slope*(X-x)-a1*X-a3))

    def multiply(self, n, P):
        if type(n) is not int:
            raise TypeError('point multiplier must be an integer')
        if n<0:
            return self.multiply(-n,self.negate(P))
        result = None
        while n:
            if n & 1:
                result = self.add(result,P)
            P = self.add(P,P);n >>= 1
        return result


@dataclass(frozen=True)
class FourTorsionCertificate:
    field: ExactGaloisField
    ainvs: tuple
    roots: tuple
    basis: tuple

    @cached_property
    def curve(self):
        return FieldCurve(self.field,self.ainvs)

    @cached_property
    def root_values(self):
        return tuple(self.field.unpack(v) for v in self.roots)

    @cached_property
    def basis_points(self):
        if len(self.basis)!=2:
            raise ValueError('expected two 4-torsion generators')
        return tuple(self.curve.point(tuple(self.field.unpack(v) for v in P)) for P in self.basis)

    @cached_property
    def points(self):
        A,B = self.basis_points
        return tuple(self.curve.add(self.curve.multiply(a,A),self.curve.multiply(b,B))
                     for b in range(4) for a in range(4))

    @cached_property
    def point_indices(self):
        if len(set(self.points))!=16:
            raise ValueError('the claimed E[4] basis does not generate sixteen points')
        return {P:i for i,P in enumerate(self.points)}

    @cached_property
    def module(self):
        F = self.field
        lookup = self.point_indices
        actions = tuple(tuple(lookup[None if P is None else tuple(F.apply(s,v) for v in P)]
                              for P in self.points) for s in range(F.degree))
        return FiniteFourModule(actions,F.generators,F.transitions)

    def verify(self):
        try:
            if not self.field.verify() or len(self.roots)!=3 or len(set(self.root_values))!=3:
                return False
            b2,b4,b6,_ = self.curve.model.b_invariants
            f = tuple(self.field.rational(c) for c in (b6/4,b4/2,b2/4,1))
            if any(polynomial_value(f,e,self.field.zero) for e in self.root_values):
                return False
            if any(self.curve.multiply(4,P) is not None for P in self.basis_points):
                return False
            return len(self.point_indices)==16 and self.module.verify()
        except (TypeError,ValueError,KeyError,ZeroDivisionError):
            return False


@dataclass(frozen=True)
class CoverTrivialization:
    """A specific isomorphism psi:D->E over the splitting field, pi=[4]psi.

    ``roots_squared`` stores three square roots in the original cubic-root
    order, with their PRODUCT fixed by the signed relative norm witness.
    This distinguishes inverses in E[4]; an unsigned square class is insufficient.
    """
    key: tuple
    f: tuple
    xi: tuple
    norm_witness: object
    jacobian_to_original: tuple
    branch_root: tuple
    roots_squared: tuple

    def _data(self, torsion):
        F = torsion.field
        f = tuple(F.rational(c) for c in self.f)
        xi = tuple(F.rational(c) for c in self.xi)
        alpha = F.unpack(self.branch_root)
        r = F.rational(self.norm_witness)
        u,shift,_,_ = tuple(F.rational(c) for c in self.jacobian_to_original)
        first = polynomial_value(tuple((i+1)*f[i+1] for i in range(4)),alpha,F.zero)
        second = polynomial_value(tuple((i+1)*(i+2)*f[i+2]/2 for i in range(3)),alpha,F.zero)
        s = (second-f[2])/3
        j = tuple(u*u*e+shift for e in torsion.root_values)
        beta = tuple(alpha+first/(v-s) for v in j)
        xia = polynomial_value(xi,alpha,F.zero)
        delta = tuple((s-v)*polynomial_value(xi,b,F.zero)/xia for v,b in zip(j,beta))
        roots = tuple(F.unpack(v) for v in self.roots_squared)
        return f,xi,alpha,r,u,first,s,beta,xia,delta,roots

    def origin(self, torsion):
        """Exact Lagrange interpolation: U(alpha)=0, U(beta_i)=1/sqrt(delta_i)."""
        F = torsion.field
        _,_,alpha,_,_,_,_,beta,_,_,roots = self._data(torsion)
        nodes = (alpha,*beta)
        values = (F.zero,*(F.one/v for v in roots))
        result = [F.zero]*4
        for i,(node,value) in enumerate(zip(nodes,values)):
            if not value:
                continue
            p,denominator = [F.one],F.one
            for j,other in enumerate(nodes):
                if j==i:
                    continue
                q = [F.zero]*(len(p)+1)
                for k,c in enumerate(p):
                    q[k] -= other*c;q[k+1] += c
                p = q;denominator *= node-other
            for k,c in enumerate(p):
                result[k] += value*c/denominator
        return tuple(result)

    def verify(self, torsion):
        try:
            if (len(self.f)!=5 or len(self.xi)!=4 or len(self.roots_squared)!=3
                    or len(self.key)!=2 or any(type(c) is not int or c<0 for c in self.key)):
                return False
            source = QuarticTwoCover(self.f)
            change = WeierstrassChange(source.ainvs,self.jacobian_to_original)
            if not change.verify() or change.target != torsion.curve.model:
                return False
            D = construct_four_cover(self.f,self.xi,norm_witness=self.norm_witness)
            if not D.check_construction():
                return False
            F = torsion.field
            f,xi,alpha,r,_,A,_,beta,xia,delta,roots = self._data(torsion)
            if polynomial_value(f,alpha,F.zero) or not A or not xia:
                return False
            if len(set((alpha,*beta)))!=4 or any(polynomial_value(f,b,F.zero) for b in beta):
                return False
            if any(v*v!=d for v,d in zip(roots,delta)) or roots[0]*roots[1]*roots[2] != -f[4]*A*r/(xia*xia):
                return False
            P = self.origin(torsion)
            # All four branch evaluations establish the two quadric equations.
            for t in (alpha,*beta):
                value = polynomial_value(xi,t,F.zero)*polynomial_value(P,t,F.zero)**2
                if value != xia*(alpha-t)/A:
                    return False
            return self.evaluate(torsion,P) is None
        except (TypeError,ValueError,ZeroDivisionError,ArithmeticError):
            return False

    def evaluate(self, torsion, point):
        """Evaluate psi, including the four points above the chosen branch root.

        At U(alpha)=0, the three completed square roots agree up to signs.
        All equal gives O; the unique opposite sign selects its actual E[2]
        point. No limit or floating-point tolerance is used at this chart.
        """
        F = torsion.field
        _,_,alpha,_,u,_,_,beta,_,_,roots = self._data(torsion)
        z = polynomial_value(point,alpha,F.zero)
        w = tuple(s*polynomial_value(point,b,F.zero) for s,b in zip(roots,beta))
        if not z:
            if not w[0] or any(v*v != w[0]*w[0] for v in w):
                raise ValueError('invalid point in the branch-origin chart')
            if w[0] == w[1] == w[2]:
                return None
            i = next((i for i in range(3) if w[(i+1)%3] == w[(i+2)%3] and w[i] != w[(i+1)%3]),None)
            if i is None:
                raise ValueError('invalid branch signs')
            e = torsion.root_values[i]
            return torsion.curve.point((e,-(torsion.curve.a[0]*e+torsion.curve.a[2])/2))
        w = tuple(v/(u*z) for v in w)
        x = torsion.root_values[0]+w[0]*w[0]
        if any(v*v != x-e for v,e in zip(w,torsion.root_values)):
            raise ValueError('the point fails the three halving equations')
        product = w[0]*w[1]*w[2]
        pairs = w[0]*w[1]+w[0]*w[2]+w[1]*w[2]
        X = x+pairs
        V = -product+(w[0]+w[1]+w[2])*pairs
        Y = V-(torsion.curve.a[0]*X+torsion.curve.a[2])/2
        result = torsion.curve.point((X,Y))
        doubled = torsion.curve.point((x,product-(torsion.curve.a[0]*x+torsion.curve.a[2])/2))
        if torsion.curve.multiply(2,result) != doubled:
            raise ArithmeticError('exact half-point formula failed duplication')
        return result

    def cocycle(self, torsion):
        F = torsion.field
        P = self.origin(torsion)
        values = tuple(torsion.point_indices[self.evaluate(torsion,tuple(F.apply(s,v) for v in P))]
                       for s in range(F.degree))
        if not torsion.module.is_cocycle(values):
            raise ArithmeticError('flex-origin differences do not satisfy the Galois cocycle law')
        return values


@dataclass(frozen=True)
class GaloisGroupEvidence:
    torsion: FourTorsionCertificate
    coverings: tuple

    @cached_property
    def class_table(self):
        return CohomologyClassTable(self.torsion.module,tuple(c.key for c in self.coverings),
                                     tuple(c.cocycle(self.torsion) for c in self.coverings))

    def verify(self):
        try:
            return (self.torsion.verify() and all(c.verify(self.torsion) for c in self.coverings)
                    and self.class_table.verify())
        except (TypeError, ValueError, KeyError, IndexError, ArithmeticError):
            return False


def evidence_record(evidence):
    """JSON-safe exact coefficients; no number-field embeddings or floats."""
    T = evidence.torsion
    F = T.field
    vector = lambda v:[str(c) for c in v]
    return {'field':{'polynomial':vector(F.polynomial),
                     'automorphism_images':[vector(v) for v in F.automorphism_images],
                     'generators':F.generators},
            'ainvs':vector(T.ainvs),'cubic_roots':[vector(v) for v in T.roots],
            'torsion_basis':[[vector(v) for v in P] for P in T.basis],
            'coverings':[{'key':c.key,'f':vector(c.f),'xi':vector(c.xi),
                          'norm_witness':str(c.norm_witness),
                          'jacobian_to_original':vector(c.jacobian_to_original),
                          'branch_root':vector(c.branch_root),
                          'square_roots':[vector(v) for v in c.roots_squared]}
                         for c in evidence.coverings]}


def evidence_from_record(record):
    """Reconstruct data, without treating a record's provenance as a proof."""
    q = lambda v:tuple(_rat(c) for c in v)
    f = record['field']
    F = ExactGaloisField(q(f['polynomial']),tuple(q(v) for v in f['automorphism_images']),tuple(f['generators']))
    T = FourTorsionCertificate(F,q(record['ainvs']),tuple(q(v) for v in record['cubic_roots']),
                              tuple(tuple(q(v) for v in P) for P in record['torsion_basis']))
    C = tuple(CoverTrivialization(tuple(c['key']),q(c['f']),q(c['xi']),_rat(c['norm_witness']),
                q(c['jacobian_to_original']),q(c['branch_root']),tuple(q(v) for v in c['square_roots']))
              for c in record['coverings'])
    return GaloisGroupEvidence(T,C)
