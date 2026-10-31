"""Reuse complete labelled Sel2 data in the zero fibre by an exact linear map.

The original rational-root relative descent remains a reference route. Here
we transfer every accepted cubic model, including local evidence, rather than
redoing a large auxiliary quartic-field/support computation. This is not a
count-only shortcut: the inverse matrix, equations and signed maps are checked.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .cubic_arithmetic import T, _poly, kummer_algebra
from .geometry import U, construct_four_cover
from .rational_root import RationalRootIdentification
from .relative import RelativeCandidate
from .marked_cover import MarkedRelativeResult
from .selmer_quartic import LabelledRelativeLifts
from .status import LocalStatus


def _rows(matrix):
    return tuple(map(tuple, matrix.tolist()))


@dataclass(frozen=True)
class ZeroBranchData:
    presentation: object

    @cached_property
    def identification(self):
        return RationalRootIdentification(self.presentation.quartic)

    @cached_property
    def algebra(self):
        return kummer_algebra(self.presentation.auxiliary.curve)

    @cached_property
    def j(self):
        u,r,_,_ = self.presentation.to_working.urst
        return sp.Poly(u*u*T+r,T,domain=sp.QQ)

    @cached_property
    def beta(self):
        I = self.identification
        h = _poly(self.algebra.f)
        return (sp.Poly(I.root,T,domain=sp.QQ)
                +I.taylor_coefficients[0]*sp.invert(self.j-I.shift,h)).rem(h)

    @cached_property
    def evaluation(self):
        """A -> Q x L: value at the selected branch root and at beta."""
        h = _poly(self.algebra.f)
        powers = [sp.Poly(1,T,domain=sp.QQ)]
        for _ in range(3):
            powers.append((powers[-1]*self.beta).rem(h))
        alpha = self.identification.root
        return sp.ImmutableMatrix([[alpha**i for i in range(4)]]
                   +[[p.nth(k) for p in powers] for k in range(3)])

    @cached_property
    def coordinate_matrix(self):
        u = self.presentation.to_working.urst[0]
        return sp.ImmutableMatrix((self.evaluation[1:,:]/u).col_join(self.evaluation[:1,:]))

    @cached_property
    def checked_identity(self):
        P,I = self.presentation,self.identification
        if (P.auxiliary.delta != P.auxiliary.algebra.one or not P.verify()
                or not I.check_marked_map() or not self.evaluation.det()):
            return False
        h = _poly(self.algebra.f)
        value = sp.Poly(0,T,domain=sp.QQ)
        for coefficient in reversed(P.quartic.f):
            value = (value*self.beta+coefficient).rem(h)
        return value.is_zero and self.coordinate_matrix.det() != 0

    def verify(self):
        return self.checked_identity


@dataclass(frozen=True)
class ZeroFibreTransport:
    branch: ZeroBranchData
    auxiliary: object
    model: object
    equation_matrix: tuple

    @property
    def coordinate_matrix(self):
        """Relative coordinates -> auxiliary cubic coordinates."""
        return self.branch.coordinate_matrix

    def verify(self, *, check_branch=True):
        try:
            B,D,R = self.branch,self.auxiliary,self.model
            P,I = B.presentation,B.identification
            if (D.curve != P.auxiliary.curve or R.f != P.quartic.f
                    or (check_branch and not B.verify())
                    or not D.verify() or not R.check_construction()):
                return False
            evaluated = B.evaluation*sp.Matrix(R.xi)
            h = _poly(B.algebra.f)
            delta = ((-B.j+I.shift)*_poly(tuple(evaluated[1:]))).rem(h)
            if evaluated[0] != 1 or tuple(delta.nth(i) for i in range(3)) != D.delta:
                return False
            if R.norm_witness != -D.norm_witness/(R.f[4]*I.taylor_coefficients[0]):
                return False
            M,Q = self.coordinate_matrix,sp.Matrix(self.equation_matrix)
            if Q.shape != (2,2) or not Q.det() or any(v.is_Rational is not True for v in Q):
                return False
            sub = dict(zip(U,M*sp.Matrix(U)))
            if any(sp.expand(q.subs(sub,simultaneous=True)-v)
                   for q,v in zip(D.quadrics,Q*sp.Matrix(R.quadrics))):
                return False
            u = P.to_working.urst[0]
            z = (M*sp.Matrix(U))[3]
            # Together with delta*V^2=X-theta*z^2 this checks the completed-y
            # sign in I o (R->C) = change^-1 o (D->working E), not just x.
            return sp.expand(R.norm_u-z*u**3*D.norm_u.subs(sub,simultaneous=True)) == 0
        except (TypeError,ValueError,ArithmeticError):
            return False

    def to_auxiliary(self, point):
        values, image = self.model.evaluate(point)
        if any(values) or not any(point):
            raise ValueError('point is not on the relative covering')
        result = tuple(self.coordinate_matrix*sp.Matrix(point))
        self.auxiliary.evaluate(result)
        return result

    def from_auxiliary(self, point):
        self.auxiliary.evaluate(point)
        result = tuple(self.coordinate_matrix.inv()*sp.Matrix(point))
        self.to_auxiliary(result)
        return result


def transport_zero_cover(branch, auxiliary):
    if not isinstance(branch,ZeroBranchData):
        branch = ZeroBranchData(branch)
    if auxiliary.curve != branch.presentation.auxiliary.curve:
        raise ValueError('zero fibre and auxiliary cover have different working targets')
    I = branch.identification
    h = _poly(branch.algebra.f)
    ratio = (_poly(auxiliary.delta)*sp.invert(-branch.j+I.shift,h)).rem(h)
    xi = tuple(branch.evaluation.inv()*sp.Matrix([1]+[ratio.nth(i) for i in range(3)]))
    f = branch.presentation.quartic.f
    R = construct_four_cover(f,xi,norm_witness=-auxiliary.norm_witness/(f[4]*I.taylor_coefficients[0]))
    sub = dict(zip(U,branch.coordinate_matrix*sp.Matrix(U)))
    monomials = [U[i]*U[j] for i in range(4) for j in range(i,4)]
    system = sp.Matrix([[sp.Poly(q,*U,domain=sp.QQ).coeff_monomial(m) for q in R.quadrics]
                         for m in monomials])
    rows = []
    for q in auxiliary.quadrics:
        transformed = sp.Poly(q.subs(sub,simultaneous=True),*U,domain=sp.QQ)
        rhs = sp.Matrix([transformed.coeff_monomial(m) for m in monomials])
        solution,parameters = system.gauss_jordan_solve(rhs)
        if parameters.rows:
            raise ArithmeticError('zero-fibre equations are not an independent pair')
        rows.append(tuple(solution))
    result = ZeroFibreTransport(branch,auxiliary,R,tuple(rows))
    if not result.verify():
        raise ArithmeticError('zero-fibre linear transport failed its exact identities')
    return result


@dataclass(frozen=True)
class TransferredLocalEvidence:
    """Whole local evidence transported by a verified Q-isomorphism.

    Do not attach the old per-prime system keys to the new equations. The
    original complete local proof and the exact global isomorphism are retained.
    """
    transport: ZeroFibreTransport
    original: object

    @property
    def status(self):
        return self.original.status

    def verify(self, model, *, max_nodes=1_000_000):
        return (model == self.transport.model and self.transport.verify()
                and self.original.verify(self.transport.auxiliary,max_nodes=max_nodes))


@dataclass(frozen=True)
class TransportedZeroRelativeResult:
    """Complete accepted relative classes obtained from a complete Sel2 torsor.

    This is deliberately not an AffineF2 norm-fibre enumeration result. Its
    completeness is the explicit equivalence with all accepted parent classes.
    """
    labelled_source: object
    candidates: tuple
    transports: tuple
    method: str = 'verified_zero_fibre_transport_from_complete_selmer2'

    @property
    def source(self):
        return self.labelled_source.quartic

    @property
    def total_model_classes(self):
        return self.labelled_source.selmer.cardinality

    @property
    def assumptions(self):
        return self.labelled_source.assumptions

    @property
    def enumeration_complete(self):
        S = self.labelled_source.selmer
        return (S.complete and self.labelled_source.ambient_label == 0
                and len(self.candidates)==S.cardinality
                and tuple(c.label for c in self.candidates)==tuple(range(S.cardinality))
                and len(self.transports)==len(self.candidates))

    @property
    def unresolved(self):
        return tuple(c for c in self.candidates if c.local.status==LocalStatus.UNKNOWN)

    @property
    def complete(self):
        return self.enumeration_complete and all(c.local.status==LocalStatus.SOLUBLE for c in self.candidates)

    @property
    def certified(self):
        return self.complete and self.labelled_source.certified and not self.assumptions

    @property
    def coverings(self):
        return tuple(c.model for c in self.candidates if c.local.status==LocalStatus.SOLUBLE)

    @property
    def unoriented_coverings(self):
        return self.coverings  # The rational branch point identifies the signs.

    def verify_evidence(self, *, max_nodes=1_000_000, geometry=True, full_geometry=False):
        if not self.complete or not self.labelled_source.verify(max_nodes=max_nodes):
            return False
        S = self.labelled_source.selmer
        for c,t in zip(self.candidates,self.transports):
            old = [d for d in S.decisions if d.label==S.ambient_label(c.label)]
            if (len(old)!=1 or old[0].model!=t.auxiliary
                    or t.branch.presentation!=self.labelled_source.presentation
                    or c.model!=t.model or c.local.transport!=t or c.local.original!=old[0].local
                    or not c.local.verify(c.model,max_nodes=max_nodes)):
                return False
            if full_geometry and not all((c.model.check_map_identity(),c.model.check_smoothness(),
                                          c.model.check_map_basepoint_free())):
                return False
        return True


def zero_fibre_lifts(cover):
    if (cover.ambient_label!=0 or not cover.selmer.complete or not cover.verify()
            or cover.marking.translation is not None):
        raise ValueError('complete labelled zero class with its exact marking is required')
    S = cover.selmer
    branch = ZeroBranchData(cover.presentation)
    if not branch.verify():
        raise ArithmeticError('zero-fibre branch identification is invalid')
    lookup = {d.label:d for d in S.decisions if d.local.status==LocalStatus.SOLUBLE}
    candidates,transports = [],[]
    for c in range(S.cardinality):
        old = lookup[S.ambient_label(c)]
        t = transport_zero_cover(branch,old.model)
        candidates.append(RelativeCandidate(c,t.model,TransferredLocalEvidence(t,old.local)))
        transports.append(t)
    relative = TransportedZeroRelativeResult(cover,tuple(candidates),tuple(transports))
    if not relative.verify_evidence():
        raise ArithmeticError('transported zero fibre failed its parent/local evidence checks')
    return LabelledRelativeLifts(cover,MarkedRelativeResult(relative,cover.marking))
