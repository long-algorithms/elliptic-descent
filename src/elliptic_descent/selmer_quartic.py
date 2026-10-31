"""Binary quartic presentations of arbitrary labelled 2-Selmer classes.

Parametrize the conic Q2=0 in the cubic-algebra auxiliary cover, and retain
both exact inverse maps and the x-theta identity. Recover the marking from
polynomial identities, not matching invariants or an imported cover list.
"""
from dataclasses import dataclass
from functools import cached_property

import sympy as sp

from .arithmetic import nonnegative_integer
from .backends.sage_conics import rational_conic_point
from .conic_parametrization import ConicParametrization, parametrize_conic
from .cubic_local import CubicCover
from .geometry import U, _rat
from .marked_cover import MarkedQuarticCover, function_equal
from .status import LocalStatus
from .two_cover import QuarticTwoCover, X, Y, Z
from .weierstrass import WeierstrassChange, x, y


@dataclass(frozen=True)
class QuarticPresentation:
    """Exact isomorphism from a quartic double cover to one auxiliary 2-cover.

    This certifies a presentation and its marking, not Selmer membership or
    class-group completeness. The LabelledQuarticCover layer supplies that link.
    """
    auxiliary: CubicCover
    parametrization: ConicParametrization

    def __post_init__(self):
        if not isinstance(self.auxiliary, CubicCover) or not isinstance(self.parametrization, ConicParametrization):
            raise TypeError('expected an auxiliary cubic cover and a conic parametrization')
        if not self.parametrization.verify(self.auxiliary.q[2]):
            raise ValueError('parametrization belongs to a different conic')
        _ = self.quartic  # Check genuine degree four and squarefreeness eagerly.

    @cached_property
    def pulled_q(self):
        sub = dict(zip(U[:3], self.parametrization.forms))
        return tuple(sp.expand(q.subs(sub, simultaneous=True)) for q in self.auxiliary.q)

    @cached_property
    def quartic(self):
        p = sp.Poly(-self.pulled_q[1], X, Z, domain=sp.QQ)
        return QuarticTwoCover(tuple(p.coeff_monomial(X**i*Z**(4-i)) for i in range(5)))

    @cached_property
    def norm_form(self):
        sub = dict(zip(U[:3], self.parametrization.forms))
        return sp.expand(self.auxiliary.norm_witness*self.auxiliary.norm_u.subs(sub, simultaneous=True))

    @property
    def map_working_xy(self):
        """Raw x-theta covering map to the working elliptic curve, not its Jacobian model."""
        a1, _, a3, _, _ = self.auxiliary.curve.ainvs
        A = self.pulled_q[0].subs({X:x, Z:1}, simultaneous=True)
        B = self.norm_form.subs({X:x, Z:1}, simultaneous=True)
        return A/y**2, B/y**3-(a1*A/y**2+a3)/2

    @cached_property
    def to_working(self):
        """Recover the unique signed (u,r,s,t) from the actual map functions.

        The canonical and raw maps both send the four branch points to O.
        Solve x_can=u^2*x_work+r and y_can=u^3*v_work. The ratio u^3/u^2
        fixes the sign without choosing a square root or dropping a translation.
        Every coefficient and both full map functions are then checked.
        """
        C = self.quartic
        px, py, _ = C.map_weighted
        F, A = C.binary_form, self.pulled_q[0]
        monomials = [X**i*Z**(4-i) for i in range(5)]
        polynomials = [sp.Poly(q, X, Z, domain=sp.QQ) for q in (A, F, px/4)]
        M = sp.Matrix([[p.coeff_monomial(m) for p in polynomials[:2]] for m in monomials])
        rhs = sp.Matrix([polynomials[2].coeff_monomial(m) for m in monomials])
        try:
            solution, parameters = M.gauss_jordan_solve(rhs)
        except ValueError as error:
            raise ArithmeticError('the quartic map has no compatible x-coordinate marking') from error
        if parameters.rows:
            raise ArithmeticError('the quartic marking is not uniquely determined')
        square, r = map(_rat, solution)
        N = sp.Poly(self.norm_form, X, Z, domain=sp.QQ)
        K = sp.Poly(py/8, X, Z, domain=sp.QQ)
        m, coefficient = N.terms()[0]
        cube = _rat(K.coeff_monomial(m)/coefficient)
        if not square or not cube or sp.expand(K.as_expr()-cube*N.as_expr()) != 0:
            raise ArithmeticError('the quartic map has an inconsistent signed norm coordinate')
        u = cube/square
        if u*u != square or u**3 != cube:
            raise ArithmeticError('marking scale is not rational')
        a1, _, a3, _, _ = self.auxiliary.curve.ainvs
        change = WeierstrassChange(C.ainvs, (u, r, a1*u/2, a3*u**3/2))
        if change.target != self.auxiliary.curve or not change.verify():
            raise ArithmeticError('recovered marking has the wrong Weierstrass target')
        return change

    @property
    def marking(self):
        return MarkedQuarticCover(self.quartic, self.to_working)

    def verify(self):
        if (not self.auxiliary.verify() or not self.parametrization.verify(self.auxiliary.q[2])
                or self.pulled_q[2] != 0 or self.quartic.binary_form != -self.pulled_q[1]):
            return False
        # Coefficients in theta of delta*p(theta)^2 = Q0-theta*F are checked
        # by the auxiliary construction and the substituted Q0,Q1,Q2 identities.
        return self.marking.verify_supplied_map(self.map_working_xy)

    def to_auxiliary(self, point):
        """A morphism in weighted coordinates: (X:Y:Z) -> (p0,p1,p2,Y)."""
        if len(point) != 3:
            raise ValueError('expected three weighted quartic coordinates')
        a, b, c = map(_rat, point)
        if not (a or c) or b*b != self.quartic.binary_form.subs({X:a, Z:c}, simultaneous=True):
            raise ValueError('point is not on this quartic')
        return tuple(p.subs({X:a, Z:c}, simultaneous=True) for p in self.parametrization.forms)+(b,)

    def from_auxiliary(self, point):
        """Inverse morphism, using both Veronese charts and scaling Y correctly."""
        if len(point) != 4:
            raise ValueError('expected four auxiliary-cover coordinates')
        p = tuple(_rat(c) for c in point)
        self.auxiliary.evaluate(p)  # Validate all equations and projective nonzero.
        a, c, scale = self.parametrization.inverse_coordinates(p[:3])
        result = (a, scale*p[3], c)
        self.to_auxiliary(result)  # Check the quartic equation after chart recovery.
        return result


def quartic_presentation(auxiliary, *, conic_point=None, parameter_change=None):
    """Parametrize Q2=0; a complete conic solver is used unless a point is supplied.

    A leading-coefficient zero is repaired by an invertible GL2 change. Among
    five integer slopes at least one is not a root of the nonzero binary quartic.
    Nothing is truncated, and a failed point computation raises rather than
    omitting the class. No local-solubility assertion follows from this function.
    """
    if not isinstance(auxiliary, CubicCover):
        raise TypeError('expected a CubicCover')
    point = rational_conic_point(auxiliary.q[2]) if conic_point is None else conic_point
    param = parametrize_conic(auxiliary.q[2], point)
    if parameter_change is not None:
        param = param.reparameterize(parameter_change)
    sub = dict(zip(U[:3], param.forms))
    form = sp.expand(-auxiliary.q[1].subs(sub, simultaneous=True))
    if form.subs({X:1, Z:0}, simultaneous=True) == 0:
        slope = next((k for k in range(5) if form.subs({X:k, Z:1}, simultaneous=True) != 0), None)
        if slope is None:
            raise ArithmeticError('the branch quartic is identically zero')
        param = param.reparameterize(((slope, 1), (1, 0)))
    result = QuarticPresentation(auxiliary, param)
    if not result.verify():
        raise ArithmeticError('quartic presentation failed exact map verification')
    return result


@dataclass(frozen=True)
class LabelledQuarticCover:
    """A quartic with explicit class provenance, basis and original target.

    Coordinates belong to the stored Selmer basis, while ambient_label and
    delta bind the presentation to its actual arithmetic class. Certification
    inherits the parent's arithmetic assumptions; local evidence is not a new
    independent class-group proof.
    """
    selmer: object
    coordinates: int
    ambient_label: int
    presentation: QuarticPresentation

    @property
    def quartic(self):
        return self.presentation.quartic

    @property
    def delta(self):
        return self.presentation.auxiliary.delta

    @property
    def target(self):
        return self.selmer.curves.original

    @cached_property
    def marking(self):
        return self.presentation.marking.retarget(self.selmer.curves.to_original)

    @property
    def assumptions(self):
        return self.selmer.assumptions

    @property
    def certified(self):
        return self.selmer.certified and self.verify()

    def verify(self, *, max_nodes=1_000_000):
        try:
            S, D = self.selmer, self.presentation.auxiliary
            if (not S.complete or type(self.coordinates) is not int or type(self.ambient_label) is not int
                    or S.ambient_label(self.coordinates) != self.ambient_label
                    or D.curve != S.curves.working
                    or S.arithmetic.coordinates(D.delta) != self.ambient_label):
                return False
            decisions = [d for d in S.decisions if d.label == self.ambient_label]
            if (len(decisions) != 1 or decisions[0].model != D
                    or decisions[0].local.status != LocalStatus.SOLUBLE
                    or not decisions[0].local.verify(D, max_nodes=max_nodes)):
                return False
            return self.presentation.verify() and self.marking.verify() and self.marking.target == self.target
        except (TypeError, ValueError, ArithmeticError):
            return False

    def evaluate(self, point):
        return self.marking.evaluate(point)

    def descend(self, **kwargs):
        """Relative lifts above this exact marked class; NOT global 4-Selmer assembly."""
        if not self.verify():
            raise ValueError('the labelled quartic evidence is invalid')
        return LabelledRelativeLifts(self, self.marking.descend(**kwargs))


@dataclass(frozen=True)
class LabelledRelativeLifts:
    representative: LabelledQuarticCover
    marked_result: object

    def __post_init__(self):
        if self.marked_result.marking != self.representative.marking:
            raise ValueError('relative result does not use the labelled marking')

    @property
    def target(self):
        return self.representative.target

    @property
    def complete(self):
        return self.marked_result.complete

    @property
    def assumptions(self):
        return tuple(dict.fromkeys((*self.representative.assumptions,
                                    *self.marked_result.relative_result.assumptions)))

    @property
    def certified(self):
        return self.representative.certified and self.marked_result.certified and not self.assumptions

    @property
    def coverings(self):
        return self.marked_result.coverings

    @property
    def unresolved(self):
        return self.marked_result.unresolved

    def verify_evidence(self, **kwargs):
        return (self.representative.verify(max_nodes=kwargs.get('max_nodes', 1_000_000))
                and self.marked_result.verify_evidence(**kwargs))


def labelled_quartic(selmer, coordinates, *, conic_point=None, parameter_change=None):
    """Construct a binary quartic for any vector of a complete labelled 2-Selmer space."""
    if not selmer.complete:
        raise ValueError('a complete labelled 2-Selmer space is required')
    c = nonnegative_integer(coordinates, 'Selmer coordinates')
    label = selmer.ambient_label(c)  # Validates the chosen basis and range.
    decisions = [d for d in selmer.decisions if d.label == label and d.local.status == LocalStatus.SOLUBLE]
    if len(decisions) != 1:
        raise ValueError('the requested Selmer class has no unique accepted model')
    presentation = quartic_presentation(decisions[0].model, conic_point=conic_point,
                                        parameter_change=parameter_change)
    result = LabelledQuarticCover(selmer, c, label, presentation)
    if not result.verify():
        raise ArithmeticError('labelled quartic failed its arithmetic, local or marking checks')
    return result
