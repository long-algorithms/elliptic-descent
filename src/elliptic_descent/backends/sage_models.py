"""Sage model normalization with explicit maps back to the user's rational curve.

Sage chooses a minimal model/isomorphism; our exact model and transformation
code independently check its coefficients and defining polynomial identity.
No Selmer completeness is implied by normalization.
"""
from dataclasses import dataclass

from .sage import _sage
from ..weierstrass import WeierstrassModel, WeierstrassChange
from ..cubic_conversion import IdentityTwoCover


@dataclass(frozen=True)
class RationalCurveContext:
    original: WeierstrassModel
    to_working: WeierstrassChange

    def __post_init__(self):
        if self.to_working.source != self.original or not self.to_working.verify():
            raise ValueError('normalization does not start at the original curve')

    @property
    def working(self):
        return self.to_working.target

    @property
    def to_original(self):
        return self.to_working.inverse()

    def identity_cover(self):
        """Quartic zero-class presentation on the working model, marked on original."""
        return IdentityTwoCover(self.working).cover.retarget(self.to_original)


def rational_curve_context(curve, *, minimal=True):
    """Accept rational ainvs, WeierstrassModel, or a Sage E/QQ.

    No nonrational base field or curve label is coerced to QQ implicitly.
    minimal=False retains the exact input presentation and identity map.
    """
    if type(minimal) is not bool:
        raise TypeError('minimal must be a boolean')
    sage = _sage()
    if hasattr(curve,'a_invariants'):
        if curve.base_ring() != sage.QQ:
            raise ValueError('expected an elliptic curve over QQ')
        original = WeierstrassModel(tuple(curve.a_invariants()))
    else:
        original = curve if isinstance(curve,WeierstrassModel) else WeierstrassModel(curve)
    if not minimal:
        return RationalCurveContext(original,WeierstrassChange(original))
    E = sage.EllipticCurve([sage.QQ(str(a)) for a in original.ainvs])
    F = E.global_minimal_model()
    iso = E.isomorphism_to(F)
    change = WeierstrassChange(original,tuple(iso.tuple()))
    if change.target != WeierstrassModel(tuple(F.a_invariants())) or not change.verify():
        raise ArithmeticError('Sage normalization failed independent exact verification')
    return RationalCurveContext(original,change)


def mark_pari_two_cover(record, change, *, translation=None, shift=None):
    """Verify an imported PARI marking and return (marked_cover, conversion).

    No class label, basis completeness, or isomorphism/translation is inferred.
    The caller supplies the explicit Jacobian change and translation. A cubic
    record is converted first; its conversion object is returned for point maps.
    """
    from ..marked_cover import MarkedQuarticCover
    from ..cubic_conversion import CubicToQuartic
    if not isinstance(change,WeierstrassChange) or change.target != WeierstrassModel(record.ainvs):
        raise ValueError('marking target differs from the imported original curve')
    if len(record.f) == 4:
        conversion = CubicToQuartic(record.f,shift)
        return conversion.marked_map(record.map_xy,change,translation),conversion
    if len(record.f) != 5 or shift is not None:
        raise ValueError('expected cubic or quartic record; shift applies only to cubics')
    result = MarkedQuarticCover(record.f,change,translation)
    if not result.verify_supplied_map(record.map_xy):
        raise ValueError('imported map differs from the declared marked covering')
    return result,None
