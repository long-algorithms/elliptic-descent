"""Choose exact fibre presentations without introducing avoidable support primes.

For split cubic algebras and a known rational 2-torsion Kummer class, use its
actual auxiliary rational point and CRT coordinate directions. The resulting
binary quartic is biquadratic. This changes only a verified presentation, not
the arithmetic class, marking, local decisions, or candidate-support theorem.
"""
import sympy as sp

from .arithmetic import quartic_support
from .geometry import _sqrt_rational


def split_torsion_parameters(auxiliary, point):
    """Return a conic point and GL2 change using labelled split coordinates.

    At the selected root the corresponding auxiliary coordinate is zero.
    The other two coordinates are rational square roots of (x(T)-theta)/delta.
    Choose one conic direction along that zero component, and one along another
    component. Diagonality then makes the branch quartic biquadratic.
    """
    E, A = auxiliary.curve, auxiliary.algebra
    point = E.point(point)
    if point is None or E.double(point) is not None or A.degrees != (1, 1, 1):
        raise ValueError('a nonzero rational 2-torsion point and a split cubic are required')
    parts = A.to_components(auxiliary.delta)
    root_index = next(i for i, g in enumerate(A.factors) if -g[0] == point[0])
    coordinates = tuple((sp.S.Zero,) if i == root_index else
                        (-_sqrt_rational((point[0]+g[0])/part[0]),)
                        for i, (g, part) in enumerate(zip(A.factors, parts)))
    P = sp.Matrix(A.from_components(coordinates))
    if auxiliary.evaluate((*P, sp.S.One)) != point:
        raise ArithmeticError('torsion Kummer point failed the auxiliary covering equations')
    pivot = next(i for i, c in enumerate(P) if c)
    old = sp.Matrix.hstack(P, *(sp.eye(3)[:, i] for i in range(3) if i != pivot))
    other = next(i for i in range(3) if i != root_index)
    directions = sp.Matrix.hstack(sp.Matrix(A.idempotents[root_index]),
                                  sp.Matrix(A.idempotents[other]))
    change = (old.inv()*directions)[1:, :]
    if not change.det():
        raise ArithmeticError('CRT conic directions are not independent modulo the point')
    return tuple(P), tuple(map(tuple, change.tolist()))


def _score(cover):
    f = cover.quartic.f
    support = quartic_support(f)
    bits = max(max(abs(int(c.p)).bit_length(), int(c.q).bit_length()) for c in f)
    return len(support), bits, support


def _prefer(candidate_score, current_score):
    """Only fewer support primes justify replacing the reference presentation."""
    return candidate_score[0] < current_score[0]


def fibre_representative(selmer, coordinates, torsion, *, initial_cover=None):
    """Select among fully verified presentations; never remove a class or a prime.

    The full conservative support is recomputed for each candidate presentation.
    Comparing its size only selects which isomorphic quartic to descend. No
    prime is omitted from the eventual relative-descent support. Equal-size
    support is NOT a reason to switch: smaller coefficients alone can hide
    much worse local reduction and cause an enormous lifting tree.
    """
    best = selmer.quartic_cover(coordinates) if initial_cover is None else initial_cover
    if initial_cover is not None and (best.selmer is not selmer or best.coordinates != coordinates
            or best.ambient_label != selmer.ambient_label(coordinates) or not best.verify()):
        raise ValueError('initial quartic must represent this exact Selmer class and basis')
    if not coordinates or selmer.arithmetic.algebra.degrees != (1, 1, 1):
        return best
    entries = tuple(e for e in torsion.entries if e.point is not None and e.coordinates == coordinates)
    if not entries:
        return best
    score = _score(best)
    for entry in entries:
        T = selmer.curves.to_working.apply(entry.point)
        point, change = split_torsion_parameters(best.presentation.auxiliary, T)
        candidate = selmer.quartic_cover(coordinates, conic_point=point, parameter_change=change)
        if candidate.ambient_label != best.ambient_label or candidate.target != best.target:
            raise ArithmeticError('presentation selection changed the labelled covering')
        candidate_score = _score(candidate)
        if _prefer(candidate_score, score):
            best, score = candidate, candidate_score
    return best
