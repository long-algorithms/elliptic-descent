"""Certified labelled 2-Selmer arithmetic over Q, using the cubic norm kernel.

This reference algorithm exhausts the finite norm kernel and screens explicit
auxiliary models locally. It is not a wrapper assigning guessed labels to a
PARI basis, nor a full 4-Selmer computation. See docs/two-selmer.md.
"""
from dataclasses import dataclass, replace
from itertools import islice
from time import perf_counter

import sympy as sp

from .accounting import valid_accounting
from .arithmetic import nonnegative_integer
from .backends.sage_models import rational_curve_context
from .backends.sage_two_selmer import cubic_square_classes
from .cubic_arithmetic import (T, _poly, kummer_algebra, kummer_value,
                               rational_two_torsion, rational_halves_of_two_torsion)
from .cubic_local import CubicCover, cubic_local_solubility
from .quotient import row_basis
from .status import LocalStatus, ProofStatus
from .weierstrass import WeierstrassChange


@dataclass(frozen=True)
class F2Labels:
    """Independent ordered ambient bit vectors; coordinates never cross bases silently."""
    ambient_dimension: int
    basis: tuple

    def __post_init__(self):
        n = nonnegative_integer(self.ambient_dimension, 'ambient dimension')
        basis = tuple(nonnegative_integer(v, 'basis label') for v in self.basis)
        if any(v >= 1 << n for v in basis) or len(row_basis(basis)) != len(basis):
            raise ValueError('basis is dependent or outside the ambient space')
        object.__setattr__(self, 'ambient_dimension', n)
        object.__setattr__(self, 'basis', basis)

    @property
    def dimension(self):
        return len(self.basis)

    def ambient(self, coordinates):
        c = nonnegative_integer(coordinates, 'coordinates')
        if c >= 1 << self.dimension:
            raise ValueError('coordinates exceed this basis dimension')
        value = 0
        for i, b in enumerate(self.basis):
            if c >> i & 1:
                value ^= b
        return value

    def coordinates(self, label):
        value = nonnegative_integer(label, 'ambient label')
        if value >= 1 << self.ambient_dimension:
            raise ValueError('label exceeds ambient dimension')
        rows = {}
        for i, b in enumerate(self.basis):
            v, c = b, 1 << i
            while v:
                pivot = v.bit_length()-1
                if pivot not in rows:
                    rows[pivot] = v, c
                    break
                r, tag = rows[pivot]
                v, c = v ^ r, c ^ tag
        result = 0
        while value:
            pivot = value.bit_length()-1
            if pivot not in rows:
                raise ValueError('class does not belong to this subspace')
            v, c = rows[pivot]
            value, result = value ^ v, result ^ c
        return result


@dataclass(frozen=True)
class TwoSelmerDecision:
    label: int
    model: CubicCover
    local: object


@dataclass(frozen=True)
class TorsionKummerEntry:
    point: object
    coordinates: int
    rational_halves: tuple


@dataclass(frozen=True)
class TorsionKummerData:
    entries: tuple
    image_basis: tuple

    @property
    def image_rank(self):
        return len(self.image_basis)

    @property
    def kernel(self):
        """E(Q)[2] intersect 2E(Q), equivalently 2E(Q)[4], as actual points."""
        return tuple(e.point for e in self.entries if e.coordinates == 0)

    @property
    def four_torsion(self):
        points = []
        for e in self.entries:
            for P in e.rational_halves:
                if P not in points:
                    points.append(P)
        return tuple(points)

    def verify(self, selmer):
        E = selmer.curves.original
        if tuple(e.point for e in self.entries) != rational_two_torsion(E):
            return False
        for e in self.entries:
            if selmer.kummer(e.point) != e.coordinates:
                return False
            if rational_halves_of_two_torsion(E, e.point) != e.rational_halves:
                return False
            if bool(e.rational_halves) != (e.coordinates == 0):
                return False
        return self.image_basis == row_basis(e.coordinates for e in self.entries)


@dataclass(frozen=True)
class TwoSelmerResult:
    curves: object
    arithmetic: object
    decisions: tuple
    timings: tuple
    selected_basis: tuple | None = None

    @property
    def candidate_space(self):
        return self.arithmetic.norm_kernel

    @property
    def total_candidates(self):
        return self.candidate_space.cardinality

    @property
    def enumeration_complete(self):
        K = self.candidate_space
        return (len(self.decisions) == K.cardinality and valid_accounting(
            K, (d.label for d in self.decisions), K.cardinality, (), len(K.basis)))

    @property
    def unresolved(self):
        return tuple(d for d in self.decisions if d.local.status == LocalStatus.UNKNOWN)

    @property
    def known_soluble_labels(self):
        return tuple(d.label for d in self.decisions if d.local.status == LocalStatus.SOLUBLE)

    @property
    def complete(self):
        if not self.enumeration_complete or any(d.local.status not in
                (LocalStatus.SOLUBLE, LocalStatus.INSOLUBLE) for d in self.decisions):
            return False
        # Never take the span and silently promote an incomplete/nonclosed set.
        labels = self.known_soluble_labels
        return len(labels) == 1 << len(row_basis(labels))

    @property
    def assumptions(self):
        return self.arithmetic.assumptions

    @property
    def certified(self):
        if not self.complete or self.arithmetic.arithmetic_proof != ProofStatus.CERTIFIED or self.assumptions:
            return False
        try:
            _ = self.labelling
        except (TypeError, ValueError):
            return False
        return True

    @property
    def labelling(self):
        if not self.complete:
            raise ValueError('partial calculation has no complete Selmer basis')
        canonical = row_basis(self.known_soluble_labels)
        chosen = canonical if self.selected_basis is None else self.selected_basis
        labels = F2Labels(self.arithmetic.dimension, chosen)
        if row_basis(labels.basis) != canonical:
            raise ValueError('selected basis does not span the computed Selmer space')
        return labels

    @property
    def basis_labels(self):
        return self.labelling.basis

    @property
    def dimension(self):
        return self.labelling.dimension

    @property
    def cardinality(self):
        return 1 << self.dimension

    def with_basis(self, ambient_labels):
        result = replace(self, selected_basis=tuple(ambient_labels))
        _ = result.labelling
        return result

    def ambient_label(self, coordinates):
        return self.labelling.ambient(coordinates)

    def coordinates_of_label(self, ambient_label):
        return self.labelling.coordinates(ambient_label)

    def coordinates(self, algebra_element):
        return self.coordinates_of_label(self.arithmetic.coordinates(algebra_element))

    def representative(self, coordinates):
        """Cubic-algebra representative, NOT a binary-quartic covering model."""
        return self.arithmetic.representative(self.ambient_label(coordinates))

    def add(self, left, right):
        labels = self.labelling
        labels.ambient(left)
        labels.ambient(right)
        return int(left) ^ int(right)

    def kummer_ambient(self, point):
        """Ambient class for a validated rational point on the ORIGINAL curve."""
        P = self.curves.to_working.apply(point)
        value = kummer_value(self.curves.working, P)
        label = self.arithmetic.coordinates(value)
        if self.arithmetic.norm_label(label):
            raise ArithmeticError('Kummer image failed the norm condition')
        return label

    def kummer(self, point):
        return self.coordinates_of_label(self.kummer_ambient(point))

    def torsion_data(self):
        E = self.curves.original
        entries = tuple(TorsionKummerEntry(P, self.kummer(P), rational_halves_of_two_torsion(E, P))
                        for P in rational_two_torsion(E))
        data = TorsionKummerData(entries, row_basis(e.coordinates for e in entries))
        if not data.verify(self):
            raise ArithmeticError('torsion Kummer/halving accounting failed')
        return data

    def transport_to(self, other, coordinates, change=None):
        """Natural class transport under an explicit original-curve isomorphism.

        Without change the original models must agree. Handles different
        working models, Sage bases and selected Selmer bases without assuming
        identical integer labels. Rational translations are NOT model changes.
        """
        if change is None:
            if self.curves.original != other.curves.original:
                raise ValueError('an explicit original-model isomorphism is required')
            change = WeierstrassChange(self.curves.original)
        if (not isinstance(change, WeierstrassChange) or change.source != self.curves.original
                or change.target != other.curves.original or not change.verify()):
            raise ValueError('wrong model-change endpoints')
        total = self.curves.to_original.then(change).then(other.curves.to_working)
        u, r, _, _ = total.urst
        f1, f2 = self.arithmetic.algebra.f, other.arithmetic.algebra.f
        if sp.expand(sum(c*(u*u*T+r)**i for i, c in enumerate(f1))
                     -u**6*sum(c*T**i for i, c in enumerate(f2))) != 0:
            raise ArithmeticError('completed-square algebra transport failed')
        value = self.representative(coordinates)
        p = sp.Poly(sum(c*(u*u*T+r)**i for i, c in enumerate(value)), T, domain=sp.QQ).rem(_poly(f2))
        return other.coordinates(tuple(p.nth(i) for i in range(3)))

    def verify_evidence(self, *, max_nodes=1_000_000):
        """Audit labels, models and local evidence, not a second class-group proof."""
        if (not self.complete or not self.curves.to_working.verify()
                or self.arithmetic.algebra != kummer_algebra(self.curves.working)
                or not self.arithmetic.verify_basis()):
            return False
        try:
            _ = self.labelling
            for d in self.decisions:
                if (d.model.curve != self.curves.working
                        or self.arithmetic.coordinates(d.model.delta) != d.label
                        or not d.model.verify()
                        or not d.local.verify(d.model, max_nodes=max_nodes)):
                    return False
        except (TypeError, ValueError):
            return False
        return True


def two_selmer(curve, *, proof=True, minimal=True, max_candidates=None,
               max_depth=None, max_nodes=None):
    """Full labelled Sel^(2)(E/Q), or explicitly partial decisions with limits.

    Reference path: certified cubic factor spaces, complete norm-kernel
    enumeration, then exact local tests. Exponential enumeration is deliberate
    and audited; no practical complexity claim or full 4-Selmer claim is made.
    """
    if type(proof) is not bool or type(minimal) is not bool:
        raise TypeError('proof and minimal must be booleans')
    limits = {}
    for name, value in (('max_candidates', max_candidates), ('max_depth', max_depth), ('max_nodes', max_nodes)):
        limits[name] = None if value is None else nonnegative_integer(value, name)
    start = perf_counter()
    curves = rational_curve_context(curve, minimal=minimal)
    arithmetic = cubic_square_classes(kummer_algebra(curves.working), proof=proof)
    arithmetic_seconds = perf_counter()-start
    selected = arithmetic.norm_kernel.elements()
    if limits['max_candidates'] is not None:
        selected = islice(selected, limits['max_candidates'])
    decisions = []
    start = perf_counter()
    for label in selected:
        delta = arithmetic.representative(label)
        model = CubicCover(curves.working, delta)
        if arithmetic.coordinates(delta) != label or not model.verify():
            raise ArithmeticError('labelled cubic construction verification failed')
        local = cubic_local_solubility(model, max_depth=limits['max_depth'], max_nodes=limits['max_nodes'])
        decisions.append(TwoSelmerDecision(label, model, local))
    result = TwoSelmerResult(curves, arithmetic, tuple(decisions),
         (('arithmetic_seconds', arithmetic_seconds), ('construction_local_seconds', perf_counter()-start)))
    if result.enumeration_complete and not result.unresolved and not result.complete:
        raise ArithmeticError('decided classes do not form a Selmer subspace')
    return result
