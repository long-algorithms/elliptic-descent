"""Labelled Cassels--Tate pairing on Sel^(2)(E/Q), via Fisher's quartic formula.

This computes the pairing and its radical, not the full 4-Selmer group. The
matrix uses the exact basis recorded by its parent TwoSelmerResult. Missing
local evaluations are unknown; a partial matrix has no reported radical.
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import count

import sympy as sp

from .arithmetic import affine_f2, nonnegative_integer
from .backends.sage_pairing import cubic_product_square_root
from .cassels_tate_forms import (FisherQuartic, FisherTriple, change_parameters,
                                 substitute_generator, invariants)
from .cassels_tate_local import evaluate_fisher


def working_invariants(curve):
    b2, b4, b6, _ = curve.b_invariants
    return (b2*b2-24*b4)/16, (-b2**3+36*b2*b4-216*b6)/32


@dataclass(frozen=True)
class NormalizedPairingQuartic:
    """A labelled quartic with common invariants and unit cubic invariant.

    Parameter changes act on X,Z only; Y is unchanged. They preserve the
    covering class. Normalization uses the retained signed Jacobian scale,
    not a fourth/sixth root inferred from j-invariants.
    """
    labelled_cover: object
    parameter_matrix: tuple

    @cached_property
    def quartic(self):
        return FisherQuartic(change_parameters(self.labelled_cover.quartic.f, self.parameter_matrix))

    @property
    def theta_invariant(self):
        S = self.labelled_cover.selmer
        b2 = S.curves.working.b_invariants[0]
        # phi=-3*theta-b2/4 identifies phi^3-3I*phi+J with the
        # monic completed-square cubic. It is not an unnamed field isomorphism.
        return substitute_generator(self.quartic.z, -3, -b2/4)

    def verify(self):
        C = self.labelled_cover
        try:
            M = sp.Matrix(self.parameter_matrix)
            u = C.presentation.to_working.urst[0]
            if M.shape != (2, 2) or M.det()**2 != 1/u**2:
                return False
            Q = self.quartic
            return (C.verify() and Q.invariants == working_invariants(C.selmer.curves.working)
                    and bool(Q.algebra.norm(Q.z)) and Q.verify_identity()
                    and C.selmer.arithmetic.coordinates(self.theta_invariant) == C.ambient_label)
        except (TypeError, ValueError, ArithmeticError):
            return False


def normalize_pairing_quartic(cover):
    if not cover.verify():
        raise ValueError('expected a verified labelled quartic')
    u = cover.presentation.to_working.urst[0]
    base = sp.diag(1/u, 1)
    # The sextic covariant has finitely many zero directions. Try the existing
    # infinity first, then determinant-one changes with infinity at (k:1).
    for k in count(-1):
        M = base if k == -1 else base*sp.Matrix([[k, -1], [1, 0]])
        f = change_parameters(cover.quartic.f, M)
        if not f[4]:
            continue
        Q = FisherQuartic(f)
        if not Q.algebra.norm(Q.z):
            continue
        result = NormalizedPairingQuartic(cover, tuple(map(tuple, M.tolist())))
        if not result.verify():
            raise ArithmeticError('quartic normalization or arithmetic-class identification failed')
        return result


@dataclass(frozen=True)
class LabelledPairingEntry:
    selmer: object
    left: int
    right: int
    normalized: tuple
    evaluation: object = None

    @property
    def is_trivial(self):
        return self.left == 0 or self.right == 0 or self.left == self.right

    @property
    def complete(self):
        if self.is_trivial:
            return self.normalized == () and self.evaluation is None
        return self.evaluation is not None and self.evaluation.complete

    @property
    def bit(self):
        if not self.complete:
            raise ValueError('this pairing entry has unresolved local evaluations')
        return 0 if self.is_trivial else self.evaluation.bit

    @property
    def assumptions(self):
        return self.selmer.assumptions

    @property
    def certified(self):
        return self.selmer.certified and self.complete

    def verify(self):
        S = self.selmer
        if not S.complete or type(self.left) is not int or type(self.right) is not int:
            return False
        try:
            S.ambient_label(self.left)
            S.ambient_label(self.right)
            if self.is_trivial:
                return self.complete  # Alternation and pairing with zero are theorems.
            if len(self.normalized) != 3 or not self.evaluation.verify():
                return False
            for N, c, Q in zip(self.normalized, (self.left, self.right, S.add(self.left, self.right)),
                               self.evaluation.triple.quartics):
                C = N.labelled_cover
                if (C.selmer is not S or C.ambient_label != S.ambient_label(c)
                        or N.quartic != Q or not N.verify()):
                    return False
            return True
        except (TypeError, ValueError, ArithmeticError):
            return False


class _PairingBuilder:
    """One computation's presentation cache; no cache shared across Selmer bases."""
    def __init__(self, selmer):
        if not selmer.complete:
            raise ValueError('a complete labelled 2-Selmer space is required')
        _ = selmer.labelling
        self.selmer, self.normalized = selmer, {}

    def cover(self, coordinates):
        if coordinates not in self.normalized:
            self.normalized[coordinates] = normalize_pairing_quartic(self.selmer.quartic_cover(coordinates))
        return self.normalized[coordinates]

    def entry(self, left, right, *, max_local_trials=None):
        S = self.selmer
        left = nonnegative_integer(left, 'left coordinates')
        right = nonnegative_integer(right, 'right coordinates')
        S.ambient_label(left)
        S.ambient_label(right)
        if not left or not right or left == right:
            return LabelledPairingEntry(S, left, right, ())
        N = tuple(self.cover(c) for c in (left, right, S.add(left, right)))
        A = S.arithmetic.algebra
        z1, z2, z3 = (n.theta_invariant for n in N)
        product = A.multiply(A.multiply(z1, z2), z3)
        mtheta = cubic_product_square_root(S.arithmetic, product)
        b2 = S.curves.working.b_invariants[0]
        mphi = substitute_generator(mtheta, -sp.Rational(1, 3), -b2/12)
        triple = FisherTriple(tuple(n.quartic for n in N), mphi)
        evaluation = evaluate_fisher(triple, max_local_trials=max_local_trials)
        result = LabelledPairingEntry(S, left, right, N, evaluation)
        if result.complete and not result.verify():
            raise ArithmeticError('pairing entry failed its exact evidence checks')
        return result


def pairing_entry(selmer, left, right, *, max_local_trials=None):
    """Compute a pair directly from quartics, independently of a stored matrix.

    Bits 0 and 1 represent 0 and 1/2 in Q/Z. An explicitly bounded local
    search may yield an unresolved entry; it never supplies a guessed bit.
    """
    if max_local_trials is not None:
        nonnegative_integer(max_local_trials, 'max_local_trials')
    return _PairingBuilder(selmer).entry(left, right, max_local_trials=max_local_trials)


def _evaluate_matrix(rows, left, right):
    """Evaluate one already materialized matrix; do not rebuild it per entry."""
    return sum(((left >> i) & 1)*((right >> j) & 1)*rows[i][j]
               for i in range(len(rows)) for j in range(len(rows))) % 2


class _MatrixOperations:
    """Exact linear algebra shared by original and explicitly rebased matrices."""
    @property
    def dimension(self):
        return self.selmer.dimension

    def evaluate(self, left, right):
        self.selmer.ambient_label(left)
        self.selmer.ambient_label(right)
        return _evaluate_matrix(self.matrix, int(left), int(right))

    @property
    def radical_basis(self):
        rows = self.matrix
        n = len(rows)
        columns = tuple(sum(rows[i][j] << i for i in range(n)) for j in range(n))
        return affine_f2(columns, 0).basis

    @property
    def radical_ambient_basis(self):
        return tuple(self.selmer.ambient_label(c) for c in self.radical_basis)

    @property
    def rank(self):
        return self.dimension-len(self.radical_basis)

    def in_radical(self, coordinates):
        self.selmer.ambient_label(coordinates)
        rows = self.matrix
        return all(_evaluate_matrix(rows, int(coordinates), 1 << j) == 0 for j in range(len(rows)))

    def with_basis(self, ambient_labels):
        return RebasedCasselsTatePairing(self, self.selmer.with_basis(ambient_labels))


@dataclass(frozen=True)
class CasselsTatePairing(_MatrixOperations):
    selmer: object
    entries: tuple

    @property
    def complete(self):
        n = self.dimension
        expected = tuple((1 << i, 1 << j) for i in range(n) for j in range(i+1, n))
        return (tuple((e.left, e.right) for e in self.entries) == expected
                and all(e.selmer is self.selmer and e.complete for e in self.entries))

    @property
    def assumptions(self):
        return self.selmer.assumptions

    @property
    def certified(self):
        return self.complete and self.selmer.certified

    @property
    def matrix(self):
        if not self.complete:
            raise ValueError('an incomplete pairing matrix has no certified entries/radical')
        n = self.dimension
        rows = [[0]*n for _ in range(n)]
        for e in self.entries:
            i, j = e.left.bit_length()-1, e.right.bit_length()-1
            rows[i][j] = rows[j][i] = e.bit
        return tuple(map(tuple, rows))

    def verify_evidence(self):
        return (self.complete and self.selmer.verify_evidence()
                and all(e.verify() for e in self.entries) and self.rank % 2 == 0)

    def direct_entry(self, left, right, **kwargs):
        return pairing_entry(self.selmer, left, right, **kwargs)


@dataclass(frozen=True)
class RebasedCasselsTatePairing(_MatrixOperations):
    original: object
    selmer: object

    def __post_init__(self):
        if (self.selmer.arithmetic is not self.original.selmer.arithmetic
                or self.selmer.decisions != self.original.selmer.decisions
                or self.selmer.curves != self.original.selmer.curves):
            raise ValueError('rebasing requires the same computed Selmer space')
        _ = self.selmer.labelling

    @property
    def complete(self):
        return self.original.complete

    @property
    def assumptions(self):
        return self.original.assumptions

    @property
    def certified(self):
        return self.original.certified and self.selmer.certified

    @property
    def matrix(self):
        # Read the original matrix once, including its completeness check.
        # Calling original.evaluate repeatedly can rebuild nested rebased
        # matrices and perform the same parent-ledger audit many times.
        rows = self.original.matrix
        images = tuple(self.original.selmer.coordinates_of_label(b) for b in self.selmer.basis_labels)
        return tuple(tuple(_evaluate_matrix(rows, a, b) for b in images) for a in images)

    def verify_evidence(self):
        return self.original.verify_evidence() and self.selmer.verify_evidence()

    def direct_entry(self, left, right, **kwargs):
        return pairing_entry(self.selmer, left, right, **kwargs)


def cassels_tate_pairing(selmer, *, max_local_trials=None):
    """Compute the labelled alternating pairing and radical on full Sel^(2).

    Only the independent upper-triangle entries require the local formula.
    Diagonal and zero entries follow from alternation. The radical is the
    liftable subspace (image of Sel^(4)); this does NOT construct those lifts
    or give Sel^(4)'s group structure. An empty or unresolved local search is
    never used to decide a pairing entry.
    """
    if max_local_trials is not None:
        nonnegative_integer(max_local_trials, 'max_local_trials')
    builder = _PairingBuilder(selmer)
    entries = tuple(builder.entry(1 << i, 1 << j, max_local_trials=max_local_trials)
                    for i in range(selmer.dimension) for j in range(i+1, selmer.dimension))
    return CasselsTatePairing(selmer, entries)
