"""Relative 4-descent for a smooth irreducible quartic over Q.

The target is a fixed quartic 2-cover, with its canonical Jacobian map.
Equivalence fixes the target C pointwise. Both orientations of D -> C are
retained. The un-oriented model quotient is available separately, for
comparison with Magma. See docs/relative-descent.md for the proof obligations.
"""
from dataclasses import dataclass, replace
from itertools import islice
from time import perf_counter

import sympy as sp

from .arithmetic import nonnegative_integer
from .backends.sage import sage_square_classes
from .geometry import _mul, _norm, _sqrt_rational, construct_four_cover
from .local_complete import complete_local_screen
from .quotient import quotient_fibre
from .status import LocalStatus, ProofStatus
from .two_cover import QuarticTwoCover


@dataclass(frozen=True)
class RelativeCandidate:
    label: int
    model: object
    local: object


@dataclass(frozen=True)
class RelativeDescentResult:
    source: QuarticTwoCover
    candidates: tuple
    total_model_classes: int
    scalar_basis: tuple
    norm_dimension: int
    arithmetic_proof: ProofStatus
    assumptions: tuple
    support: tuple
    timings: tuple

    @property
    def enumeration_complete(self):
        return len(self.candidates) == self.total_model_classes

    @property
    def unresolved(self):
        return tuple(c for c in self.candidates if c.local.status == LocalStatus.UNKNOWN)

    @property
    def complete(self):
        return self.enumeration_complete and not self.unresolved

    @property
    def certified(self):
        return self.complete and self.arithmetic_proof == ProofStatus.CERTIFIED and not self.assumptions

    @property
    def unoriented_coverings(self):
        """Models modulo the target involution Y -> -Y; an explicit quotient."""
        return tuple(c.model for c in self.candidates if c.local.status == LocalStatus.SOLUBLE)

    @property
    def coverings(self):
        """All maps to the fixed C, including both inequivalent orientations.

        Irreducibility implies Q[t]/f is a field. An element v with v^2 in Q
        has positive norm, so no scalar-square stabilizer reverses orientation.
        """
        return tuple(oriented for model in self.unoriented_coverings
                     for oriented in (model, replace(model, norm_witness=-model.norm_witness)))

    @property
    def scope(self):
        return "relative 2-coverings of the fixed quartic; no target automorphisms quotiented"

    def verify_evidence(self, *, max_nodes=1_000_000, geometry=True, full_geometry=False):
        """Check model/map/local evidence, not independently redo class groups.

        Arithmetic completeness is a trusted certified Sage dependency. False
        can mean the independent local verifier's resource budget was too small.
        """
        if not self.complete or not self.source.check_map():
            return False
        for candidate in self.candidates:
            D = candidate.model
            if D.f != self.source.f or not candidate.local.verify(D, max_nodes=max_nodes):
                return False
            if geometry and not D.check_construction():
                return False
            if full_geometry and not all((D.check_map_identity(), D.check_smoothness(),
                                          D.check_map_basepoint_free())):
                return False
        return True


def relative_four_descent(f, *, proof=True, max_candidates=None,
                         max_depth=None, max_nodes=None, verify_geometry=True):
    """Compute every locally soluble relative 4-covering in the supported domain.

    f: five rational coefficients, in ascending order, or QuarticTwoCover.
    No artificial limits are applied by default. Finite limits give partial
    results and NEVER turn unresolved candidates into local obstructions.
    Reducible/singular quartics are rejected rather than labelled complete.

    This does not compute Sel^(4)(E) or a Mordell--Weil rank. The fixed
    source may itself be locally insoluble; an empty result alone makes no
    Sha claim without separate source local-solubility information.
    """
    if type(proof) is not bool or type(verify_geometry) is not bool:
        raise TypeError("proof and verify_geometry must be booleans")
    for name, value in (("max_candidates", max_candidates), ("max_depth", max_depth),
                        ("max_nodes", max_nodes)):
        if value is not None:
            nonnegative_integer(value, name)
    source = f if isinstance(f, QuarticTwoCover) else QuarticTwoCover(tuple(f))
    if not source.check_map():
        raise ArithmeticError("canonical quartic map verification failed")
    start = perf_counter()
    context = sage_square_classes(source.f, proof=proof)
    space = context.space
    fibre = space.norm_fibre()
    # Include the sign class: omitting -1 would retain duplicate models.
    scalars = tuple(context.coordinates((c, 0, 0, 0)) for c in (-1, *space.support))
    quotient, scalar_basis = quotient_fibre(fibre, scalars)
    arithmetic_seconds = perf_counter()-start
    candidates = []
    selected = quotient.elements()
    if max_candidates is not None:
        selected = islice(selected, int(max_candidates))
    start = perf_counter()
    for label in selected:
        xi = (sp.S.One, sp.S.Zero, sp.S.Zero, sp.S.Zero)
        for i, generator in enumerate(space.generators):
            if label >> i & 1:
                xi = _mul(xi, generator, space.f)
        witness = _sqrt_rational(_norm(xi, space.f)/space.f[4])
        if context.coordinates(xi) != label:
            raise ArithmeticError("number-field coordinate verification failed")
        model = construct_four_cover(space.f, xi, norm_witness=witness)
        if verify_geometry and not model.check_construction():
            raise ArithmeticError("covering geometry verification failed")
        local = complete_local_screen(model, max_depth=max_depth, max_nodes=max_nodes)
        candidates.append(RelativeCandidate(label, model, local))
    local_seconds = perf_counter()-start
    return RelativeDescentResult(source, tuple(candidates), quotient.cardinality,
                                 scalar_basis, len(fibre.basis), space.arithmetic_proof,
                                 space.assumptions, space.support,
                                 (("arithmetic_seconds", arithmetic_seconds),
                                  ("construction_and_local_seconds", local_seconds)))
