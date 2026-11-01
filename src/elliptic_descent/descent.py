"""Curve-level descent with explicit stage boundaries and per-session reuse.

All mathematics is delegated to the existing exact constructors. A session
reuses arithmetic objects and witnesses, never an evidence-verification verdict.
Its proof setting, original model and chosen Selmer basis are not interchangeable
with another session's. Resource limits retain partial ledgers; arithmetic and
backend errors still raise. See docs/public-api.md.
"""
from collections import Counter
from dataclasses import dataclass, fields
from enum import Enum
from itertools import islice
from time import perf_counter

from .arithmetic import nonnegative_integer
from .cassels_tate import (
    CasselsTatePairing, _PairingBuilder, normalize_pairing_quartic,
)
from .cubic_local import CubicCover, cubic_local_solubility
from .fibre_presentations import fibre_representative
from .four_selmer_fibres import FourSelmerFibres, PairingExclusion, _assemble_one
from .four_selmer_group import four_selmer_group
from .limits import FieldDegreeLimit
from .status import CompletionStatus, LocalStatus
from .two_selmer import TwoSelmerDecision, TwoSelmerResult, two_selmer


class DescentStage(str, Enum):
    TWO_SELMER = 'two_selmer'
    PAIRING = 'pairing'
    FIBRES = 'fibres'
    GROUP = 'group'


@dataclass(frozen=True)
class DescentLimits:
    """Independent stage budgets; None is unbounded, zero is a real limit.

    Local node/depth budgets apply per candidate/place as in the low-level
    routines, not to the pipeline as a whole. field_degree bounds the accepted
    normal field, not discovery runtime. No wall-clock guarantee is implied.
    Already completed cached computations need no new discovery budget.
    """
    two_selmer_candidates: int | None = None
    two_selmer_depth: int | None = None
    two_selmer_nodes: int | None = None
    pairing_trials: int | None = None
    fibres: int | None = None
    relative_candidates: int | None = None
    relative_depth: int | None = None
    relative_nodes: int | None = None
    field_degree: int | None = None

    def __post_init__(self):
        for f in fields(self):
            value = getattr(self, f.name)
            if value is not None:
                if isinstance(value, str):
                    raise TypeError(f'{f.name} must be an exact integer, not text')
                object.__setattr__(self, f.name, nonnegative_integer(value, f.name))

    @property
    def relative(self):
        return dict(max_candidates=self.relative_candidates, max_depth=self.relative_depth,
                    max_nodes=self.relative_nodes)


@dataclass(frozen=True)
class FourDescentResult:
    """An immutable pipeline snapshot, not a replacement for its evidence.

    finished concerns the requested stage; complete ALWAYS means that the full
    group has been constructed. Thus a verified staged 2-Selmer calculation is
    finished but is not a complete 4-descent. Models/class keys use this exact
    result's basis and marking. Parent assumptions can never be dropped.
    """
    requested_stage: DescentStage
    limits: DescentLimits
    selmer: object = None
    pairing: object = None
    fibres: object = None
    group: object = None
    stop_reason: str | None = None
    limit_detail: tuple = ()
    work: tuple = ()
    timings: tuple = ()

    def __post_init__(self):
        object.__setattr__(self, 'requested_stage', DescentStage(self.requested_stage))
        if not isinstance(self.limits, DescentLimits):
            raise TypeError('limits must be a DescentLimits instance')

    def _linked(self):
        return (self.selmer is not None
                and (self.pairing is None or self.pairing.selmer is self.selmer)
                and (self.fibres is None or (self.pairing is not None
                     and self.fibres.selmer is self.selmer and self.fibres.pairing is self.pairing))
                and (self.group is None or (self.fibres is not None and self.group.fibres is self.fibres)))

    @property
    def stage(self):
        for stage in reversed(tuple(DescentStage)):
            if self._object(stage) is not None:
                return stage
        return None

    def _object(self, stage):
        return self.selmer if stage == DescentStage.TWO_SELMER else getattr(self, stage.value)

    @property
    def curve(self):
        if self.selmer is None:
            raise ValueError('the original curve has not yet been normalized')
        return self.selmer.curves.original

    @property
    def finished(self):
        obj = self._object(self.requested_stage)
        return self._linked() and obj is not None and obj.complete and self.stop_reason is None

    @property
    def complete(self):
        return self._linked() and self.group is not None and self.group.complete

    @property
    def assumptions(self):
        return tuple(dict.fromkeys(a for obj in (self.selmer, self.pairing, self.fibres, self.group)
                                   if obj is not None for a in obj.assumptions))

    def verify_evidence(self, *, max_nodes=1_000_000):
        """Recheck the latest computed stage, including its parents.

        True for a complete intermediate/conditional stage does NOT certify a
        full unconditional group. No result of this method is memoized.
        """
        if not self._linked() or self.stage is None:
            return False
        obj = self._object(self.stage)
        if not obj.complete:
            return False
        if self.stage == DescentStage.PAIRING:
            # Pairing has its own exact local certificates and verifier API.
            return obj.verify_evidence()
        return obj.verify_evidence(max_nodes=max_nodes)

    @property
    def certified(self):
        return (self.complete and not self.assumptions and self.fibres.certified
                and self.verify_evidence())

    @property
    def completion_status(self):
        if not self.complete or not self.verify_evidence():
            return CompletionStatus.INCOMPLETE
        if self.assumptions or not self.fibres.certified:
            return CompletionStatus.CONDITIONAL
        return CompletionStatus.CERTIFIED

    def require_group(self, *, certified=False):
        if type(certified) is not bool:
            raise TypeError('certified must be a boolean')
        if not self.complete:
            raise ValueError('no complete 4-Selmer group; inspect the retained stages and stop_reason')
        if certified and not self.certified:
            raise ValueError('the complete group is not unconditionally certified')
        return self.group

    @property
    def cardinality(self):
        return self.require_group().cardinality

    @property
    def invariant_factors(self):
        return self.require_group().invariant_factors

    def covering(self, key):
        """Retrieve an actual known marked class; never launch another descent.

        A partial fibre assembly may have some fully assembled fibres. Their
        known representatives are usable, without asserting exhaustiveness.
        """
        if not self._linked() or self.fibres is None:
            raise ValueError('no global covering representatives are available')
        self.fibres.projection(key)  # Validates both key coordinates.
        return self.fibres.fibre(key[0]).classes[key[1]].covering

    def summary(self):
        """JSON-compatible metadata; certification performs a fresh evidence check."""
        return {
            'requested_stage': self.requested_stage.value,
            'stage': None if self.stage is None else self.stage.value,
            'finished': self.finished, 'complete': self.complete,
            'certified': self.certified, 'assumptions': list(self.assumptions),
            'stop_reason': self.stop_reason, 'limit_detail': dict(self.limit_detail),
            'original_ainvs': None if self.selmer is None else [str(a) for a in self.curve.ainvs],
            'selmer_dimension': self.selmer.dimension if self.selmer is not None and self.selmer.complete else None,
            'known_global_classes': len(self.fibres.known_classes) if self.fibres is not None else 0,
            'cardinality': self.cardinality if self.complete else None,
            'invariant_factors': list(self.invariant_factors) if self.complete else None,
            'work': dict(self.work), 'timings': dict(self.timings),
        }


class _SessionPairingBuilder(_PairingBuilder):
    def __init__(self, session):
        super().__init__(session._selmer)
        self.session = session

    def cover(self, coordinates):
        if coordinates not in self.normalized:
            self.normalized[coordinates] = normalize_pairing_quartic(self.session.quartic_cover(coordinates))
            self.session._work['quartic_normalizations'] += 1
        return self.normalized[coordinates]


class DescentSession:
    """Synchronous reusable arithmetic for one original model and proof policy.

    Reuse is local to this session, not a global cache or persisted proof store.
    run() can be called again with larger limits. Completed candidate decisions,
    pairing entries and fibres are reused; unresolved attempts are retried when
    their budgets change. This mutable workspace is not thread-safe.
    """
    def __init__(self, curve, *, proof=True, minimal=True):
        if type(proof) is not bool or type(minimal) is not bool:
            raise TypeError('proof and minimal must be booleans')
        self._curve, self._proof, self._minimal = curve, proof, minimal
        self._seed = self._selmer = self._pairing = self._builder = self._torsion = None
        self._two_attempts, self._two_decided = {}, {}
        self._entries, self._entry_attempts = {}, {}
        self._quartics, self._presentations, self._fibre_attempts, self._fibre_complete = {}, {}, {}, {}
        self._assemblies, self._groups = {}, {}
        self._work, self._timings = Counter(), Counter()
        self._last_result = None

    @classmethod
    def from_two_selmer(cls, selmer, *, proof=True):
        """Reuse an existing complete, evidence-checked basis without relabelling.

        proof controls NEW relative arithmetic only. A conditional parent keeps
        every assumption even when proof=True is requested for later stages.
        Changing the Selmer basis requires a separate session.
        """
        if not isinstance(selmer, TwoSelmerResult) or not selmer.complete or not selmer.verify_evidence():
            raise ValueError('expected a complete, evidence-verified labelled 2-Selmer result')
        obj = cls(selmer.curves.original, proof=proof, minimal=False)
        obj._seed = obj._selmer = selmer
        return obj

    @property
    def proof(self):
        return self._proof

    @property
    def last_result(self):
        """Last successfully retained snapshot, also accessible after an exception."""
        return self._last_result

    def cache_info(self):
        """A copy of actual construction counters; verification is never counted as a hit."""
        return dict(self._work)

    def _timed(self, name, function):
        start = perf_counter()
        try:
            return function()
        finally:
            self._timings[name] += perf_counter()-start

    def _two_selmer(self, limits):
        if self._selmer is not None and self._selmer.complete:
            return self._selmer
        if self._seed is None:
            self._seed = two_selmer(self._curve, proof=self.proof, minimal=self._minimal, max_candidates=0)
            self._work['arithmetic_contexts'] += 1
        base = self._seed
        labels = base.arithmetic.norm_kernel.elements()
        if limits.two_selmer_candidates is not None:
            labels = islice(labels, limits.two_selmer_candidates)
        decisions = []
        start = perf_counter()
        for label in labels:
            key = (label, limits.two_selmer_depth, limits.two_selmer_nodes)
            decision = self._two_decided.get(label, self._two_attempts.get(key))
            if decision is None:
                model = CubicCover(base.curves.working, base.arithmetic.representative(label))
                if base.arithmetic.coordinates(model.delta) != label or not model.verify():
                    raise ArithmeticError('labelled cubic construction verification failed')
                local = cubic_local_solubility(model, max_depth=limits.two_selmer_depth,
                                                max_nodes=limits.two_selmer_nodes)
                decision = TwoSelmerDecision(label, model, local)
                self._two_attempts[key] = decision
                self._work['two_selmer_candidates'] += 1
                if local.status != LocalStatus.UNKNOWN:
                    self._two_decided[label] = decision
            decisions.append(decision)
        result = TwoSelmerResult(base.curves, base.arithmetic, tuple(decisions),
            (('arithmetic_seconds', dict(base.timings).get('arithmetic_seconds', 0.0)),
             ('construction_local_seconds', perf_counter()-start)))
        if result.enumeration_complete and not result.unresolved and not result.complete:
            raise ArithmeticError('decided classes do not form a Selmer subspace')
        self._selmer = result
        return result

    def quartic_cover(self, coordinates):
        """On-demand canonical 2-cover, only after a complete 2-Selmer stage.

        This returns the exact same witness object on reuse. Certification and
        evidence verification remain live methods on that object.
        """
        S = self._selmer
        if S is None or not S.complete:
            raise ValueError('run the complete two_selmer stage before requesting a quartic')
        c = nonnegative_integer(coordinates, 'Selmer coordinates')
        S.ambient_label(c)
        if c not in self._quartics:
            self._quartics[c] = S.quartic_cover(c)
            self._work['quartic_presentations'] += 1
        return self._quartics[c]

    def _compute_pairing(self, limits):
        if self._pairing is not None and self._pairing.complete:
            return self._pairing
        if self._builder is None:
            self._builder = _SessionPairingBuilder(self)
        entries = []
        for i in range(self._selmer.dimension):
            for j in range(i+1, self._selmer.dimension):
                pair = (1 << i, 1 << j)
                key = (*pair, limits.pairing_trials)
                e = self._entries.get(pair, self._entry_attempts.get(key))
                if e is None:
                    e = self._builder.entry(*pair, max_local_trials=limits.pairing_trials)
                    self._entry_attempts[key] = e
                    self._work['pairing_entries'] += 1
                    if e.complete:
                        self._entries[pair] = e
                entries.append(e)
        self._pairing = CasselsTatePairing(self._selmer, tuple(entries))
        return self._pairing

    def _compute_fibres(self, limits, prune, reuse_zero_fibre):
        S, B = self._selmer, self._pairing
        if not B.verify_evidence():
            raise ValueError('pairing evidence is invalid')
        if self._torsion is None:
            self._torsion = S.torsion_data()
            self._work['torsion_tables'] += 1
        torsion = self._torsion
        image = tuple(c for c in range(S.cardinality) if B.in_radical(c))
        exclusions = tuple(PairingExclusion(c, next(1 << j for j in range(S.dimension)
                                     if B.evaluate(c, 1 << j)))
                           for c in range(S.cardinality) if c not in image) if prune else ()
        chosen = image if prune else tuple(range(S.cardinality))
        if limits.fibres is not None:
            chosen = chosen[:limits.fibres]
        assembled = []
        for c in chosen:
            # Reference-vs-zero-transport remains a distinct requested route.
            route = (c, reuse_zero_fibre)
            key = (*route, limits.relative_candidates, limits.relative_depth, limits.relative_nodes)
            f = self._fibre_complete.get(route, self._fibre_attempts.get(key))
            if f is None:
                if c not in self._presentations:
                    self._presentations[c] = fibre_representative(
                        S, c, torsion, initial_cover=self.quartic_cover(c))
                    self._work['fibre_presentations'] += 1
                f = _assemble_one(self._presentations[c], torsion, c in image, proof=self.proof,
                                  relative_limits=limits.relative, reuse_zero_fibre=reuse_zero_fibre)
                self._fibre_attempts[key] = f
                self._work['fibre_assemblies'] += 1
                if f.verify_structure(S, torsion):
                    self._fibre_complete[route] = f
            assembled.append(f)
        # Reuse the immutable assembly itself, preserving its group-element identity.
        key = (prune, tuple(id(f) for f in assembled))
        if key not in self._assemblies:
            result = FourSelmerFibres(S, B, torsion, tuple(assembled), exclusions, prune, ())
            if result.complete:
                _ = result.cardinality
                _ = result.kernel_map
            self._assemblies[key] = result
        return self._assemblies[key]

    def _compute_group(self, fibres, backend, limit):
        key = (id(fibres), backend)
        if key in self._groups:
            result = self._groups[key]
            evidence = result.galois_evidence
            if evidence is not None and limit is not None and evidence.torsion.field.degree > limit:
                raise FieldDegreeLimit(evidence.torsion.field.degree, limit)
            return result
        self._work['group_constructions'] += 1
        result = four_selmer_group(fibres, backend=backend, max_field_degree=limit)
        self._groups[key] = result
        return result

    def run(self, *, until='group', limits=None, prune=True, reuse_zero_fibre=True, backend='auto'):
        """Advance synchronously; retain partial results without guessing later stages.

        Unexpected arithmetic/backend failures propagate. last_result retains
        earlier successful stages. Only the typed field-degree resource limit
        is converted to an explicit partial result, never a false obstruction.
        """
        target = DescentStage(until)
        limits = DescentLimits() if limits is None else limits
        if not isinstance(limits, DescentLimits):
            raise TypeError('limits must be a DescentLimits instance')
        if type(prune) is not bool or type(reuse_zero_fibre) is not bool:
            raise TypeError('prune and reuse_zero_fibre must be booleans')
        if backend not in ('auto', 'galois'):
            raise ValueError("backend must be 'auto' or 'galois'")
        current = dict(selmer=None, pairing=None, fibres=None, group=None)

        def snapshot(reason=None, detail=()):
            result = FourDescentResult(target, limits, **current, stop_reason=reason,
                         limit_detail=detail, work=tuple(sorted(self._work.items())),
                         timings=tuple(sorted(self._timings.items())))
            self._last_result = result
            return result

        S = current['selmer'] = self._timed('two_selmer_seconds', lambda: self._two_selmer(limits))
        snapshot()
        if not S.complete:
            return snapshot('two_selmer_incomplete')
        if target == DescentStage.TWO_SELMER:
            return snapshot()
        B = current['pairing'] = self._timed('pairing_seconds', lambda: self._compute_pairing(limits))
        snapshot()
        if not B.complete:
            return snapshot('pairing_incomplete')
        if target == DescentStage.PAIRING:
            return snapshot()
        F = current['fibres'] = self._timed('fibres_seconds',
            lambda: self._compute_fibres(limits, prune, reuse_zero_fibre))
        snapshot()
        if not F.complete:
            return snapshot('fibres_incomplete')
        if target == DescentStage.FIBRES:
            return snapshot()
        try:
            current['group'] = self._timed('group_seconds',
                lambda: self._compute_group(F, backend, limits.field_degree))
        except FieldDegreeLimit as error:
            return snapshot('field_degree_limit', (('degree', error.degree), ('limit', error.limit)))
        return snapshot()


def four_descent(curve, *, proof=True, minimal=True, until='group', limits=None,
                 prune=True, reuse_zero_fibre=True, backend='auto'):
    """Compute a curve-level 4-descent, or an explicitly staged/limited result.

    Use DescentSession directly to retain arithmetic across repeated calls.
    Complete default mode returns the full group and marked representatives;
    it does not assert the Mordell--Weil rank or that any covering has a point.
    """
    return DescentSession(curve, proof=proof, minimal=minimal).run(
        until=until, limits=limits, prune=prune, reuse_zero_fibre=reuse_zero_fibre, backend=backend)
