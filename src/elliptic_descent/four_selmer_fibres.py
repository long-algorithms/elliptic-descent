"""Global 4-Selmer fibres as verified orbits of relative covering diagrams.

This implements the projection, global equivalence and kernel injection in
E(Q)[2] -> Sel2 -> Sel4 -> Sel2. It does NOT yet implement addition in Sel4.
Class keys are presentation-dependent; the actual marked covering is retained.
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import islice
from time import perf_counter

import sympy as sp

from .arithmetic import AffineF2, nonnegative_integer
from .backends.sage_fibres import oriented_class_matcher
from .cassels_tate import cassels_tate_pairing
from .equivalence import ScalarSquareEquivalence
from .fibre_actions import quartic_deck_action, FibreEquivalence, zero_fibre_kummer
from .quotient import normal_form, quotient_fibre
from .two_selmer import F2Labels


def action_orbits(permutations, size):
    """Exact orbit partition. Reject nonpermutations; never infer an orbit size."""
    size = nonnegative_integer(size,'action size')
    permutations = tuple(tuple(p) for p in permutations)
    for p in permutations:
        if len(p)!=size or any(type(i) is not int for i in p) or set(p)!=set(range(size)):
            raise ValueError('an action row is not a permutation of the diagrams')
    unseen,orbits = set(range(size)),[]
    while unseen:
        todo,orbit = [min(unseen)],set()
        while todo:
            i = todo.pop()
            if i in orbit:
                continue
            orbit.add(i)
            todo.extend(p[i] for p in permutations if p[i] not in orbit)
        unseen.difference_update(orbit)
        orbits.append(tuple(sorted(orbit)))
    return tuple(orbits)


@dataclass(frozen=True)
class TorsionActionRecord:
    point: object  # on the original E
    action: object
    permutation: tuple
    witnesses: tuple


@dataclass(frozen=True)
class GlobalFibreClass:
    """One orbit's marked 4-cover; not a Sel4 group-law coordinate."""
    image_coordinates: int
    orbit_index: int
    relative_indices: tuple
    covering: object

    @property
    def key(self):
        return self.image_coordinates,self.orbit_index


@dataclass(frozen=True)
class PairingExclusion:
    coordinates: int
    test_vector: int

    def verify(self, pairing):
        return pairing.evaluate(self.coordinates,self.test_vector) == 1


@dataclass(frozen=True)
class AssembledFibre:
    lifts: object
    expected_liftable: bool
    actions: tuple
    orbits: tuple
    matching_assumptions: tuple = ()

    @property
    def representative(self):
        return self.lifts.representative

    @property
    def coordinates(self):
        return self.representative.coordinates

    @property
    def relative(self):
        return self.lifts.marked_result.relative_result

    @property
    def models(self):
        return self.relative.coverings

    @property
    def assumptions(self):
        return tuple(dict.fromkeys((*self.lifts.assumptions,*self.matching_assumptions)))

    def verify_structure(self, selmer, torsion):
        if (self.representative.selmer is not selmer or not self.lifts.complete
                or type(self.expected_liftable) is not bool):
            return False
        models = self.models
        if not self.expected_liftable:
            return not models and not self.actions and not self.orbits
        if len(models) != selmer.cardinality:
            return False
        entries = torsion.entries
        if tuple(a.point for a in self.actions) != tuple(e.point for e in entries):
            return False
        try:
            if self.orbits != action_orbits((a.permutation for a in self.actions),len(models)):
                return False
        except (TypeError,ValueError):
            return False
        if any(len(orbit) != 1<<torsion.image_rank for orbit in self.orbits):
            return False
        lookup = {e.point:i for i,e in enumerate(entries)}
        identity = tuple(range(len(models)))
        E = selmer.curves.original
        for a,e in zip(self.actions,entries):
            if (a.action.presentation != self.representative.presentation
                    or a.action.point != selmer.curves.to_working.apply(e.point)
                    or len(a.witnesses) != len(models)):
                return False
            if (a.permutation == identity) != (e.coordinates == 0):
                return False
            if e.coordinates and any(i==j for i,j in enumerate(a.permutation)):
                return False  # Nontrivial translations of the Sel2 torsor are free.
            for b in self.actions:
                c = self.actions[lookup[E.add(a.point,b.point)]]
                if tuple(a.permutation[b.permutation[i]] for i in range(len(models))) != c.permutation:
                    return False
        return True

    def verify_evidence(self, selmer, torsion, *, max_nodes=1_000_000):
        if not self.verify_structure(selmer,torsion) or not self.lifts.verify_evidence(max_nodes=max_nodes):
            return False
        models = self.models
        for a in self.actions:
            if not a.action.verify():
                return False
            for i,(j,w) in enumerate(zip(a.permutation,a.witnesses)):
                if (w.transported.action != a.action or w.transported.source != models[i]
                        or w.target != models[j] or not w.verify()):
                    return False
        return True

    @property
    def classes(self):
        return tuple(GlobalFibreClass(self.coordinates,i,indices,
                      self.representative.marking.compose(self.models[indices[0]]))
                     for i,indices in enumerate(self.orbits))

    def equivalence(self, source_index, target_index):
        """Actual target -> source matrix when the diagrams define the same global class."""
        for value in (source_index,target_index):
            if type(value) is not int or not 0<=value<len(self.models):
                raise ValueError('relative diagram index is outside this fibre')
        for action in self.actions:
            if action.permutation[source_index] == target_index:
                return action.witnesses[source_index]
        return None


@dataclass(frozen=True)
class KernelInjection:
    """Explicit Sel2 -> zero fibre; kernel is the computed torsion Kummer image."""
    selmer: object
    fibre: AssembledFibre
    torsion: object
    labels: tuple  # Sel2 coordinates for each actual relative diagram

    @cached_property
    def quotient(self):
        S = self.selmer
        return quotient_fibre(AffineF2(S.dimension,0,tuple(1<<i for i in range(S.dimension))),
                              self.torsion.image_basis)[0]

    def apply(self, coordinates):
        self.selmer.ambient_label(coordinates)
        index = self.labels.index(int(coordinates))
        orbit = next(i for i,indices in enumerate(self.fibre.orbits) if index in indices)
        return self.fibre.classes[orbit]

    def quotient_coordinates(self, coordinates):
        self.selmer.ambient_label(coordinates)
        value = normal_form(int(coordinates),self.torsion.image_basis)
        return F2Labels(self.selmer.dimension,self.quotient.basis).coordinates(value)

    def from_quotient(self, coordinates):
        value = F2Labels(self.selmer.dimension,self.quotient.basis).ambient(coordinates)
        return self.apply(value)

    @property
    def kernel_basis(self):
        return self.torsion.image_basis

    def verify(self):
        S,F,D = self.selmer,self.fibre,self.torsion
        if (F.coordinates != 0 or not F.verify_structure(S,D)
                or len(self.labels)!=len(F.models) or any(type(c) is not int for c in self.labels)
                or set(self.labels)!=set(range(S.cardinality))):
            return False
        for label,model in zip(self.labels,F.models):
            if S.coordinates(zero_fibre_kummer(F.representative,model,S.arithmetic)) != label:
                return False
        # Verify every action on ACTUAL classes, not just its expected image size.
        for a,e in zip(F.actions,D.entries):
            if any(self.labels[j] != self.labels[i]^e.coordinates for i,j in enumerate(a.permutation)):
                return False
        for i in range(len(self.labels)):
            for j in range(len(self.labels)):
                same_orbit = any(i in orbit and j in orbit for orbit in F.orbits)
                same_coset = normal_form(self.labels[i]^self.labels[j],D.image_basis) == 0
                if same_orbit != same_coset:
                    return False
        return len(F.orbits) == self.quotient.cardinality


@dataclass(frozen=True)
class FourSelmerFibres:
    """Set of global classes with projection and exact kernel accounting.

    A class key (Sel2 coordinate, orbit index) is not an additive coordinate.
    Addition and invariant-factor generators belong to the subsequent milestone.
    """
    selmer: object
    pairing: object
    torsion: object
    fibres: tuple
    exclusions: tuple
    pruned: bool
    timings: tuple

    @property
    def assumptions(self):
        return tuple(dict.fromkeys((*self.selmer.assumptions,*self.pairing.assumptions,
                                   *(a for f in self.fibres for a in f.assumptions))))

    @property
    def image_basis(self):
        return self.pairing.radical_basis

    @property
    def image_classes(self):
        labels = F2Labels(self.selmer.dimension,self.image_basis)
        return tuple(sorted(labels.ambient(i) for i in range(1<<labels.dimension)))

    @property
    def kernel_dimension(self):
        return self.selmer.dimension-self.torsion.image_rank

    @property
    def predicted_cardinality(self):
        """Exact-sequence prediction, distinct from a completed explicit assembly."""
        return 1 << (self.kernel_dimension+len(self.image_basis))

    @property
    def complete(self):
        S,B,D = self.selmer,self.pairing,self.torsion
        if not S.complete or not B.complete or B.selmer is not S:
            return False
        expected = self.image_classes if self.pruned else tuple(range(S.cardinality))
        if tuple(f.coordinates for f in self.fibres) != expected:
            return False
        excluded = tuple(c for c in range(S.cardinality) if c not in self.image_classes) if self.pruned else ()
        if tuple(e.coordinates for e in self.exclusions) != excluded or not all(e.verify(B) for e in self.exclusions):
            return False
        return all(f.expected_liftable == B.in_radical(f.coordinates) and f.verify_structure(S,D)
                   for f in self.fibres)

    @property
    def certified(self):
        return (self.complete and self.selmer.certified and self.pairing.certified
                and not self.assumptions and all(f.lifts.certified for f in self.fibres))

    @property
    def cardinality(self):
        if not self.complete:
            raise ValueError('partial assembly has no computed global cardinality')
        count = sum(len(f.orbits) for f in self.fibres)
        if count != self.predicted_cardinality:
            raise ArithmeticError('explicit orbits disagree with the exact sequence')
        return count

    @property
    def known_classes(self):
        return tuple(c for f in self.fibres for c in f.classes)

    def fibre(self, coordinates):
        self.selmer.ambient_label(coordinates)
        found = [f for f in self.fibres if f.coordinates == int(coordinates)]
        if len(found)!=1:
            raise ValueError('this fibre was excluded by the pairing or has not been computed')
        return found[0]

    def projection(self, key):
        if len(key)!=2 or any(type(v) is not int for v in key):
            raise ValueError('expected a (Sel2 coordinate, orbit index) key')
        f = self.fibre(key[0])
        if not 0<=key[1]<len(f.orbits):
            raise ValueError('unknown global class key')
        return f.coordinates

    @cached_property
    def kernel_map(self):
        f = self.fibre(0)
        if not f.verify_structure(self.selmer,self.torsion):
            raise ValueError('zero fibre is incomplete or has invalid action accounting')
        labels = tuple(self.selmer.coordinates(zero_fibre_kummer(f.representative,D,self.selmer.arithmetic))
                       for D in f.models)
        result = KernelInjection(self.selmer,f,self.torsion,labels)
        if not result.verify():
            raise ArithmeticError('actual zero-fibre labels disagree with the kernel exact sequence')
        return result

    @property
    def zero_class(self):
        return self.kernel_map.apply(0)

    def verify_evidence(self, *, max_nodes=1_000_000):
        return (self.complete and self.selmer.verify_evidence(max_nodes=max_nodes)
                and self.pairing.verify_evidence() and self.torsion.verify(self.selmer)
                and all(f.verify_evidence(self.selmer,self.torsion,max_nodes=max_nodes) for f in self.fibres)
                and self.kernel_map.verify() and len(self.known_classes)==self.predicted_cardinality)


def _assemble_one(cover, torsion, expected_liftable, *, proof, relative_limits):
    lifts = cover.descend(proof=proof,**relative_limits)
    relative = lifts.marked_result.relative_result
    if not lifts.complete:
        return AssembledFibre(lifts,expected_liftable,(),())
    models = relative.coverings
    if not models:
        if expected_liftable:
            raise ArithmeticError('pairing predicts a lift, but complete relative descent is empty')
        return AssembledFibre(lifts,False,(),())
    S = cover.selmer
    if not expected_liftable:
        raise ArithmeticError('relative descent found a lift outside the pairing radical')
    if len(models) != S.cardinality:
        raise ArithmeticError('the relative Selmer torsor does not have the expected size')
    actions = tuple(quartic_deck_action(cover.presentation,S.curves.to_working.apply(e.point))
                    for e in torsion.entries)
    transported = tuple(tuple(g.transport(D) for D in models) for g in actions)
    matcher = None
    if len(actions)>1:
        matcher = oriented_class_matcher(cover.quartic.f,models,
                  tuple(t.raw for row in transported for t in row),proof=proof)
    records = []
    for entry,g,row in zip(torsion.entries,actions,transported):
        matches = []
        for i,t in enumerate(row):
            if entry.point is None:
                e = FibreEquivalence(t,models[i],ScalarSquareEquivalence(sp.S.One,(1,0,0,0)))
                if not e.verify():
                    raise ArithmeticError('identity action failed exact verification')
                matches.append((i,e))
            else:
                matches.append(matcher.match(t))
        records.append(TorsionActionRecord(entry.point,g,tuple(i for i,_ in matches),tuple(w for _,w in matches)))
    orbits = action_orbits((r.permutation for r in records),len(models))
    assumptions = () if matcher is None else matcher.arithmetic.space.assumptions
    result = AssembledFibre(lifts,True,tuple(records),orbits,assumptions)
    if not result.verify_structure(S,torsion):
        raise ArithmeticError('torsion action, stabilizers or orbit accounting are inconsistent')
    return result


def assemble_four_selmer_fibres(selmer, pairing=None, *, proof=True, prune=True,
                               max_fibres=None, max_candidates=None, max_depth=None, max_nodes=None):
    """Assemble explicit global fibres with equivalence witnesses, not group addition.

    prune=False independently descends ALL Sel2 classes and checks liftability
    against the pairing. Any limit retains partial status; no incomplete fibre
    is turned into a global quotient or an empty-fibre assertion.
    """
    if type(proof) is not bool or type(prune) is not bool:
        raise TypeError('proof and prune must be booleans')
    if not selmer.complete:
        raise ValueError('a complete labelled 2-Selmer space is required')
    if max_fibres is not None:
        max_fibres = nonnegative_integer(max_fibres,'max_fibres')
    limits = {name:None if value is None else nonnegative_integer(value,name)
              for name,value in (('max_candidates',max_candidates),('max_depth',max_depth),('max_nodes',max_nodes))}
    start = perf_counter()
    B = cassels_tate_pairing(selmer) if pairing is None else pairing
    if B.selmer is not selmer or not B.complete or not B.verify_evidence():
        raise ValueError('pairing must be complete, verified and in this exact Selmer basis')
    torsion = selmer.torsion_data()
    preparation_seconds = perf_counter()-start
    image = tuple(c for c in range(selmer.cardinality) if B.in_radical(c))
    exclusions = tuple(PairingExclusion(c,next(1<<j for j in range(selmer.dimension) if B.evaluate(c,1<<j)))
                       for c in range(selmer.cardinality) if c not in image) if prune else ()
    chosen = iter(image if prune else range(selmer.cardinality))
    if max_fibres is not None:
        chosen = islice(chosen,max_fibres)
    start = perf_counter()
    fibres = tuple(_assemble_one(selmer.quartic_cover(c),torsion,c in image,
                                proof=proof,relative_limits=limits) for c in chosen)
    result = FourSelmerFibres(selmer,B,torsion,fibres,exclusions,prune,
             (('pairing_and_torsion_seconds',preparation_seconds),('fibre_assembly_seconds',perf_counter()-start)))
    if result.complete:
        _ = result.cardinality
        _ = result.kernel_map
    return result
