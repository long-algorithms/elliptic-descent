"""Explicit group arithmetic on the actual global covering classes of Selmer-4.

A normal form is accepted only after geometric theta-cochain comparison binds
it to a constructed covering and verifies generator relations. The two exact-
sequence isomorphism cases need no algebraic-closure computation. An abstract
exponent-four presentation by itself never supplies geometric class labels.
"""
from dataclasses import dataclass, replace
from functools import cached_property
from time import perf_counter

from .exponent_four import ExponentFourPresentation
from .four_selmer_fibres import GlobalFibreClass
from .two_selmer import F2Labels


def _divide(left, right):
    return tuple(a/b for a, b in zip(left, right))


def _product(generators, a, b, k):
    result = [1]*16
    for i, values in enumerate(generators):
        selected = (a >> i) & 1 if i < k else (b >> (i-k)) & 1
        if selected:
            result = [x*y for x, y in zip(result, values)]
    return tuple(result)


@dataclass(frozen=True)
class ThetaGroupLink:
    """Evidence linking extension normal forms and relations to actual covers."""
    method: str
    normal_keys: tuple
    lift_keys: tuple
    context: object = None
    trivializations: tuple = ()
    cochains: tuple = ()
    normal_witnesses: tuple = ()
    relation_witnesses: tuple = ()
    inverse_witness: object = None

    def verify(self, fibres, presentation):
        F, P = fibres, presentation
        if not F.complete:
            return False
        classes = F.known_classes
        keys = tuple(c.key for c in classes)
        k, r = P.kernel_dimension, P.image_dimension
        if (k != F.kernel_dimension or r != len(F.image_basis)
                or P.doubles != tuple(F.kernel_map.quotient_coordinates(c) for c in F.image_basis)
                or len(self.normal_keys) != P.cardinality or len(set(self.normal_keys)) != P.cardinality
                or set(self.normal_keys) != set(keys) or len(self.lift_keys) != r):
            return False
        image = F2Labels(F.selmer.dimension, F.image_basis)
        for b in range(1 << r):
            for a in range(1 << k):
                key = self.normal_keys[a+(b << k)]
                if F.projection(key) != image.ambient(b):
                    return False
                if b == 0 and key != F.kernel_map.from_quotient(a).key:
                    return False
        if tuple(self.normal_keys[(1 << j) << k] for j in range(r)) != self.lift_keys:
            return False
        if self.method == 'kernel_inclusion_isomorphism':
            return r == 0 and self.context is None and not any((self.trivializations, self.cochains,
                                                               self.normal_witnesses, self.relation_witnesses))
        if self.method == 'projection_isomorphism':
            return k == 0 and self.context is None and not any((self.trivializations, self.cochains,
                                                               self.normal_witnesses, self.relation_witnesses))
        if self.method == 'cyclic_signed_inverse':
            return (k == r == 1 and P.doubles == (1,) and self.context is None
                    and not any((self.trivializations, self.cochains, self.normal_witnesses, self.relation_witnesses))
                    and self.inverse_witness is not None
                    and self.inverse_witness.source_key == self.normal_keys[2]
                    and self.inverse_witness.target_key == self.normal_keys[3]
                    and self.inverse_witness.verify(F))
        if self.method != 'exact_theta_cochains' or self.context is None:
            return False
        C = self.context
        if not C.verify():
            return False
        if (C.curve != F.selmer.curves.original or len(self.trivializations) != len(keys)
                or len(self.cochains) != len(keys) or len(self.normal_witnesses) != P.cardinality
                or len(self.relation_witnesses) != k+r):
            return False
        for element, trivialization, cochain in zip(classes, self.trivializations, self.cochains):
            if (trivialization.covering != element.covering or not trivialization.verify()
                    or C.cochain(trivialization) != cochain):
                return False
        lookup = dict(zip(keys, self.cochains))
        generator_keys = tuple(F.kernel_map.from_quotient(1 << i).key for i in range(k))+self.lift_keys
        generators = tuple(lookup[key] for key in generator_keys)
        for b in range(1 << r):
            for a in range(1 << k):
                index = a+(b << k)
                quotient = _divide(_product(generators, a, b, k), lookup[self.normal_keys[index]])
                if not self.normal_witnesses[index].verify(C, quotient):
                    return False
        for j, (values, witness) in enumerate(zip(generators, self.relation_witnesses)):
            expected = F.zero_class.key if j < k else F.kernel_map.apply(F.image_basis[j-k]).key
            quotient = _divide(tuple(v*v for v in values), lookup[expected])
            if not witness.verify(C, quotient):
                return False
        return True


@dataclass(frozen=True)
class SignedInverseWitness:
    """Changing the relative norm sign negates the full marked map to E."""
    source_key: tuple
    target_key: tuple
    opposite_index: int

    def verify(self, fibres):
        from .marked_cover import function_equal
        from .weierstrass import y
        try:
            f = fibres.fibre(self.source_key[0])
            source = f.classes[self.source_key[1]].covering.model
            if (type(self.opposite_index) is not int or not 0 <= self.opposite_index < len(f.models)
                    or f.models[self.opposite_index] != replace(source, norm_witness=-source.norm_witness)):
                return False
            target = next(c for c in f.classes if self.opposite_index in c.relative_indices)
            if target.key != self.target_key or f.representative.marking.translation is not None:
                return False
            marking = f.representative.marking
            X, Y = marking.map_xy
            a1, _, a3, _, _ = marking.target.ainvs
            return (function_equal(X.subs(y, -y), X, source.f)
                    and function_equal(Y.subs(y, -y), -Y-a1*X-a3, source.f))
        except (ValueError, IndexError, StopIteration):
            return False


def _signed_inverse(fibres, key):
    f = fibres.fibre(key[0])
    model = f.classes[key[1]].covering.model
    opposite = replace(model, norm_witness=-model.norm_witness)
    index = f.models.index(opposite)
    target = next(c for c in f.classes if index in c.relative_indices)
    result = SignedInverseWitness(key, target.key, index)
    if not result.verify(fibres):
        raise ArithmeticError('negative marked covering did not verify')
    return result


def _theta_link(fibres, presentation, lift_keys):
    from .backends.sage_theta_group import ThetaContext, flex_trivialization
    F, P = fibres, presentation
    C = ThetaContext(F.selmer.curves.original)
    classes = F.known_classes
    trivializations = tuple(flex_trivialization(c.covering) for c in classes)
    cochains = tuple(C.cochain(t) for t in trivializations)
    lookup = {c.key: values for c, values in zip(classes, cochains)}
    k, r = P.kernel_dimension, P.image_dimension
    generator_keys = tuple(F.kernel_map.from_quotient(1 << i).key for i in range(k))+lift_keys
    generators = tuple(lookup[key] for key in generator_keys)
    keys, witnesses, relations = [], [], []
    image = F2Labels(F.selmer.dimension, F.image_basis)
    for b in range(1 << r):
        candidates = F.fibre(image.ambient(b)).classes
        for a in range(1 << k):
            values = _product(generators, a, b, k)
            matches = []
            for candidate in candidates:
                witness = C.coboundary(_divide(values, lookup[candidate.key]))
                if witness is not None:
                    matches.append((candidate.key, witness))
            if len(matches) != 1:
                raise ArithmeticError('geometric theta sum has no unique assembled covering class')
            key, witness = matches[0]
            keys.append(key)
            witnesses.append(witness)
    for j, values in enumerate(generators):
        expected = F.zero_class.key if j < k else F.kernel_map.apply(F.image_basis[j-k]).key
        witness = C.coboundary(_divide(tuple(v*v for v in values), lookup[expected]))
        if witness is None:
            raise ArithmeticError('geometric doubling disagrees with the exact Selmer sequence')
        relations.append(witness)
    return ThetaGroupLink('exact_theta_cochains', tuple(keys), lift_keys, C,
                          trivializations, cochains, tuple(witnesses), tuple(relations))


@dataclass(frozen=True)
class FourSelmerGroup:
    """Finite abelian group with explicit marked covering representatives.

    Public elements are the existing GlobalFibreClass objects; tuple keys refer
    only to this assembly. Invariant coordinates and extension coordinates are
    separate, explicit presentations, never a re-interpretation of orbit keys.
    """
    fibres: object
    presentation: ExponentFourPresentation
    link: ThetaGroupLink
    timings: tuple

    @cached_property
    def _lookup(self):
        return {c.key: c for c in self.fibres.known_classes}

    @cached_property
    def _normal_lookup(self):
        k = self.presentation.kernel_dimension
        return {key: (i % (1 << k), i >> k) for i, key in enumerate(self.link.normal_keys)}

    @cached_property
    def complete(self):
        return self.fibres.complete and self.link.verify(self.fibres, self.presentation)

    @property
    def assumptions(self):
        return self.fibres.assumptions

    @property
    def certified(self):
        return self.complete and self.fibres.certified and not self.assumptions

    @property
    def cardinality(self):
        if not self.complete:
            raise ValueError('group-to-covering identification is incomplete or invalid')
        return self.presentation.cardinality

    @property
    def invariant_factors(self):
        if not self.complete:
            raise ValueError('invariant coordinates are not bound to valid covering classes')
        return self.presentation.invariant_factors

    @property
    def elements(self):
        if not self.complete:
            raise ValueError('incomplete group identification')
        return tuple(self._lookup[key] for key in self.link.normal_keys)

    @property
    def zero(self):
        return self._from_normal((0, 0))

    @property
    def generators(self):
        return tuple(self._from_normal(v) for v in self.presentation.invariant_generators)

    @property
    def relations(self):
        """Rows on (kernel generators, chosen quotient lifts), with exact carries."""
        k, r = self.presentation.kernel_dimension, self.presentation.image_dimension
        rows = []
        for i in range(k):
            rows.append(tuple(2*int(i == j) for j in range(k+r)))
        for j, d in enumerate(self.presentation.doubles):
            rows.append(tuple(-((d >> i) & 1) for i in range(k))
                        +tuple(2*int(j == i) for i in range(r)))
        return tuple(rows)

    def _key(self, value):
        if isinstance(value, GlobalFibreClass):
            if value.key not in self._lookup or value != self._lookup[value.key]:
                raise ValueError('covering element belongs to a different assembly')
            return value.key
        key = tuple(value)
        if len(key) != 2 or any(type(v) is not int for v in key) or key not in self._lookup:
            raise ValueError('unknown global covering key in this group')
        return key

    def normal_coordinates(self, value):
        if not self.complete:
            raise ValueError('incomplete group identification')
        return self._normal_lookup[self._key(value)]

    def _from_normal(self, value):
        if not self.complete:
            raise ValueError('incomplete group identification')
        a, b = self.presentation.point(value)
        return self._lookup[self.link.normal_keys[a+(b << self.presentation.kernel_dimension)]]

    def coordinates(self, value):
        return self.presentation.to_invariants(self.normal_coordinates(value))

    def element(self, coordinates):
        return self._from_normal(self.presentation.from_invariants(coordinates))

    def add(self, left, right):
        return self._from_normal(self.presentation.add(self.normal_coordinates(left), self.normal_coordinates(right)))

    def negate(self, value):
        return self._from_normal(self.presentation.negate(self.normal_coordinates(value)))

    def multiply(self, scalar, value):
        return self._from_normal(self.presentation.multiply(scalar, self.normal_coordinates(value)))

    def order(self, value):
        return self.presentation.order(self.normal_coordinates(value))

    def projection(self, value):
        return self.fibres.projection(self._key(value))

    def include_two_selmer(self, coordinates):
        if not self.complete:
            raise ValueError('incomplete group identification')
        return self.fibres.kernel_map.apply(coordinates)

    def verify_evidence(self, *, max_nodes=1_000_000):
        # Do not use a cached status to validate altered witnesses.
        return (self.fibres.verify_evidence(max_nodes=max_nodes)
                and self.link.verify(self.fibres, self.presentation))


def four_selmer_group(fibres, *, lift_keys=None, force_theta=False):
    """Give the completed global fibre assembly its actual cohomological group law.

    Default exact algebraic-closure calculations have no numerical tolerance or
    artificial class limit. Failure to compare a class raises, never assigns an
    arbitrary group structure. Conditional parent assumptions remain attached.
    """
    if type(force_theta) is not bool:
        raise TypeError('force_theta must be a boolean')
    F = fibres
    if not F.complete or not F.verify_evidence():
        raise ValueError('complete, evidence-verified global fibres are required')
    start = perf_counter()
    k = F.kernel_dimension
    P = ExponentFourPresentation(k, tuple(F.kernel_map.quotient_coordinates(c) for c in F.image_basis))
    if lift_keys is None:
        lift_keys = tuple(F.fibre(c).classes[0].key for c in F.image_basis)
    else:
        lift_keys = tuple(tuple(key) for key in lift_keys)
    if len(lift_keys) != len(F.image_basis) or any(F.projection(key) != c for key, c in zip(lift_keys, F.image_basis)):
        raise ValueError('one actual covering lift is required for each recorded radical basis vector')
    if not F.image_basis and not force_theta:
        link = ThetaGroupLink('kernel_inclusion_isomorphism',
                              tuple(F.kernel_map.from_quotient(a).key for a in range(1 << k)), ())
    elif k == 0 and not force_theta:
        image = F2Labels(F.selmer.dimension, F.image_basis)
        link = ThetaGroupLink('projection_isomorphism',
                              tuple(F.fibre(image.ambient(b)).classes[0].key for b in range(P.cardinality)), lift_keys)
    elif k == 1 and P.doubles == (1,) and not force_theta:
        inverse = _signed_inverse(F, lift_keys[0])
        keys = (F.zero_class.key, F.kernel_map.from_quotient(1).key,
                lift_keys[0], inverse.target_key)
        link = ThetaGroupLink('cyclic_signed_inverse', keys, lift_keys, inverse_witness=inverse)
    else:
        link = _theta_link(F, P, lift_keys)
    result = FourSelmerGroup(F, P, link, (('geometric_group_identification_seconds', perf_counter()-start),))
    if not result.complete:
        raise ArithmeticError('group normal forms or geometric relation witnesses failed verification')
    return result
