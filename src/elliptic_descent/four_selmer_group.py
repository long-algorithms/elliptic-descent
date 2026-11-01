"""Group arithmetic attached to the assembled global covering classes.

General addition is computed in actual Galois cocycles of the marked covers.
Only when the exact sequence uniquely determines every operation do we use
its kernel/projection directly. No arbitrary labelling by an abstract group
of the expected cardinality is permitted.
"""
from dataclasses import dataclass
from functools import cached_property
from itertools import product
from time import perf_counter

from .arithmetic import nonnegative_integer
from .galois_covering import GaloisGroupEvidence
from .limits import FieldDegreeLimit


def _key(key):
    if not isinstance(key,tuple) or len(key)!=2 or any(type(v) is not int or v<0 for v in key):
        raise ValueError('expected an exact (image coordinate, orbit index) key')
    return key


@dataclass(frozen=True)
class FiniteClassLaw:
    """An explicit operation on actual class keys, with exhaustive finite checks."""
    keys: tuple
    table: tuple
    zero: int

    def verify(self):
        n = len(self.keys)
        try:
            if not n or len(set(_key(k) for k in self.keys))!=n or len(self.table)!=n:
                return False
            if type(self.zero) is not int or not 0<=self.zero<n:
                return False
            for i,row in enumerate(self.table):
                if len(row)!=n or any(type(v) is not int or not 0<=v<n for v in row):
                    return False
                if set(row)!=set(range(n)) or row[self.zero]!=i:
                    return False
            if any(self.table[i][j]!=self.table[j][i] for i in range(n) for j in range(n)):
                return False
            if any(self.table[self.table[i][j]][k]!=self.table[i][self.table[j][k]]
                   for i in range(n) for j in range(n) for k in range(n)):
                return False
            return all(self.table[self.table[i][i]][self.table[i][i]]==self.zero for i in range(n))
        except (TypeError,ValueError,IndexError):
            return False

    def multiply(self, scalar, index):
        if type(scalar) is not int or type(index) is not int or not 0<=index<len(self.keys):
            raise ValueError('invalid scalar or class index')
        result = self.zero
        for _ in range(scalar % 4):
            result = self.table[result][index]
        return result

    def order(self, index):
        if type(index) is not int or not 0<=index<len(self.keys):
            raise ValueError('class index outside the group')
        return 1 if index==self.zero else 2 if self.table[index][index]==self.zero else 4

    @cached_property
    def decomposition(self):
        """Invariant factors WITH actual generators and a bijection to all covers."""
        if not self.verify():
            raise ValueError('invalid group law')
        def extend(coordinates,g,order):
            result = {}
            for c,i in coordinates.items():
                for a in range(order):
                    j = self.table[i][self.multiply(a,g)]
                    result[c+(a,)] = j
            if len(set(result.values()))!=len(result):
                raise ArithmeticError('invariant-factor generators are dependent')
            return result
        four,two = [],[]
        doubles = {self.zero}
        for i in range(len(self.keys)):
            d = self.table[i][i]
            if d not in doubles:
                four.append(i)
                doubles |= {self.table[x][d] for x in tuple(doubles)}
        coordinates = {():self.zero}
        for g in four:
            coordinates = extend(coordinates,g,4)
        for i in range(len(self.keys)):
            if self.order(i)==2 and i not in coordinates.values():
                two.append(i)
                coordinates = extend(coordinates,i,2)
        if len(coordinates)!=len(self.keys):
            raise ArithmeticError('generators do not cover every global class')
        generators = tuple(two+four)
        factors = (2,)*len(two)+(4,)*len(four)
        coordinates = {():self.zero}
        for g,order in zip(generators,factors):
            coordinates = extend(coordinates,g,order)
        return factors,generators,coordinates


def forced_class_law(keys, zero_key, projections, inclusion_keys):
    """The group law forced by specified exact-sequence maps, when unique.

    This checks finite data only. Their geometric meaning is supplied by the
    parent fibre evidence. It never extends this shortcut to a general fibre.
    """
    keys = tuple(keys)
    index = {k:i for i,k in enumerate(keys)}
    if not keys or len(index)!=len(keys) or zero_key not in index or len(projections)!=len(keys):
        raise ValueError('invalid exact-sequence class data')
    if (not inclusion_keys or len(inclusion_keys) & (len(inclusion_keys)-1)
            or any(type(p) is not int or not 0 <= p < len(inclusion_keys) for p in projections)):
        raise ValueError('invalid Selmer coordinate range')
    if any(k not in index or projections[index[k]]!=0 for k in inclusion_keys):
        raise ValueError('kernel inclusion does not land in the actual zero fibre')
    z = index[zero_key]
    if set(inclusion_keys)!={k for k,p in zip(keys,projections) if p==0} or inclusion_keys[0]!=zero_key:
        raise ValueError('kernel inclusion is not onto the zero fibre')
    if all(p==0 for p in projections):
        labels = tuple(inclusion_keys.index(k) for k in keys)
        table = tuple(tuple(index[inclusion_keys[a^b]] for b in labels) for a in labels)
        return 'kernel_isomorphism',FiniteClassLaw(keys,table,z)
    if len(set(inclusion_keys))==1:
        by_image = {p:i for i,p in enumerate(projections)}
        if len(by_image)!=len(keys):
            raise ValueError('projection with trivial kernel is not injective')
        table = tuple(tuple(by_image[a^b] for b in projections) for a in projections)
        return 'projection_isomorphism',FiniteClassLaw(keys,table,z)
    if len(keys)==4:
        # On these four specific elements, zero and doubling determine the law.
        doubles = tuple(index[inclusion_keys[p]] for p in projections)
        cyclic = next((i for i in range(4) if doubles[i]!=z),None)
        if cyclic is not None:
            coordinates = (z,cyclic,doubles[cyclic],next(i for i in range(4) if i not in (z,cyclic,doubles[cyclic])))
            inverse = {v:i for i,v in enumerate(coordinates)}
            table = tuple(tuple(coordinates[(inverse[i]+inverse[j])%4] for j in range(4)) for i in range(4))
        else:
            table = tuple(tuple(j if i==z else i if j==z else z if i==j else
                                next(k for k in range(4) if k not in (z,i,j))
                                for j in range(4)) for i in range(4))
        law = FiniteClassLaw(keys,table,z)
        if any(table[i][i]!=d for i,d in enumerate(doubles)) or not law.verify():
            raise ValueError('inconsistent order-four exact-sequence maps')
        return 'unique_order_four_extension',law
    return None


def _sequence_law(fibres):
    return forced_class_law(tuple(C.key for C in fibres.known_classes),fibres.zero_class.key,
                           tuple(fibres.projection(C.key) for C in fibres.known_classes),
                           tuple(fibres.kernel_map.apply(c).key for c in range(fibres.selmer.cardinality)))


def _check_sequence(fibres, law):
    keys = tuple(C.key for C in fibres.known_classes)
    if law.keys!=keys or law.keys[law.zero]!=fibres.zero_class.key or not law.verify():
        return False
    index = {k:i for i,k in enumerate(keys)}
    S = fibres.selmer
    inclusion = tuple(index[fibres.kernel_map.apply(c).key] for c in range(S.cardinality))
    for i,a in enumerate(keys):
        if law.table[i][i]!=inclusion[fibres.projection(a)]:
            return False  # [2] in S4 equals inclusion after projection, not projection itself.
        for j,b in enumerate(keys):
            if fibres.projection(keys[law.table[i][j]]) != fibres.projection(a)^fibres.projection(b):
                return False
    return all(law.table[inclusion[a]][inclusion[b]]==inclusion[a^b]
               for a in range(S.cardinality) for b in range(S.cardinality))


def _bound_to_fibres(evidence, fibres):
    classes = fibres.known_classes
    if (evidence.torsion.ainvs!=fibres.selmer.curves.original.ainvs
            or tuple(c.key for c in evidence.coverings)!=tuple(C.key for C in classes)):
        return False
    for witness,C in zip(evidence.coverings,classes):
        D = C.covering.model
        M = C.covering.marking
        if (M.translation is not None or witness.f!=D.f or witness.xi!=D.xi
                or witness.norm_witness!=D.norm_witness or witness.jacobian_to_original!=M.change.urst):
            return False
    return True


@dataclass(frozen=True,eq=False)
class FourSelmerElement:
    group: object
    index: int

    def __post_init__(self):
        if type(self.index) is not int or not 0<=self.index<len(self.group.law.keys):
            raise ValueError('invalid global Selmer element index')

    @property
    def key(self):
        return self.group.law.keys[self.index]

    @property
    def covering(self):
        return self.group.fibres.known_classes[self.index].covering

    @property
    def coordinates(self):
        return self.group.coordinates(self)

    def order(self):
        return self.group.law.order(self.index)

    def __eq__(self,other):
        return isinstance(other,FourSelmerElement) and self.group is other.group and self.index==other.index

    def __hash__(self):
        return hash((id(self.group),self.index))

    def __add__(self,other):
        return self.group.add(self,other)

    def __neg__(self):
        return self.group.negate(self)

    def __sub__(self,other):
        return self+(-other)

    def __rmul__(self,scalar):
        return self.group.multiply(scalar,self)


@dataclass(frozen=True)
class FourSelmerGroup:
    fibres: object
    law: FiniteClassLaw
    method: str
    galois_evidence: object
    timings: tuple

    @property
    def complete(self):
        return self.fibres.complete and len(self.law.keys)==self.fibres.cardinality

    @property
    def assumptions(self):
        return self.fibres.assumptions

    @property
    def certified(self):
        return self.fibres.certified and not self.assumptions and self.verify_evidence(parent=False)

    @property
    def cardinality(self):
        if not self.complete:
            raise ValueError('an incomplete assembly has no completed group cardinality')
        return len(self.law.keys)

    @property
    def zero(self):
        return FourSelmerElement(self,self.law.zero)

    @property
    def elements(self):
        return tuple(FourSelmerElement(self,i) for i in range(self.cardinality))

    def element(self, key):
        key = _key(key)
        try:
            return FourSelmerElement(self,self.law.keys.index(key))
        except ValueError as error:
            raise ValueError('no such global class in this group') from error

    def _index(self, element):
        if not isinstance(element,FourSelmerElement) or element.group is not self:
            raise ValueError('elements from different group presentations cannot be mixed')
        return element.index

    def add(self, left, right):
        return FourSelmerElement(self,self.law.table[self._index(left)][self._index(right)])

    def multiply(self, scalar, element):
        if type(scalar) is not int:
            raise TypeError('Selmer scalar must be an exact integer')
        return FourSelmerElement(self,self.law.multiply(scalar,self._index(element)))

    def negate(self, element):
        return self.multiply(-1,element)

    def projection(self, element):
        return self.fibres.projection(self.law.keys[self._index(element)])

    def include_two_selmer(self, coordinates):
        return self.element(self.fibres.kernel_map.apply(coordinates).key)

    @property
    def invariant_factors(self):
        return self.law.decomposition[0]

    @property
    def generators(self):
        return tuple(FourSelmerElement(self,i) for i in self.law.decomposition[1])

    @cached_property
    def _coordinates(self):
        return {i:c for c,i in self.law.decomposition[2].items()}

    def coordinates(self, element):
        return self._coordinates[self._index(element)]

    def from_coordinates(self, coordinates):
        coordinates = tuple(coordinates)
        if len(coordinates)!=len(self.invariant_factors) or any(type(c) is not int for c in coordinates):
            raise ValueError('wrong invariant-factor coordinate tuple')
        reduced = tuple(c%d for c,d in zip(coordinates,self.invariant_factors))
        return FourSelmerElement(self,self.law.decomposition[2][reduced])

    @property
    def relations(self):
        return tuple(tuple(n if i==j else 0 for j in range(len(self.generators)))
                     for i,n in enumerate(self.invariant_factors))

    def verify_evidence(self, *, parent=True, max_nodes=1_000_000):
        if not self.complete or not _check_sequence(self.fibres,self.law):
            return False
        if parent and not self.fibres.verify_evidence(max_nodes=max_nodes):
            return False
        if self.method=='geometric_galois_cocycles':
            E = self.galois_evidence
            if not isinstance(E,GaloisGroupEvidence) or not _bound_to_fibres(E,self.fibres) or not E.verify():
                return False
            T = E.class_table
            if self.law.keys!=T.keys or self.law.zero!=T.zero_index:
                return False
            if any(self.law.table[i][j]!=T.add_indices(i,j) for i in range(self.cardinality) for j in range(self.cardinality)):
                return False
        else:
            expected = _sequence_law(self.fibres)
            if expected is None or (self.method,self.law)!=expected or self.galois_evidence is not None:
                return False
        try:
            factors,gens,coords = self.law.decomposition
            return (len(coords)==self.cardinality and len(set(coords.values()))==self.cardinality
                    and all(self.law.order(g)==d for g,d in zip(gens,factors)))
        except (TypeError,ValueError,ArithmeticError):
            return False


def four_selmer_group(fibres, *, backend='auto', max_field_degree=None, evidence=None):
    """Group operations on completed global fibres, with actual covering witnesses.

    ``backend='galois'`` forces the algebraic cocycle reference path even when
    a uniquely determined exact-sequence shortcut is available. Supplying exact
    Galois evidence permits replay without Sage. Discovery failures raise; they
    never produce a fabricated group with the predicted invariant factors.
    """
    if backend not in ('auto','galois'):
        raise ValueError("backend must be 'auto' or 'galois'")
    if max_field_degree is not None:
        nonnegative_integer(max_field_degree,'maximum field degree')
    if not fibres.complete or not fibres.verify_evidence():
        raise ValueError('a complete, evidence-verified global fibre assembly is required')
    start = perf_counter()
    shortcut = _sequence_law(fibres) if backend=='auto' and evidence is None else None
    if shortcut is not None:
        method,law = shortcut
    else:
        from .backends.sage_four_group import discover_group_evidence
        E = discover_group_evidence(fibres,max_field_degree=max_field_degree) if evidence is None else evidence
        if not isinstance(E,GaloisGroupEvidence) or not _bound_to_fibres(E,fibres) or not E.verify():
            raise ValueError('the Galois data do not certify these specific marked coverings')
        if max_field_degree is not None and E.torsion.field.degree>max_field_degree:
            raise FieldDegreeLimit(E.torsion.field.degree, max_field_degree)
        table = E.class_table
        law = FiniteClassLaw(table.keys,tuple(tuple(table.add_indices(i,j) for j in range(len(table.keys)))
                                             for i in range(len(table.keys))),table.zero_index)
        method,evidence = 'geometric_galois_cocycles',E
    result = FourSelmerGroup(fibres,law,method,evidence,(('group_arithmetic_seconds',perf_counter()-start),))
    if not result.verify_evidence(parent=False):
        raise ArithmeticError('group arithmetic, exact-sequence maps or covering association failed verification')
    return result
