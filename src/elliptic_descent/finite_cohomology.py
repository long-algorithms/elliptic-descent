"""Exact finite H^1 arithmetic for an E[4] module, without Sage.

Elements of (Z/4)^2 are encoded a+4*b. A supplied finite Galois action is
checked against its generator transitions. Cocycles are quotiented by all
sixteen actual coboundaries, not identified from group cardinalities.
"""
from dataclasses import dataclass
from functools import cached_property


def m4_add(a, b):
    return ((a+b) & 3) + 4*(((a//4+b//4) & 3))


def m4_neg(a):
    return (-a & 3) + 4*(-(a//4) & 3)


def m4_scale(n, a):
    return (n*a & 3) + 4*(n*(a//4) & 3)


def _label(a):
    return type(a) is int and 0 <= a < 16


@dataclass(frozen=True)
class FiniteFourModule:
    """All G elements, their action on E[4], and right-generator products.

    Group identity is row 0. ``transitions[s][j]`` is s*g_j. The field
    certificate supplies the proof that these transitions are actual Galois
    automorphism composition. This class alone checks finite consistency.
    """
    actions: tuple
    generators: tuple
    transitions: tuple

    def verify(self):
        n = len(self.actions)
        if (not n or len(self.transitions) != n
                or any(len(row) != len(self.generators) for row in self.transitions)):
            return False
        if (len(set(self.generators)) != len(self.generators)
                or any(type(g) is not int or not 0 < g < n for g in self.generators)):
            return False
        if self.actions[0] != tuple(range(16)):
            return False
        for row in self.actions:
            if len(row) != 16 or any(not _label(v) for v in row) or set(row) != set(range(16)):
                return False
            # Checking the two module generators determines every linear map.
            if any(row[a+4*b] != m4_add(m4_scale(a,row[1]),m4_scale(b,row[4]))
                   for a in range(4) for b in range(4)):
                return False
        for j,g in enumerate(self.generators):
            column = tuple(row[j] for row in self.transitions)
            if any(type(t) is not int or not 0 <= t < n for t in column) or set(column) != set(range(n)):
                return False
            if self.transitions[0][j] != g:
                return False
            for s in range(n):
                if self.actions[column[s]] != tuple(self.actions[s][v] for v in self.actions[g]):
                    return False
        seen,todo = {0},[0]
        while todo:
            for t in self.transitions[todo.pop()]:
                if t not in seen:
                    seen.add(t);todo.append(t)
        return len(seen) == n

    def is_cocycle(self, values):
        if len(values) != len(self.actions) or any(not _label(v) for v in values) or values[0] != 0:
            return False
        return all(values[self.transitions[s][j]] == m4_add(values[s],self.actions[s][values[g]])
                   for s in range(len(values)) for j,g in enumerate(self.generators))

    def coboundary(self, point):
        if not _label(point):
            raise ValueError('expected an E[4] coordinate in 0,...,15')
        return tuple(m4_add(row[point],m4_neg(point)) for row in self.actions)

    @cached_property
    def coboundaries(self):
        return tuple(sorted(set(self.coboundary(t) for t in range(16))))

    def canonical(self, cocycle):
        c = tuple(cocycle)
        if not self.is_cocycle(c):
            raise ValueError('the values do not form a 1-cocycle for this action')
        return min(tuple(m4_add(a,b) for a,b in zip(c,d)) for d in self.coboundaries)

    def sum(self, left, right):
        if not self.is_cocycle(left) or not self.is_cocycle(right):
            raise ValueError('addition requires two cocycles on the same module')
        return self.canonical(tuple(m4_add(a,b) for a,b in zip(left,right)))

    def equivalence_witness(self, left, right):
        """Return t with right-left = sigma(t)-t, or None if inequivalent."""
        if not self.is_cocycle(left) or not self.is_cocycle(right):
            raise ValueError('equivalence requires two cocycles')
        difference = tuple(m4_add(b,m4_neg(a)) for a,b in zip(left,right))
        return next((t for t in range(16) if self.coboundary(t) == difference), None)


@dataclass(frozen=True)
class CohomologyClassTable:
    """Actual covering keys mapped to distinct H^1 classes, with checked sums.

    The association with a covering MUST be proved by a separate geometric
    certificate. This object does not promote arbitrary labels to Selmer data.
    """
    module: FiniteFourModule
    keys: tuple
    cocycles: tuple

    @cached_property
    def signatures(self):
        return tuple(self.module.canonical(c) for c in self.cocycles)

    @cached_property
    def signature_index(self):
        if len(set(self.signatures)) != len(self.signatures):
            raise ValueError('two allegedly distinct coverings have cohomologous cocycles')
        return {s:i for i,s in enumerate(self.signatures)}

    def verify(self):
        try:
            if (not self.module.verify() or not self.keys or len(set(self.keys)) != len(self.keys)
                    or len(self.cocycles) != len(self.keys)):
                return False
            index = self.signature_index
            if (0,)*len(self.module.actions) not in index:
                return False
            return all(self.module.sum(a,b) in index for a in self.cocycles for b in self.cocycles)
        except (TypeError,ValueError):
            return False

    def add_indices(self, left, right):
        for v in (left,right):
            if type(v) is not int or not 0 <= v < len(self.keys):
                raise ValueError('class index outside the computed group')
        signature = self.module.sum(self.cocycles[left],self.cocycles[right])
        try:
            return self.signature_index[signature]
        except KeyError as error:
            raise ArithmeticError('sum is absent from the supposedly complete covering list') from error

    @property
    def zero_index(self):
        return self.signature_index[(0,)*len(self.module.actions)]
