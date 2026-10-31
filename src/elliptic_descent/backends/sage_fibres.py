"""Exact square-class matching for a finite set of oriented relative covers.

Enlarge support explicitly, solve the rational-scalar image equation over F2,
then return an independently checked square-root and signed norm witness. This
is equivalence discovery, not a new local-solubility/completeness computation.
"""
from dataclasses import dataclass

import sympy as sp

from .sage_etale import sage_etale_square_classes
from .sage_pairing import cubic_product_square_root
from ..arithmetic import prime_support, quartic_support, square_class, FiniteNormSpace
from ..status import ProofStatus
from ..equivalence import ScalarSquareEquivalence, orientation_certificate
from ..etale import EtaleQuarticAlgebra
from ..fibre_actions import FibreEquivalence
from ..quotient import row_basis, normal_form
from ..two_selmer import F2Labels


@dataclass(frozen=True)
class OrientedClassMatcher:
    arithmetic: object
    models: tuple
    scalar_values: tuple
    scalar_labels: tuple
    keys: tuple

    @property
    def algebra(self):
        return self.arithmetic.algebra

    @property
    def scalar_rows(self):
        return row_basis(self.scalar_labels)

    def key(self, xi):
        return normal_form(self.arithmetic.coordinates(xi),self.scalar_rows)

    def scalar_square(self, source, target):
        """Return target.xi=c*v^2*source.xi with its prescribed norm sign, or None."""
        if source.f != self.algebra.f or target.f != source.f:
            raise ValueError('equivalence requires the same intermediate quartic')
        labels = F2Labels(self.arithmetic.dimension,self.scalar_labels)
        difference = self.arithmetic.coordinates(source.xi)^self.arithmetic.coordinates(target.xi)
        try:
            bits = labels.coordinates(difference)
        except ValueError:
            return None
        c = sp.prod(v for i,v in enumerate(self.scalar_values) if bits>>i & 1)
        ratio = self.algebra.multiply(target.xi,self.algebra.inverse(source.xi))
        ratio = tuple(v/c for v in ratio)
        # The helper works with any labelled product algebra and checks v^2
        # in rational power-basis arithmetic, despite its historical cubic name.
        v = cubic_product_square_root(self.arithmetic,ratio)
        if c*c*self.algebra.norm(v)*source.norm_witness != target.norm_witness:
            flip = orientation_certificate(self.algebra).reversing_multiplier
            if flip is None:
                return None
            v = self.algebra.multiply(v,flip)
        witness = ScalarSquareEquivalence(c,v)
        if not witness.verify(source,target):
            raise ArithmeticError('discovered oriented equivalence failed exact verification')
        return witness

    def match(self, transported):
        key = self.key(transported.raw.xi)
        found = []
        for i,(model,k) in enumerate(zip(self.models,self.keys)):
            if k != key:
                continue
            witness = self.scalar_square(transported.raw,model)
            if witness is not None:
                evidence = FibreEquivalence(transported,model,witness)
                if not evidence.verify():
                    raise ArithmeticError('global diagram equivalence failed verification')
                found.append((i,evidence))
        if len(found) != 1:
            raise ArithmeticError(f'torsion action has {len(found)} matches, expected exactly one')
        return found[0]


def oriented_class_matcher(f, models, additional=(), *, proof=True):
    """Create a shared coordinate context for all action images and target models.

    Outside support the monic order is integral etale, coordinates are integral,
    and each element has unit norm, so it is a unit at every prime. This ensures
    scalar-image testing sees every possible rational scalar identification.
    """
    A = EtaleQuarticAlgebra(f)
    values = []
    for D in (*models,*additional):
        if D.f != A.f:
            raise ValueError('matcher models have different intermediate quartics')
        values.extend(c.q for c in D.xi)
        values.append(A.norm(D.xi))
    support = tuple(sorted(set(quartic_support(A.f))|set(prime_support(values))))
    arithmetic = (rational_product_square_classes(A,support) if A.degrees == (1,1,1,1)
                  else sage_etale_square_classes(A,support,proof=proof))
    scalar_values,scalar_labels = [],[]
    for c in (-1,*support):
        label = arithmetic.coordinates(A.scalar(c))
        if len(row_basis((*scalar_labels,label))) > len(scalar_labels):
            scalar_values.append(sp.Integer(c))
            scalar_labels.append(label)
    rows = row_basis(scalar_labels)
    keys = tuple(normal_form(arithmetic.coordinates(D.xi),rows) for D in models)
    result = OrientedClassMatcher(arithmetic,tuple(models),tuple(scalar_values),tuple(scalar_labels),keys)
    # Distinct norm-square classes may have the same coordinates only through
    # a scalar. Check oriented nonredundancy rather than assuming it from counts.
    for i,D in enumerate(models):
        for j in range(i):
            if keys[i] == keys[j] and result.scalar_square(D,models[j]) is not None:
                raise ArithmeticError('duplicate fixed-target relative covering in matcher input')
    return result


@dataclass(frozen=True)
class _RationalFactor:
    degree: int = 1


@dataclass(frozen=True)
class RationalProductSquareClasses:
    """Elementary Q^4 square classes: sign and supported prime valuation parity.

    No Sage arithmetic or assumed basis is used. This is also an independent
    reference for split matching. It makes no local-solubility claim.
    """
    algebra: object
    space: object

    @property
    def factors(self):
        return (_RationalFactor(),)*4

    @property
    def dimension(self):
        return 4*(len(self.space.support)+1)

    def coordinates(self, element):
        width = len(self.space.support)+1
        parts = self.algebra.to_components(element)
        return sum(square_class(part[0],self.space.support) << (width*i) for i,part in enumerate(parts))


def rational_product_square_classes(algebra, support):
    if algebra.degrees != (1,1,1,1):
        raise ValueError('elementary rational matching requires four linear factors')
    generators = tuple(algebra.embed_factor_unit(i,(c,)) for i in range(4) for c in (-1,*support))
    space = FiniteNormSpace(algebra.f,generators,tuple(support),ProofStatus.CERTIFIED,(),
                            'exact rational sign and prime valuation parity')
    result = RationalProductSquareClasses(algebra,space)
    if any(result.coordinates(g) != 1<<i for i,g in enumerate(generators)):
        raise ArithmeticError('rational product basis failed exact round trips')
    return result
