"""Release-wide audit, including a separate fully split candidate enumeration."""
import os
from dataclasses import replace
from itertools import product

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import EtaleQuarticAlgebra, relative_four_descent, construct_four_cover
from elliptic_descent.arithmetic import square_class
from elliptic_descent.local_complete import complete_local_screen
from elliptic_descent.status import LocalStatus


@pytest.mark.parametrize('f,multiplicity', [
    ((1,-1,0,0,1), 2), ((6,0,-5,0,1), 2),
    ((0,1,-1,0,1), 1), ((2,0,-3,0,1), 1), ((4,0,-5,0,1), 1),
])
def test_all_types_retain_an_auditable_candidate_ledger(f, multiplicity):
    r = relative_four_descent(f)
    assert r.accounting_consistent and r.enumeration_complete and r.certified
    assert r.total_model_classes == r.candidate_space.cardinality
    assert {c.label for c in r.candidates} == set(r.candidate_space.elements())
    assert len(r.coverings) == multiplicity*len(r.unoriented_coverings)
    assert r.verify_evidence()
    assert len(r.candidates) >= 2
    # This keeps the length unchanged, unlike ordinary truncation tests.
    duplicate = replace(r, candidates=(r.candidates[0],)*len(r.candidates))
    assert not duplicate.enumeration_complete and not duplicate.certified
    assert not duplicate.verify_evidence()
    outside = replace(r.candidates[0], label=1 << r.candidate_space.nvariables)
    forged = replace(r, candidates=(outside, *r.candidates[1:]))
    assert not forged.complete and not forged.verify_evidence()
    wrong_count = replace(r, total_model_classes=r.total_model_classes+1)
    assert not wrong_count.complete


def test_empty_descent_has_a_consistent_empty_arithmetic_fibre():
    r = relative_four_descent((-3,52,-232,60,-4))
    assert r.candidate_space.particular is None and r.total_model_classes == 0
    assert r.accounting_consistent and r.certified and r.verify_evidence()


def test_split_classes_against_independent_rational_enumeration():
    f = (4,0,-5,0,1)
    algebra = EtaleQuarticAlgebra(f)
    result = relative_four_descent(f)
    # Fix first component to 1 using the diagonal rational scalar. Choose the
    # next two in Q({2,3},2); the norm condition determines the last. No Sage
    # Selmer space or affine quotient implementation is used in this enumeration.
    reps = tuple(sp.Rational(sign*2**i*3**j) for sign,i,j in product((-1,1),range(2),range(2)))
    rep_by_class = {square_class(a,(2,3)): a for a in reps}

    def key(xi):
        parts = [c[0] for c in algebra.to_components(xi)]
        return tuple(square_class(c/parts[0],(2,3)) for c in parts)

    seen, accepted = set(), set()
    for a,b in product(reps, repeat=2):
        c = rep_by_class[square_class(a*b,(2,3))]
        xi = algebra.from_components(((1,),(a,),(b,),(c,)))
        label = key(xi)
        assert label not in seen
        seen.add(label)
        model = construct_four_cover(f,xi)
        local = complete_local_screen(model)
        assert local.status != LocalStatus.UNKNOWN and local.verify(model)
        if local.status == LocalStatus.SOLUBLE:
            accepted.add(label)
    assert len(seen) == result.total_model_classes == 64
    assert accepted == {key(D.xi) for D in result.coverings}
    assert len(accepted) == 4
