"""Labelled pairing matrices, independent PARI ranks and exact local evidence."""
import os
from dataclasses import replace
from functools import lru_cache
from itertools import product, permutations

import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import two_selmer, WeierstrassModel, WeierstrassChange
from elliptic_descent.cassels_tate import cassels_tate_pairing, pairing_entry, normalize_pairing_quartic
from elliptic_descent.hilbert import hilbert_bit, local_square

FIXTURES = ((0,0,0,-1,1), (0,0,0,-1,0), (0,0,0,1,0), (0,0,0,4,0),
            (0,-5,0,-16,80), (0,-1,1,-10,-20), (0,-1,1,-929,-10595))


@lru_cache(maxsize=16)
def selmer(a):
    return two_selmer(a)


@lru_cache(maxsize=16)
def pairing(a):
    return cassels_tate_pairing(selmer(a), max_local_trials=10000)


@pytest.mark.parametrize('a', FIXTURES)
def test_pairing_matrix_and_radical_against_pari_rank(a):
    S=selmer(a); P=pairing(a)
    assert P.complete and P.certified and P.verify_evidence()
    independent=sage.pari.ellinit(list(a)).ellrank()
    assert P.rank==int(independent[2])
    assert len(P.radical_basis)==S.dimension-P.rank
    assert len(P.radical_ambient_basis)==len(P.radical_basis)
    for c,d,e in product(range(S.cardinality),repeat=3):
        assert P.evaluate(c,c)==0
        assert P.evaluate(c,d)==P.evaluate(d,c)
        assert P.evaluate(c,d^e)==P.evaluate(c,d)^P.evaluate(c,e)
    for entry in S.torsion_data().entries:
        assert P.in_radical(entry.coordinates)
    if a==FIXTURES[-1]:
        assert P.matrix==((0,1),(1,0)) and P.rank==2 and not P.radical_basis
    else:
        assert P.rank==0


def test_direct_formula_agrees_for_each_nonzero_pair_not_just_matrix_extension():
    a=FIXTURES[-1]; P=pairing(a)
    for c,d in permutations((1,2,3),2):
        direct=P.direct_entry(c,d,max_local_trials=10000)
        assert direct.complete and direct.verify() and direct.bit==P.evaluate(c,d)==1
    for c in range(4):
        assert P.direct_entry(c,c).bit==0 and P.direct_entry(0,c).verify()


def test_rebasing_and_independent_recomputation_preserve_the_form():
    S=selmer(FIXTURES[-1]); P=pairing(FIXTURES[-1])
    b0,b1=S.basis_labels
    new=P.with_basis((b0^b1,b0))
    assert new.certified and new.verify_evidence() and new.rank==P.rank
    direct=cassels_tate_pairing(new.selmer,max_local_trials=10000)
    assert direct.verify_evidence() and direct.matrix==new.matrix
    for c,d in product(range(4),repeat=2):
        cc=new.selmer.coordinates_of_label(S.ambient_label(c))
        dd=new.selmer.coordinates_of_label(S.ambient_label(d))
        assert new.evaluate(cc,dd)==P.evaluate(c,d)


def test_limits_conditional_arithmetic_and_missing_entries_fail_closed():
    S=selmer(FIXTURES[-1])
    partial=cassels_tate_pairing(S,max_local_trials=0)
    assert not partial.complete and not partial.certified and not partial.verify_evidence()
    for prop in ('matrix','rank','radical_basis'):
        with pytest.raises(ValueError): getattr(partial,prop)
    with pytest.raises(ValueError): cassels_tate_pairing(two_selmer(FIXTURES[0],max_candidates=0))
    conditional=cassels_tate_pairing(two_selmer(FIXTURES[-1],proof=False),max_local_trials=10000)
    assert conditional.complete and conditional.assumptions and not conditional.certified
    assert conditional.verify_evidence()
    P=pairing(FIXTURES[-1])
    assert not replace(P,entries=()).complete
    assert not replace(P,entries=P.entries*2).verify_evidence()
    record=P.entries[0]
    assert not replace(record,left=0).verify()
    bad=replace(record.evaluation,local=record.evaluation.local[:-1])
    assert not replace(record,evaluation=bad).verify()
    for value in (-1,S.cardinality):
        with pytest.raises(ValueError): P.evaluate(value,0)
    with pytest.raises(TypeError): P.evaluate(1.0,1)


def test_normalization_recovers_classes_and_handles_zero_divisor_invariants():
    for a in (FIXTURES[1],FIXTURES[4],FIXTURES[-1]):
        S=selmer(a)
        for c in range(S.cardinality):
            C=S.quartic_cover(c)
            N=normalize_pairing_quartic(C)
            assert N.verify()
            assert S.coordinates(N.theta_invariant)==c
            assert N.quartic.algebra.norm(N.quartic.z)
            assert not replace(N,parameter_matrix=((2,0),(0,2))).verify()


def test_hilbert_symbols_against_pari_including_two_and_infinity():
    values=[-11,-7,-1,2,3,5,sp.Rational(1,8),sp.Rational(-9,5)]
    for a,b,p in product(values,values,('real',2,3,5,7,11)):
        expected=int(sage.pari.hilbert(sage.QQ(str(a)),sage.QQ(str(b)),0 if p=='real' else p))
        assert (-1)**hilbert_bit(a,b,p)==expected


def test_rational_kummer_classes_lie_in_radical_and_model_changes_commute():
    base=WeierstrassModel((-1,1)); S=selmer(base.ainvs); P=pairing(base.ainvs)
    assert P.in_radical(S.kummer((0,1)))
    change=WeierstrassChange(base,(sp.Rational(1,2),1,-1,3))
    other=two_selmer(change.target)
    Q=cassels_tate_pairing(other)
    assert Q.certified and Q.verify_evidence()
    for c,d in product(range(S.cardinality),repeat=2):
        assert Q.evaluate(S.transport_to(other,c,change),S.transport_to(other,d,change))==P.evaluate(c,d)


def test_parameter_changed_representatives_give_the_same_nonzero_value():
    from elliptic_descent.cassels_tate import _PairingBuilder
    S=selmer(FIXTURES[-1])
    builder=_PairingBuilder(S)
    for c in (1,2,3):
        C=S.quartic_cover(c,parameter_change=((1,1),(0,1)))
        builder.normalized[c]=normalize_pairing_quartic(C)
    entry=builder.entry(1,2,max_local_trials=20000)
    assert entry.complete and entry.verify() and entry.bit==pairing(FIXTURES[-1]).evaluate(1,2)==1
