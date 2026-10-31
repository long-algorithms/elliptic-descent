"""Exact published pairing, local symbols, support and corrupted-evidence checks."""
from dataclasses import replace
from itertools import permutations, product
from random import Random

import pytest
import sympy as sp

from elliptic_descent.arithmetic import prime_support
from elliptic_descent.cassels_tate_forms import FisherQuartic, FisherTriple, change_parameters, substitute_generator
from elliptic_descent.cassels_tate_local import evaluate_fisher, pairing_support, local_pairing
from elliptic_descent.cubic_arithmetic import kummer_algebra
from elliptic_descent.geometry import _sqrt_rational
from elliptic_descent.hilbert import local_square, hilbert_bit
from elliptic_descent.weierstrass import WeierstrassModel

# Fisher (2022), Example 3.4: fixed coefficients and product square root.
PAPER = ((-64,-164,-52,68,-11), (-3,-52,-232,-60,-4), (-53,102,32,-78,-31))
M = tuple(sp.Rational(c, 9) for c in (936032,-8656,20))


def test_published_nonzero_571a1_pairing_and_local_contributions():
    F = FisherTriple(PAPER, M)
    assert F.quartics[0].invariants == (44608,18842960)
    assert F.gamma == (-sp.Rational(16,3),-sp.Rational(64,9),sp.Rational(20,9))
    assert F.verify() and pairing_support(F) == (2,3,5,7,571)
    report = evaluate_fisher(F)
    assert report.complete and report.bit == 1 and report.verify()
    assert report.local[0].bit == 1
    assert all(r.bit == 0 for r in report.local[1:])


@pytest.mark.parametrize('order', list(permutations(range(3))))
def test_direct_pairings_on_all_nonzero_571a1_classes(order):
    F = FisherTriple(tuple(PAPER[i] for i in order), M)
    report = evaluate_fisher(F)
    assert report.verify() and report.bit == 1


def test_product_square_root_sign_changes_only_reciprocal_local_factors():
    plus = evaluate_fisher(FisherTriple(PAPER, M))
    minus = evaluate_fisher(FisherTriple(PAPER, tuple(-c for c in M)))
    assert plus.verify() and minus.verify() and plus.bit == minus.bit == 1
    assert tuple(-c for c in plus.triple.gamma) == minus.triple.gamma


@pytest.mark.parametrize('p', ['real',2,3,5,7,11,571])
def test_exact_hilbert_symmetry_bilinearity_and_square_factors(p):
    for a,b,c in ((-1,-1,2),(2,3,5),(-7,11,13),(sp.Rational(1,8),sp.Rational(-9,5),7)):
        assert hilbert_bit(a,b,p) == hilbert_bit(b,a,p)
        assert hilbert_bit(a,b*c,p) == hilbert_bit(a,b,p)^hilbert_bit(a,c,p)
        assert hilbert_bit(a, -a, p) == 0
        assert hilbert_bit(a, b*sp.Rational(16,25), p) == hilbert_bit(a,b,p)
        assert local_square(sp.Rational(49,36),p)
    if p == 2:
        assert hilbert_bit(-1,-1,2) == 1 and hilbert_bit(2,5,2) == 1
        assert local_square(17,2) and not local_square(5,2)
    if p == 'real':
        assert hilbert_bit(-1,-1,p) == 1 and not local_square(-1,p)


def test_hilbert_reciprocity_for_seeded_rational_pairs():
    rng=Random(504)
    for _ in range(35):
        a=sp.Rational(rng.choice([i for i in range(-31,32) if i]),rng.randrange(1,30))
        b=sp.Rational(rng.choice([i for i in range(-31,32) if i]),rng.randrange(1,30))
        assert sum(hilbert_bit(a,b,p) for p in ('real',*prime_support((2,a,b)))) % 2 == 0


def test_no_local_budget_can_be_misread_as_a_zero_pairing():
    F=FisherTriple(PAPER,M)
    limited=evaluate_fisher(F,max_local_trials=0)
    assert not limited.complete and not limited.verify()
    assert all(r.bit is None for r in limited.local)
    with pytest.raises(ValueError): _=limited.bit
    with pytest.raises(ValueError): evaluate_fisher(F,max_local_trials=-1)
    with pytest.raises(TypeError): evaluate_fisher(F,max_local_trials=1.0)
    report=local_pairing(F,2,max_trials=1)
    assert report.bit is None and not report.verify(F)


def test_corrupted_product_identity_support_and_local_evidence_fail():
    with pytest.raises(ValueError): FisherTriple(PAPER,(1,0,0))
    with pytest.raises(ValueError): FisherTriple((PAPER[0],PAPER[1],tuple(4*c for c in PAPER[2])),M)
    report=evaluate_fisher(FisherTriple(PAPER,M))
    assert not replace(report,support=report.support[:-1]).verify()
    assert not replace(report,local=report.local[:-1]).verify()
    first=report.local[0]
    assert not replace(first,bit=1-first.bit).verify(report.triple)
    assert not replace(first,point=(0,0)).verify(report.triple)
    assert not replace(first,point=(1.0,0)).verify(report.triple)
    for args in ((0,1,2),(1,0,2),(1,1,4)):
        with pytest.raises(ValueError): hilbert_bit(*args)
    with pytest.raises(TypeError): hilbert_bit(1.0,1,2)


def split_triple():
    # Three nonzero torsion Kummer classes on y^2=x^3-x, expressed with I=3,J=0.
    data = (((-4,0,24,0,-4),4), ((4,16,24,16,5),-2), ((32,128,192,128,34),-4))
    Q=[]
    for f,u in data:
        q=FisherQuartic(change_parameters(f,((sp.Rational(1,u),0),(0,1))))
        k=0
        while not q.algebra.norm(q.z):
            q=FisherQuartic(change_parameters(change_parameters(f,((sp.Rational(1,u),0),(0,1))),((k,-1),(1,0))))
            k+=1
        Q.append(q)
    A=Q[0].algebra
    value=A.multiply(A.multiply(Q[0].z,Q[1].z),Q[2].z)
    parts=A.to_components(value)
    root=A.from_components(tuple((_sqrt_rational(p[0]),) for p in parts))
    return Q,A,root


def test_split_algebra_all_independent_square_root_sign_choices():
    Q,A,root=split_triple()
    assert A.degrees == (1,1,1)
    for signs in product((-1,1),repeat=3):
        flip=A.from_components(tuple((s,) for s in signs))
        triple=FisherTriple(tuple(Q),A.multiply(root,flip))
        assert triple.verify()
        report=evaluate_fisher(triple)
        assert report.verify() and report.bit == 0


def test_zero_divisor_invariant_is_detected_and_changed_by_exact_gl2():
    f=FisherQuartic((-4,0,24,0,-4))
    assert f.algebra.norm(f.z) == 0
    q=FisherQuartic(change_parameters(f.f,((2,-1),(1,0))))
    assert q.invariants == f.invariants and q.algebra.norm(q.z) != 0
    assert q.verify_identity()
