"""Radical/liftability checks using the relative solver, not pairing-based pruning."""
import os
from functools import lru_cache
from itertools import product
import pytest
import sympy as sp

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage = pytest.importorskip('sage.all', reason='requires SageMath runtime')

from elliptic_descent import two_selmer, relative_four_descent, WeierstrassChange
from elliptic_descent.cassels_tate import cassels_tate_pairing, pairing_entry
from elliptic_descent.cassels_tate_forms import FisherQuartic, change_parameters, substitute_generator


@lru_cache(maxsize=8)
def computed(a):
    S=two_selmer(a)
    return S,cassels_tate_pairing(S,max_local_trials=20000)


def test_nonzero_pairing_with_nonzero_radical_and_arbitrary_basis_change():
    a=(0,0,0,-289,0)  # Full rational 2-torsion, with nontrivial Sha[2].
    S,P=computed(a)
    assert S.dimension==4 and P.complete and P.certified and P.verify_evidence()
    assert P.rank==int(sage.pari.ellinit(list(a)).ellrank()[2])==2
    assert len(P.radical_basis)==2
    torsion={S.kummer(T.point) for T in S.torsion_data().entries}
    radical={c for c in range(S.cardinality) if P.in_radical(c)}
    assert radical==torsion
    basis=list(S.basis_labels)
    basis[0]^=basis[1]; basis[2]^=basis[3]; basis.reverse()
    Q=P.with_basis(basis)
    assert Q.certified and Q.verify_evidence()
    for c,d in product(range(S.cardinality),repeat=2):
        cc=Q.selmer.coordinates_of_label(S.ambient_label(c))
        dd=Q.selmer.coordinates_of_label(S.ambient_label(d))
        assert P.evaluate(c,d)==Q.evaluate(cc,dd)
    # Independently evaluate a sum in each slot, not just a basis matrix entry.
    direct=pairing_entry(S,3,12,max_local_trials=20000)
    assert direct.complete and direct.verify() and direct.bit==P.evaluate(3,12)


def test_three_dimensional_pairing_and_rational_point_subspace():
    a=(0,0,1,-7,6)
    S,P=computed(a)
    assert S.dimension==3 and P.certified and P.verify_evidence()
    assert P.rank==0==int(sage.pari.ellinit(list(a)).ellrank()[2])
    from elliptic_descent.quotient import row_basis
    points=((0,2),(1,0),(2,0))
    assert len(row_basis(S.kummer(Q) for Q in points))==3
    for Q in points:
        assert P.in_radical(S.kummer(Q))
    direct=pairing_entry(S,3,5,max_local_trials=20000)
    assert direct.complete and direct.verify() and direct.bit==0


@pytest.mark.parametrize('a', [(0,0,0,-1,1),(0,0,0,-1,0)])
def test_every_small_class_liftability_matches_radical_without_pruning(a):
    S,P=computed(a)
    for c in range(S.cardinality):
        C=S.quartic_cover(c)
        lifts=C.descend(proof=True)
        assert lifts.complete and lifts.certified and lifts.verify_evidence()
        assert bool(lifts.coverings)==P.in_radical(c)


def test_published_nonzero_classes_have_no_lifts_and_exhaust_nonradical_labels():
    a=(0,-1,1,-929,-10595)
    S,P=computed(a)
    quartics=((-64,-164,-52,68,-11),(-3,-52,-232,-60,-4),(-53,102,32,-78,-31))
    labels=set()
    b2=S.curves.working.b_invariants[0]
    for f in quartics:
        Q=FisherQuartic(change_parameters(f,((sp.Rational(1,2),0),(0,1))))
        value=substitute_generator(Q.z,-3,-b2/4)
        c=S.coordinates(value)
        labels.add(c)
        assert c and not P.in_radical(c)
        # Direct relative descent on the published small presentation. It does
        # not consult the pairing; labels are bound by the cubic invariant.
        result=relative_four_descent(f,proof=True)
        assert result.certified and result.verify_evidence() and not result.coverings
    assert labels=={c for c in range(S.cardinality) if not P.in_radical(c)}


def test_nonzero_pairing_is_natural_under_nonminimal_original_model_change():
    a=(0,-1,1,-929,-10595)
    S,P=computed(a)
    change=WeierstrassChange(S.curves.original,(sp.Rational(1,2),1,-1,3))
    other=two_selmer(change.target)
    Q=cassels_tate_pairing(other,max_local_trials=20000)
    assert Q.certified and Q.verify_evidence()
    for c,d in product(range(S.cardinality),repeat=2):
        assert P.evaluate(c,d)==Q.evaluate(S.transport_to(other,c,change),S.transport_to(other,d,change))
