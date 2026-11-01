"""Rational origin algebras and exact field-component flex maps."""
import os
from dataclasses import replace
import pytest

if os.environ.get('ELLIPTIC_DESCENT_REQUIRE_SAGE') == '1':
    import sage.all as sage
else:
    sage=pytest.importorskip('sage.all',reason='requires SageMath runtime')

from elliptic_descent import construct_four_cover, MarkedQuarticCover
from elliptic_descent.backends.sage_origin_torsors import origin_geometry, verify_field_flex

@pytest.mark.parametrize('f,xi', [
    ((4,0,-5,0,1),(1,0,0,0)),
    ((6,0,-5,0,1),(1,0,0,0)),
    ((1,-1,0,0,1),(0,-1,0,0)),
], ids=['split','two_quadratic','irreducible'])
def test_rational_origin_algebra_and_component_markings(f,xi):
    D=construct_four_cover(f,xi)
    A=origin_geometry(MarkedQuarticCover.canonical(f).compose(D))
    assert A.verify() and A.algebra.degree==16
    assert sum(g.degree() for g in A.algebra.factors)==16
    for g,K,root,point,B in A.flexes:
        assert verify_field_flex(A.covering,point,K,B)
        assert not verify_field_flex(A.covering,point,K,2*B+sage.identity_matrix(K,4))
    inverse=A.inverse_flex_matrix
    for g,K,root,point,B in A.flexes:
        actual=sage.matrix(K,4,4,lambda i,j:A.algebra.evaluate(inverse[i,j],root))
        assert actual*B==sage.identity_matrix(K,4)
