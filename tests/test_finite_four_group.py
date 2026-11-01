"""Finite arithmetic checks; synthetic modules here are not Selmer computations."""
from dataclasses import replace
from itertools import product

import pytest

from elliptic_descent.finite_cohomology import (
    m4_add, m4_neg, m4_scale, FiniteFourModule, CohomologyClassTable,
)
from elliptic_descent.four_selmer_group import FiniteClassLaw, forced_class_law


def cyclic_module(n, generator_action=None):
    identity = tuple(range(16))
    action = identity if generator_action is None else generator_action
    actions = [identity]
    for _ in range(1,n):
        actions.append(tuple(actions[-1][v] for v in action))
    assert tuple(actions[-1][v] for v in action) == identity
    return FiniteFourModule(tuple(actions), (1,) if n>1 else (),
                            tuple(((i+1)%n,) if n>1 else () for i in range(n)))


@pytest.mark.parametrize('a,b', [(0,0),(3,1),(4,12),(7,10),(15,15)])
def test_mod_four_coordinates_do_not_carry_between_factors(a,b):
    assert m4_add(a,b) == (a%4+b%4)%4+4*((a//4+b//4)%4)
    assert m4_add(a,m4_neg(a)) == 0
    assert m4_scale(4,a) == 0 and m4_scale(-1,a) == m4_neg(a)


def test_cyclic_four_and_klein_four_have_distinct_laws_not_just_counts():
    M=cyclic_module(4)
    cocycles=tuple(tuple(m4_scale(s,k) for s in range(4)) for k in range(4))
    C=CohomologyClassTable(M,tuple((k,0) for k in range(4)),cocycles)
    assert M.verify() and C.verify()
    law=FiniteClassLaw(C.keys,tuple(tuple(C.add_indices(i,j) for j in range(4)) for i in range(4)),C.zero_index)
    assert law.verify() and law.decomposition[0] == (4,)
    assert law.order(1)==4 and law.table[1][1]==2
    klein=FiniteClassLaw(C.keys,tuple(tuple(i^j for j in range(4)) for i in range(4)),0)
    assert klein.verify() and klein.decomposition[0] == (2,2)
    assert not replace(C,cocycles=(cocycles[0],)*4).verify()
    assert not replace(C,cocycles=cocycles[:3],keys=C.keys[:3]).verify()


def test_nontrivial_galois_action_and_explicit_coboundary_witness():
    M=cyclic_module(2,tuple(m4_neg(i) for i in range(16)))
    assert M.verify()
    for t in range(16):
        d=M.coboundary(t)
        assert M.is_cocycle(d) and M.canonical(d)==(0,0)
    # All values at complex conjugation are cocycles, modulo the even vectors.
    C=CohomologyClassTable(M,((0,0),(0,1),(0,2),(0,3)),
                          ((0,0),(0,1),(0,4),(0,5)))
    assert C.verify()
    assert C.add_indices(1,1)==0 and C.add_indices(1,2)==3
    assert M.equivalence_witness((0,1),(0,3)) is not None
    assert M.equivalence_witness((0,1),(0,4)) is None
    assert not M.is_cocycle((1,0))
    with pytest.raises(ValueError): M.canonical((1,0))
    with pytest.raises(ValueError): M.coboundary(16)
    with pytest.raises(ValueError): C.add_indices(-1,0)


@pytest.mark.parametrize('mutation', ['action','transition','generators','identity'])
def test_finite_module_corruption_fails(mutation):
    M=cyclic_module(4)
    if mutation=='action':
        rows=list(M.actions); row=list(rows[1]);row[1],row[2]=row[2],row[1];rows[1]=tuple(row)
        M=replace(M,actions=tuple(rows))
    elif mutation=='transition':
        M=replace(M,transitions=((0,),)*4)
    elif mutation=='generators':
        M=replace(M,generators=())
    else:
        M=replace(M,actions=(tuple(reversed(range(16))),*M.actions[1:]))
    assert not M.verify()
    assert not replace(cyclic_module(1),transitions=((0,),)).verify()


@pytest.mark.parametrize('orders', [(),(2,),(4,),(2,2),(2,4),(4,4),(2,2,4)])
def test_actual_generators_and_coordinates_exhaust_every_element(orders):
    coordinates=tuple(product(*(range(n) for n in orders)))
    index={c:i for i,c in enumerate(coordinates)}
    table=tuple(tuple(index[tuple((a+b)%n for a,b,n in zip(c,d,orders))]
                      for d in coordinates) for c in coordinates)
    # Deliberately scramble the presentation and move its identity.
    permutation=tuple(reversed(range(len(coordinates))))
    inverse={v:i for i,v in enumerate(permutation)}
    moved=tuple(tuple(inverse[table[i][j]] for j in permutation) for i in permutation)
    law=FiniteClassLaw(tuple((0,i) for i in permutation),moved,inverse[0])
    assert law.verify()
    factors,gens,labels=law.decomposition
    assert factors==orders and len(labels)==len(coordinates)
    assert set(labels.values())==set(range(len(coordinates)))
    for c,i in labels.items():
        total=law.zero
        for coefficient,generator in zip(c,gens):
            total=law.table[total][law.multiply(coefficient,generator)]
        assert total==i
        assert law.multiply(4,i)==law.zero
        assert law.table[i][law.multiply(-1,i)]==law.zero
    for generator,order in zip(gens,factors):
        assert law.order(generator)==order
        assert law.multiply(order,generator)==law.zero


def test_unique_exact_sequence_shortcuts_and_nonunique_refusal():
    keys=((0,0),(0,1),(1,0),(1,1))
    method,law=forced_class_law(keys,keys[0],(0,0,1,1),(keys[0],keys[1]))
    assert method=='unique_order_four_extension' and law.decomposition[0]==(4,)
    # Entire group is the kernel: actual inclusion labels, deliberately permuted.
    method,law=forced_class_law(keys,keys[2],(0,0,0,0),(keys[2],keys[0],keys[3],keys[1]))
    assert method=='kernel_isomorphism' and law.decomposition[0]==(2,2)
    # Trivial kernel determines law via the actual projection.
    method,law=forced_class_law(keys,keys[1],(3,0,1,2),(keys[1],)*4)
    assert method=='projection_isomorphism' and law.decomposition[0]==(2,2)
    # For size 8 with nonzero kernel and image, count and doubling do NOT
    # label all cross-fibre sums. Require geometric evidence instead.
    keys=tuple((c,i) for c in range(4) for i in range(2))
    assert forced_class_law(keys,keys[0],tuple(c for c,i in keys),
                            (keys[0],keys[1],keys[0],keys[1])) is None


def test_invalid_finite_laws_and_sequence_data():
    law=FiniteClassLaw(((0,0),(0,1)),((0,1),(1,0)),0)
    assert not replace(law,table=((0,0),(1,0))).verify()
    assert not replace(law,zero=2).verify()
    assert not replace(law,keys=((0,0),(0,0))).verify()
    with pytest.raises(ValueError): law.multiply(1.0,0)
    with pytest.raises(ValueError): law.order(2)
    with pytest.raises(ValueError): forced_class_law(law.keys,law.keys[0],(0,0),(law.keys[0],))
    with pytest.raises(ValueError): forced_class_law(law.keys,law.keys[0],(0,1),(law.keys[0],)*3)
