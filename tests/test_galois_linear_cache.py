"""Exact fast automorphisms agree with independent Horner substitution.

Only rational linear algebra and repeated immutable field data are cached.
No field, covering, or cocycle verification verdict is cached.
"""
from dataclasses import replace
import json
from pathlib import Path
from random import Random

import pytest
import sympy as sp

from elliptic_descent.galois_covering import ExactGaloisField, evidence_from_record, polynomial_value

DATA = Path(__file__).parent / 'data'


def certificate():
    return evidence_from_record(json.loads((DATA / 'group_split_full.json').read_text())['galois_evidence'])


def test_direct_power_basis_and_cached_automorphisms_match_horner():
    field = certificate().torsion.field
    rng = Random(5060601)
    vectors = [(sp.S.Zero,)*field.degree,
               (sp.Rational(2, 7),)+(sp.S.Zero,)*(field.degree-1)]
    vectors += [tuple(sp.Rational(rng.randrange(-4, 5), rng.randrange(1, 6))
                      for _ in range(field.degree)) for _ in range(5)]
    values = []
    for vector in vectors:
        old = polynomial_value(tuple(field.rational(c) for c in vector),
                               field.generator, field.zero)
        value = field.unpack(vector)
        assert value == old and field.pack(value) == vector
        values.append(value)
        for i, image in enumerate(field.images):
            expected = polynomial_value(tuple(field.rational(c) for c in vector), image, field.zero)
            assert field.apply(i, value) == expected
            # Repeated use exercises cached columns, not just their construction.
            assert field.apply(i, value) == expected
    for i in range(field.degree):
        a, b = values[-2:]
        assert field.apply(i, a+b) == field.apply(i, a)+field.apply(i, b)
        assert field.apply(i, a*b) == field.apply(i, a)*field.apply(i, b)
    assert field.verify()
    for invalid in (-1, field.degree, True, 1.0):
        with pytest.raises(ValueError):
            field.apply(invalid, field.one)
    with pytest.raises(ValueError):
        field.unpack((1,))
    with pytest.raises(TypeError):
        field.unpack((1.0,)+(0,)*(field.degree-1))


def test_caches_do_not_turn_corrupt_field_or_torsion_data_into_evidence():
    evidence = certificate()
    assert evidence.verify()
    original = evidence.torsion
    field = original.field
    # Warm caches on the real field, then replace the public evidence.
    field.apply(1, field.generator**3)
    bad_images = (field.automorphism_images[0],)*field.degree
    bad_field = replace(field, automorphism_images=bad_images)
    assert not bad_field.verify()
    cover = evidence.coverings[0]
    old_data = cover._data(original)
    old_origin = cover.origin(original)
    assert cover._data(original) is old_data
    assert cover.origin(original) is old_origin
    equivalent = replace(original)
    assert cover._data(equivalent) == old_data
    assert cover._data(equivalent) is not old_data
    assert cover.origin(equivalent) == old_origin
    assert cover.origin(equivalent) is not old_origin
    wrong_curve = replace(original, ainvs=(0, 0, 0, -2, 0))
    assert not cover.verify(wrong_curve)
    # A fresh signed model cannot reuse the original's warmed origin/data.
    bad_cover = replace(cover, norm_witness=-cover.norm_witness)
    assert not bad_cover.verify(original)
    assert cover.verify(original)


def test_fast_path_retains_every_original_cocycle_and_sum():
    evidence = certificate()
    assert evidence.verify()
    table = evidence.class_table
    assert table.verify()
    # Independently recompute the action with scalar Horner substitution and
    # compare each complete cocycle, not just an abstract group cardinality.
    field = evidence.torsion.field
    for cover, cocycle in zip(evidence.coverings, table.cocycles):
        P = cover.origin(evidence.torsion)
        direct = []
        for image in field.images:
            point = tuple(polynomial_value(tuple(field.rational(c) for c in field.pack(v)),
                                           image, field.zero) for v in P)
            direct.append(evidence.torsion.point_indices[cover.evaluate(evidence.torsion, point)])
        assert tuple(direct) == cocycle


def test_generator_words_match_all_direct_field_and_torsion_actions():
    evidence = certificate()
    T, F = evidence.torsion, evidence.torsion.field
    assert F.verify()
    lookup = {v:i for i,v in enumerate(F.images)}
    # Both multiplication order and every conjugate are checked, independently
    # of the cached generator permutations used by the optimized implementation.
    expected = tuple(tuple(lookup[polynomial_value(
        tuple(F.rational(c) for c in F.automorphism_images[g]), image, F.zero)]
        for g in F.generators) for image in F.images)
    assert F.transitions == expected
    for s,image in enumerate(F.images):
        for i,P in enumerate(T.points):
            actual = None if P is None else tuple(polynomial_value(
                tuple(F.rational(c) for c in F.pack(v)),image,F.zero) for v in P)
            assert T.module.actions[s][i] == T.point_indices[actual]
    for C in evidence.coverings:
        assert C.cocycle(T) == C.cocycle_direct(T)


def test_generator_closure_rejects_a_distinct_false_nongenerator_image():
    evidence = certificate()
    F = evidence.torsion.field
    assert F.verify()
    nongenerator = next(i for i in range(1,F.degree) if i not in F.generators)
    images = list(F.automorphism_images)
    images[nongenerator] = F.pack(F.rational(123))
    bad = replace(F, automorphism_images=tuple(images))
    assert not bad.verify()
    assert not replace(F,generators=()).verify()
    assert not replace(F,generators=F.generators[:1]).verify()
    # A different generating order changes traversal words, not the named
    # automorphisms, their action or any complete covering cocycle.
    other = replace(F,generators=tuple(reversed(F.generators)))
    assert other.verify()
    T = replace(evidence.torsion,field=other)
    assert T.verify() and T.module.actions == evidence.torsion.module.actions
    for C in evidence.coverings:
        assert C.cocycle(T) == C.cocycle(evidence.torsion)


def test_nonabelian_generator_composition_has_the_correct_direction():
    # Construct Q(r,i), r^4=2, i^2=-1 in the independent eight-dimensional
    # tensor basis r^k*i^l. Its Galois group is dihedral, so left/right
    # multiplication cannot accidentally be interchanged as in an abelian test.
    def multiply(a,b):
        result = [sp.S.Zero]*8
        for j,x in enumerate(a):
            for k,y in enumerate(b):
                r,i = (j%4+k%4), (j//4+k//4)
                result[r%4+4*(i%2)] += x*y*2**(r//4)*(-1)**(i//2)
        return sp.Matrix(result)
    one,r,i = (sp.eye(8)[:,j] for j in (0,1,4))
    alpha = r+i
    columns = [one]
    for _ in range(8):
        columns.append(multiply(columns[-1],alpha))
    matrix = sp.Matrix.hstack(*columns[:8])
    inverse = matrix.inv()
    polynomial = tuple(-v for v in inverse*columns[8])+(sp.S.One,)
    ir = multiply(i,r)
    conjugates = (r+i, ir+i, -r+i, -ir+i, r-i, ir-i, -r-i, -ir-i)
    F = ExactGaloisField(polynomial,tuple(tuple(inverse*v) for v in conjugates),(1,4))
    assert F.verify()
    lookup = {v:j for j,v in enumerate(F.images)}
    direct = tuple(tuple(lookup[polynomial_value(
        tuple(F.rational(c) for c in F.automorphism_images[g]),image,F.zero)]
        for g in F.generators) for image in F.images)
    assert F.transitions == direct
    assert F.apply(1,F.images[4]) != F.apply(4,F.images[1])
    for s,word in enumerate(F._generator_action_data[1]):
        value = F.generator
        for j in word:
            value = F.apply(F.generators[j],value)
        assert value == F.images[s]


def test_subfield_inversion_matches_dense_field_arithmetic_and_rejects_zero():
    F = certificate().torsion.field
    values = [F.rational(1),F.rational(-7),F.generator,F.generator+1,
              F.generator**2-F.generator+2]
    values += [v for P in certificate().torsion.basis for v in (F.unpack(P[0]),F.unpack(P[1])) if v]
    for value in values:
        expected = F.one/value
        actual = F._subfield_inverse(value)
        assert actual == expected and value*actual == F.one
        assert F._subfield_inverse(value) is actual
    with pytest.raises(ZeroDivisionError):F.inverse(F.zero)
    with pytest.raises(ZeroDivisionError):F._subfield_inverse(F.zero)


@pytest.mark.parametrize('name',['group_split_full.json','group_zero_fibre.json'])
def test_denominator_free_origins_and_torsion_values_match_the_reference(name):
    record = json.loads((DATA/name).read_text())
    E = evidence_from_record(record.get('galois_evidence',record))
    T,F = E.torsion,E.torsion.field
    assert T.verify()
    for C in E.coverings:
        f,xi,alpha,r,u,A,s,beta,xia,delta,roots = C._data(T)
        j = tuple(u*u*e+F.rational(C.jacobian_to_original[1]) for e in T.root_values)
        assert beta == tuple(alpha+A/(v-s) for v in j)
        assert delta == tuple((s-v)*polynomial_value(xi,b,F.zero)/xia for v,b in zip(j,beta))
        nodes = (alpha,*beta)
        values = (F.zero,*(F.one/v for v in roots))
        original = [F.zero]*4
        for i,(node,value) in enumerate(zip(nodes,values)):
            if not value:continue
            p,denominator = [F.one],F.one
            for j,other in enumerate(nodes):
                if j==i:continue
                q = [F.zero]*(len(p)+1)
                for k,c in enumerate(p):q[k]-=other*c;q[k+1]+=c
                p=q;denominator*=node-other
            for k,c in enumerate(p):original[k]+=value*c/denominator
        P = C.origin(T)
        assert P == tuple(original)
        assert C.verify(T)
        for index in range(F.degree):
            point = tuple(F.apply(index,v) for v in P)
            assert C._torsion_value(T,point) == T.point_indices[C.evaluate(T,point)]
        assert C.cocycle(T) == C.cocycle_direct(T)
