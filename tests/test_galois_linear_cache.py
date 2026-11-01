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

from elliptic_descent.galois_covering import evidence_from_record, polynomial_value

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
