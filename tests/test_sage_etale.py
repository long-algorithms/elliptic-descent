"""Real Sage integration for labelled product square classes; never mock arithmetic."""
import os
from random import Random

import pytest
import sympy as sp

if os.environ.get("ELLIPTIC_DESCENT_REQUIRE_SAGE") == "1":
    import sage.all as sage
else:
    sage = pytest.importorskip("sage.all", reason="requires SageMath runtime")

from elliptic_descent import EtaleQuarticAlgebra, relative_four_descent
from elliptic_descent.arithmetic import square_class
from elliptic_descent.backends.sage import sage_square_classes
from elliptic_descent.backends.sage_etale import sage_etale_square_classes
from elliptic_descent.geometry import _mul, _norm, construct_four_cover
from elliptic_descent.quotient import row_basis
from elliptic_descent.status import ProofStatus

# Includes isomorphic but differently presented quadratic factors and class-group 2-torsion.
FIXTURES = (
    (1, -1, 0, 0, 1),        # 4
    (-3, 4, -1, -3, 1),     # 3+1
    (6, 0, -5, 0, 1),       # 2+2
    (2, 0, -3, 0, 1),       # 2+1+1
    (4, 0, -5, 0, 1),       # 1+1+1+1
    (16, 0, -10, 0, 1),     # two copies of Q(sqrt(2))
    (-5, 0, 4, 0, 1),       # Q x Q x Q(sqrt(-5))
)


@pytest.fixture(scope="module", params=FIXTURES)
def context(request):
    return sage_etale_square_classes(request.param, support=(2,), proof=True)


def test_product_basis_norms_and_coordinate_roundtrips(context):
    A = context.algebra
    assert context.space.arithmetic_proof == ProofStatus.CERTIFIED
    assert not context.space.assumptions
    assert context.dimension == sum(context.factor_dimensions)
    assert len(context.offsets) == len(context.factors)+1
    assert context.space.f == A.f
    assert context.norm_fibre() == context.space.norm_fibre()
    for i, g in enumerate(context.space.generators):
        assert context.coordinates(g) == 1 << i and A.is_unit(g)
        assert context.coordinates(context.representative(1 << i)) == 1 << i
        assert square_class(_norm(g, A.f), context.space.support) == context.norm_columns[i]
        native = sp.S.One
        for factor, coefficients in zip(context.factors, A.to_components(g)):
            element = factor.as_field_element(coefficients)
            value = element if factor.degree == 1 else element.norm()
            native *= sp.Rational(str(value))
        assert native == A.norm(g) == _norm(g, A.f)
    rng = Random(882)
    for _ in range(10):
        a, b = (rng.randrange(1 << context.dimension) for _ in range(2))
        x, y = context.representative(a), context.representative(b)
        assert context.coordinates(x) == a
        assert context.coordinates(A.multiply(x, y)) == a ^ b
        assert context.norm_coordinates(a) == square_class(_norm(x, A.f), context.space.support)
        assert context.norm_coordinates(a ^ b) == context.norm_coordinates(a) ^ context.norm_coordinates(b)
    # Arithmetic norm constraints integrate with the existing bounded ledger;
    # no local solubility or relative-descent completeness is inferred here.
    batch = context.space.enumerate(limit=3)
    for candidate in batch.candidates:
        assert context.coordinates(candidate.xi) == candidate.label
        assert _norm(candidate.xi, A.f) == A.f[4]*candidate.norm_witness**2
        assert construct_four_cover(A.f, candidate.xi, norm_witness=candidate.norm_witness).check_construction()


def test_scalar_image_and_validation(context):
    A = context.algebra
    for c in (-1, 2, -8, sp.Rational(1, 4)):
        assert context.norm_coordinates(context.coordinates(A.scalar(c))) == 0
    assert context.representative(0) == A.one
    with pytest.raises(ValueError):
        context.coordinates(A.zero)
    with pytest.raises(ValueError):
        context.coordinates((1, 2, 3, 4, 5))
    for label in (-1, 1 << context.dimension, sp.Rational(1, 2)):
        with pytest.raises(ValueError):
            context.representative(label)
        with pytest.raises(ValueError):
            context.norm_coordinates(label)
    with pytest.raises(TypeError):
        context.representative(1.0)
    if len(A.factors) > 1:
        with pytest.raises(ValueError):
            context.coordinates(A.idempotents[0])
    # A prime unramified in the global polynomial gives odd factor valuations.
    with pytest.raises(ValueError):
        context.coordinates(A.scalar(101))


def test_fully_split_against_explicit_rational_square_classes():
    context = sage_etale_square_classes((4, 0, -5, 0, 1), support=())
    assert context.dimension == 4 and context.factor_dimensions == (1, 1, 1, 1)
    assert all(factor.field is sage.QQ for factor in context.factors)
    assert context.norm_fibre().cardinality == 8
    assert set(context.norm_fibre().elements()) == {mask for mask in range(16) if mask.bit_count() % 2 == 0}
    for label in range(16):
        parts = context.algebra.to_components(context.representative(label))
        assert parts == tuple(((-1)**((label >> i) & 1),) for i in range(4))
    batch = context.space.enumerate(limit=100)
    assert batch.norm_fibre_exhausted and batch.total_count == 8
    assert context.coordinates(context.algebra.scalar(-1)) == 15


def test_class_group_two_torsion_survives_in_a_product():
    context = sage_etale_square_classes((-5, 0, 4, 0, 1), support=(), proof=True)
    factor = next(factor for factor in context.factors if factor.degree == 2)
    assert factor.field.class_number(proof=True) == 2
    assert factor.dimension == 2
    assert context.dimension == 4
    assert any(abs(g.norm()) != 1 for g in factor.generators)
    embedded = context.algebra.embed_factor_unit(factor.index, (2, 0))
    assert context.coordinates(embedded) != 0
    assert context.algebra.norm(embedded) == 4


def test_isomorphic_quadratic_factors_keep_independent_bases():
    context = sage_etale_square_classes((16, 0, -10, 0, 1), support=())
    assert len(context.factors) == 2
    assert context.factors[0].polynomial != context.factors[1].polynomial
    assert context.factors[0].field.is_isomorphic(context.factors[1].field)
    labels = []
    for factor, offset in zip(context.factors, context.offsets):
        minus_one = context.algebra.embed_factor_unit(factor.index, (-1, 0))
        label = context.coordinates(minus_one)
        assert label != 0
        assert label >> offset < 1 << factor.dimension
        assert label % (1 << offset) == 0
        labels.append(label)
    assert labels[0] & labels[1] == 0
    assert context.coordinates(context.algebra.scalar(-1)) == labels[0] ^ labels[1]


def test_irreducible_adapter_agreement_without_assuming_same_basis():
    f = (1, -1, 0, 0, 1)
    old = sage_square_classes(f, proof=True)
    new = sage_etale_square_classes(f, proof=True)
    assert new.dimension == len(old.space.generators)
    images = tuple(new.coordinates(g) for g in old.space.generators)
    assert len(row_basis(images)) == new.dimension
    assert len(row_basis(old.coordinates(g) for g in new.space.generators)) == new.dimension
    def convert(label):
        result = 0
        for j, image in enumerate(images):
            if label >> j & 1:
                result ^= image
        return result
    assert {convert(v) for v in old.space.norm_fibre().elements()} == set(new.norm_fibre().elements())


def test_nonmonic_default_support_and_explicit_empty_support():
    context = sage_etale_square_classes((12, 0, -10, 0, 2))
    assert context.algebra.f[4] == 2 and context.space.f[4] == 2
    assert context.space.support == (2, 3)
    batch = context.space.enumerate(limit=4)
    assert batch.candidates
    for c in batch.candidates:
        assert _norm(c.xi, context.space.f) == 2*c.norm_witness**2
    small = sage_etale_square_classes((4, 0, -5, 0, 1), support=())
    default = sage_etale_square_classes((4, 0, -5, 0, 1))
    assert not small.space.support and default.space.support == (2, 3)
    assert default.dimension == 12


def test_proof_assumptions_and_explicit_settings():
    context = sage_etale_square_classes((2, 0, -3, 0, 1), support=(2,), proof=False)
    assert context.space.arithmetic_proof == ProofStatus.CONDITIONAL
    assert context.space.assumptions
    assert all(f.arithmetic_proof == (ProofStatus.CERTIFIED if f.degree == 1 else ProofStatus.CONDITIONAL)
               for f in context.factors)
    assert context.space.enumerate(limit=1).assumptions == context.space.assumptions
    split = sage_etale_square_classes((4, 0, -5, 0, 1), support=(), proof=False)
    assert split.space.arithmetic_proof == ProofStatus.CERTIFIED and not split.space.assumptions
    for invalid in (None, 0, 1, "false"):
        with pytest.raises(TypeError):
            sage_etale_square_classes((4, 0, -5, 0, 1), proof=invalid)
    for support in ((2, 2), (4,), (0,)):
        with pytest.raises(ValueError):
            sage_etale_square_classes((4, 0, -5, 0, 1), support)
    with sage.proof.WithProof('number_field', False):
        certified = sage_etale_square_classes((6, 0, -5, 0, 1), (), proof=True)
        assert certified.space.arithmetic_proof == ProofStatus.CERTIFIED


def test_rational_root_types_still_require_separate_descent_equivalences():
    # PR2 supports two quadratic factors, not the odd-degree-factor cases.
    for f in (FIXTURES[1], FIXTURES[3], FIXTURES[4], FIXTURES[6]):
        with pytest.raises(NotImplementedError):
            relative_four_descent(f)
