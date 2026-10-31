# Quartic etale arithmetic: v0.4 preparation, PR 1

**Update:** the subsequent [two-quadratic descent step](two-quadratic-descent.md)
now integrates this backend for type `(2,2)`. The original PR1 scope and tests
are described below; rational-root types still lack a complete descent API.

This change supplies the arithmetic foundation for all five squarefree quartic
factorization patterns. It does **not** extend `relative_four_descent` to
reducible quartics, decide orientation equivalences, or complete v0.4.
The package version remains 0.3.0. The existing irreducible descent and its
regression suite are unchanged.

## Exact decomposition and representations

For a squarefree quartic `f`, write its monic normalization as a product of
pairwise coprime monic irreducibles `g_i`. Then

\[
 A=\mathbb Q[t]/(f)\simeq\prod_i K_i,\qquad K_i=\mathbb Q[t]/(g_i).
\]

`EtaleQuarticAlgebra` retains the **original** geometric polynomial `f`,
including its leading coefficient. Normalizing the algebra presentation never
changes the curve to `y^2=f(x)/lc(f)`.

Global elements are four exact rational coefficients in ascending powers of
`t`. Component elements have exactly `degree(g_i)` ascending coefficients.
For a linear factor, that single coefficient is the rational value at its
root, not a formal field generator. Components are ordered by increasing
degree, then lexicographically by their monic ascending coefficients.

Distinct polynomial factors stay distinct even if the fields they define are
isomorphic. In particular `(t^2-2)(t^2-8)` has two separately labelled factors;
merging them would corrupt the algebra and its square-class space.

The Chinese remainder maps use idempotents

\[
 e_i=(f_{\mathrm{monic}}/g_i)\,
       ((f_{\mathrm{monic}}/g_i)^{-1}\bmod g_i)\pmod {f_{\mathrm{monic}}}.
\]

They satisfy `e_i^2=e_i`, `e_i e_j=0` for `i != j`, and `sum(e_i)=1`.
`to_components` reduces modulo each `g_i`; `from_components` reconstructs with
the idempotents. `verify_decomposition` checks factorization, idempotents,
and both compositions on rational bases. Multiplication uses polynomial
remainders, independently of the geometric kernel's fixed-degree reduction.
Norms are products of component resultants, independently of the geometric
kernel's global multiplication-matrix determinant.

Zero divisors are legitimate algebra elements, but they cannot be inverted or
used as multiplicative square classes. A component generator is embedded as
`(1,...,g,...,1)`, equivalently `1+e_i*(g-1)`, **not** as `e_i*g`.

```python
from elliptic_descent import EtaleQuarticAlgebra

A = EtaleQuarticAlgebra((6, 0, -5, 0, 1))  # (t^2-2)(t^2-3)
assert A.factorization_type == (2, 2)
assert A.verify_decomposition()
u = (1, 2, 3, 4)
assert A.from_components(A.to_components(u)) == u
assert A.multiply(u, A.inverse(u)) == A.one
```

## Factorwise square classes and norms

Run the following with `sage -python`:

```python
from elliptic_descent.backends.sage_etale import sage_etale_square_classes

context = sage_etale_square_classes((6, 0, -5, 0, 1), proof=True)
print(context.space.support, context.factor_dimensions)
for j, generator in enumerate(context.space.generators):
    assert context.coordinates(generator) == 1 << j
    assert context.coordinates(context.representative(1 << j)) == 1 << j
batch = context.space.enumerate(limit=16)
print(batch.total_count, batch.norm_fibre_exhausted)
```

For rational prime support `S`, the new backend computes

\[
 A(S,2)=\bigoplus_i K_i(S_i,2),\qquad
 S_i=\{\mathfrak p:\mathfrak p\mid p,\ p\in S\}.
\]

Each number-field factor uses Sage's `selmer_space` with every prime over S,
including the class-group contribution. Rational factors use
`QQ.selmer_space`; they are not artificial degree-one number fields.

Bit-vector coordinates concatenate the individual factor bases in the same
order as the polynomial factors. `offsets` records every block boundary.
Each embedded generator is checked using three norms: the native Sage factor
norm, the CRT/resultant norm, and the existing global determinant norm.

The total norm is the product of factor norms, or XOR of their rational
square-class coordinate vectors. `norm_columns`, `norm_coordinates(label)`
and `norm_fibre()` expose this map and its affine fibre over `[lc(f)]`.
`representative(label)` returns an actual unit in the original quartic power
basis. Only the coordinate round trip is an equality of labels; replacing an
arbitrary element by the chosen representative is equality **modulo squares**.

`support=None` uses the existing conservative quartic support, including primes
of inter-factor resultants through the full quartic discriminant. `support=()`
means empty support and is not replaced by the default. Supplied support is
validated as a set of distinct rational primes. An odd valuation outside it is
rejected, not silently accepted or silently added to S.

## Proof status and boundary of the result

`proof=True` is passed explicitly to every number-field calculation. The
certified Sage class-group/unit computation remains a trusted dependency;
round-trip checks do not independently certify its exhaustiveness.

`proof=False` propagates a GRH assumption when a nonrational factor is present.
A fully split algebra has only rational arithmetic, so it remains unconditional
with either flag. Factor-level proof statuses are retained and aggregated.

The existing `FiniteNormSpace` and `CandidateBatch` types are reused. Their
exhaustion flag concerns only the specified finite norm fibre. No local
solubility, rational-scalar quotient, stabilizer, or orientation identification
is computed by this new backend. In particular, a finite norm fibre is not a
complete relative Selmer set.

The irreducible-only `sage_square_classes` and `relative_four_descent` functions
retain their previous contracts. Tests ensure that reducible input still raises
`NotImplementedError` at the complete-descent boundary. Subsequent PRs will
integrate this arithmetic with the necessary equivalence and local arguments.

## Verification and benchmarks

The new tests cover all five factorization patterns, nonmonic and rational
coefficients, a rational root at zero, isomorphic quadratic factors, component
norms, zero divisors, invalid shapes and inexact inputs. Sage tests additionally
cover class-group 2-torsion in `Q(sqrt(-5))`, explicit fully split enumeration,
proof-setting propagation, out-of-support valuations, and agreement with the
old irreducible backend **without assuming the same choice of basis**.

```sh
PYTHONPATH=src python -m pytest -q
ELLIPTIC_DESCENT_REQUIRE_SAGE=1 sage -python -m pytest -q
sage -python benchmarks/etale.py --sage --verify --repeats 3
```

The dedicated Sage CI job requires these tests: absence of Sage is an error,
not a skipped integration check. Python-only jobs explicitly skip that module.
The JSONL benchmark corpus covers all splitting types and extra arithmetic
edge cases. It separates CRT construction, factorwise arithmetic, and audit
costs, and records support, dimensions, norms, proof assumptions, environment
and cumulative process memory. Repetitions share a process; they are not
cold-cache timings or comparisons against Magma.

## References

- [Sage number-field Selmer spaces](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html#sage.rings.number_field.number_field.NumberField_generic.selmer_space).
- [Sage rational-field Selmer spaces](https://doc.sagemath.org/html/en/reference/rings_standard/sage/rings/rational_field.html#sage.rings.rational_field.RationalField.selmer_space).
- [Sage's p-Selmer group implementation and class-group contribution](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/selmer_group.html).
- Schaefer--Stoll, *How to do a p-descent on an elliptic curve* (2004), and the
  [general descent bibliography](descent-references.md).
