# Rational-root relative descent: v0.4 development step 3

This step adds the remaining degree-four factorization patterns `3+1`,
`2+1+1` and `1+1+1+1`. The same explicit norm-pair construction now supports
**every squarefree quartic of degree four over Q**, with its canonical marking.
The package version remains **0.3.0** until the separate v0.4 release/audit
step. This is not the full curve-level 4-Selmer group or an exact rank algorithm.

## The one-or-two orientation rule

Let `A=Q[t]/(f)` and let `e` be the CRT idempotent of a rational linear factor.
Then `v=1-2e` has `v^2=1` and `N_A(v)=-1`: its component at that factor is -1
and every other component is 1. The norm-pair equivalence

\[
 (\xi,r)\sim(c v^2\xi,c^2N_A(v)r)
\]

with `c=1` therefore identifies `(xi,r)` with `(xi,-r)` **over the identity
of the target quartic**. Multiplication by v gives a rational invertible
coordinate transformation of the quadric intersection. It fixes the quadratic
map coordinates and reverses the product-norm coordinate. The corresponding
diagram commutes; we do not identify different target points by fiat.

`OrientationCertificate` records the exact reversing multiplier. Its verifier
checks the labelled factorization, `v^2=1`, and `N_A(v)=-1`. The
`identification(model)` method returns a `ScalarSquareEquivalence`, whose
verifier checks the actual quadratic and weighted-map identities against the
model with the opposite norm witness. `verify_evidence()` audits this
identification for each returned locally soluble representative by default.

For types `4` and `2+2`, no such stabilizer exists, by the lemma in
[two-quadratic-descent.md](two-quadratic-descent.md). The output rule is:

| Factorization | Fixed-target orientations per accepted scalar-square class |
|---|---:|
| `4`, `2+2` | 2 |
| `3+1`, `2+1+1`, `1+1+1+1` | 1 |

Thus the two output lists `coverings` and `unoriented_coverings` have the same
length for rational-root inputs. There is no extra factor-permutation quotient,
no division by the number of rational roots, and no quotient by the rational
2- or 4-torsion subgroup. Flipping any one rational component supplies an
identification; flipping two such components does not reverse the orientation.
Choosing the first labelled rational factor makes the evidence reproducible.

## A unified construction, not an incomplete shortcut

The factorwise backend computes the full `A(S,2)`, including class-group
contributions. We impose the **single total norm condition** and quotient by
**diagonal** rational scalars, as in the preceding two-quadratic step. The
support uses the entire quartic discriminant, including inter-factor resultant
primes. The finite-support and scalar-image arguments work in this product
algebra just as in the other squarefree cases.

After the exact orientation identification, every quotient representative is
constructed and tested locally. The geometry and finite-place code are shared
with the earlier implementation. Complete mode has no artificial enumeration,
node, or depth limit; explicitly limited runs can remain partial. The
norm-pair parametrization, finite-support bound, quotient and local decision
arguments therefore give exhaustiveness in the same fixed-target equivalence.

The rational branch point identifies the source C with its Jacobian E, so the
problem can also be interpreted as ordinary 2-descent on E. This PR **does not**
substitute a returned 2-Selmer basis for all its elements or silently change the
marked map. It keeps the unified construction and adds an independently
verified identification for testing and future optimized dispatch. A faster
ordinary-2-descent backend can later use this identification without changing
the mathematical contract.

## Exact identification with the marked Jacobian

For a chosen rational root alpha write

\[
 f(\alpha+h)=A h+B h^2+C h^3+D h^4,
\]

where `A=f'(alpha) != 0` and `D=lc(f)`. The reciprocal substitution

\[
 u=\frac{A}{x-\alpha},\qquad v=\frac{A y}{(x-\alpha)^2}
\]

gives the derivative-based Weierstrass model

\[
 v^2=u^3+B u^2+ACu+A^2D.
\]

For the canonical Jacobian convention used by `QuarticTwoCover`, set
`s=(B-f_2)/3`. The isomorphism is

\[
 \iota(X:Y:Z)=
 \left(s+\frac{AZ}{X-\alpha Z},
       -\frac{AY}{(X-\alpha Z)^2}\right).
\]

The minus sign is significant. With that sign the existing canonical map
satisfies **`pi_C = [2] o iota`**, not its negative. The selected branch point
`(alpha,0,1)` maps to the elliptic origin. On affine Jacobian points `(u,v)`,
the inverse in weighted quartic coordinates is

\[
 (X:Y:Z)=(\alpha(u-s)+A:\;-Av:\;u-s).
\]

The inverse includes the two possible quartic points at infinity. Its affine
form has no basepoint because `Z=0` forces `X=A != 0`; its extension at the
elliptic origin is the chosen branch point. Smoothness gives the unique global
extensions of these inverse rational maps.

`RationalRootIdentification.check_marked_map()` checks the derivative model,
its translation to the canonical Jacobian, the inverse equation, and the
identity with the duplication formula in the exact function field. The test
suite covers every rational root of each fixture, not just the chosen root,
as well as point round trips, infinity, weighted scaling, and rational
4-torsion. These checks are stronger than equality of j-invariants.

```python
from elliptic_descent import relative_four_descent, RationalRootIdentification

result = relative_four_descent((4, 0, -5, 0, 1), proof=True)
assert result.certified and result.orientation.multiplicity == 1
assert result.verify_evidence()
identification = RationalRootIdentification(result.source, root=2)
assert identification.check_marked_map()
assert identification.to_jacobian((2, 0, 1)) is None
P = (8, 12)  # On the canonical Jacobian y^2=(x-5)(x-4)(x+4).
assert identification.to_jacobian(identification.from_jacobian(P)) == P
assert result.source.evaluate(identification.from_jacobian(P)) == (5, 0)
```

## Real-root isolation repair

SymPy's isolating intervals can be open and touch another, rational root. For
example, `(t^2-2)(t^2-1)` can have the interval `(-2,-1)` for `-sqrt(2)` directly
before the singleton `(-1,-1)`. Previously the real solver attempted to refine
the singleton away from the preceding endpoint and raised an error.

The solver now refines each non-singleton interval until both endpoints are
nonroots of the full polynomial, in addition to checking the required sign.
The stored closed intervals are disjoint; the next rational singleton needs
no impossible refinement. The independent verifier checks exact rational
endpoints and root counts, and rejects an interval that includes two roots.
Regression tests exercise both signs, mixed rational/irrational roots,
nonintegral roots and a deliberately malformed interval.

## Independent checks and limitations

Every source here has a rational point, namely a branch point. Its locally
soluble relative 2-coverings are consequently counted by `Sel^(2)(E/Q)`.
The Sage integration suite compares the **complete fixed-target count** to
`2**E.mwrank_curve().selmer_rank()`, using the derivative-based cubic and an
integral minimal model for the independent eclib calculation. This is a
cardinality comparison, not a class-by-class comparison of covering models.

The benchmark corpus covers a zero branch point, a translated branch point,
one-real-root and three-real-root cubic factors, real and imaginary quadratic
factors, fully split curves with rational 4-torsion, negative leading
coefficients, and rational coefficients. No live Magma execution is claimed.
All old irreducible and two-quadratic regression suites remain enabled.

A certified empty descent cannot occur for a rational-root source: the trivial
relative 2-covering exists. Empty *partial* output is still possible when no
candidates have been processed. Tests explicitly distinguish those situations.
No new empty Sha example is asserted; the existing irreducible 571a1 fixtures
continue to exercise that separate case.

`proof=False` remains conditional when a nonrational factor relies on
uncertified class-group arithmetic. Fully split rational arithmetic carries
no such dependency and can still be unconditional with that flag. Manual
result construction is not a certificate. The evidence verifier does not redo
class-group computations independently.

The input remains a genuine quartic: five exact coefficients with nonzero
leading coefficient and nonzero discriminant. A degree-three polynomial with
a branch point at infinity must not be silently interpreted as a degree-four
algebra. Such presentation conversion, arbitrary externally supplied markings,
full 4-Selmer-group assembly, general minimization/reduction and large-height
point searches are not added by this PR.

## Reproduction and performance records

```sh
sage -pip install -e '.[test]'
ELLIPTIC_DESCENT_REQUIRE_SAGE=1 sage -python -m pytest -q
sage -python benchmarks/rational_root.py --verify --repeats 3
```

The benchmark emits JSONL with proof assumptions, orientation evidence,
root and derivative-model data, support, candidate counts, local work,
model changes, separate solver/verification times and cumulative process
memory. Verification includes the marked-map identity and the independent
eclib count. Repetitions share a process and are not cold-cache timings.
No performance advantage over other systems is inferred from CI measurements.

## References

- Fisher, *Some improvements to 4-descent on an elliptic curve* (2008), and
  Fisher--Sills, *Local solubility and height bounds for coverings of elliptic
  curves* (2012): [core references](../README.md#references).
- [Magma descent documentation](https://magma.maths.usyd.edu.au/magma/handbook/elliptic_curves_over_q_and_number_fields): rational-root quartics reduce to ordinary 2-descent; this implementation retains explicit maps to the original source.
- [Sage/eclib Selmer-rank interface](https://doc.sagemath.org/html/en/reference/libs/sage/libs/eclib/interface.html#sage.libs.eclib.interface.mwrank_EllipticCurve.selmer_rank): independent cardinality check, not a Mordell--Weil-rank substitution.
- [Broader descent bibliography](descent-references.md).
