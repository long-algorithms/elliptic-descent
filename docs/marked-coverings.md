# Marked curves and coverings — v0.5 development step 1

This component is included in [v0.5.0](v0.5.md). These interfaces preserve the user's elliptic
curve and exact covering maps while using the existing relative-descent engine.
They do not compute a labelled 2-Selmer group, a pairing matrix, or a full
4-Selmer group. See [the eight-step plan](v0.5-plan.md).

## Exact Weierstrass changes and directions

`WeierstrassModel` stores the five rational coefficients of

    y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6.

Two coefficients abbreviate (a4,a6). Inexact and singular inputs are rejected.
Points are affine rational pairs; `None` denotes the origin. Point validation,
negation and chord/tangent addition handle infinity, inverse pairs and vertical
tangents. These are interface utilities, not a replacement for Sage arithmetic.

`WeierstrassChange(source, (u,r,s,t))` has direction **source -> target**:

    x' = (x-r)/u^2,
    y' = (y-s*(x-r)-t)/u^3,       u != 0.

The target coefficients are derived explicitly, not inferred from invariants.
This is the convention of Sage's `WeierstrassIsomorphism`, not necessarily the
direction suggested by every method named 'change model'. `verify()` checks
the polynomial identity and discriminant scaling. Inverse and composition
(`first.then(second)` means `second o first`) keep their endpoints explicit.
Weighted elliptic coordinates use (2,3,1); quartic double-cover coordinates
use (1,2,1). The API never treats these two triples as interchangeable.

The optional Sage `rational_curve_context(E)` returns the exact original model,
a minimal working model, and maps in both directions. It accepts rational
coefficients, `WeierstrassModel`, or an actual Sage curve over QQ. An extension
or finite field is rejected, not coerced. The chosen Sage transformation is
independently verified using the polynomial formulas. `minimal=False` preserves
the original presentation with identity maps.

## What a marking means

A `MarkedQuarticCover` consists of a `QuarticTwoCover`, an explicit
`WeierstrassChange` from its canonical Jacobian to the target E, and an optional
rational point T on E. Its map is

    C --pi_can--> Jac(C) --isomorphism--> E --translation by T--> E.

Every stage is retained. Translation changes the covering class in general;
it cannot be dropped during normalization. The construction is a 2-covering
by construction, not because a random rational function happens to land on E.
`verify_supplied_map` compares both coordinate functions with this complete
marking. A constant map, the degree-one identity on a rational-root model, or
a map with the wrong translation is rejected even if its target equation holds.

Use `retarget(change)` to postcompose a target-model isomorphism. It also
transforms T. Use `compose(D)` to retain the degree-sixteen map of a relative
4-covering D to the original E. `evaluate` follows the proper maps and exact
group law at exceptional points. Generic `map_xy` expressions have poles and
removable indeterminacies and should not be substituted naively at such points.

Functions on the source are reduced exactly to A(x)+B(x)*y using y^2=f(x).
This bounds expression growth when composing translations. The multivariate
4-cover composite is kept factored, rather than forcibly expanding a very large
fraction. Its verifier checks the exact identities at each stage; the retained
geometric construction verifier supplies its degree and map-to-source proof.
It is not a new independent class-group or local-solubility test.

## Cubic presentations, including the branch point at infinity

`CubicToQuartic` accepts exactly four ascending coefficients of a squarefree
cubic f. Choose a rational s with f(s)!=0. The default tries 0,1,2,3, one of
which must work. The substitution

    x = s + 1/t,       y = w/t^2

produces the genuine degree-four polynomial g(t)=t^4*f(s+1/t), whose constant
term is zero and whose leading coefficient is f(s). This is **not** padding a
cubic with a zero quartic coefficient and using an incorrect degree-four algebra.
On weighted coordinates its inverse and forward transformations are

    quartic -> cubic: (U:W:V) -> (s*U+V : W : U),
    cubic -> quartic: (X:Y:Z) -> (Z : Y : X-s*Z).

They are invertible projective transformations. In particular, the cubic branch
point (1:0:0) goes to (0:0:1), and points with affine x=s go to quartic points
at infinity. Both are handled without dividing by zero. `verify()` checks the
binary-form identity and inverse matrices exactly.

`transport_map` only precomposes rational functions; it does not certify an
arbitrary supplied covering. `marked_map` additionally requires an explicit
canonical-Jacobian isomorphism and translation, and checks equality with their
known 2-covering. `mark_pari_two_cover(record, change, translation=..., shift=...)`
uses the same checks and returns both the marking and any cubic-conversion
object. It also checks the exact target against the imported ainvs. It does
not certify PARI's basis, assign Selmer coordinates, or guess the marking.
The tests use synthetic imported records with known maps; no claim is made that
automatic marking recovery for all `ell2cover` outputs is implemented.

## The zero class is multiplication by two

`IdentityTwoCover(E)` constructs a quartic presentation of [2]:E->E. Complete
square with Y=2y+a1*x+a3 to obtain

    Y^2 = 4*x^3 + b2*x^2 + 2*b4*x + b6.

Convert this cubic using s. If J is the canonical Jacobian of the resulting
quartic, the explicit change J -> original E has parameters

    (u,r,s_parameter,t) = (-2, -4*s, -a1, -4*a3).

The sign u=-2 is intentional. `verify()` checks the original equation and
compares the complete map with the independent duplication formula. The
identification map from the quartic to E has degree one; the covering map is
its composition with [2]. These must not be confused when representing the
zero element of the 2-Selmer group.

```python
from elliptic_descent import WeierstrassModel, IdentityTwoCover

E = WeierstrassModel((-1, 1))
zero = IdentityTwoCover(E)
assert zero.verify()
assert zero.from_quartic(zero.to_quartic(None)) is None
P = (0, 1)
assert zero.cover.evaluate(zero.to_quartic(P)) == E.double(P)
```

With Sage available, a marked relative calculation is explicit:

```python
result = zero.cover.descend(proof=True)
assert result.target == E
assert result.certified and result.verify_evidence()
for D in result.coverings:
    assert D.target == E and D.verify()
```

This computes a relative descent above one marked 2-cover, **not** the entire
4-Selmer group of E. The wrapper retains the raw result, its proof assumptions,
partial decisions and unresolved cases. Its existence does not implement v0.5
steps 2--7. It does not identify global duplicates under translations or torsion.

## Checks and provenance

The Python suite covers nonzero a1/a3, rational/nonminimal changes, inverse
composition diagrams, weighted scaling, zero-class doubling, branch points,
points at infinity, rational 2-/4-torsion, translations and wrong-map rejection.
Seeded exact examples exercise multiple rational presentations. A Sage suite
compares transformation tuples, discriminants, minimal models, and point
addition/doubling with Sage, and checks marked relative outputs against an
independent eclib count for one rational-point source. It is not a separate
independent 4-descent implementation.

`benchmarks/marked_models.py` times model/marking construction separately from
exact verification. Its scope is interface arithmetic, not full Selmer descent.
Optional Sage checks test normalization and point images. All older release
benchmarks remain in CI, and distribution tests ensure new modules ship in
both wheel and sdist. Release-wide acceptance is recorded in [v0.5.0](v0.5.md).

## Primary references

- [Sage Weierstrass isomorphism conventions and formulas](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/weierstrass_morphism.html).
- [Sage rational elliptic-curve arithmetic](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_rational_field.html).
- [PARI elliptic-curve covering interface](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html).
- [Existing canonical quartic construction](algorithm.md) and [core sources](../README.md#references), including Fisher--Sills for the canonical covering map.

The cubic conversion and identity-cover formulas above are verified directly
in code; they do not rely on an unrecorded heuristic isomorphism computation.
