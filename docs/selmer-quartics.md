# Marked quartics for labelled classes — v0.5 development step 3

`S.quartic_cover(c)` now constructs a binary-quartic 2-covering for **any vector**
of a complete labelled `two_selmer` result: zero, a basis vector, or any sum.
The arithmetic label, inverse presentation maps and exact map to the user's
original elliptic curve are retained. The package remains **0.4.0**; this step
does not compute the Cassels–Tate pairing or assemble the global 4-Selmer group.

## From the labelled auxiliary model to a quartic

Keep the conventions of [cubic 2-Selmer arithmetic](two-selmer.md). For
`L=Q[theta]/h`, a unit delta with norm r squared gives

\[
 \delta(u_0+u_1\theta+u_2\theta^2)^2=Q_0+Q_1\theta+Q_2\theta^2,
 \qquad D_\delta:\ Q_2=0,\quad Q_1+z^2=0.
\]

The conic `Q2=0` is nonsingular. Indeed the coefficient-of-theta-squared
functional on L is `a -> Tr(a/h'(theta))`, so Q2 is the trace quadratic form
with invertible multiplier `delta/h'(theta)`. Separability and invertibility
make this pairing nondegenerate, for fields and products alike.

A local point on D projects to a nonzero point on this conic: if all three u
coordinates vanished then the second equation would force z=0. Thus a Selmer
class gives a conic soluble at every completion. Hasse–Minkowski supplies a
rational conic point, even when the genus-one covering has no rational point.
No search for a rational point on that covering is required.

The adapter uses an immediate coordinate point when available, otherwise the
Sage interface to PARI's complete ternary quadratic-form solver `qfsolve` [1,2].
The returned point is checked exactly. A backend failure raises an exception;
it does not silently omit a class or certify an empty descent. An explicit
rational conic point may also be supplied for reproducibility or testing.

For a point P on `q=Q2`, let H be its Hessian, choose coordinate vectors V,W
completing P to a basis, and put `d=XV+ZW`. The three binary quadratics

\[
 p(X,Z)=q(d)P-(P^T H d)d
\]

satisfy `q(p)=0`. Their coefficient matrix M in `(X^2,XZ,Z^2)` is checked to be
invertible over Q. Hence this is an isomorphism of P^1 with the conic, not a
possibly multiple parametrization or a rational map with an omitted basepoint.
Common rational content is removed from the three forms together; the resulting
quartic and all maps are derived from those actual forms, so no scalar twist
is silently discarded.

Set

\[
 F(X,Z)=-Q_1(p(X,Z)),\qquad C_\delta:Y^2=F(X,Z).
\]

The explicit isomorphism to the auxiliary curve is

\[
 C_\delta\longrightarrow D_\delta,\qquad
 (X:Y:Z)\longmapsto(p_0:p_1:p_2:Y).
\]

All four coordinates on the right have weight two on the source. There is no
common zero because M is invertible. For the inverse, write
`M^-1*(u0,u1,u2)^T=(a,b,c)^T`, so `ac=b^2`. Use

\[
 (X:Y:Z)=(a:az:b)\quad\text{when }a\ne0,
 \qquad (b:cz:c)\quad\text{otherwise}.
\]

The second chart includes the tangent parameter and points at infinity; when
a=0, necessarily b=0 and c is nonzero. The Y scaling is essential for inverse
weighted-projective coordinates. Both charts and their scaling laws are tested.
Smoothness follows from the explicit isomorphism to the smooth auxiliary curve.
In particular the binary branch quartic is squarefree.

The public quartic constructor requires a nonzero degree-four leading
coefficient. If `F(1,0)=0`, move infinity by the invertible substitution
`(X,Z) -> (kX+Z,X)`. At least one integer k among 0,1,2,3,4 is not a zero of
the nonzero polynomial `F(k,1)`. This finite choice is an exact presentation
change, not a candidate cutoff or a deletion of a branch point.

## The arithmetic class is retained, not inferred from invariants

Let `A=Q0(p)` and `N=r*N_L(p0+p1*theta+p2*theta^2)`. On the quartic the raw map
to the working elliptic curve is

\[
 x_E=A/Y^2,\qquad
 y_E=N/Y^3-(a_1x_E+a_3)/2.
\]

The three exact coefficient identities give

\[
 \delta\left(\frac{p_0+p_1\theta+p_2\theta^2}{Y}\right)^2=x_E-\theta.
\]

More strongly, the explicit inverse presentation maps identify this entire
marked covering with the already labelled auxiliary model. This is the class
certificate; having the same Jacobian or merely the correct number of models
would not establish it.

`LabelledQuarticCover` records the parent Selmer space, the requested coordinates
in its chosen basis, the ambient label, and the auxiliary arithmetic element.
The verifier checks the label against the backend coordinate map and the unique
accepted local decision for that class. It checks the conic and covering maps,
not merely a copied label field. `S.representative(c)` retains its previous
meaning: a cubic-algebra representative. It is not replaced by this new model
interface.

## Recovering the exact signed marking

The canonical map of `QuarticTwoCover` has coordinates
`x_can=P4/(4Y^2)`, `y_can=P6/(8Y^3)`. Solve the full coefficient identities

\[
 P_4/4=k A+r_0F,\qquad P_6/8=\ell N.
\]

They determine `k=u^2` and `ell=u^3`. Taking `u=ell/k` and checking both powers
fixes the sign without a square-root choice. The resulting Weierstrass change
from the canonical Jacobian to the working E is

\[
 (u,r_0,s,t)=(u,r_0,a_1u/2,a_3u^3/2).
\]

Its target coefficients and polynomial identity are checked, followed by equality
of both complete covering-map functions. Both canonical and constructed maps
send every branch point to O, so their target translation is zero. This is a
property of this construction, not a licence to discard translations from
arbitrary imported maps. The final marking composes with the exact working-to-
original model change from step 1. Nonminimal models and nonzero a1/a3 are kept.

For the zero class the auxiliary point `(1:0:0:0)` lies over O. Its quartic
presentation represents multiplication by two after identification with E,
not the identity degree-one morphism. Tests compare it with doubling at finite,
infinite, 2-torsion and 4-torsion points. No extra torsion quotient is applied.

## API and evidence

```python
from elliptic_descent import two_selmer

S = two_selmer((0, -5, 0, -16, 80), proof=True)
assert S.certified
C = S.quartic_cover(S.add(1, 2))  # The sum, not just a selected basis element.
assert C.certified and C.verify()
assert S.coordinates(C.delta) == C.coordinates
assert C.ambient_label == S.ambient_label(C.coordinates)
assert C.marking.target == S.curves.original
print(C.quartic.f)               # Five ascending rational coefficients.
print(C.marking.map_xy)          # Map to the user's original E.
```

The parent space must be complete. Partial computations cannot manufacture a
full basis by constructing models for a few discovered classes. A complete
conditional parent may produce a fully checked exact model, but its assumptions
remain attached and `C.certified` remains false. Fully rational arithmetic does
not acquire a spurious GRH assumption.

`C.presentation.to_auxiliary(point)` and `from_auxiliary(point)` are the exact
inverse projective maps; `C.evaluate(point)` uses the original marked target.
`conic_point=` and `parameter_change=` on `quartic_cover` permit alternative
verified presentations. A different parameterization changes equations but not
the ambient class or target marking. Changing the Selmer basis changes coordinate
integers according to the explicitly recorded basis, not the underlying class.

`C.descend(...)` connects directly to the existing marked relative solver and
returns `LabelledRelativeLifts`. This wrapper retains the intermediate class and
combines arithmetic assumptions from both stages. In particular, an unconditional
relative computation does not retroactively certify a conditional parent Selmer
calculation. Explicit relative budgets still return partial results. The returned
coverings are classes of diagrams over this fixed intermediate quartic; they
are **not** deduplicated global 4-Selmer classes.

## Tests, benchmarks and trust boundary

The new Sage corpus constructs **all 19 classes** across the seven existing
2-Selmer fixtures, not just their bases. It includes the zero-dimensional group,
no/cyclic/full rational 2-torsion, rational 4-torsion, and all three nonzero
classes of 571a1. The latter exercise conic parametrization without a rational
point on the genus-one curve. Class sums, alternate bases, alternate conic
points, GL2 parameter changes and original-model transport are checked.

Python tests verify the inverse Veronese charts, arithmetic x-theta identity,
zero-class doubling, branch and infinite points, rational scaling, and malformed
inputs. Sage tests retain conditional assumptions, reject altered labels and
models, and exercise the path to relative descent. Existing local, arithmetic,
model, descent, packaging and minimum-dependency tests remain enabled.

```sh
sage -python -m pytest -q
sage -python benchmarks/selmer_quartics.py --verify > selmer-quartics-results.jsonl
```

Each benchmark record contains the arithmetic label and delta, all quartic
coefficients, the conic matrix and the exact Jacobian/target changes. These
permit the geometric construction to be replayed without Sage. The record is
not an independent class-group certificate. Selmer arithmetic is timed once per
curve; quartic construction (including internal checks) and a subsequent evidence
audit have separate timings. There is no fixed performance threshold or live
Magma comparison.

Number-field completeness remains a certified Sage dependency. The conic solver's
output is checked as a rational point, while the ensuing parametrization, inverse
maps and marking are checked independently in exact arithmetic. This does not
constitute formal verification or an external mathematical review. Parameters
can yield large coefficients; general minimization/reduction is not implemented
here. Failure of a backend, memory or runtime resources raises or leaves a parent
computation partial; it never justifies dropping a class.

The next step is the labelled Cassels–Tate pairing and its radical. Global fibre
equivalence, 4-Selmer group operations and the curve-level API remain later steps.

## Primary references and conventions

1. [Sage quadratic-equation interface](https://doc.sagemath.org/html/en/reference/quadratic_forms/sage/quadratic_forms/qfsolve.html): `qfsolve` and the underlying PARI/Denis Simon algorithms. Only point solving is delegated here; parametrizations and their inverses are constructed and verified directly.
2. [PARI quadratic-form functions](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Linear_algebra__Vectors__matrices_and_sets.html#qfsolve): exact rational quadratic-form solving and its return conventions.
3. [Cubic 2-Selmer arithmetic and its references](two-selmer.md): norm-kernel labels, auxiliary curves and local evidence.
4. [Marked-covering interfaces](marked-coverings.md) and the [core descent sources](../README.md#references), including Fisher–Sills for the canonical quartic map.
