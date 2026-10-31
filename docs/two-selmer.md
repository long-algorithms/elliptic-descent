# Labelled 2-Selmer arithmetic — v0.5 development step 2

The package remains **0.4.0**. `two_selmer(E)` computes the full 2-Selmer group
of a rational elliptic curve in explicit cubic-algebra coordinates, together
with Kummer maps and rational 2-/4-torsion data. It does **not** compute the
4-Selmer group, its fibres, a Cassels--Tate pairing, or a Mordell--Weil basis.
[Binary-quartic representatives with their marked maps](selmer-quartics.md)
are now available as the subsequent development step.

This is a direct reference implementation, not an adapter which attaches
arbitrary bit labels to PARI's `ell2cover` output. It uses Sage for certified
number-field arithmetic and exact local-test routines shared with the existing
package. Exhaustive norm-kernel enumeration is intentional; it is not claimed
to outperform local-image linear algebra or established descent packages.

## Exact arithmetic model

For the retained working Weierstrass model put

\[
 v=y+(a_1x+a_3)/2,\qquad
 h(x)=x^3+(b_2/4)x^2+(b_4/2)x+b_6/4,
 \qquad L=\mathbb Q[\theta]=\mathbb Q[t]/(h).
\]

The cubic is monic and squarefree; its factors can have types `3`, `2+1`, or
`1+1+1`. All components remain labelled. The new cubic CRT arithmetic checks
factor reconstruction, idempotents and basis round trips; component resultants
and multiplication-matrix determinants independently check norms.

The standard split Galois-module sequence

\[
 0\longrightarrow E[2]\longrightarrow \operatorname{Res}_{L/\mathbb Q}\mu_2
 \xrightarrow{N}\mu_2\longrightarrow0
\]

splits by diagonal signs because the degree is odd. Therefore

\[
 H^1(\mathbb Q,E[2])\simeq
 \ker\left(L^*/L^{*2}\xrightarrow{N}\mathbb Q^*/\mathbb Q^{*2}\right).
\]

There is **no rational-scalar quotient**. Indeed, a rational scalar c has norm
c cubed, so a scalar in the norm kernel is already a rational square. Reusing
the quartic relative-descent scalar quotient here would lose actual classes.
The description and Kummer convention are explained by Stoll [1, pp. 9--10].

## Finite support and arithmetic completeness

The support S contains 2, all primes in the numerator/denominator of disc(h),
and all coefficient denominators. Outside S, h is integral with unit
discriminant and L is unramified. This is deliberately conservative, not a
claim to be the smallest possible support.

The factorwise backend computes every `K_i(S_i,2)`, using every prime ideal
above S and including the class-group 2-torsion contribution, not just S-units.
Rational factors use `QQ.selmer_space`. The product coordinates concatenate
these independent bases. The norm map is a matrix of rational square classes,
and its full kernel is enumerated. Basis independence/exhaustiveness depends
on the certified Sage computation; coordinate and norm checks alone would not
establish it. `proof=True` is passed explicitly, independently of global flags.

Why is the support exhaustive? Outside S a local Kummer class has even
valuations in every factor. If x is nonintegral, monicity gives
`v(h(x))=3*v(x)`, so `v(x)` is even; all factors `x-theta` have that valuation.
If x is integral, at most one factor can vanish modulo the unramified prime;
it has residue degree one, and its valuation is even since h(x) is a square.
Rational 2-torsion values follow by the derivative formula or continuity.
Thus every locally soluble cohomology class belongs to this finite norm kernel.

Every kernel class is still tested locally. A norm-compatible space is not
mistaken for the Selmer group, and no equality of its cardinality with an
external program is used as a completeness proof.

## Auxiliary local-test models

For an element delta of norm r squared, expand

\[
 \delta(u_0+u_1\theta+u_2\theta^2)^2=Q_0+Q_1\theta+Q_2\theta^2.
\]

The degree-four genus-one model in P^3 is

\[
 Q_2=0,\qquad Q_1+z^2=0.
\]

It maps to the working E by weighted coordinates of weights `(2,3,1)`:

\[
 (X:Y:Z)=
 \left(Q_0:\ rN(u)-(a_1Q_0z+a_3z^3)/2:\ z\right).
\]

The norm identity gives the exact Weierstrass equation. It is checked both by
an independently calculated norm form and by polynomial reduction modulo the
two quadrics. These are 2-coverings of E of map degree four, despite their
embedding as intersections of quadrics; embedding degree is not descent degree.

Over an algebraic closure, evaluate at the three distinct roots theta_i and
rescale by square roots of delta_i. The equations become

\[
 w_i^2=X-\theta_i z^2\quad (i=1,2,3).
\]

Eliminating X gives two diagonal quadrics whose pencil has four distinct roots.
Hence the intersection is smooth of genus one. The map has no basepoint:
`z=Q0=0` would force all three w_i and hence every coordinate to vanish.
A generic point has four lifts after the norm-sign condition; even sign changes
are its E[2] deck transformations. Choosing the square roots with product r
shows that its twisting class is exactly delta in the norm-kernel description.
For delta=1 the point `(1,0,0,0)` lies over the origin and the cover is trivial.
This supplies an explicit local-membership test for each cohomology class.

At the real place, the signs of delta at the ordered real roots must be those
of `x-theta`: a prefix of plus signs followed by minus signs. Unlike the even
quartic construction, a common sign change is **not** allowed. Complex factors
impose no obstruction. Exact disjoint root intervals and sign certificates give
a complete real decision. Necessity and sufficiency can be checked away from
the finite exceptional fibres because these models are smooth.

At finite places the existing complete smooth-intersection routines are reused,
including exact integral model changes and independently checked modular/Hensel
certificates. The determinant-pencil good-reduction bound covers all remaining
primes. Unbounded mode terminates by the same finite-branching/compactness and
smooth Hensel argument as the existing solver. Explicit search limits retain
unknown outcomes. The auxiliary models are not yet binary-quartic presentations
suitable for automatically feeding the relative 4-descent interface.

### A compact singular-reduction obstruction

The 571a1 fixture exposed a large residue tree at p=571. The optional
`radical_plane_mod_p2` certificate avoids enumerating it. At an odd prime,
suppose the reductions of both Hessians have a common two-dimensional radical
W. Record an integral basis matrix N, invertible modulo p, with W in its last
two columns. Verify that the two binary forms on the complementary coordinates
have no common projective zero over the algebraic closure of F_p. Hence every
primitive residue solution has its first two coordinates divisible by p.

The radical condition says that all Hessian entries involving W are divisible
by p. Write a possible solution in the new coordinates as `(p*a,b)`. Modulo
p squared its values therefore reduce to the forms restricted to W. Divide
those restricted forms by p; if this second pair of binary forms has no common
projective zero, then b is also divisible by p. This contradicts primitivity,
so the original system has no Q_p point.

The verifier checks the matrix identity, determinant modulo p, divisibility,
and two small polynomial gcds, including the point at infinity. It does not
trust the detector's nullspace computation or enumerate P^3(F_p). This is a
sufficient test only: failure falls back to the unchanged complete solver.
The certificate is not used at p=2 or to override a zero work budget or a
requested maximum depth below two. Small-prime exhaustive residue tests and
changed-coordinate/tampered-certificate tests check the implementation.
For the four rejected 571a1 norm classes, this certificate proves the required
local obstruction directly, rather than raising an arbitrary search limit.

## Labels, basis changes and partial results

There are two different coordinate systems. **Ambient labels** are bit vectors
in the certified L(S,2) basis. **Selmer coordinates** are bit vectors in the
chosen basis of its locally soluble subspace. They are not interchangeable.

`basis_labels` records the ambient vector of each Selmer basis element.
`ambient_label`, `coordinates_of_label`, `coordinates`, and `representative`
translate explicitly. `representative(c)` returns three rational coefficients
of an element of L, not a covering curve or an assertion of a rational point.
Addition is XOR in one specified basis. A basis change is validated with
`with_basis`; it cannot silently change the underlying Selmer subspace.

After complete enumeration, the accepted labels must contain exactly all
vectors of their span. The implementation checks this instead of taking a
span and silently admitting untested or rejected classes. The zero class must
be present. The ledger checks duplicate, out-of-kernel and missing labels.
A partial result exposes decisions and known soluble labels but refuses to
report a dimension or complete basis. Incomplete candidate enumeration and
unresolved local tests are independent conditions.

## Kummer and torsion bookkeeping

A point on the original user curve is first transported by the retained
origin-preserving model change. Its Kummer value is `x(P)-theta`. The origin
maps to 1. At a nonzero rational 2-torsion point with x-coordinate alpha,
replace the zero component by `h'(alpha)`; every other component stays
`alpha-theta`. The norm is then `h'(alpha)^2`, so the element is invertible
and belongs to the norm kernel. This is the usual extension of the Kummer map,
not a convention making all torsion classes zero.

`torsion_data()` records every rational 2-torsion point, its actual Selmer
coordinates, the image basis, the kernel, and **all rational halves**. To find
halves of nonzero T, solve the exact rational-root polynomial

\[
 h'(x)^2-4h(x)(2x+h_2+x(T))=0.
\]

For each rational root with h(x) nonzero, test both rational square roots and
check doubling on the original curve. This finite calculation exhausts the
possible rational halves; it is not a bounded point search. Halves of the
origin are exactly E(Q)[2]. Their union gives E(Q)[4]. The verifier checks

\[
 \ker(E(\mathbb Q)[2]\to\operatorname{Sel}^{(2)}(E))
 =E(\mathbb Q)[2]\cap2E(\mathbb Q)=2E(\mathbb Q)[4].
\]

Consequently rational 4-torsion is accounted for rather than assuming that the
map from E(Q)[2] is injective. The tests compare halves with Sage division points.
They also verify the Kummer homomorphism on rational point sums and exact
chord identities, including nonzero a1/a3 and nonminimal presentations.

## Original models and natural transport

The `RationalCurveContext` from step 1 retains the original and working models
and exact maps in both directions. `minimal=False` deliberately retains the
input presentation. Kummer inputs and torsion outputs use the **original** E.

`transport_to(other, c, change)` transports a Selmer class under a specified
isomorphism of original models, without assuming identical field generators
or Sage bases. For a working-model change with `x_old=u^2*x_new+r`, the algebra
map substitutes `theta_old=u^2*theta_new+r`; completed-square cubics satisfy
`h_old(u^2*t+r)=u^6*h_new(t)`. The extra u squared in Kummer values is itself a
square. The polynomial identity and endpoints are checked. Rational translations
are not silently treated as origin-preserving model changes.

## Usage

```python
from elliptic_descent import two_selmer

S = two_selmer((0, 0, 0, 4, 0), proof=True)  # y^2=x^3+4x
assert S.certified and S.verify_evidence()
print(S.dimension, S.basis_labels)
assert S.coordinates(S.representative(0)) == 0
assert S.kummer((0, 0)) == 0       # twice the rational point (2,4)
assert S.kummer((2, 4)) != 0
T = S.torsion_data()
assert T.verify(S)
assert len(T.four_torsion) == 4
```

`proof=False` propagates GRH assumptions for nonrational factors. A wholly split
cubic uses only rational arithmetic and need not acquire such an assumption.
A complete conditional calculation can still have verified exact local evidence;
that does not make its number-field arithmetic unconditional.

## Verification and scope

The Sage tests compare complete dimensions with independent eclib computations,
including zero Selmer dimension, no/cyclic/full rational 2-torsion, rational
4-torsion and the 571a1 nontrivial-Sha fixture. Thus the algorithm is not merely
returning the subgroup generated by found rational points. Tests also cover
basis changes, model transport, partial/conditional outcomes and invalid labels.
The benchmark records dimensions, support, candidate counts, basis representatives,
torsion image/kernel data and separate solver/evidence timings.

The independently checked norms, maps and local certificates do not provide a
second independent proof of class-group completeness. Certified Sage arithmetic
remains a trusted dependency. Manual dataclass construction is not a certificate,
and the implementation is not formally verified or externally mathematically
reviewed. eclib tests compare cardinalities, not a class-by-class correspondence.
No live Magma comparison or performance superiority is claimed.

The subsequent [step 3 interface](selmer-quartics.md) now converts arbitrary
labelled classes to marked binary quartics via the auxiliary conic, preserving
the arithmetic label and original target. This does not turn a class list into
a global 4-Selmer group. Pairings, global assembly and exact Mordell--Weil
rank/saturation remain later work.

## Primary references

1. Michael Stoll, [*Selmer groups and descent*](https://people.maths.bris.ac.uk/~matyd/Trieste2017/Stoll.pdf), 2017, Section 3, pp. 9--10: cubic norm kernel, finite support and the x-theta Kummer map.
2. Schaefer--Stoll, [*How to do a p-descent on an elliptic curve*](https://www.mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf), *Transactions AMS* 356 (2004), 1209--1231: reduction to certified number-field computations.
3. [Sage number-field Selmer spaces](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html), with proof flags and class-group contributions.
4. [Sage/eclib `selmer_rank`](https://doc.sagemath.org/html/en/reference/libs/sage/libs/eclib/interface.html#sage.libs.eclib.interface.mwrank_EllipticCurve.selmer_rank), including the 571a1 regression example. A Selmer dimension is not a Mordell--Weil rank.
