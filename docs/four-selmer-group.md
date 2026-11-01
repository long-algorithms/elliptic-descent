# Exact 4-Selmer group arithmetic — v0.5 development step 6

This proposed step adds operations on the **actual global covering classes**
constructed by [step 5](four-selmer-fibres.md). The package version remains
**0.4.0**. Its general algorithm is a deliberately expensive reference route:
trivialize each marked covering over an exact normal number field, compute its
Galois cocycle in E[4], and add cocycles modulo actual coboundaries. It is not
Fisher's efficient obstruction-algebra implementation, and no performance
advantage is claimed. The accompanying Sage discovery/integration tests and all
eight group benchmarks are acceptance gates, not results inferred from local
Python replay.

## Why abstract invariant factors do not assign an operation to our models

Step 5 supplies the set S4, a zero class, and the exact maps

    S2 --i--> S4 --rho--> R,
    ker(i) = delta(E(Q)[2]),        ker(rho) = im(i).

Multiplication by 2 inside S4 is **i composed with rho**; it is not rho itself.
These data determine the abstract invariant factors, but in general do not
identify which representative in a third fibre is the sum of two specified
coverings. Assigning coordinates arbitrarily within fibres would pass size and
order checks while mislabelling the actual sum.

The reference backend therefore computes a cohomology class from each model's
retained equations and signed covering map. A missing sum or two cohomologous
supposedly distinct global classes is an error, not permission to renumber them.

## Exact trivialization of a marked 4-covering

Let D -> C -> E be a recorded diagram. Write `C: Y^2=F(X,Z)`, with affine
quartic f, relative norm pair `(xi,r)`, `N(xi)=f4*r^2`, and the existing map

    (X:Y:Z) = (-q0 : f4*r*N(U) : q1),
    xi*U^2 = q0+q1*theta+q2*theta^2+q3*theta^3.

The intermediate marking is a specified origin-preserving Weierstrass change
from the canonical Jacobian J of C to the user's original E. The assembly
currently supplies zero target translation; this backend checks that condition
instead of silently dropping a translation from a different external input.
Let its parameters be `(u,r0,s0,t0)`, in the convention
`x_J=u^2*x_E+r0`.

Complete the square on E: `v=y+(a1*x+a3)/2`, so
`v^2=h(x)=product(x-e_i)`. Choose a root alpha of f over the algebraic closure,
and put

    A = f'(alpha),       b = (f''(alpha)/2-f2)/3,
    j_i = u^2*e_i+r0,    beta_i = alpha+A/(j_i-b),
    d_i = (b-j_i)*xi(beta_i)/xi(alpha).

The beta_i are the other three branch roots. Choose square roots s_i of d_i,
with the **signed** product

    s_1*s_2*s_3 = -f4*A*r / xi(alpha)^2.

The implementation fixes the first two roots and obtains the third from this
product, checking all three squares. On D set

    w_i = s_i*U(beta_i)/(u*U(alpha)).

The norm-pair equations give

    w_i^2 = x_base-e_i,
    w_1*w_2*w_3 = v_base,

where `(x_base,v_base)` is the branch-root identification C -> E applied to
D -> C. This identification has `C -> E = [2]` after the identification,
including the sign of its y-coordinate. The signed product above is essential
for that equality; replacing it with an arbitrary norm square root can invert
a class of order four.

Write `P=w1*w2+w1*w3+w2*w3`. The classical half-point formula is

    x_psi = x_base+P,
    v_psi = -v_base+(w1+w2+w3)*P,
    y_psi = v_psi-(a1*x_psi+a3)/2.

It gives `[2]psi=(x_base,y_base)`, so the original diagram is `[4]psi`.
The formula is the usual division-by-two formula with the product-positive
square-root convention; see Bekker–Zarhin [2], equations (3)–(6). Its curve
and duplication identities are also checked symbolically in the tests.
The rational map psi extends to a morphism between smooth proper genus-one
curves. Since `[2]psi` has degree four, psi has degree one and is an isomorphism.

### The origin and exceptional chart are explicit

The point Q_D mapping to O is constructed by exact Lagrange interpolation:

    U(alpha)=0,         U(beta_i)=1/s_i.

At all four branch roots one checks

    xi(t)*U(t)^2 = xi(alpha)*(alpha-t)/A.

Thus Q_D lies on the two defining quadrics. No algebraic point search or
numerical limit is involved.

For evaluation when `U(alpha)=0`, the three nonzero numbers `s_i*U(beta_i)`
agree up to signs. If all agree, psi is O. Otherwise the unique opposite sign
identifies the actual point `(e_i, -(a1*e_i+a3)/2)` of E[2]. This covers the
four exceptional points without division by zero. Away from this chart the
half-point formula is evaluated and its duplication identity is checked.

## Galois cocycles and the genuine sum of two coverings

Take a finite normal field L containing all trivialization data and E[4].
For every automorphism sigma of L/Q, compute

    z_D(sigma) = psi_D(sigma(Q_D)) in E[4].

Because `[4]psi_D` is defined over Q, the conjugate trivializations differ by
translations by E[4]. With the origin convention above, z_D satisfies

    z_D(sigma*tau) = z_D(sigma) + sigma(z_D(tau)).

The code checks this law for every sigma and every recorded group generator;
the generators are checked to reach all automorphisms. Induction on words
then gives the full cocycle law. This constructs the actual class of the
marked 4-covering, not a cocycle chosen to fit a predicted group.

E[4] is represented using two actual algebraic points A4,B4. Their sixteen
linear combinations must be distinct and killed by four. Galois images are
computed as point images, giving the actual module action. Coordinates are
`a+4*b`, with a,b modulo four; there is no carry between the two coordinates.

Two cocycles represent the same class precisely when their difference is
`sigma(T)-T` for one of those sixteen points. All sixteen coboundaries are
considered. Inflation to absolute Galois cohomology is injective here, so
passing to this finite normal field loses no class identifications. For every
pair of retained coverings, add the cocycles and locate the unique existing
class modulo coboundaries. If it is missing or duplicated, verification fails.
A change of trivializing origin changes a cocycle by a coboundary, not the
underlying class or its addition law.

## Discovery and independent exact replay

Sage's QQbar finds **exact** roots, not floating-point approximations. Its
`number_field_elements_from_algebraics` combines the data; a normal closure is
taken when needed. Every embedding conversion is checked. The discovery
backend exports only rational power-basis data:

- A monic irreducible field polynomial and degree-many distinct conjugates
  lying in that field; these prove normality in characteristic zero.
- Automorphism generator indices, with composition checked by exact substitution.
- Cubic roots, an actual E[4] basis, and each model's branch root and signed
  square roots.

The independent Python verifier checks irreducibility, every conjugate,
generator reachability, torsion equations/actions, all covering data, the
signed identities, origins and cocycles. It works without Sage or an algebraic
root-discovery routine. It binds each witness to the **same** f, xi, norm sign,
marking and class key from the parent fibre assembly. No claim is made that
this replay independently recomputes class groups or the entire Selmer set.

`max_field_degree` can reject an excessive field degree, raising an exception
without returning a completed group. It is not a runtime/memory bound on field
discovery. A common normal closure can be large; this is an exact reference
algorithm, not an efficient algorithm for every input. Discovery or resource
failure is never treated as a proof that a class is absent.

## Exact-map shortcuts that really determine the law

The default `backend='auto'` uses the existing exact-sequence maps only in
three situations where they uniquely determine operations on the named models:

1. R=0: the recorded kernel injection identifies every class with its actual
   coordinate in S2/delta(E(Q)[2]); addition is inherited from that quotient.
2. The kernel is zero: the recorded projection is an isomorphism onto R.
3. There are four classes: the recorded zero and **actual doubling map**
   distinguish the cyclic and Klein cases and uniquely fix their operations.

The last rule uses the values of `i(rho(x))` on the specific four classes,
not cardinality alone. It is intentionally **not** extended to arbitrary
nontrivial extensions of order eight or more. `backend='galois'` forces the
geometric reference route even when a shortcut is available; tests require
agreement class-by-class on a complete example.

## Operations, coordinates and exact-sequence compatibility

```python
from elliptic_descent import (
    two_selmer, cassels_tate_pairing, assemble_four_selmer_fibres,
    four_selmer_group,
)

S = two_selmer((0, 0, 0, 4, 0), proof=True)
F = assemble_four_selmer_fibres(S, cassels_tate_pairing(S))
G = four_selmer_group(F)  # Equivalent: F.group()
assert G.certified and G.verify_evidence()
assert G.invariant_factors == (4,)
g = G.generators[0]
assert g.order() == 4 and 4*g == G.zero
assert G.from_coordinates(g.coordinates) == g
assert 2*g == G.include_two_selmer(G.projection(g))
print(g.key, g.covering.target)
```

Addition, subtraction, negation and integer multiplication return elements of
the same explicitly constructed group. Each element retains its global key
and actual marked covering. Elements from different presentations cannot be
mixed silently. The API does not invent a natural identification by matching
invariant factors; rebase the parent arithmetic or use explicit covering data.

After the operation is known, actual order-four generators are chosen with
independent doubles; then order-two generators extend their subgroup. Every
invariant-factor coordinate combination is enumerated and checked to give a
bijection with the covering classes. The diagonal relations and generator
orders are checked. The group law is also checked for closure, identity,
inverses, commutativity, associativity and exponent dividing four.

For every class and every pair of classes, verify

    2*x = i(rho(x)),    rho(x+y)=rho(x)+rho(y),
    i(c+d)=i(c)+i(d).

These checks supplement geometric addition; they do not substitute for it.
An incomplete fibre assembly is rejected. Conditional parent assumptions are
retained, so an exact operation on a conditional class list is not certified
unconditionally. Manually constructing result dataclasses is not a proof.

## Tests and acceptance gates

The Python suite includes synthetic module/group tests explicitly labelled as
such, exact signed trivializations on real PR5 models, and a degree-16 normal
field certificate for **all four global classes** on `y^2=x^3-x`. The latter
computes cross-fibre sums geometrically and agrees with the separately forced
projection law. It replays the parent's 64 geometric equivalence witnesses;
it is not a new computation of that parent's Selmer arithmetic. A second real
record distinguishes cyclic order four using the actual kernel/projection.
An eight-dimensional field example checks zero-fibre cocycles against their
recorded torsion orbits. Altered roots, norm signs, torsion bases, conjugates,
keys, models, generators and addition tables must fail.

Genuine Sage tests separately exercise field discovery, every one of the eight
previous fibre fixtures, forced-versus-geometric comparison, the 2-by-4 case,
the nontrivial-radical case, original models, conditional arithmetic and partial
inputs. These tests are **not** replaced by the Python evidence fixtures.
The workflow keeps all previous regression and distribution gates and adds
separate eight-curve group benchmarks with exact replay:

```sh
sage -python -m pytest -q tests/test_sage_four_selmer_group.py
sage -python benchmarks/four_selmer_group.py --verify > groups.jsonl
python benchmarks/replay_four_group.py groups.jsonl
```

The benchmark exports complete parent geometry and exact group evidence.
`replay_four_group.py` verifies both, rather than accepting a table because its
size and invariant factors match. Construction and subsequent verification
times are separate. Shared CI timings are not performance guarantees.

This step is a proposal until the exact final revision passes those genuine
Sage and CI gates. Local geometric replay does not establish that every Sage
API call or every normal-closure calculation has succeeded. Neither the code
nor its tests constitute formal verification or external mathematical review.

## Exact generator-based certificate verification

The portable checker stores all field conjugates and the complete Galois
module, but need not repeat dense field arithmetic for every automorphism.
For a degree-n irreducible defining polynomial f, it checks f(sigma(a))=0
for each proposed generating automorphism. Such a substitution is an actual
field automorphism. It then checks exactly that these substitutions permute
the n distinct listed images and reach all of them from a. This proves that
every image is a root and that the list gives the entire Galois group. The
retained words recover both left and right multiplication in their correct
order; this is checked on a nonabelian dihedral example as well as the
covering fixtures.

Actual images of all sixteen E[4] points are computed for each generator.
Composition along the verified field words supplies their action for every
automorphism. For each covering, the checker evaluates the flex-origin
cocycle geometrically on the generators. The identity pi=[4]psi implies that
conjugate trivializations differ by an E[4] translation, so these values
uniquely determine the full cocycle. Every edge relation

    c(s*g)=c(s)+s(c(g))

is checked, including every edge to an already visited group element. No
value is chosen to force an expected cardinality, and inconsistency fails.
`cocycle_direct` retains the previous all-automorphism geometric evaluation;
regressions compare the complete cocycles and every group sum against it.

Substitution columns, evaluated immutable field data and flex origins are
cached as exact rational algebra. They do not cache a verification verdict.
The cache is tied to its actual field/torsion object, and replaced evidence
must be checked again. Corrupted nongenerator conjugates, inadequate generating
sets, altered signed models and alternate generator orders are tested.

### Avoiding large-field inversions without changing the maps

The same branch-root formulas can be evaluated without repeated division in
the full normal field. Invert f'(theta) and xi(theta) modulo the quartic f
over Q, then evaluate these degree-at-most-three Bezout polynomials at alpha.
Their products with f'(alpha) and xi(alpha) are checked to equal one. Writing
j_i=u^2*e_i+r and s as above, the identity

    product_i (j_i-s) = -lc(f)*f'(alpha)^2

then gives every beta_i using multiplication and rational division only. The
three cubic-algebra square classes remain exactly the original expressions.

For interpolation, put d_i=sqrt(delta_i)*f'(beta_i). Their product is

    product_i d_i = -disc(f)*r_relative / (lc(f)*xi(alpha)^2).

Its reciprocal therefore uses rational division only. Synthetic division of
f by t-beta_i and multiplication by the other two d_j give the original
Lagrange origin, with the same normalization, not merely an isomorphic point.
The reciprocal product identity and every branch evaluation are still checked.
For cocycle values, homogeneous half-point coordinates are compared directly
with the sixteen actual E[4] points, and doubling is checked against their
verified labels. The reference affine evaluator and generic interpolation are
retained in differential tests, including exceptional branch charts and signs.

Torsion-point additions in a large normal field use another exact reduction:
the first linear dependence among 1,v,...,v^n finds the minimal relation of
a denominator v. Its nonzero constant term yields v inverse as a polynomial
in v, and multiplication by v verifies that inverse. This works in the small
subfield generated by v, rather than requiring a dense extended GCD in the
whole degree-n field. The power-space search is bounded by n+1 for mathematical
reasons, not by a heuristic resource limit. Rational linear elimination and
inverse multiplication are exact, and test fixtures compare with ordinary
field inversion. No class or automorphism is omitted by these optimizations.

## References

1. T. A. Fisher, [*Some improvements to 4-descent on an elliptic curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf), ANTS VIII (2008), especially the exact sequence and group-law sections. The present reference route uses direct finite cocycles rather than claiming to implement Fisher's optimized obstruction-algebra algorithm.
2. B. M. Bekker and Y. G. Zarhin, [*The divisibility by 2 of rational points on elliptic curves*](https://arxiv.org/abs/1702.02255), Section 2, equations (3)–(6). The product sign is translated explicitly above.
3. [Sage exact algebraic fields](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/qqbar.html) and [number fields](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html): exact common-field conversion, normal closures and automorphisms.
4. [Global fibre construction](four-selmer-fibres.md), [labelled Selmer arithmetic](two-selmer.md), and [marked covering conventions](marked-coverings.md): the parent computations and trust boundary.
