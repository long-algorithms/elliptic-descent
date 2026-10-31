# Global 4-Selmer fibres — v0.5 development step 5

`assemble_four_selmer_fibres(S, B)` constructs the global classes above the
image of Selmer-4 in a complete labelled Selmer-2 space. It retains explicit
marked 4-coverings, the projection, the kernel injection and exact global
covering-equivalence witnesses. Here B is the verified pairing on the **same
chosen basis** of S; it is computed when omitted.

The package remains **0.4.0**. This step assembles the set of classes and the
maps in the exact sequence, not addition on the whole 4-Selmer group. A key
`(image_coordinates, orbit_index)` is a presentation-dependent identifier,
not an additive coordinate or an invariant-factor generator. Group arithmetic
is the next step. Inputs remain rational elliptic curves with the explicit
markings furnished by the preceding steps.

## Why concatenating relative lists is wrong

Write S2 and S4 for the 2- and 4-Selmer groups and R for the radical of the
Cassels–Tate pairing on S2. The cohomological exact sequence is

    E(Q)[2] --delta--> S2 --i--> S4 --rho--> S2,
    im(rho) = R.

The first map is the ordinary Kummer map restricted to rational 2-torsion.
Its kernel is `2 E(Q)[4]`, not necessarily zero. Its image D is available as
actual S2 coordinate vectors from the torsion interface of step 2. Hence

    ker(rho) = S2 / D,
    |S4| = |R| * |S2| / |D|.

These identities are consistency checks, not an algorithm for assigning
unlabelled coverings to global classes. In particular, dividing each list by
`|E(Q)[2]|` would be wrong in the presence of rational 4-torsion.

Fix a marked 2-cover C -> E representing c in R. The everywhere locally soluble
relative covers D4 -> C, up to isomorphism over the identity of C, form a torsor
under S2. To see the local condition, trivialize one locally soluble relative
cover over each completion: a twist is locally soluble precisely when its
class is a local Kummer image. Thus, once nonempty, this relative Selmer set
has |S2| elements. Its **global** fibre over c is different: the intermediate
marking can change by an automorphism of C over E.

Those automorphisms are exactly translations by E[2](Q). Indeed geometrically
C -> E is multiplication by two. An automorphism over E can only be a deck
translation; a nontranslation origin-preserving automorphism cannot compose
with [2] to give [2]. Descent to Q requires the translating point to be rational.
Changing the intermediate diagram by translation T acts on the relative torsor
by the Kummer class delta(T). This is also the connecting-map interpretation
of the exact sequence above. Consequently the action kernel is `2 E(Q)[4]`,
every orbit has |D| elements, and the orbits are exactly the global fibre.

Conversely, an isomorphism of two global 4-coverings induces an isomorphism of
their canonical intermediate 2-coverings, since quotienting the geometric E[4]
deck group by its characteristic E[2] subgroup is intrinsic. For two diagrams
above our fixed C, that induced map is one of these rational E[2] translations.
Different S2 labels cannot give equivalent global coverings. This proves the
exhaustiveness and nonredundancy of the orbit construction, not merely its count.
See Fisher [1], equation (1) and Sections 4–5 for the exact sequence and the
fixed-target norm-pair equivalence underlying this construction.

## Explicit intermediate translations

Each labelled quartic retains its conic presentation of the cubic-algebra
auxiliary cover. For a nonzero rational 2-torsion point T, let e be the idempotent
of its rational root in the monic cubic algebra. The unit

    w = 2e - 1

has square one and norm one. Its signs are + at the root corresponding to T
and - at the other two roots. This is the even-sign realization of that actual
E[2] point, not an arbitrary numbering of three involutions. Multiplication by w
on the three u coordinates, with the fourth coordinate z fixed, preserves the
auxiliary quadrics and their complete map to E.

If p(X,Z)=P*(X²,XZ,Z²) is the retained conic parametrization and W is the
multiplication matrix of w, then `P^-1 W P` preserves the Veronese conic.
Images of the parameters `(1:0)`, `(0:1)` and `(1:1)` recover a rational matrix M.
All coefficients of its symmetric-square action are then checked; three sample
points alone are not the verification. The resulting weighted map is

    g_T(X:Y:Z) = (aX+bZ : kY : cX+dZ),
    M = [[a,b],[c,d]],
    p(M(X,Z)) = k W p(X,Z),     F(M(X,Z)) = k² F(X,Z).

The second identity verifies the quartic equation; the first ties the map to
the specified deck translation of the labelled auxiliary curve. It therefore
commutes with the map to the **original** E, not just its Jacobian up to sign.
We retain both M and k, since an unrecorded sign on Y would change the action.
The zero torsion point gives the identity transformation.

## Action on relative norm pairs

Put A=Q[theta]/f for the intermediate quartic. Its branch points are permuted
by g_T, and none is moved to infinity in a genuine degree-four presentation.
Define the invertible algebra map and unit

    sigma(theta) = (d*theta-b)/(a-c*theta),
    t = a-c*theta,            Norm(t) = k².

For a norm pair `(xi,r)` the transported pair is

    xi_raw = t * sigma(xi),    r_raw = k*r.

The coordinate change is `V=sigma(U)`. If q is the vector of four coefficients
of `xi U²`, multiplication by t after sigma is represented by a rational matrix L.
Its first two columns are

    (a,-c,0,0)^T,  (-b,d,0,0)^T.

The lower-right 2-by-2 block is invertible. Thus the defining quadric ideals
are carried into each other and the induced map on C is exactly g_T. The
norm coordinate transforms by k, including its sign. The verifier checks the
algebra root identity, norm, all four quadratic identities and the Y identity.
These stagewise identities avoid expanding the full degree-sixteen composite.

## Matching oriented classes and returning witnesses

The transported pair is matched to the finite list of accepted relative models.
An explicit support includes the full quartic support, every coordinate
denominator and every element norm. Outside this support the order is integral
and etale, and each element is integral with unit norm, hence a unit in every
factor. Any rational scalar identification therefore has even valuations there.
The scalar image is generated by -1 and the supported rational primes.

Factorwise square-class coordinates solve the diagonal-scalar equation. For a
matched target we find exact c and v satisfying

    xi_target = c*v²*xi_raw,
    r_target = c²*Norm(v)*r_raw.

Square roots never extend the number fields. They are checked in rational
power-basis arithmetic. A norm-sign reversal is allowed only when the existing
orientation certificate supplies a reversing unit. Thus opposite orientations
of types 4 and 2+2 are not accidentally merged. In the fully split case matching
uses elementary rational sign/valuation coordinates rather than class groups.

The existing scalar-square verifier checks all quadratic and weighted-map
identities. Combining its multiplication matrix with sigma inverse gives an
explicit rational 4-by-4 matrix from the selected target model to the original
source model. This is the `FibreEquivalence` returned by `fibre.equivalence(i,j)`.
It is an isomorphism over E whose intermediate map is the recorded T translation.
A false return between two computed diagrams means they are in different orbits;
this is not an arbitrary-external-cover equivalence API.

Every torsion point's entire permutation is retained. Checks require bijections,
the E[2] group law, the correct stabilizers and actual orbit partitions. No orbit
is inferred from the expected size, and a missing or ambiguous match raises an
error rather than deleting a class.

## The actual kernel injection, not an arbitrary labelling

Over the zero intermediate class C, choose a rational branch root alpha. Put
`A0=f'(alpha)`, `s=(f''(alpha)/2-f2)/3`. The earlier rational-root identification
sends C to its canonical Jacobian J with

    x_J = s + A0/(x_C-alpha),       pi_C = [2] composed with this identification.

For the retained change J -> working E, write `j=u²*theta+r` in the working
cubic algebra. Then the other three branch points have expression

    beta = alpha + A0/(j-s).

For a relative norm pair xi, its ordinary 2-covering class after this
identification is represented by

    delta_xi = (s-j) * xi(beta)/xi(alpha).

All denominators are units, and `f(beta)=0` is checked exactly. On the relative
curve the defining norm-pair equation gives the function-field identity

    x_work-theta = delta_xi * (U(beta)/(u*U(alpha)))².

This identifies the actual cohomology class of that diagram. Each computed
class is converted to S2 coordinates using the certified arithmetic context.
The labels must give a bijection with **all** S2 vectors. For each rational T,
the computed permutation must send a label eta to `eta + delta(T)`.

`kernel_map.apply(c)` returns the global class corresponding to i(c). Its
kernel basis is the actual torsion Kummer image. `quotient_coordinates` and
`from_quotient` give explicit coordinates on S2/D, and `zero_class` is i(0).
Pairwise verification checks that two zero-fibre diagrams are in the same orbit
if and only if their labels differ by D. This verifies the exact-sequence map,
not just the equation for its order. Changing the selected rational branch
point shifts labels by a torsion Kummer class and leaves the quotient unchanged.

## Projection, exclusions and partial computations

A global class records its S2 image, orbit index, all relative indices in the
orbit, and a chosen marked 4-cover of the user's original E. `projection(key)`
returns its labelled S2 image. For a class outside the pairing radical, an
explicit vector with pairing value one certifies exclusion.

Default assembly descends only radical classes. `prune=False` invokes relative
descent on every S2 class and checks the resulting positive/empty answer against
the pairing. The relative solver itself never uses the pairing to discard a
candidate. This is a reference cross-check, not a separate independent entire
4-descent implementation.

`max_fibres`, `max_candidates`, `max_depth` and `max_nodes` are explicit limits.
An incomplete relative list is not quotiented as though it were exhaustive.
Such a result has `complete=False`; `cardinality` refuses to give a completed
answer. `predicted_cardinality` is separately labelled as the exact-sequence
prediction from the complete parent data. No limit turns unknown into empty.
An incomplete pairing has no radical and is rejected by this assembly interface.

Arithmetic assumptions are retained from the parent Selmer space, pairing,
relative descents and any matching context. A complete conditional result can
have exactly checked geometric evidence but is not certified. The package
continues to trust certified Sage number-field arithmetic for exhaustiveness;
recorded matrices and norm identities do not constitute a second class-group
proof. Manual result construction is not a certificate.

## Usage and reproducibility

```python
from elliptic_descent import two_selmer, cassels_tate_pairing, assemble_four_selmer_fibres

S = two_selmer((0, 0, 0, 4, 0), proof=True)
B = cassels_tate_pairing(S)
F = assemble_four_selmer_fibres(S, B)
assert F.certified and F.verify_evidence()
assert F.cardinality == 4
assert F.projection(F.zero_class.key) == 0
for C in F.known_classes:
    assert C.covering.target == S.curves.original
    print(C.key, C.relative_indices)
```

This example has a nonzero rational 2-torsion point with a rational half. Its
Kummer image is zero, so its action fixes every relative class; dividing by two
would be incorrect. Tests contrast it with cyclic 2-torsion without rational
4-torsion and full rational 2-torsion with both possible stabilizer behaviours.
They also cover actual coordinate witnesses, original models, basis changes,
corrupted actions/labels/orbits, and partial and conditional calculations.

The benchmark exports every intermediate conic, relative norm pair, torsion
permutation and explicit equivalence matrix. For zero-fibre labels it additionally
exports cubic square-root witnesses comparing the derived class with its stored
S2 representative. `replay_fibres.py` rechecks this geometry and kernel accounting
without Sage or class discovery. It does not redo the parent Selmer computation
or the complete local enumeration. Existing arithmetic, pairing, relative-descent,
packaging and dependency-version regressions remain enabled.

```sh
sage -python benchmarks/four_selmer_fibres.py --verify > four-selmer-fibres.jsonl
python benchmarks/replay_fibres.py four-selmer-fibres.jsonl
```

General minimization/reduction and faster group-law-based twisting are not added
here; repeatedly descending all image classes can be expensive. No fixed timing
or performance superiority is claimed. Group addition, invariant-factor generators
and the fully integrated curve-level interface remain later steps. This is not
formal verification or an external mathematical review.

## References

1. T. A. Fisher, [*Some improvements to 4-descent on an elliptic curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf), ANTS VIII (2008): equation (1), Sections 4–5, and the later group-law constructions. The explicit action and zero-fibre formulas above are derived in the package's retained coordinate conventions and checked directly.
2. [Labelled 2-Selmer arithmetic](two-selmer.md): the Kummer map, rational halves and the finite-support norm-kernel computation.
3. [Marked quartic representatives](selmer-quartics.md): conic presentations, exact inverse maps and original-target markings.
4. [Relative descent](relative-descent.md) and its factorization-specific extensions: fixed-target norm-pair equivalence and completeness.
5. [Cassels–Tate pairing](cassels-tate.md): the labelled radical and its checked local evidence.
