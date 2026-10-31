# Relative 4-descent: baseline and factorization-specific extensions

This document specifies the irreducible baseline of `relative_four_descent`.
The API also supports two irreducible quadratic factors and all rational-root
patterns. [Two-quadratic descent](two-quadratic-descent.md) and
[rational-root descent](rational-root-descent.md) supply the product-algebra
extensions and their different orientation rules. The irreducible argument
below is retained as the original case, not applied unchanged to linear factors.
It does not claim a full curve-level 4-Selmer-group computation, a rank
algorithm, minimization/reduction, or a practical algorithm for all input sizes.

## Input, output, and the equivalence being computed

Let

\[
 C: Y^2=F(X,Z),\qquad F(X,Z)=\sum_{i=0}^4 f_iX^iZ^{4-i},
\]

where `f=(f_0,...,f_4)` consists of exact rationals, `f_4 != 0`, and
`f(t)=F(t,1)` is squarefree and irreducible over Q. The input is this **fixed
quartic**, with the canonical degree-four map to its Jacobian described below.
It can be passed as coefficients or as a `QuarticTwoCover` object. A different,
arbitrarily supplied map to a different elliptic-curve model is not silently
substituted for this marking.

The output consists of all everywhere locally soluble relative 2-coverings
`D -> C`, up to isomorphisms **over the identity of C**, with explicit equations
and maps. Composing with `C -> E` makes these 4-coverings of E. Thus the degrees
of the maps are four and sixteen, respectively, not two and four.

An entry in `result.coverings` includes its orientation. The separate
`result.unoriented_coverings` quotient identifies a map with its composition
with `(X:Y:Z) -> (X:-Y:Z)`. It does not quotient by other target automorphisms,
rational torsion translations, or arbitrary isomorphisms of abstract curves.
The source's Jacobian alone would not determine the covering problem.

The source need not be locally soluble for the routine to execute. In that
case it still returns the set specified above, possibly empty. An empty result
is not by itself a statement about a nonzero element of Sha: that interpretation
also requires independent information that C is everywhere locally soluble.

## The norm-pair parametrization

Put `K=Q[theta]=Q[t]/(f)` and `a=f_4`. For

\[
 \xi\in K^*,\qquad r\in\mathbb Q^*,\qquad N_{K/\mathbb Q}(\xi)=ar^2,
\]

expand

\[
 \xi U^2=q_0+q_1\theta+q_2\theta^2+q_3\theta^3,
 \quad U=u_0+u_1\theta+u_2\theta^2+u_3\theta^3.
\]

Set `D_xi={q_2=q_3=0}` in projective three-space. Its map to the fixed C is

\[
 (X:Y:Z)=(-q_0:arN(U):q_1).
\]

The polynomial norm identity proves that the image lies on C. Over an
algebraic closure, evaluation at the four distinct roots diagonalizes the
algebra. The resulting two diagonal quadrics have four distinct pencil roots.
Consequently their intersection is smooth of genus one. Its common map
basepoints would satisfy all four `q_i=0`, hence `xi*U^2=0`; invertibility of xi
and the reduced algebra force `U=0`, which is not a projective point.

For a generic point of C, four choices of square roots modulo simultaneous
sign satisfy the prescribed product sign. This proves that the map has degree
four. It is a relative 2-cover: over the algebraic closure its deck group is the
even sign changes modulo simultaneous sign, identified with E[2].

The standard descent construction is exhaustive for the everywhere locally
soluble relative covers in question. Here is the obstruction issue that must
not be omitted. The pullbacks to D of the four ramification points of `C -> P^1`
are linearly equivalent divisors of degree four. They determine a rational
geometric divisor class and a possible Brauer obstruction to its realization.
A point on D over each completion kills that obstruction locally. The
local-to-global injectivity for `Br(Q)` kills it globally, giving the required
projective-three-space model. Factoring the conjugate rank-one quadratic forms
recovers xi and the product sign r. This is the construction in Fisher,
Section 3, with equivalence and its cohomological interpretation in Sections
4--5; see also Womack, Chapter 2.

### Equivalence and the two orientations

For a fixed target C the norm pairs have the equivalence

\[
 (\xi,r)\sim(c v^2\xi,\ c^2N(v)r),\qquad c\in\mathbb Q^*,\ v\in K^*.
\]

Multiplication by v is a rational invertible linear change of the four
coordinates. With `U_old=v*U_new`, the quadrics scale by c and the map coordinates
scale by `(c,c^2,c)`, exactly the weighted-projective equivalence on C.
Conversely, an isomorphism over C preserves its pulled-back hyperplane class
and the four labelled rank-one forms. It therefore gives this scalar-square
relation. Fisher's Lemma 5.1 and Theorem 5.2 explain why allowing determinant
sign minus one changes the target map and loses an orientation.

After forgetting r, equivalence is precisely equality in
`K*/(Q* K*2)`. For an **irreducible** quartic the two signs of r give distinct
fixed-target coverings. Indeed a scalar-square stabilizer has `v^2 in Q*`.
Then v is rational or lies in a quadratic subfield of the quartic field, so
`N_K/Q(v)` is positive: it is respectively a fourth power or the square of a
quadratic norm. It cannot reverse r. Both signs are therefore returned, with
no duplicate oriented classes. This argument is one reason not to extend the
same implementation silently to reducible etale algebras.

## Finite support and an exhaustive quotient

The default rational support S contains 2, the primes in the numerator and
denominator of a and `disc(f)`, and all coefficient denominators. Outside S the
power-basis order is etale and a is a unit.

The support argument in [arithmetic.md](arithmetic.md) proves that every
locally soluble construction class has a representative in `K(S,2)`, after
multiplication by a rational scalar. Briefly, at a good prime a generic local
point gives `xi*U^2=c(X-theta*Z)`. Normalize X,Z to be primitive. In the
unramified algebra, the valuations of `X-theta*Z` are even in every factor:
either it is a unit, or only a single degree-one factor is divisible and its
valuation is even by `F(X,Z)=Y^2`. Thus the valuations of xi at primes over p
all have the same parity, removable by multiplication by a rational p.
Only finitely many primes need this adjustment.

Sage computes the **entire** `K(S,2)`, including the class-group contribution,
not merely S-units. Exact F2 linear algebra then imposes

\[
 [N(\xi)]=[a]\quad\text{in }\mathbb Q^*/\mathbb Q^{*2}.
\]

This is an affine fibre, which can be empty. The image of rational scalars in
this fibre is generated by `-1` and the primes in S. To justify that no other
scalar identification is missing, suppose that a rational scalar has image
in `K(S,2)`. At every prime outside S the ramification index is one, so its
rational valuation must be even. It therefore lies in `Q(S,2)` modulo rational
squares. Row reduction computes the image, including dependencies caused by
quadratic subfields, and gives a canonical transversal. All its elements are
processed unless the user explicitly requests a candidate limit.

The norm equation by itself is not a Selmer condition. **Every resulting model
is tested locally.** Combining the norm-pair parametrization, the support
argument, the full rational-scalar quotient, both orientations, and the local
decisions proves exhaustiveness and nonredundancy. There is no additional
unverified “fake Selmer equals Selmer” inference.

## Complete local decisions

### Real place

At the ordered real roots `alpha_i` of f, the equation `xi*U^2=c(X-theta*Z)`
requires the signs of `xi(alpha_i)` to be a threshold pattern, up to a common
sign. Conversely, such a pattern permits a real x in the corresponding interval
and square roots in every real factor; complex factors impose no obstruction.
This constructs a real point. Necessity can be tested away from the finitely
many exceptional fibres because a smooth real curve with a point has an open
real analytic arc.

Root isolation and interval evaluation are exact. Since f is squarefree and xi
is a unit, refinement eventually separates every sign. `RealSquareClassResult`
records isolating rational intervals; its verifier checks their root counts,
ordering, and signs rather than trusting a floating-point root approximation.

### Finite places

Before enumeration, exact integral model changes saturate the equation rows
and remove common linear factors from reductions. A square factor in a single
row also permits a coordinate scaling. Every strict step decreases the
nonzero determinant-pencil discriminant valuation; at most four neutral steps
are allowed between strict decreases. This finite preprocessing is not a full
minimization algorithm and never limits the ensuing decision search. Its
rational coordinate and equation matrices are retained and checked exactly.
They are invertible over Q_p, so no local branch is discarded.

The p-adic routine then enumerates all normalized primitive residue solutions in
projective charts. An exact finite-field eliminant solves two quadratics,
including every degenerate fibre; every reconstructed root is checked against
both original equations. Higher digits satisfy a 2-by-3 linear system over
F_p and are all enumerated. These are accelerations of the reference residue
enumeration, not heuristic filters.

A strong Hensel witness certifies solubility when the equation valuations exceed
twice the valuation of a nonzero Jacobian minor. An empty residue level proves
insolubility. There is no default depth or work cutoff. Termination on the
smooth supported domain follows from compactness: an infinite finitely
branching tree of surviving residues has a p-adic limit point. At such a point
smoothness gives a nonzero minor, and sufficiently accurate approximations meet
the Hensel condition. If there is no p-adic point, some finite level is empty.
A requested resource limit returns `unknown`, not an obstruction.

The required finite places are bounded by the regular determinant-pencil
argument in [local.md](local.md): outside the primes dividing twice its leading
coefficient and discriminant the model has smooth genus-one reduction. The
Hasse bound gives a residue point and Hensel lifts it. This bound may be much
larger than a minimal bad-prime set, but it is mathematically sufficient.

## Verification and trust boundary

`check_construction()` verifies the defining coefficients, norm witness, and
norm form. The last is computed independently by a resultant rather than the
multiplication-matrix determinant used by the constructor. The construction
theorem then supplies smoothness, map degree and basepoint freeness. The slower
generic Groebner checks remain available and are cross-tested on small models.

Positive finite-place certificates and real signs have inexpensive exact
verifiers. A negative modular certificate is independently checked by direct
enumeration modulo its recorded prime power when affordable. Otherwise a
separate verifier enumerates all initial F_p residues without the solver's
eliminants, and checks higher levels using rank-0/rank-1 compatibility instead
of its RREF implementation. Exact change-of-model matrices transport these
certificates back to the original input. The separate verification budget can
be insufficient even though the solver's negative decision is correct.
`verify_evidence()` is an audit of the recorded geometric/local evidence, not a
second independent class-group computation. It returns false for unresolved
local decisions or an insufficient certificate-verification budget.

`proof=True` explicitly requests Sage's unconditional arithmetic. That backend
is part of the trusted computation. `proof=False` records a GRH assumption and
cannot produce `certified=True`. Large integer primes use Sage's
`is_prime(proof=True)`; no probable-prime result is labelled certified.

The result fields distinguish enumeration, decisions, and proof:

- `enumeration_complete`: every quotient class was processed;
- `complete`: enumeration is complete and all local decisions are known;
- `certified`: complete, unconditional arithmetic, and no unproved assumptions.

Manually constructing or modifying a result dataclass does not constitute a
proof. A result with no returned coverings is an empty descent only when the
completeness and certification requirements are satisfied.

## References

- T. Fisher, [*Some improvements to 4-descent on an elliptic curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf) (2008), Sections 3--5, particularly Lemma 5.1 and Theorem 5.2.
- T. Womack, [*Explicit Descent on Elliptic Curves*](https://johncremona.github.io/theses/womack.pdf) (2003), Chapter 2.
- T. Fisher and G. Sills, [*Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944) (2011), Lemma 3.1 for the quartic-to-Jacobian map.
- Sage, [number-field Selmer spaces](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html) and [certified integer primality](https://doc.sagemath.org/html/en/reference/rings_standard/sage/rings/integer.html).
- Magma, [higher descent handbook examples](https://magma.maths.usyd.edu.au/magma/handbook/text/1570), used as published regression fixtures, not as an execution backend.
