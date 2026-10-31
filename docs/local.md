# Local solubility and certificates (v0.2)

The reference routines accept two nonzero homogeneous quadrics in four
variables with exact rational coefficients. They clear denominators, divide
each equation by its content, and fix its sign. This does not change its zero
locus over Q_p or R. All decisions use exact integer/rational arithmetic.

The routines return `soluble`, `insoluble`, or `unknown`. An unknown result has
no certificate. Singular input models may be passed to the finite-place
routine for testing, but the good-reduction bound requires a regular squarefree
pencil; unsupported models raise an exception rather than returning an empty
answer. None of these routines enumerates Selmer classes.

## Projective enumeration at finite primes

Every point of P^3(Q_p) has an integral representative whose first unit
coordinate is 1. The four disjoint charts enumerate all residues modulo p:
coordinates before the pivot are zero, and coordinates after it are arbitrary.
Each surviving residue is lifted by all p^3 choices of digits on the other
three coordinates. Earlier coordinates are allowed to become nonzero multiples
of p at higher precision. Omitting these lifts would incorrectly exclude
points.

The solver counts every tested residue, including rejected residues, against
`max_nodes`. It also respects `max_depth`. Either exhausted budget returns
`unknown`, retaining no negative conclusion. Large primes use lazy digit
iteration, so a small node budget does not first allocate p^3 objects.

### Positive certificates

An exact nonzero common rational zero is a global point certificate.
Otherwise let `x` be a tested integral point with pivot 1, and let `J` be a
2-by-2 Jacobian minor using two non-pivot coordinates. With the remaining
coordinate held fixed, the strong multivariate Hensel criterion is

\[
 \min_i v_p(Q_i(x)) > 2v_p(\det J(x)).
\]

A nonzero determinant satisfying this inequality certifies a Q_p point.
The inequality follows from the usual adjugate/Newton Hensel argument and
works also at p=2 and with singular reduction modulo p. A nonsingular residue
is the special case where the determinant is a unit.

The certificate records the point, pivot, minor, prime, and normalized-model
fingerprint. The verifier recomputes the values and Jacobian independently of
the residue search. No floating-point approximation enters this check.

### Negative certificates

If every primitive residue at some modulus p^k has been excluded, there is no
Q_p point. The negative certificate records that modulus and model identity.
Its independent verifier directly enumerates all normalized projective
residues modulo p^k rather than calling the solver's recursive lifting code.
There are

\[
 p^{3(k-1)}(p^3+p^2+p+1)
\]

such residues. The verifier has its own budget and returns `False` when it
cannot finish; this means the supplied evidence was not verified, not that
the model is soluble. Lack of a smooth point modulo p alone is never used as
an insolubility criterion.

## Real certificates

A supplied rational point, or one found in a small fixed search, proves real
solubility. A definite rational combination of the two quadratic forms proves
real insolubility: a common real zero would have to be zero. Definiteness is
checked exactly by Sylvester's criterion on leading principal minors.
The bounded pencil search is not a complete real-solubility algorithm. Failure
to find either type of certificate returns `unknown`.

For a complete real algorithm in a future release, see Fisher–Sills, Section
2.4. The current reference path is intentionally simpler and conservative.

## Why only finitely many primes need testing

Let A and B be integral Hessian matrices, so `2Q_1(u)=u^T A u` and similarly
for B over Z. No division by 2 is used in the finite-place polynomial tests.
Choose M=A+jB with nonzero determinant, trying j=0,...,4, and put

\[
 g(t)=\det(tM+B).
\]

For a regular squarefree pencil, at least one such j works. Include all primes
dividing `2*lc(g)*disc(g)`. Outside this set M remains invertible and the
self-adjoint operator `M^{-1}B` has distinct eigenvalues over the algebraic
closure. Simultaneous diagonalization then shows that the two quadrics form
a smooth complete intersection of genus one. It is geometrically connected.
The Hasse bound guarantees a point over F_p, and smooth Hensel lifting gives
a Q_p point. Thus the finite set is a valid, possibly redundant, local bound.

This argument requires the squarefree regular determinant pencil; the code
refuses to infer a bound if its discriminant vanishes. Primality/factorization
limits are described in [arithmetic.md](arithmetic.md).

`everywhere_local_solubility` combines the real test, these finite-place tests,
and the good-reduction bound. One insoluble place proves failure everywhere
local. One rational point proves solubility at every place. Otherwise all
required places must be certified soluble; any missing decision leaves the
answer unknown.

## Verification scope

The mathematical certificates prove local claims about a fixed model, not
that the input is the desired elliptic-curve covering. Geometry, maps, finite
candidate enumeration, and equivalence remain separate checks. Directly
constructing or modifying a `LocalResult` or `LocalScreen` is not a proof;
consumers of external records must verify the evidence against the model.

This is an independently auditable baseline, not an implementation of the
optimized Fisher–Sills transformations. The tests include mod-4 obstructions,
2-adic singular lifting, all projective charts, rational equation scalings,
forged certificates, budget exhaustion, and random comparisons with an
independent primitive-vector enumeration.

## References

- T. Fisher and G. Sills, [*Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944), 2011, Section 2.
- J. E. Cremona, T. Fisher and M. Stoll, [*Minimisation and reduction of 2-, 3- and 4-coverings*](https://arxiv.org/abs/0908.1741), 2010.
