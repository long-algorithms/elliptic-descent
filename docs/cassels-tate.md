# Labelled Cassels–Tate pairing — v0.5 development step 4

`cassels_tate_pairing(S)` computes the alternating pairing on a complete labelled
2-Selmer space over Q, in the **same basis as S**, together with local evidence
and its radical. The package version remains **0.4.0**. This step identifies the
subspace that lifts to the 4-Selmer group; it does not construct that group,
identify global fibre equivalences, or determine a Mordell–Weil basis.

## Formula, not inference from a rank

The implementation uses Fisher's Theorem 3.1 [1] on three binary quartics
representing `c`, `d` and `c+d`. These come from the labelled representative
interface of [step 3](selmer-quartics.md). The matrix is not guessed from a rank
bound, the radical, a list of soluble classes, or an external program's output.
PARI's independent implementation supplies a **comparison rank** in tests only.

All coefficient vectors are ascending. For

    g(X,Z)=a X^4+b X^3 Z+c X^2 Z^2+d X Z^3+e Z^4,

use the classical invariants `I=12ae-3bd+c^2` and
`J=72ace-27ad^2-27b^2e+9bcd-2c^3`, and the cubic algebra

    L_phi = Q[phi]/(phi^3-3I*phi+J).

Its cubic invariant is `z(g)=(4a*phi+3b^2-8ac)/3`. With Fisher's Hessian h,
set `G=(4*phi*g+h)/3` and

    H = (d^2 G/dX^2)/12 + (2/9)*(I-phi^2)*Z^2.

The exact identity `z(g)*G=H^2` is independently checked by cubic-algebra
multiplication. For three equal-invariant quartics, retain an actual unit m
satisfying `m^2=z(g1)*z(g2)*z(g3)`. If gamma is the coefficient of phi squared
in `(z(g2)*z(g3)/m)*H1`, the pairing is the product, over all places, of

    (g2(1,0), gamma(X_v,Z_v))_v,

where `g1(X_v,Z_v)` is a square and gamma is nonzero. The implementation stores
bits: **0 means +1, hence 0 in Q/Z; 1 means -1, hence 1/2 in Q/Z**.

The actual cubic invariants must be used, not representatives reduced modulo
squares. Such a reduction would change gamma. A product square root is found
componentwise in the parent's existing rational/number-field factors, without
extending those fields. Squaring in the rational product algebra verifies it.
All independent sign choices in a split product are legitimate by Fisher's
Remark 3.2(vi); tests exercise them, rather than assuming only an overall sign.

## Common invariants and arithmetic labels

A stored quartic has an explicit signed Weierstrass change to the working E.
If its scale parameter is u, substitute `(X,Z) -> (X/u,Z)`. The normalized
invariants are

    I0 = c4(E)/16,       J0 = c6(E)/32.

A subsequent determinant-one parameter change can avoid a zero divisor in
`z(g)` and a zero leading coefficient. The deterministic sequence first tries
the existing infinity, then `(X,Z)->(kX-Z,X)` for k=0,1,... . Only finitely many
directions are excluded by the branch quartic and its nonzero sextic covariant.
The transformation and its invariants are checked; no isomorphism is inferred
merely from equality of j-invariants.

The completed-square cubic used by the parent Selmer calculation has generator
theta. Its identification with this invariant algebra is explicit:

    phi = -3*theta-b2(E)/4,       theta = -phi/3-b2(E)/12.

For each normalized quartic, convert `z(g)` by this formula and check its exact
square-class coordinate against the stored **ambient Selmer label**. The three
labels must be those of c, d and c+d in the parent's chosen basis. This check
binds the invariant-theory formula to the arithmetic classes, including all
rational-torsion factorization cases. The retained quartic's marking and local
membership evidence are also checked by the representative interface.

The product square root is computed in the theta algebra and transported back
to phi. No field-generator identification, target translation, or basis change
is implicit in this conversion.

## Finite local support

Following Fisher's Remark 3.3, the finite set includes:

- 2, 3, 5 and 7;
- primes in the numerator or denominator of the discriminant of g1 and in
  the denominators of its coefficients;
- primes in a common denominator and the numerator content of gamma;
- primes in the numerator or denominator of `g2(1,0)`.

Outside it, the residue characteristic is odd and at least 11, g1 has smooth
reduction, gamma has nonzero reduction, and the first Hilbert argument is a
unit. Hasse's bound gives a smooth residue point outside the at most four
geometric points above the zero directions of gamma. At the chosen point both
Hilbert arguments are units, so the omitted local factor is +1. The local
factor is independent of the allowed point choice by the pairing theorem.

Support is recomputed by the verifier. Missing or extra substituted support
metadata is not silently accepted. Factorization and primality certification
use the existing exact arithmetic infrastructure.

## Exact local witnesses, including p=2 and the real place

No floating-point or approximate p-adic value is exported. A **rational**
parameter `(X,Z)` where g1 is a nonzero square in Q_p proves the existence of
a point over Q_p with that parameter. Nonzero rational square classes are
decided from valuation parity and unit residues: modulo p for odd p, modulo 8
at p=2, and sign over R. The same data give the Hilbert symbol exactly.

If the first Hilbert argument is already a local square, the local factor is
trivial without selecting a point. Otherwise the finite-place search traverses
a dense sequence in both disjoint projective charts `(t:1)` and `(1:p*t)`,
including infinity, and accepts only a point where g1 is a nonzero local square
and gamma is nonzero. For a smooth locally soluble source there is a point off
the finite excluded divisor; square classes are locally constant there, so
this exhaustive search terminates. This is a reference algorithm, not a claim
of efficient precision growth on difficult models.

At the real place, exact root isolation divides the line into sign intervals.
Three rational samples in each interval avoid the at most two roots of gamma.
Every recorded interval/sample uses exact rationals. The verifier recomputes
the square test, gamma value and Hilbert symbol from the parameter, independently
of the point-search traversal.

An optional `max_local_trials` bounds each local search. Exhaustion returns an
unknown entry, never a guessed zero. A partial pairing refuses to provide a
matrix, rank or radical. In dimension zero or one there are no independent
entries: alternation determines the matrix without a local search.

## Labelled matrix and basis changes

The upper triangle is computed directly; the diagonal is zero and symmetry
follows from alternation over F2. These mathematical identities are used to
assemble the matrix, not advertised as independent experimental discoveries.
`pairing_entry(S,c,d)` and `P.direct_entry(c,d)` independently run the quartic
formula for arbitrary inputs and are available for tests on sums and reversed
arguments. They do not just return the corresponding matrix expression.

`P.evaluate(c,d)` is bilinear evaluation in **S's chosen coordinates**.
`radical_basis` uses those coordinates; `radical_ambient_basis` gives the same
vectors in the parent's finite arithmetic space. The kernel is calculated by
exact F2 elimination. By the Cassels–Tate theorem this is the image of
`Sel^(4)(E) -> Sel^(2)(E)`, not the kernel of multiplication by two inside an
already constructed 4-Selmer group.

`P.with_basis(ambient_labels)` validates a new basis of the same computed space
and performs the explicit congruence transformation. It does not assume that
integer labels from unrelated computations are interchangeable. Natural model
transport uses the parent's exact `transport_to` interface. Tests recompute
pairings on changed models and bases, including a nonzero form with a nonzero
radical.

```python
from elliptic_descent import two_selmer, cassels_tate_pairing

S = two_selmer((0, -1, 1, -929, -10595), proof=True)  # 571a1
P = cassels_tate_pairing(S)
assert P.certified and P.verify_evidence()
assert P.matrix == ((0, 1), (1, 0))
assert P.rank == 2 and P.radical_basis == ()
assert P.evaluate(1, 2) == 1  # Value 1/2 in Q/Z, not the rational number 1.
```

A complete conditional parent yields a pairing with the same assumptions, not
an unconditional result. Exact local evidence can still verify independently
of the class-group proof mode. An incomplete parent has no full labelled basis
and is rejected. Manually constructing result dataclasses is not a certificate.

## Validation and independent paths

The test corpus includes Fisher's published nonzero 571a1 example with exact
quartics and product root, all six ordered pairs of its three nonzero classes,
all split-product root-sign choices, Hilbert reciprocity and independent PARI
Hilbert-symbol comparisons. Further tests cover changed bases and Weierstrass
models, alternate class presentations, zero and low-dimensional groups,
rational Kummer images in the radical, and partial/corrupted evidence.

Higher-dimensional cases include a three-dimensional zero pairing and a
four-dimensional 2-Selmer group with pairing rank two and a two-dimensional
radical. Independent ordinary PARI `ellrank` computations compare the **pairing
rank**; they do not supply any entries of our matrix. PARI implements Fisher's
method too, so this is an independent implementation, not an unrelated formula.

Separate small-case tests run the existing relative solver **without using the
pairing to prune classes**, and compare lift existence with radical membership.
The small published 571a1 forms are assigned to our labels through their cubic
invariants, and their empty relative descents test the nonradical direction.
These are not merely tests of the matrix's internal linear algebra.

## Reproducible evidence and limits

```sh
sage -python -m pytest -q
sage -python benchmarks/cassels_tate.py --verify > cassels-tate-results.jsonl
python benchmarks/replay_pairing.py cassels-tate-results.jsonl
```

The benchmark records all matrices, basis labels, radicals, normalized quartics,
product roots, gamma forms, support primes, local rational evaluation points,
and parameter changes. Parent Selmer computation, pairing construction and
subsequent verification have separate timings. The replay verifies every formula
and local factor and reconstructs the matrix/radical accounting **without Sage**.
It checks label metadata against the recorded basis but does not re-establish
that basis or the underlying class-group computations. Benchmark records never
claim to contain a complete 4-Selmer group.

Trusted dependencies include the cited formula, certified Sage class-group/unit
arithmetic for the parent space, and the exact arithmetic runtime. Square-root
discovery is followed by rational-algebra verification; local search is followed
by independent evaluation checks. This is not formal verification or an external
mathematical review. No live Magma run or cross-system speed superiority is
claimed. Large coefficients and slow exhaustive local searches remain possible.

The next step is **global fibre assembly and equivalence**, with the actual
rational-torsion action and Selmer maps. Knowing the radical alone is insufficient
to distinguish different extensions or determine the 4-Selmer group law.

## Primary references

1. T. Fisher, [*On binary quartics and the Cassels–Tate pairing*](https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf), *Research in Number Theory* 8 (2022), article 74; [arXiv:2208.14977](https://arxiv.org/abs/2208.14977). Theorem 3.1, Remarks 3.2–3.3, and Example 3.4 are the implemented formula, choice/support arguments, and published regression.
2. [PARI elliptic-curve documentation](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html#ellrank): the third `ellrank` output is the pairing rank. It is not a Mordell–Weil rank or a labelled pairing matrix.
3. [Sage number-field element square roots](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field_element.html): the discovery adapter explicitly disallows extensions; every returned root is squared and verified in the rational product algebra.
4. [Labelled 2-Selmer arithmetic](two-selmer.md), [marked quartics](selmer-quartics.md), and the [general bibliography](descent-references.md).
