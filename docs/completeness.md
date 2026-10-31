# What “complete 4-descent” means

This is a specification and release-gate checklist, **not a claim about the
current implementation**. In v0.1.x only the supplied-square-class geometric
constructor and its exact checks are implemented. None of the Selmer,
local-solubility, equivalence, or exhaustive-enumeration stages exists yet.

## Objects, levels, and equivalence

**Geometric prototype.** For squarefree `f` and a supplied `xi` of square norm,
construct the two quadrics `q2=q3=0` and polynomial map to `y²=f(x)`.
Exact algebraic identities establish that the map lands on the quartic;
projective smoothness and basepoint-freeness are *separate* checks. This is not
a Selmer or local-solubility calculation.

**Relative 4-descent (first mathematical goal).** Input a smooth, everywhere
locally soluble **2-covering with its specified map to E**, not merely a quartic
with the same Jacobian. Enumerate all relevant everywhere locally soluble
4-coverings above it, with their covering maps, up to a **declared equivalence
of diagrams of coverings**. The chosen quotient (including any torsion action,
known-point identifications, or map automorphisms) must be specified; an
unlabelled list of quadric intersections is insufficient.

**Full 4-Selmer descent.** Determine `Sel^(4)(E/Q)` as a finite abelian group,
with labelled elements, representatives, maps, the doubling map to
`Sel^(2)(E/Q)`, and exact accounting for rational 2- and 4-torsion. A collection
of descents for arbitrarily chosen 2-covering *basis elements* need not give a
full Selmer computation: combinations must be accounted for.

**Mature software.** Minimize/reduce models, search for rational points, and
provide robust performance on difficult inputs; these do not change the
mathematical definition of a complete descent. A locally soluble 4-covering
may have no rational point, so failure of point search is not a failed
4-descent calculation.

For a complete relative descent on its declared domain, an implementation must
prove the following for each supported valid input:

1. **Soundness.** Every reported class defines the claimed smooth genus-one
   covering with correct maps, and satisfies all local Selmer conditions.
2. **Exhaustiveness.** Every relevant locally soluble lift is represented,
   with a complete finite-support and square-class argument, certified arithmetic
   (units, class groups, norm conditions), and justified finite local sieve.
3. **Nonredundancy.** Equivalence testing uses the declared quotient correctly.
4. **Evidence.** Each positive/negative local decision, discarded candidate,
   arithmetic dependency, and global enumeration/exhaustiveness claim is
   backed by recorded computation or a cited proven reduction. An inability to
   decide is **unknown**, never locally insoluble.
5. **Termination in the supported domain.** Complete mode cannot silently
   truncate enumeration or a local search. A resource limit gives a partial
   result with explicitly unresolved work.

A *fake Selmer set*, an approximate norm-compatible set, or a list of valid
coverings does not meet the exhaustiveness requirement without proof that it
coincides with the desired relative Selmer set.

## Results must distinguish three independent issues

- **Local verdict:** `unknown`, `soluble`, `insoluble`.
- **Proof strength:** `unverified`, `conditional`, `certified`; document exact
  dependencies (including GRH-conditional class-group computations).
- **Exhaustiveness and equivalence:** whether all candidate classes have been
  enumerated and whether the declared quotient is completely accounted for.

`DescentResult` in `status.py` is an initial *decision ledger*, not an
algorithm or proof-checking kernel. Its `complete` and `certified` properties
are derived conservatively from the above flags. For now evidence fields are
human-readable references; future backends must attach verifiable certificates
and must set a flag only after the required mathematical proof has been
established. Manual construction of a `DescentResult` is **not** a certificate.

In particular:

| Candidate list | Enumeration | Decisions | Equivalence | Status |
|---|---|---|---|---|
| empty | partial/unknown | — | any | incomplete |
| empty | exhaustive, certified | — | checked, certified | complete, certified |
| nonempty | exhaustive | some unknown | checked | incomplete |
| decided | exhaustive, conditional | decided, conditional | checked | complete, conditional |
| decided | exhaustive, certified | decided, certified | checked, certified | complete, certified (if no assumptions) |

`coverings=[]` with `complete=False` must never be described as an empty
4-descent. Inputs not supported (e.g. a missing quartic factorization type)
should report `unsupported` at the algorithm boundary, not fabricate a
complete empty list.

## Mathematical consistency checks (global stage)

Let `B` denote the Cassels pairing on the 2-Selmer group and `R` its radical.
For the standard isogeny cohomology sequence, the image of
`Sel^(4)(E) -> Sel^(2)(E)` is `R`, and there is an exact segment

    E(Q)[2] --delta--> Sel^(2)(E) --> Sel^(4)(E) --> Sel^(2)(E).

Consequently, with the usual map conventions,

    |Sel^(4)(E)| = (|Sel^(2)(E)| / |im(delta)|) * |R|.

Care is required: `delta` need not be injective when E has rational 4-torsion.
This count is a necessary consistency check, not an independent proof of the
correct covering representatives or group law. The exact maps, pairing,
quotients, and group structure must also be checked. The Selmer group may be
nontrivial even if particular nontrivial 2-coverings have empty descents.

Complete 4-descent does **not** imply a proof of the rank of `E(Q)`, a full
Mordell–Weil basis, or triviality of `Sha(E)[4]`.

## Release gates

- **v0.1.x — geometric baseline (current):** quartic algebra and norm witness,
  quadric equations, map identity, exact model-smoothness and basepoint checks,
  negative inputs, preliminary timing harness. Not complete descent.
- **v0.2 — arithmetic and local foundation:** finitely supported square classes
  and norm constraints, local solubility with `unknown` outcomes, Sage/PARI
  adapter, and explicit arithmetic proof assumptions. No completeness claim
  outside families with a documented proven bound.
- **v0.3 — complete irreducible-quartic relative 4-descent:** soundness,
  exhaustiveness, nonredundancy, and unconditional certification for supported
  nonsingular irreducible quartics; positive and empty descent tests.
- **v0.4 — every nonsingular quartic factorization type:** `4`, `3+1`, `2+2`,
  `2+1+1`, `1+1+1+1` with consistent covering equivalence and torsion cases.
- **v0.5 — curve-level full 4-Selmer:** group, doubling, Cassels pairing,
  exact torsion bookkeeping, corresponding coverings for all relevant classes.
- **v1.0 — practical full 4-descent over Q:** all previous gates, plus robust
  minimization/reduction, point search, and reproducible certified benchmarks.

These are *acceptance criteria*, not target dates. No release may be described
as achieving one of these gates unless it has corresponding proofs, passing
regressions, and a documented supported-input domain.

## Regression and performance requirements

1. Test norm identities and maps against independent exact calculations;
   verify nonsingularity and absence of map basepoints over an algebraic closure.
2. Generate coordinate changes, alternate field generators, rational-square
   rescalings, and generalized quartic models. Ensure maps commute. Inject
   deliberately invalid witnesses and degenerate quadrics.
3. For global arithmetic, include nontrivial class-group 2-torsion, large units,
   multiple quartic splitting types, rational 2-/4-torsion, and unconditional
   versus GRH-conditional computation.
4. Include positive descents, genuinely empty descents (e.g. published
   `571a1` covering fixtures), nontrivial `Sha[2]` and `Sha[4]`, prime 2
   singular reductions, and nonterminating/budget-exhausted local searches.
   Every published fixture must specify *its actual covering map*, not just E.
5. Compare against Magma on mathematically equivalent tasks, not identical
   coefficient lists. Compare per-covering equivalence, maps, known counts,
   and group-level invariants; use independent arithmetic where possible.
6. Benchmark kernels, fixed-covering descent, and curve-level descent
   separately, with proof assumptions, cold/warm caches, time, memory,
   candidate counts, coefficient growth, random seeds, and environment.
   Never compare a partial result with an unconditional complete one as if
   they represented the same task.

In v0.1.x the performance runner measures geometry only. Its `--verify` flag
checks the algebra, map identity, smoothness, and basepoint freeness
separately from the construction timing; none of these certify local
solubility or complete 4-descent. Later versions should retain the slow exact reference path so every
optimized implementation can be differential-tested.

## References

- Merriman–Siksek–Smart, *Explicit 4-descents on an elliptic curve* (1996).
- T. Womack, *Explicit Descent on Elliptic Curves* (2003),
  <https://johncremona.github.io/theses/womack.pdf>.
- T. Fisher, *Some improvements to 4-descent on an elliptic curve* (2008),
  <https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf>.
- Cremona–Fisher–Stoll, *Minimisation and reduction of 2-, 3- and 4-coverings*
  (2010), <https://arxiv.org/abs/0908.1741>.
- Fisher–Sills, *Local solubility and height bounds for coverings of elliptic
  curves* (2011), <https://arxiv.org/abs/1103.4944>.
- T. Fisher, *On binary quartics and the Cassels–Tate pairing* (2022),
  <https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf>.
