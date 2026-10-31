# Changelog

## Unreleased — v0.4 development

- Add exact quartic etale arithmetic and factorwise square-class spaces (PR #5).
- Extend complete relative descent to two irreducible quadratic factors, with
  explicit orientation hypotheses and verified scalar-square transformations.
- Add independent eclib Selmer-cardinality comparisons, locally obstructed
  sources, additional exact tests and type-2+2 benchmark artifacts.
- Extend relative descent to all rational-root patterns with an explicit
  norm-minus-one involution, avoiding double-counted orientations.
- Add exact rational-root/Jacobian identifications, inverse point maps and
  verification of the canonical map against doubling.
- Fix touching rational-root isolation intervals and extend independent eclib
  count checks, regression tests and benchmark artifacts.
- Retain version 0.3.0 pending the final v0.4 release/audit step. Full curve-level
  4-Selmer assembly, degree-three input conversion and point search remain separate.

## 0.3.0 — Complete relative descent for irreducible quartics

- Assemble an exhaustive norm-pair descent with the rational-scalar quotient
  and both fixed-target orientations; keep the involution quotient explicit.
- Add canonical quartic-to-Jacobian maps and composition of covering maps.
- Add complete real sign decisions and unbounded p-adic decisions on smooth
  models, with exact finite-field elimination and linearized lifting.
- Resolve difficult reductions through verified integral coordinate/equation
  changes; retain independent positive and negative certificate verification.
- Add a fast exact construction audit backed by the published construction
  theorem; retain generic Groebner audits as an independent reference.
- Extend certified primality beyond 64 bits through Sage/PARI.
- Add published positive/empty descents, torsion and quotient-dependency cases,
  tamper/budget tests, end-to-end benchmarks, and reproducible arithmetic exports.
- Document the supported domain, completeness argument, equivalence, and trusted
  arithmetic dependencies. Full curve-level 4-Selmer computation remains future work.

## 0.2.0 — Arithmetic and local foundations

- Add finite-support square classes, exact affine norm equations, explicit norm
  witnesses, and bounded candidate enumeration with a separate exhaustion flag.
- Add optional Sage `K(S,2)` and PARI `ell2cover` adapters; preserve original
  quartic models, maps, proof settings, and class-group contributions.
- Add exact p-adic and real certificates, independent verification, explicit
  unknown outcomes, and a conservative finite good-reduction prime bound.
- Add a candidate-to-local-screening pipeline that never claims complete
  relative descent and never discards unresolved candidates.
- Add arithmetic/local/backend regression tests and dedicated Sage CI.
- Add JSONL foundation benchmarks, proof/scope metadata, and CI artifacts.
- Document support bounds, local correctness arguments, dependencies, and
  outstanding v0.3 obligations.

## 0.1.x — Geometric reference implementation

- Construct candidate quadric intersections from supplied quartic square classes.
- Retain explicit maps and verify exact algebraic identities, projective
  smoothness, and absence of map basepoints.
- Establish conservative status bookkeeping and formal completeness gates.
