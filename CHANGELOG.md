# Changelog

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
