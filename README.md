# elliptic-descent

Exact, testable algorithms for elliptic-curve descent over the rational numbers.

**v0.3.0 implements certified relative 4-descent for smooth irreducible quartics.**
Given a quartic 2-covering with its canonical map to its Jacobian, it enumerates
all everywhere locally soluble lifts, with explicit quadric equations and maps,
up to the fixed-target equivalence specified below. Complete mode has no
artificial candidate or local-search cutoff.

This is **not** a full curve-level 4-Selmer-group computation, a Mordell–Weil
rank algorithm, or a performance match for Magma. Reducible quartic arithmetic,
full Selmer-group assembly, minimization/reduction, and practical point search
remain separate milestones.

## Install and run

The complete descent requires SageMath. Library files use ordinary Python
syntax; the geometry and local reference kernels also run with Python 3.10+
and SymPy without Sage.

```sh
sage -pip install -e '.[test]'
ELLIPTIC_DESCENT_REQUIRE_SAGE=1 sage -python -m pytest -q
sage -python benchmarks/relative.py --verify --repeats 1
```

```python
from elliptic_descent import relative_four_descent

# Ascending coefficients: -7*x^4 + 12*x^3 + 20*x^2 - 120*x + 172.
result = relative_four_descent((172, -120, 20, 12, -7), proof=True)
assert result.complete and result.certified
assert result.verify_evidence()

print(result.source.ainvs)                 # Canonical Jacobian model
print(len(result.unoriented_coverings))    # 1, as in the published Magma example
print(len(result.coverings))               # 2 maps to the fixed quartic
D = result.coverings[0]
print(D.quadrics)
print(D.map_xyz)                           # Degree-4 map D -> quartic
print(result.source.compose(D))            # Degree-16 map D -> Jacobian
```

### Precisely what is counted

`result.coverings` retains **both signs of the map to the fixed quartic**.
Isomorphisms must commute with the identity of that quartic.
`result.unoriented_coverings` additionally identifies the target involution
`(X:Y:Z) -> (X:-Y:Z)`. For an irreducible quartic these counts differ by a factor
of two. Neither list silently quotients by other target automorphisms or
rational torsion translations. A quartic with merely the same Jacobian is not
a substitute for the specified target.

Inputs are five exact rational coefficients or a `QuarticTwoCover` object.
The canonical quartic-to-Jacobian map is part of that object. An independently
supplied elliptic-curve model requires an explicit isomorphism before its
coordinates can be compared; it is not silently inferred from a label.

### Complete versus bounded calculations

```python
partial = relative_four_descent((1, -1, 0, 0, 1), max_candidates=1,
                                max_nodes=10_000, max_depth=8)
print(partial.enumeration_complete, partial.complete, partial.certified)
print(len(partial.unresolved))
```

With no limits, every finite quotient class and required local place is decided
on the supported smooth domain. Explicit limits can return a partial result;
`unknown` is never converted to insolubility. A returned empty list certifies
an empty relative descent only when `certified` is true. The input quartic
need not itself be locally soluble, so interpreting an empty result in terms
of Sha requires that additional information.

`proof=True` explicitly requests Sage's unconditional number-field arithmetic.
`proof=False` retains a GRH assumption and cannot return `certified=True`.
`verify_evidence()` checks exact geometric/local evidence, not a second
independent class-group computation. Its separate verification budget may be
insufficient; failure to verify is not a proof of solubility.

## Implemented components

| Component | Capability |
|---|---|
| Geometry | Exact norm-pair construction, maps, independent norm-form and polynomial audits |
| Arithmetic | Sage `K(S,2)`, class-group contribution, affine norm fibre, rational-scalar quotient |
| Equivalence | Fixed-target norm-pair equivalence and both orientations; scalar dependencies retained |
| Real place | Complete exact sign test using isolated roots, including zero, two, or four real roots |
| Finite places | Terminating smooth-model search, exact eliminants and linear lifts, verified integral model changes |
| Certificates | Hensel witnesses, independent modular obstruction checks, real signs and model transformations |
| Integration | Canonical Jacobian maps and optional PARI `ell2cover` model/map import |

The earlier bounded routines remain available as independent reference paths.
They are not silently relabelled complete. The geometric constructor handles
squarefree reducible quartics, but the **relative descent** currently requires
irreducibility. Large-prime certification uses Sage/PARI rather than accepting
probable primes. Runtime and memory can still be substantial on difficult inputs.

## Tests and benchmarks

```sh
# Python-only kernels; Sage test modules are explicitly skipped.
python -m pip install -e '.[test]'
python -m pytest -q
python benchmarks/geometry.py --repeats 3 --verify
python benchmarks/foundation.py --repeats 3 --verify

# Complete end-to-end benchmark corpus, including independently checked obstructions.
sage -python benchmarks/relative.py --fixture all --verify --repeats 3
```

CI tests Python 3.10 and 3.13 plus a digest-pinned Sage 10.6 environment.
The Sage job requires its integration tests and exports JSONL benchmarks.
Records separate solver and verification times, arithmetic and local stages,
proof assumptions, quotient counts, unresolved cases, coefficient sizes,
local work, model changes, environment and process peak memory. These are
in-process repetitions, not cold-cache or cross-system speed comparisons.

## Documentation

- [Relative descent: contract and completeness argument](docs/relative-descent.md)
- [v0.3 acceptance and regression corpus](docs/v0.3.md)
- [Construction and conventions](docs/algorithm.md)
- [Finite-support arithmetic](docs/arithmetic.md)
- [Reference local algorithms](docs/local.md)
- [Full project completeness specification](docs/completeness.md)
- [Changelog](CHANGELOG.md)

The next milestone is v0.4, covering every nonsingular quartic factorization
type. The full 4-Selmer group and its group operations are v0.5 work.

## References

The main sources are Womack's
[*Explicit Descent on Elliptic Curves*](https://johncremona.github.io/theses/womack.pdf)
(2003), Fisher's
[*Some improvements to 4-descent*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf)
(2008), Cremona–Fisher–Stoll's
[minimization and reduction algorithms](https://arxiv.org/abs/0908.1741)
(2010), and Fisher–Sills's
[local-solubility and height algorithms](https://arxiv.org/abs/1103.4944)
(2011). Specific theorems and implementation conventions are cited in the docs.

## License

MIT. Copyright (c) 2026 Christopher D. Long.
Sage/PARI remain separate dependencies under their own licenses.
