# elliptic-descent

Exact, testable algorithms for elliptic-curve descent over the rational numbers.

**Complete relative 4-descent is available for smooth quartics of type 4 or 2+2.**
The package version remains **v0.3.0**; the two-quadratic extension is
the second v0.4 development step, not the completed v0.4 release.
Given a quartic 2-covering with its canonical map to its Jacobian, it enumerates
all everywhere locally soluble lifts, with explicit quadric equations and maps,
up to the fixed-target equivalence specified below. Complete mode has no
artificial candidate or local-search cutoff.

This is **not** a full curve-level 4-Selmer-group computation, a Mordell–Weil
rank algorithm, or a performance match for Magma. Rational-root descent cases,
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
`(X:Y:Z) -> (X:-Y:Z)`. For either supported factorization type these counts
differ by a factor of two. The [two-quadratic proof and regression
tests](docs/two-quadratic-descent.md) explain why this does not extend to
rational-root cases. Neither list silently quotients by other target automorphisms or
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
| Arithmetic | Sage field/product square classes, class-group contribution, norm fibre, diagonal scalar quotient |
| Equivalence | Fixed-target norm-pair equivalence and both orientations; scalar dependencies retained |
| Real place | Complete exact sign test using isolated roots, including zero, two, or four real roots |
| Finite places | Terminating smooth-model search, exact eliminants and linear lifts, verified integral model changes |
| Certificates | Hensel witnesses, independent modular obstruction checks, real signs and model transformations |
| Integration | Canonical Jacobian maps and optional PARI `ell2cover` model/map import |

The earlier bounded routines remain available as independent reference paths.
They are not silently relabelled complete. The geometric constructor handles
squarefree reducible quartics, but the **relative descent** currently supports
only irreducible quartics and products of two irreducible quadratics. Large-prime
certification uses Sage/PARI rather than accepting probable primes. Runtime and memory can still be substantial on difficult inputs.

## Tests and benchmarks

```sh
# Python-only kernels; Sage test modules are explicitly skipped.
python -m pip install -e '.[test]'
python -m pytest -q
python benchmarks/geometry.py --repeats 3 --verify
python benchmarks/foundation.py --repeats 3 --verify

# Complete end-to-end benchmark corpus, including independently checked obstructions.
sage -python benchmarks/relative.py --fixture all --verify --repeats 3
sage -python benchmarks/two_quadratic.py --verify --repeats 3
```

CI tests Python 3.10 and 3.13 plus a digest-pinned Sage 10.6 environment.
The Sage job requires its integration tests and exports JSONL benchmarks.
Records separate solver and verification times, arithmetic and local stages,
proof assumptions, quotient counts, unresolved cases, coefficient sizes,
local work, model changes, environment and process peak memory. These are
in-process repetitions, not cold-cache or cross-system speed comparisons.

## Documentation

- [Relative descent: contract and completeness argument](docs/relative-descent.md)
- [Two-quadratic descent: equivalence and independent checks](docs/two-quadratic-descent.md)
- [v0.3 acceptance and regression corpus](docs/v0.3.md)
- [Construction and conventions](docs/algorithm.md)
- [Finite-support arithmetic](docs/arithmetic.md)
- [Reference local algorithms](docs/local.md)
- [Full project completeness specification](docs/completeness.md)
- [Changelog](CHANGELOG.md)

The next milestone is v0.4, covering every nonsingular quartic factorization
type; rational-root cases remain after the completed two-quadratic step.
The full 4-Selmer group and its group operations are v0.5 work.

## References

These are the **core sources for our 4-descent implementation**, rather than a
comprehensive literature survey. See the separate
[**descent-algorithms bibliography**](docs/descent-references.md) for 2-, 3-,
5-, 6-, 7-, 8-, 9- and 12-descent, isogeny descents, arithmetic and software.
The proof and coordinate conventions are in
[the relative-descent documentation](docs/relative-descent.md).

- **S. Siksek (1995)**,
  [*Descents on Curves of Genus 1* (PhD thesis)](https://samirsiksek.github.io/siksek.github.io/papers/phdnew.pdf).
  Geometric foundations of genus-one higher descent.
- **J. R. Merriman, S. Siksek and N. P. Smart (1996)**,
  [*Explicit 4-descents on an elliptic curve*](https://samirsiksek.github.io/siksek.github.io/papers/fourdesc.pdf),
  *Acta Arithmetica* **77**, 385–404. Original 4-descent construction.
- **T. Womack (2003)**,
  [*Explicit Descent on Elliptic Curves* (PhD thesis)](https://johncremona.github.io/theses/womack.pdf).
  Detailed quartic-algebra constructions and algorithms.
- **E. F. Schaefer and M. Stoll (2004)**,
  [*How to do a p-descent on an elliptic curve*](https://mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf),
  *Trans. Amer. Math. Soc.* **356**, 1209–1231.
  Selmer groups, finite support, units and class groups.
- **T. A. Fisher (2008)**,
  [*Some improvements to 4-descent on an elliptic curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf),
  *ANTS VIII*, LNCS **5011**, 125–138. Less redundant global arithmetic.
- **J. E. Cremona, T. A. Fisher and M. Stoll (2010)**,
  [*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic curves*](https://arxiv.org/abs/0908.1741),
  *Algebra & Number Theory* **4**, 763–820. Small integral models.
- **T. A. Fisher and G. F. Sills (2012)**,
  [*Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944),
  *Math. Comp.* **81**, 1635–1662.
  Local conditions, maps and point searches (2011 preprint).

## License

MIT. Copyright (c) 2026 Christopher D. Long.
Sage/PARI remain separate dependencies under their own licenses.
