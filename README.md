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

The following are the principal mathematical sources, along with related
algorithms and software references. The descriptions indicate their relevance
to this project; citing a work does not imply that all of its algorithms are
implemented. For theorem-level citations and our coordinate conventions, see
[relative descent](docs/relative-descent.md),
[arithmetic](docs/arithmetic.md), and
[local solubility](docs/local.md).

### Foundational 4-descent

- **S. Siksek (1995),**
  [*Descents on Curves of Genus 1* (PhD thesis)](https://samirsiksek.github.io/siksek.github.io/papers/phdnew.pdf).
  Early development of explicit higher descent on genus-one curves.
- **J. R. Merriman, S. Siksek and N. P. Smart (1996),**
  [*Explicit 4-descents on an elliptic curve*](https://samirsiksek.github.io/siksek.github.io/papers/fourdesc.pdf),
  *Acta Arithmetica* **77**(4), 385–404.
  Original published explicit 4-descent construction
  ([DOI](https://doi.org/10.4064/aa-77-4-385-404)).
- **J. E. Cremona (1997),**
  [*Higher Descents on Elliptic Curves*](https://johncremona.github.io/papers/d2.pdf).
  Expository account of second descent, including 2-isogeny cases.
- **T. Womack (2003),**
  [*Explicit Descent on Elliptic Curves* (PhD thesis)](https://johncremona.github.io/theses/womack.pdf).
  Detailed practical 4-descent, quartic algebras, local tests, quadric
  intersections, reduction and point searches.
- **S. Stamminger (2005),**
  [*Explicit 8-Descent on Elliptic Curves* (PhD thesis)](https://www.mathe2.uni-bayreuth.de/stoll/papers/Stamminger-8descent.pdf).
  Further descent building on 4-coverings and explicit higher Selmer groups.
- **T. A. Fisher (2008),**
  [*Some improvements to 4-descent on an elliptic curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf),
  *Algorithmic Number Theory* (ANTS VIII), LNCS **5011**, 125–138.
  Reducing repeated class-group and unit-group calculations.

### Algorithms, arithmetic and genus-one models

- **J. E. Cremona (1997),**
  [*Algorithms for Modular Elliptic Curves*, 2nd ed.](https://johncremona.github.io/book/fulltext/index.html).
  Elliptic-curve algorithms, 2-descent and computational foundations.
- **J. E. Cremona (2001),**
  [*Classical invariants and 2-descent on elliptic curves*](https://doi.org/10.1006/jsco.1998.1004),
  *Journal of Symbolic Computation* **31**, 71–87.
  Binary-quartic invariants and the input 2-covering theory.
- **E. F. Schaefer and M. Stoll (2004),**
  [*How to do a p-descent on an elliptic curve*](https://www.mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf),
  *Transactions of the AMS* **356**, 1209–1231.
  Selmer groups via étale algebras, units and class groups.
- **J. E. Cremona and T. A. Fisher (2009),**
  [*On the equivalence of binary quartics*](https://johncremona.github.io/papers/quartequiv.pdf),
  *Journal of Symbolic Computation* **44**, 673–682.
  Equivalence testing, including exceptional and reducible cases.
- **J. E. Cremona, T. A. Fisher and M. Stoll (2010),**
  [*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic curves*](https://arxiv.org/abs/0908.1741),
  *Algebra & Number Theory* **4**(6), 763–820.
  Exact model minimization and reduction, including residue characteristic two.
- **T. A. Fisher and G. F. Sills (2012),**
  [*Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944),
  *Mathematics of Computation* **81**, 1635–1662.
  Local-solubility algorithms, explicit covering maps and height bounds
  (preprint issued in 2011).
- **J. E. Cremona, T. A. Fisher, C. O'Neil, D. Simon and M. Stoll (2015),**
  [*Explicit n-descent on elliptic curves. III. Algorithms*](https://arxiv.org/abs/1107.3516),
  *Mathematics of Computation* **84**, 895–922.
  General explicit-descent infrastructure, particularly relevant to future
  descent degrees.
- **T. A. Fisher (2015),**
  [*Higher descents on an elliptic curve with a rational 2-torsion point*](https://arxiv.org/abs/1509.03234).
  Specialized higher-descent methods when rational 2-torsion is available.
- **T. A. Fisher (2022),**
  [*On binary quartics and the Cassels–Tate pairing*](https://doi.org/10.1007/s40993-022-00376-z),
  *Research in Number Theory* **8**, article 74.
  Explicit 2-Selmer pairing formula useful for higher-descent rank bounds.

### Software, documentation and further bibliographies

- [**SageMath:** number-field arithmetic and `selmer_space`](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html) —
  certified square-class computations and the Sage integration layer.
- [**PARI/GP:** elliptic-curve functions](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html) —
  especially `ell2cover` and `ellrank`.
- [**eclib / mwrank**](https://github.com/JohnCremona/eclib) —
  John Cremona's open-source 2-descent and Mordell–Weil software.
- [**Magma handbook:** elliptic-curve descent](https://magma.maths.usyd.edu.au/magma/handbook/elliptic_curves_over_q_and_number_fields) —
  4-descent APIs and published regression examples; not an execution
  dependency.
- [**John Cremona:** papers](https://johncremona.github.io/papers/index.html)
  and [theses](https://johncremona.github.io/theses/index.html) —
  original papers, dissertations and computational references.
- [**Tom Fisher:** research publications](https://www.dpmms.cam.ac.uk/~taf1000/research.html) —
  papers on higher descents, genus-one models, invariant theory and pairings.

## License

MIT. Copyright (c) 2026 Christopher D. Long.
Sage/PARI remain separate dependencies under their own licenses.
