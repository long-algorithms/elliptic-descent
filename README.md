# elliptic-descent

**v0.5.0: explicit 4-Selmer groups and marked coverings over the rationals.**

The curve-level API computes the labelled 2-Selmer group, Cassels–Tate pairing,
global 4-Selmer classes and their group law. It returns generators, invariant
factors, relations, coordinates and explicit covering maps to the original
elliptic curve. Exact evidence, arithmetic assumptions and partial results are
retained throughout. See the [release contract and audit](docs/v0.5.md).

This is an exact reference implementation, not an exact Mordell–Weil rank or
saturation algorithm, a general minimization/point-search package, or a claimed
performance match for Magma. Difficult number fields and local calculations can
require substantial time and memory.

## Install and run

Full descent requires SageMath. From a checkout of the release:

```sh
sage -pip install '.[test]'
ELLIPTIC_DESCENT_REQUIRE_SAGE=1 sage -python -m pytest -q
```

A downloaded release wheel can instead be installed with
`sage -pip install ./elliptic_descent-0.5.0-py3-none-any.whl`.
The optional `fast` extra enables exact FLINT acceleration; it does not change
proof requirements. Pure geometry, evidence replay and many tests also run with
Python 3.10+ and SymPy without Sage. The package is not being published to PyPI
by this release process.

Run examples under Sage's Python, not its preparser:

```python
from elliptic_descent import four_descent

R = four_descent((0, 0, 0, 4, 0), proof=True)  # y^2 = x^3 + 4*x
G = R.require_group(certified=True)
assert G.invariant_factors == (4,)
g = G.generators[0]
assert g.order() == 4 and 4*g == G.zero
assert 2*g == G.include_two_selmer(G.projection(g))
assert R.covering(g.key) is g.covering
print(g.covering.target)  # The user's original elliptic-curve model.
```

## Stages, reuse and certification

`DescentSession` reuses arithmetic contexts, completed decisions and constructed
coverings. `finished` means the requested stage completed; `complete` means the
full group was constructed; `certified` additionally requires unconditional
arithmetic and a fresh evidence check. Verification verdicts are not cached.

```python
from elliptic_descent import DescentSession, DescentLimits

session = DescentSession((0, 0, 0, -1, 1), proof=True)
partial = session.run(limits=DescentLimits(two_selmer_candidates=0))
assert not partial.complete
result = session.run()  # Resume with no artificial candidate limit.
```

Explicit limits preserve partial ledgers; unknown is never silently replaced
by insolubility. `proof=False` retains any number-field assumptions. Later
unconditional work cannot erase an earlier conditional dependency. Sessions
are in-memory workspaces, not durable checkpoints or thread-safe global caches.
See [the public API](docs/public-api.md) for the full status and budget contract.

## Relative descent remains available

The lower-level solver handles every squarefree **degree-four** quartic over Q:

```python
from elliptic_descent import relative_four_descent

R = relative_four_descent((172, -120, 20, 12, -7), proof=True)
assert R.complete and R.certified and R.verify_evidence()
D = R.coverings[0]
print(D.quadrics)
print(D.map_xyz)  # Degree-four map to the specified intermediate quartic.
```

Here `coverings` counts diagrams over a fixed quartic, not global classes over
an elliptic curve. For factor types `4` and `2+2`, both map orientations are
retained; a rational-root involution identifies them in the other three types.
`unoriented_coverings` is the further quotient by the quartic involution.
The curve-level API separately performs the verified rational-torsion orbit
quotient and exact-sequence accounting. See [global fibres](docs/four-selmer-fibres.md)
and the retained [v0.4 relative contract](docs/v0.4.md).

Cubic presentation conversion is explicit through `CubicToQuartic`; a cubic is
not silently treated as a quartic algebra. Singular curves, inexact inputs and
unsupported base fields are rejected.

## Validation and benchmarks

```sh
sage -python benchmarks/four_selmer_group.py --verify > groups.jsonl
python benchmarks/replay_four_group.py groups.jsonl
sage -python benchmarks/independent_mw.py --fixture all
sage -python benchmarks/public_api.py --verify
```

CI retains the complete relative, arithmetic, local, pairing, fibre, group and
API regression suites, Python-version and minimum-SymPy checks, and wheel/sdist
audits. Six independently bounded rational quotients check group structure;
four also compare every named covering and sum using rational points and exact
Sage division tests. This is distinct from replaying our own geometric evidence.
Published Magma examples and eclib/PARI comparisons remain in the corpus; no
new live Magma comparison or external mathematical review is claimed.

Release assets include the wheel, source distribution, checksums and a
single-revision verification bundle. Publication requires every named workflow,
all expected benchmark fixtures and the installed-distribution audit to pass.
The [release guide](docs/v0.5.md) distinguishes independent checks, shared-kernel
comparisons and portable evidence replay. Passing tests are not formal
verification or an independent computation of every number-field class group.

## Documentation

- [v0.5 release contract, independent checks and audit](docs/v0.5.md)
- [Public API and certification](docs/public-api.md)
- [Labelled 2-Selmer arithmetic](docs/two-selmer.md), [quartic representatives](docs/selmer-quartics.md) and [pairing](docs/cassels-tate.md)
- [Global fibres](docs/four-selmer-fibres.md) and [exact group arithmetic](docs/four-selmer-group.md)
- [Completeness specification](docs/completeness.md), [development milestones](docs/v0.5-plan.md) and [changelog](CHANGELOG.md)

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
