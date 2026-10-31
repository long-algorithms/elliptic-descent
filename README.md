# elliptic-descent

Exact, testable algorithms for elliptic-curve descent, beginning with explicit
4-descent over the rational numbers.

**v0.2.0: arithmetic and local-solubility foundations.** The package constructs
candidate 4-coverings, computes finite-support square classes through Sage,
solves their norm conditions, and screens the resulting models locally with
checkable certificates. It does **not yet implement complete relative
4-descent**, compute the full 4-Selmer group, or determine Mordell–Weil ranks.

## Install and test

The geometry, bit-vector arithmetic, and local routines require Python 3.10+
and SymPy. Number-field square classes and PARI imports additionally require
SageMath. All library files use ordinary Python syntax.

```sh
python -m pip install -e '.[test]'
python -m pytest -q
python benchmarks/geometry.py --repeats 3 --verify
python benchmarks/foundation.py --repeats 3 --verify
```

For the complete suite, including the genuine Sage/PARI integration tests:

```sh
sage -pip install -e '.[test]'
ELLIPTIC_DESCENT_REQUIRE_SAGE=1 sage -python -m pytest -q
sage -python benchmarks/foundation.py --sage --repeats 1 --verify
```

The Python-only suite explicitly skips the Sage module. The dedicated CI job
requires Sage and cannot silently skip these tests. CI stores benchmark JSONL
artifacts; experiments and proof details belong in `docs/`, not the root README.

## Construct and check a covering

```python
from elliptic_descent import construct_four_cover, everywhere_local_solubility

# Ascending coefficients: f(t)=t^4-t+1, xi=-theta.
D = construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0))
assert D.check_algebra() and D.check_map_identity()
assert D.check_smoothness() and D.check_map_basepoint_free()
local = everywhere_local_solubility(D, point=(1, 0, 0, 0))
assert local.status.value == 'soluble'
assert all(report.verify(D) for report in local.reports)
```

## Enumerate a finite norm fibre and screen it

Run this example with `sage -python`:

```python
from elliptic_descent.backends.sage import sage_square_classes
from elliptic_descent import screen_norm_candidates

context = sage_square_classes((1, -1, 0, 0, 1), proof=True)
result = screen_norm_candidates(context.space, limit=128,
                                max_depth=6, max_nodes=100_000)
print(result.arithmetic.total_count)
print(result.arithmetic.norm_fibre_exhausted)
print(len(result.retained), len(result.unresolved))
assert not result.relative_descent_complete
```

The last assertion is intentional: exhausting a norm fibre is **not** the same
as exhausting the relative Selmer set with the correct equivalences and
orientation data. A smaller user-supplied support is allowed for experiments
but does not inherit a descent-completeness claim.

## Implemented scope

| Component | v0.2 capability |
|---|---|
| Geometry | Exact quartic-algebra construction, maps, smoothness/basepoint checks |
| Arithmetic | Sage `K(S,2)`, including class-group 2-torsion; affine norm equations; witnesses |
| PARI | Import `ell2cover` basis models/maps; independently check map identities |
| Finite places | Bounded exhaustive residue lifting, strong Hensel witnesses, modular obstructions |
| Real place | Rational-point witnesses and definite-pencil obstructions; otherwise `unknown` |
| All other primes | A conservative good-reduction bound for regular squarefree pencils |
| Verification | Independent certificate checks, negative tests, Sage integration, JSONL benchmarks |

The number-field adapter currently accepts irreducible quartics. The geometric
kernel also handles squarefree reducible quartics. Local searches are bounded;
resource exhaustion is `unknown`, never insolubility. The reference prime
validation/factorization layer supports prime factors below `2**64`.

`proof=True` explicitly requests Sage's certified arithmetic. `proof=False`
retains a GRH assumption. Hand-constructed status objects are bookkeeping,
not mathematical certificates. Returned local witnesses can be checked using
`report.verify(model)` independently of the class-group assumptions.

## Documentation and development

- [Construction and conventions](docs/algorithm.md)
- [Square classes, finite support, and proof dependencies](docs/arithmetic.md)
- [Local algorithms and certificate verification](docs/local.md)
- [Completeness specification and release gates](docs/completeness.md)
- [v0.2 acceptance and reproducibility checklist](docs/v0.2.md)
- [Changelog](CHANGELOG.md)

The next mathematical milestone is **v0.3: complete certified relative
4-descent for irreducible quartics**. Later milestones add all factorization
types, curve-level 4-Selmer groups, and practical minimization/reduction/search.
Performance work retains the slow exact reference paths for differential tests.

## References

The main mathematical sources are Womack's
[*Explicit Descent on Elliptic Curves*](https://johncremona.github.io/theses/womack.pdf)
(2003), Fisher's
[*Some improvements to 4-descent*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf)
(2008), Cremona–Fisher–Stoll's
[minimization and reduction algorithms](https://arxiv.org/abs/0908.1741)
(2010), and Fisher–Sills's
[local-solubility and height algorithms](https://arxiv.org/abs/1103.4944)
(2011). More specific implementation references appear in the documentation.

## License

MIT. Copyright (c) 2026 Christopher D. Long.
Sage/PARI remain separate dependencies under their own licenses.
