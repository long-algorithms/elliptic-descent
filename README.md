# elliptic-descent

Exact, testable computational tools for elliptic-curve descent, beginning with
explicit 4-descent over the rational numbers. The long-term goal is certified
and efficient 4-descent and reusable algorithms for higher descent, with Sage,
PARI/GP, and eclib interoperability.

**Status: v0.1.0 is an exact geometric prototype, not an implementation of
complete 4-descent.** Given a squarefree quartic `f` and an element `xi` of
`Q[t]/(f)` satisfying a square-norm condition, it computes two quadrics in four
variables and an explicit rational map to `y^2=f(x)`. It does **not** compute
Selmer groups, locally soluble coverings, rank bounds, or Mordell–Weil points.

## Install and run

Requires Python 3.10+ and SymPy; no Magma or Sage installation is necessary for
this first reference kernel. Sage/PARI adapters are planned. This small
backend-neutral implementation is deliberately optimized for auditability,
not speed.

```sh
python -m pip install -e '.[test]'
python -m pytest -q
python benchmarks/geometry.py --repeats 5 --verify
```

```python
from elliptic_descent import construct_four_cover

# f(t)=t^4-t+1, xi=-theta, Norm(xi)=1
D = construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0))
print(D.quadrics)       # homogeneous quadrics q2, q3
print(D.map_xyz)        # (X,Y,Z) for the quartic model
assert D.check_algebra()
assert D.evaluate((1, 0, 0, 0))[0] == (0, 0)
```

See [construction and conventions](docs/algorithm.md) for the norm identity,
map, limitations, and references. The project now has a
[completeness specification and release gates](docs/completeness.md), including
explicit distinction between partial, conditional, and certified results.
The `status.py` bookkeeping types reserve these distinctions for future descent
backends; they do **not** calculate or certify Selmer groups. The geometric
prototype additionally provides optional exact `check_smoothness()` and
`check_map_basepoint_free()` checks over an algebraic closure.

## Development plan

1. **v0.1.x** — exact quartic-algebra geometry, candidate quadrics, map
   identities, projective smoothness/basepoint checks, negative tests,
   status bookkeeping, and JSON microbenchmarks (current prototype).
2. **v0.2** — square-class enumeration over finite supports, norm-space
   computation, certified class-group/units dependencies, PARI/Sage adapters,
   and exact local solubility with honest unknown outcomes.
3. **v0.3** — complete certified relative 4-descent for irreducible quartics.
4. **v0.4** — all nonsingular quartic factorization types.
5. **v0.5** — full curve-level 4-Selmer group, maps, and Cassels pairing.
6. **v1.0** — practical minimization, reduction, point search, and benchmarks.

See the [formal release gates](docs/completeness.md) for precise acceptance
criteria and required independent tests. These milestones are not timelines.

Correctness, local solubility, completeness, and arithmetic proof assumptions
will be represented separately. Benchmark fixtures will distinguish timings of
geometry, global arithmetic, local filtering, and point searches.

## References

- Womack, *Explicit Descent on Elliptic Curves* (2003),
  https://johncremona.github.io/theses/womack.pdf
- Fisher, *Some improvements to 4-descent on an elliptic curve* (2008),
  https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf
- Cremona–Fisher–Stoll, *Minimisation and reduction of 2-, 3- and 4-coverings*
  (2010), https://arxiv.org/abs/0908.1741
- Fisher–Sills, *Local solubility and height bounds for coverings* (2011),
  https://arxiv.org/abs/1103.4944

## License

MIT. Copyright (c) 2026 Christopher D. Long.
