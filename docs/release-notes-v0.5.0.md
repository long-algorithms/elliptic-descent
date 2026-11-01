# elliptic-descent v0.5.0

Explicit 4-Selmer groups over Q, with marked covering representatives and exact
evidence. This completes the eight-stage v0.5 implementation program.

## Included

- A curve-level `four_descent` API and reusable `DescentSession`.
- Labelled 2-Selmer classes, arbitrary marked quartic representatives and the
  Cassels–Tate pairing with exact local contributions.
- Global 4-Selmer fibre assembly with rational-torsion equivalence witnesses.
- Addition, inverses, orders, generators, relations, invariant factors and
  coordinates on the actual marked covering classes.
- Original-model preservation, exact-sequence maps, explicit conditional and
  partial states, and optional exact FLINT acceleration.

## Release validation

All required workflows run on the tagged commit. The attached audit report
contains the actual Sage test totals and complete benchmark inventory. Six
independent rational-quotient audits check group structure; four identify every
named class and sum using rational covering points and exhaustive rational
4-division. Existing eclib/PARI and published Magma comparisons, exact evidence
replays, API reuse tests and installed-distribution audits remain required.

The wheel, sdist, checksums and verification bundle all refer to this release
revision. See `docs/v0.5.md` for the supported domain, independent-check scope,
proof dependencies and reproducibility instructions.

## Installation and limitations

Install the wheel with Sage's Python, for example:

```sh
sage -pip install ./elliptic_descent-0.5.0-py3-none-any.whl
```

Sage is required for full descent; ordinary Python supports the exact reference
geometry and portable replay. This release is not a general point-search,
minimization or exact rank/saturation package. It can be expensive on difficult
inputs. Certified number-field arithmetic remains a Sage dependency; no formal
verification, external mathematical review, live Magma comparison or performance
superiority is claimed. This GitHub release does not publish to PyPI.
