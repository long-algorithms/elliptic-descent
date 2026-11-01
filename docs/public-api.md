# Curve-level public API — v0.5 development step 7

The development API accepts an exact rational elliptic curve and coordinates
the labelled 2-Selmer computation, Cassels–Tate pairing, global fibre assembly,
and explicit 4-Selmer group law. The package remains **0.4.0** pending the final
v0.5 release audit. The lower-level interfaces remain available unchanged.

## One call, with all intermediate evidence retained

Run under Sage's Python, not the Sage preparser:

```python
from elliptic_descent import four_descent

R = four_descent((0, 0, 0, 4, 0), proof=True)  # y^2 = x^3 + 4*x
G = R.require_group(certified=True)
assert R.complete and R.verify_evidence()
assert G.invariant_factors == (4,)
g = G.generators[0]
assert g.order() == 4 and 4*g == G.zero
assert 2*g == G.include_two_selmer(G.projection(g))
assert R.covering(g.key) is g.covering
assert g.covering.target == R.curve
```

Inputs are exact two/five Weierstrass coefficients, `WeierstrassModel`, or a
Sage elliptic curve over QQ. Rational coefficients and nonminimal models are
supported. The original model is retained; normalization is linked by explicit
origin-preserving maps. Nonrational base fields, singular curves and inexact
coefficients are rejected by the existing model interface. No curve label is
silently resolved through an external database.

`R.selmer`, `R.pairing`, `R.fibres` and `R.group` expose the existing mathematical
objects. Before a stage is computed its field is `None`. `R.curve` is the
original model. `R.cardinality` and `R.invariant_factors` require a complete group;
a fibre-stage count is deliberately accessed as `R.fibres.cardinality` instead.
The pairing radical is an image, not a substitute for the constructed group.

`R.covering(key)` retrieves a known marked global class without another solve.
A key is `(Sel2 coordinate, orbit index)` in this exact presentation, **not**
an additive coordinate. Use `G.from_coordinates(...)` for invariant-factor
coordinates. Do not mix elements from different group objects. A partial
assembly may contain fully constructed fibres whose representatives are usable,
but their list is not asserted exhaustive.

## Completion, assumptions and verification are different

- `finished` means the **requested stage** completed. `until='two_selmer'`
  can finish without having computed a pairing or any 4-coverings.
- `complete` always means the **full 4-Selmer group** has been constructed.
  An intermediate result never reports full completion.
- `certified` additionally requires unconditional parent arithmetic and a
  fresh successful evidence check. A complete conditional group is not certified.

`verify_evidence()` rechecks the latest available complete stage and its
parents. It can return true for an intermediate or conditional computation;
that does not make it a complete unconditional 4-descent. The existing
finite-place certificate verification budget can itself be insufficient.
Verification failure is not evidence of insolubility.

`completion_status` uses the existing `incomplete`, `conditional`, and
`certified` values for the whole computation. `assumptions` is the union of
all retained stage assumptions. In particular, later `proof=True` arithmetic
cannot remove a conditional parent's class-group assumptions. Rational-only
arithmetic is not automatically marked conditional merely because `proof=False`
was requested. `summary()` supplies JSON-compatible metadata and performs a
fresh certification check; it is not an independent certificate format.

A group object exists only after its geometric identification is verified.
No rank, saturation, rational point on a covering, or exact Tate–Shafarevich
order is inferred from completing this API.

## Staged computation and reuse

```python
from elliptic_descent import DescentSession

session = DescentSession((0, 0, 0, -1, 1), proof=True)
S2 = session.run(until='two_selmer')
assert S2.finished and not S2.complete
C = session.quartic_cover(1)       # Construct only when requested.
assert session.quartic_cover(1) is C
R = session.run()                 # Continue through pairing, fibres and group.
again = session.run()
assert again.group is R.group
assert again.work == R.work       # No new mathematical construction.
```

The four `until` values are `two_selmer`, `pairing`, `fibres`, and `group`.
The default is `group`. `DescentStage` exposes the same values as an enum.
`DescentSession.from_two_selmer(S, proof=...)` accepts an existing complete,
evidence-verified labelled space. It uses that **exact basis**, original model
and arithmetic context. To change the basis, first call `S.with_basis(...)`
and create a separate session. There is no global cache indexed only by a
curve label or j-invariant.

A session reuses the finite cubic arithmetic context, completed candidate
local decisions, canonical quartic presentations, normalized pairing quartics,
completed pairing entries, torsion data, selected fibre presentations,
completed fibre actions/equivalences, and group objects. New budgets retry
unknown decisions; known decisions are not recomputed. The ordinary zero-fibre
reference route and the optimized transported route are cached separately,
as are forced-Galois and automatic group-law calculations.

Only exact data and constructed objects are reused—not the return value of
`verify_evidence()` or an unconditional-certification verdict. Frozen assembly
objects also reuse their immutable class wrappers rather than rebuilding all
covering wrappers on every element lookup. Replacing a dataclass constructs
new caches; modified witnesses are still checked by the existing verifiers.

This is a synchronous in-memory workspace, not a background task, durable
checkpoint, process-wide cache, or thread-safe object. It can retain substantial
memory. The full solver still enumerates the relevant global classes and, in
the general case, constructs their cocycles. On-demand access does **not** mean
that full group computation avoids those necessary representatives.

## Explicit limits and resuming

```python
from elliptic_descent import DescentLimits, DescentSession

session = DescentSession((0, 0, 0, -1, 1))
partial = session.run(limits=DescentLimits(two_selmer_candidates=0))
assert not partial.complete and partial.pairing is None
assert partial.stop_reason == 'two_selmer_incomplete'
R = session.run()                 # Reuse arithmetic, remove the explicit limit.
assert R.complete
```

Each field accepts an exact nonnegative integer or `None` (unbounded):

| Limit | Scope |
|---|---|
| `two_selmer_candidates` | Prefix of cubic norm-kernel candidates |
| `two_selmer_depth`, `two_selmer_nodes` | Existing local search budgets per candidate/place |
| `pairing_trials` | Local rational-witness trials in pairing entries |
| `fibres` | Prefix of intermediate classes selected for fibre assembly |
| `relative_candidates`, `relative_depth`, `relative_nodes` | Existing relative-descent budgets per fibre/place |
| `field_degree` | Maximum accepted common normal-field degree for group arithmetic |

Zero is a real limit, not an alias for unbounded. Budgets are not global CPU,
memory, or wall-clock limits. Completed cached stages/decisions may be returned
without fresh discovery even with a smaller subsequent search budget. Until a
stage completes, a different prefix limit selects that prefix for the returned
snapshot; earlier larger ledgers remain unchanged. A cold session is the way
to reproduce an exact bounded run from scratch.

`field_degree` is an **acceptance bound**, even for reused group evidence. It
is checked after field discovery and may not save discovery time. The automatic
unique exact-sequence shortcuts need no field and can succeed with degree cap
zero. `backend='galois'` explicitly requests the normal-field route.

A stopped calculation retains its earlier data and reports one of
`two_selmer_incomplete`, `pairing_incomplete`, `fibres_incomplete`, or
`field_degree_limit`. A field-degree stop includes the actual degree and cap
in `limit_detail`; it keeps complete fibres but does not fabricate group data.
The low-level `FieldDegreeLimit` exception remains a subclass of `RuntimeError`
for compatibility.

Only this typed resource exception is converted to a partial group-stage
result. Unexpected backend errors, invalid evidence and arithmetic
inconsistencies still raise. `session.last_result` retains the earlier
successfully returned stages for diagnosis. A failed solver is never an empty
Selmer group or an insolubility proof. `prune=False` and
`reuse_zero_fibre=False` retain the independently useful reference routes.

## Performance measurements and validation

`cache_info()` and each result's `work` record actual constructor calls.
Canonical-quartic counts do not include alternative conic parametrizations
tried internally by fibre selection. Cumulative stage times are in `timings`.
They include checks performed by those constructors; warm API calls still
perform consistency checks and explicit verification remains separate work.
No universal speedup, large-input bound, or Magma comparison is claimed.

```sh
sage -python -m pytest -q tests/test_sage_descent_api.py
sage -python benchmarks/public_api.py --verify
```

The new integration corpus compares the API with the existing explicit
pipeline, including actual class-key addition tables, and checks partial
resumption, conditional parents, changed models/bases, field caps and separate
reference routes. These paths share kernels; this is not a second independent
descent algorithm. Unit tests use explicit stubs only for orchestration/error
contracts, never as a substitute for genuine Sage arithmetic.

Four API benchmark fixtures record cold staged time, warm time, construction
counts, original markings, full group invariants and exact reference equality.
They require a zero warm-construction delta, not a fragile timing threshold.
The existing eight group benchmarks, fibre/pairing/descent regressions,
minimum-dependency checks and distribution audits remain enabled.

Certified number-field completeness still relies on Sage. Exact geometric
and cohomological evidence is checked by the existing verifiers. The API adds
no new arithmetic theorem, formal verification, independent class-group proof,
or claim of an externally reviewed release. The final release audit is step 8.

## Mathematical and lower-level interfaces

See [the development plan](v0.5-plan.md), [labelled 2-Selmer arithmetic](two-selmer.md),
[marked quartics](selmer-quartics.md), [the pairing](cassels-tate.md),
[global fibre equivalence](four-selmer-fibres.md), and
[the exact group law](four-selmer-group.md). Those documents retain the
mathematical references and verification obligations. The root
[core references](../README.md#references) and
[broader bibliography](descent-references.md) are unchanged.
