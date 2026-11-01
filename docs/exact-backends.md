# Exact arithmetic backends

The portable checker works with SymPy's standard rational polynomial domains.
Install `.[fast]` to let SymPy use python-flint for the same exact rational
polynomial arithmetic. This does not replace any test by floating-point
arithmetic, change the Galois certificate, omit conjugates or cache a proof
verdict. The nonaccelerated path remains supported and tested.

```sh
sage -pip install -e '.[test,fast]'
SYMPY_GROUND_TYPES=flint sage -python benchmarks/four_selmer_group.py --verify
python benchmarks/compare_exact_backends.py
```

The differential runner starts separate Python and FLINT processes, verifies
an actual degree-16 covering certificate with both, compares every labelled
cocycle and all 16 class sums, and rejects the same altered signed model. Group
CI jobs use python-flint 0.8.0; the ordinary Python, minimum-SymPy and earlier
Sage regression jobs retain their existing environments. The additional FLINT
Python job runs the entire Python suite as well. Benchmark records identify
the actual ground domain and python-flint version, so timing comparisons do
not conceal that change. Both fresh discovery and independent evidence replay
remain required, with the same job timeouts and mathematical assertions.

References for the accelerator: [SymPy polynomial domains](https://docs.sympy.org/latest/modules/polys/domainsintro.html)
and [python-flint](https://python-flint.readthedocs.io/en/latest/).

The signed-product and branch-evaluation checks also clear their nonzero
field denominators rather than computing dense inverses. Nonvanishing is
verified first, so these are equivalent polynomial identities. A regression
forbids nonrational division during covering verification and separately
compares both forms of each identity on actual degree-16 covering data.
