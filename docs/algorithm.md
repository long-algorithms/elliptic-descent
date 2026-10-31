# Initial exact construction

The present implementation is **not full 4-descent**: it constructs the
algebraic pair of quadrics associated to a supplied square-class representative.
It does not impose Selmer conditions, compute class groups, test local
solubility, minimize, reduce, or enumerate all coverings.

Let `f(t) = a*t^4 + ...` be squarefree over Q and `L = Q[t]/(f)`;
write `theta` for the residue class of `t`. Given `xi in L*` with
`N(xi) = a*r^2` for `r in Q*`, put

```
U = u0 + u1*theta + u2*theta^2 + u3*theta^3
xi * U^2 = q0 + q1*theta + q2*theta^2 + q3*theta^3.
```

Then `q2=q3=0` defines a candidate degree-4 genus-one model in `P^3`.
On this intersection of quadrics the map

```
(X,Y,Z) = (-q0, a*r*N(U), q1)
```

satisfies `Y^2 = F(X,Z)`, where `F(X,Z) = Z^4*f(X/Z)`.
Here `(X,Y,Z)` has weights `(1,2,1)`; the polynomial coordinate functions
have degrees `(2,4,2)` in the `u` variables. A point with all three map
coordinates zero is outside the domain of the displayed rational map and
requires separate treatment before any stronger claim about covering maps.

The identity follows from multiplicativity of the norm and
`N(q0+q1*theta)=F(-q0,q1)/a`. The `check_algebra()` method verifies the
quartic norm identity, the square witness, and the polynomial multiplication
identities. It is *not* a local or global correctness/completeness certificate.

## Conventions and scope

- Arrays are in ascending coefficient order: `f=(c0,...,c4)` and
  `xi=(b0,...,b3)`. No floating-point inputs are admitted.
- The quotient algebra may be a product of fields; squarefreeness ensures it
  is etale. There is no irreducibility requirement for this geometric stage.
- A norm witness is chosen as the nonnegative rational square root by default,
  or supplied explicitly. Choosing the opposite sign reverses `Y`.
- No curve-level `E` or map to an elliptic curve is stored in v0.1.0.
- No assertion is made that a candidate is smooth or everywhere locally soluble.

## Sources

- T. Womack, *Explicit Descent on Elliptic Curves* (PhD thesis, 2003),
  https://johncremona.github.io/theses/womack.pdf
- T. Fisher, *Some improvements to 4-descent on an elliptic curve* (2008),
  https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf
- J. E. Cremona, T. A. Fisher, M. Stoll, *Minimisation and reduction of
  2-, 3- and 4-coverings of elliptic curves* (2010),
  https://arxiv.org/abs/0908.1741
- T. Fisher, G. Sills, *Local solubility and height bounds for coverings
  of elliptic curves* (2011), https://arxiv.org/abs/1103.4944
