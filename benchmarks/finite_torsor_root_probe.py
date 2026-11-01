"""Exact root-discovery probe for the rational finite-torsor group backend.

These polynomials were exported from an actual contracted product of covering
origin algebras for y^2=x^3-x. This probe checks root discovery only, not the
geometric group law; all resulting isomorphisms require rational verification.
"""
import json
from time import perf_counter
from sage.all import QQ, PolynomialRing, NumberField

COEFFICIENTS = (
440694610019959252569146555782136956578732251075042246957076483405633089,
887006634561692488396403055888126095606281733718456921741946758800,
-1971458719618126144737552321968171644740250330828598967431704536,
-2329245699335486583571705458711161866964620402453040418224,
4237712661559359150467837814251381069349416112369869596,
-100057391482200294278622231000506388482034844912,
-4283433103290910966686968696800552793115317096,
4928963315519854254738471828762586424400,
2765471145518108074420974425008675718,
-5169987758499548877753027504336,
-2041661436762140510487408872,
-869411160184765114000,
2883858450691031964,
-3373517031632,
99712168,-16,1)

start = perf_counter()
R = PolynomialRing(QQ, 't')
t = R.gen()
g = t**16+88*t**14+1180*t**12+3560*t**10+108166*t**8+752232*t**6-3057124*t**4+52469464*t**2+68674369
assert g.is_irreducible()
K = NumberField(g, 'a')
S = PolynomialRing(K, 'x')
f = S(list(COEFFICIENTS))
roots = f.roots(multiplicities=False)
assert all(f(r) == 0 for r in roots)
print(json.dumps({'seconds': perf_counter()-start, 'root_count': len(roots),
                  'source_polynomial': [str(c) for c in COEFFICIENTS],
                  'field_polynomial': [str(c) for c in g.list()],
                  'roots': [[str(c) for c in r.list()] for r in roots],
                  'scope': 'exact field-root discovery only; not a Selmer certificate'}, sort_keys=True))
