# Elliptic-curve descent: references and resources

An annotated bibliography covering algorithms **beyond our implemented
4-descent**. This is a research map, not a feature list. The root
[README](../README.md#references) contains only the principal references
directly used by our current 4-descent implementation.

Publication years below refer to final journal publications where known.
Some linked preprints appeared earlier. Each item distinguishes the sort of
descent or auxiliary computation it actually describes.

## General Selmer groups and explicit n-descent

- **Schaefer–Stoll (2004)**,
  [*How to do a p-descent on an elliptic curve*](https://mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf),
  *Trans. AMS* **356**, 1209–1231.
  Reduction of p-Selmer calculations to number-field units and class groups.
- **Cremona–Fisher–O'Neil–Simon–Stoll (2008)**,
  [*Explicit n-descent on elliptic curves, I. Algebra*](https://arxiv.org/abs/math/0606580).
  Cohomological and algebraic framework, particularly for 3-descent.
- **The same authors (2009)**,
  [*Explicit n-descent on elliptic curves, II. Geometry*](https://arxiv.org/abs/math/0611606).
  Genus-one normal curves representing Selmer classes.
- **The same authors (2015)**,
  [*Explicit n-descent on elliptic curves, III. Algorithms*](https://arxiv.org/abs/1107.3516),
  *Math. Comp.* **84**, 895–922. Practical explicit-descent algorithms.

## 2-descent and 2-isogeny descent

- **J. E. Cremona (1997)**,
  [*Higher Descents on Elliptic Curves*](https://johncremona.github.io/papers/d2.pdf).
  Second descent after a 2-isogeny descent.
- **J. E. Cremona (1997)**,
  [*Algorithms for Modular Elliptic Curves*, second edition](https://johncremona.github.io/book/fulltext/index.html).
  Algorithms for Mordell–Weil computation and 2-descent.
- **J. E. Cremona (2001)**, *Classical invariants and 2-descent on elliptic curves*,
  *J. Symbolic Comput.* **31**, 71–87;
  [Cremona's paper archive](https://johncremona.github.io/papers/index.html).
  Binary quartics and their invariants.
- **D. Simon (2002)**,
  [*Computing the rank of elliptic curves over number fields*](https://doi.org/10.1112/S1461157000000668),
  *LMS J. Comput. Math.* **5**, 7–17.
  Practical 2-descent over general number fields.

## 3-descent and 3-isogenies

- **Cremona–Fisher–O'Neil–Simon–Stoll (2008–2015)**,
  the [three explicit n-descent papers](https://arxiv.org/abs/1107.3516)
  above. Full 3-Selmer calculations and explicit plane cubics.
- **H. Cohen and F. Pazuki (2009)**,
  [*Elementary 3-descent with a 3-isogeny*](https://arxiv.org/abs/0903.4963),
  *Acta Arith.* **140**, 369–404. Descent for rational 3-torsion families.
- **M. van Beek and T. A. Fisher (2018)**,
  [*Computing the Cassels–Tate pairing on 3-isogeny Selmer groups via
  cubic norm equations*](https://arxiv.org/abs/1711.02432),
  *Acta Arith.* **185**, 367–396.
  Refined rank bounds without S-unit computations in the pairing algorithm.

## 4-descent and 8-descent

- **S. Siksek (1995)**,
  [*Descents on Curves of Genus 1* (thesis)](https://samirsiksek.github.io/siksek.github.io/papers/phdnew.pdf).
  Early higher-descent foundations.
- **J. R. Merriman, S. Siksek and N. P. Smart (1996)**,
  [*Explicit 4-descents on an elliptic curve*](https://samirsiksek.github.io/siksek.github.io/papers/fourdesc.pdf),
  *Acta Arith.* **77**, 385–404. Original explicit 4-descent.
- **T. Womack (2003)**,
  [*Explicit Descent on Elliptic Curves* (thesis)](https://johncremona.github.io/theses/womack.pdf).
  Extensive computational descriptions.
- **S. Stamminger (2005)**,
  [*Explicit 8-Descent on Elliptic Curves* (thesis)](https://opus.constructor.university/frontdoor/index/index/docId/424).
  Further descent from 4-coverings and 8-Selmer examples.
- **T. A. Fisher (2008)**,
  [*Some improvements to 4-descent on an elliptic curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/fourdesc.pdf),
  *ANTS VIII*, LNCS **5011**, 125–138. Improved global computation.
- **T. A. Fisher (2017)**,
  [*Higher descents on an elliptic curve with a rational 2-torsion point*](https://arxiv.org/abs/1509.03234),
  *Math. Comp.* **86**, 2493–2518. Specialized higher 2-power descent.

## 5-descent, 7-descent and isogenies

- **T. A. Fisher (2000)**,
  [*On 5 and 7 descents for elliptic curves* (thesis)](https://www.dpmms.cam.ac.uk/~taf1000/thesis.html).
  Rational torsion/isogeny families; PhD awarded in 2001.
- **T. A. Fisher (2001)**,
  [*Some examples of 5 and 7 descent for elliptic curves over Q*](https://doi.org/10.1007/s100970100030),
  *J. Eur. Math. Soc.* **3**, 169–201. Specialized descents.
- **T. A. Fisher (2013)**,
  [*Explicit 5-descent on elliptic curves*](https://antsmath.org/ANTSX/fisher/paper.pdf),
  ANTS X proceedings (conference 2012).
  Degree-5 genus-one normal curves as Pfaffian models.
- **T. A. Fisher (2011 preprint)**,
  [*Minimisation and reduction of 5-coverings of elliptic curves*](https://arxiv.org/abs/1112.5131).
  Integral models and lattice reduction for quintic coverings.
- **B. Creutz and R. L. Miller (2012)**,
  [*Second isogeny descents and the Birch and Swinnerton-Dyer conjectural formula*](https://arxiv.org/abs/1105.4018),
  *J. Algebra* **372**, 673–701. Secondary descent for prime-degree isogenies.

## 6-descent and 12-descent

- **T. A. Fisher (2008)**,
  [*Finding rational points on elliptic curves using 6-descent and 12-descent*](https://arxiv.org/abs/0711.3774),
  *J. Algebra* **320**, 853–884.
  Combines 2- and 3-coverings, or 3- and 4-coverings, to search for large points.
- **T. A. Fisher (2014)**,
  [*Minimal models for 6-coverings of elliptic curves*](https://www.dpmms.cam.ac.uk/~taf1000/papers/minmodels6.pdf),
  *LMS J. Comput. Math.* **17** (special issue A), 112–127.
  Combination without field extensions, via pairs of cubics.

## 9-descent and second p-descent

- **B. Creutz (2010)**,
  [*Explicit Second p-Descent on Elliptic Curves* (thesis)](https://opus.constructor.university/frontdoor/index/index/docId/224).
  Second descent on genus-one torsors.
- **B. Creutz (2014)**,
  [*Second p-descents on elliptic curves*](https://arxiv.org/abs/1209.3085),
  *Math. Comp.* **83**, 365–409.
  Practical 9-descent as second 3-descent.

## Models, local solubility, pairings and point searches

- **J. E. Cremona and T. A. Fisher (2009)**,
  *On the equivalence of binary quartics*, *J. Symbolic Comput.* **44**, 673–682;
  [Cremona's archive](https://johncremona.github.io/papers/index.html).
  Identifying equivalent 2-covering models.
- **J. E. Cremona, T. A. Fisher and M. Stoll (2010)**,
  [*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic curves*](https://arxiv.org/abs/0908.1741),
  *Algebra & Number Theory* **4**, 763–820. Efficient integral models.
- **T. A. Fisher and G. F. Sills (2012)**,
  [*Local solubility and height bounds for coverings of elliptic curves*](https://arxiv.org/abs/1103.4944),
  *Math. Comp.* **81**, 1635–1662. Rigorous local tests and covering maps.
- **T. A. Fisher (2012)**,
  [*The Hessian of a genus one curve*](https://www.dpmms.cam.ac.uk/~taf1000/papers/g1hessians.pdf),
  *Proc. London Math. Soc.* **104**, 613–648. Covariants and invariants.
- **T. A. Fisher (2022)**,
  [*On binary quartics and the Cassels–Tate pairing*](https://arxiv.org/abs/2208.14977),
  *Research in Number Theory* **8**, article 74.
  Explicit 2-Selmer pairing, useful for higher-descent rank bounds.

## Software, code and living bibliographies

The following are sources, benchmarks or potential dependencies, **not** claims
that their respective algorithms are implemented in this repository.

- [SageMath: number-field Selmer-space documentation](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html).
- [PARI/GP: elliptic curves and `ell2cover`](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html).
- [eclib/mwrank: John Cremona's C++ code](https://github.com/JohnCremona/eclib).
- [Magma: elliptic-curve descent documentation](https://magma.maths.usyd.edu.au/magma/handbook/elliptic_curves_over_q_and_number_fields).
- [John Cremona's papers](https://johncremona.github.io/papers/index.html)
  and [theses](https://johncremona.github.io/theses/index.html).
- [Tom Fisher's papers](https://www.dpmms.cam.ac.uk/~taf1000/research.html)
  and [thesis collection](https://www.dpmms.cam.ac.uk/~taf1000/phd_theses.html).
- [ANTS X: explicit 5-descent paper and examples](https://mathweb.ucsd.edu/~kedlaya/ants10/paper-fisher.html).

## Editorial conventions

Prefer an author-hosted manuscript, arXiv, DOI or official software manual.
Distinguish a *full n-Selmer computation* from isogeny descent, descent on a
specified covering, and algorithms constructing covering models. New items
should include year, authors and an explanation of algorithmic relevance.
Bibliographic inclusion is not an implementation-status claim.
