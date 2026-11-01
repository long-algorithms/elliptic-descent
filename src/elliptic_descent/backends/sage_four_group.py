"""Discover an exact common splitting field and flex-origin certificates.

QQbar is used for exact algebraic roots, never approximate root matching.
All output is exported as rational power-basis data and independently checked
by galois_covering.py. This deliberately expensive reference backend does
not claim to be Fisher's efficient obstruction-algebra implementation.
"""
from itertools import product

import sympy as sp

from .sage import _sage
from ..arithmetic import nonnegative_integer
from ..galois_covering import (ExactGaloisField, FourTorsionCertificate,
                               CoverTrivialization, GaloisGroupEvidence)


def _sage_rational(sage, value):
    return sage.QQ(str(value))


def discover_group_evidence(fibres, *, max_field_degree=None):
    """Return certificates tied to every actual global class in ``fibres``.

    No class is dropped if the optional degree cap is exceeded: the method
    raises, and no completed group object is returned. The cap bounds the
    accepted field degree, not the time spent discovering that field.
    """
    if not fibres.complete:
        raise ValueError('group arithmetic requires a complete fibre assembly')
    cap = None if max_field_degree is None else nonnegative_integer(max_field_degree,'maximum field degree')
    sage = _sage()
    from sage.rings.qqbar import number_field_elements_from_algebraics
    QQbar = sage.QQbar
    E0 = fibres.selmer.curves.original
    ainvs = tuple(_sage_rational(sage,c) for c in E0.ainvs)
    E = sage.EllipticCurve(QQbar,ainvs)
    b2,b4,b6,_ = E0.b_invariants
    R = sage.PolynomialRing(sage.QQ,'t')
    h = R([_sage_rational(sage,c) for c in (b6/4,b4/2,b2/4,1)])
    roots = tuple(h.roots(ring=QQbar,multiplicities=False))
    if len(roots)!=3:
        raise ArithmeticError('the completed-square cubic did not have three distinct roots')
    torsion = [E(0)]
    for e in roots:
        torsion.append(E(e,-(ainvs[0]*e+ainvs[2])/2))
        others = tuple(i for i in range(3) if roots[i]!=e)
        for signs in product((-1,1),repeat=2):
            w = [QQbar(0)]*3
            for i,sign in zip(others,signs):
                w[i] = sign*(e-roots[i]).sqrt()
            pairs = w[0]*w[1]+w[0]*w[2]+w[1]*w[2]
            X = e+pairs
            V = sum(w)*pairs
            P = E(X,V-(ainvs[0]*X+ainvs[2])/2)
            if 2*P != E(e,-(ainvs[0]*e+ainvs[2])/2):
                raise ArithmeticError('algebraic torsion halving failed')
            if P not in torsion:
                torsion.append(P)
    if len(torsion)!=16 or any(4*P!=E(0) for P in torsion):
        raise ArithmeticError('the full E[4] module was not constructed')
    basis = next(((A,B) for A in torsion for B in torsion
                  if 2*A!=E(0) and 2*B!=E(0)
                  and len({a*A+b*B for a in range(4) for b in range(4)})==16),None)
    if basis is None:
        raise ArithmeticError('failed to choose a basis of E[4]')

    numbers = list(roots)+[c for P in basis for c in P.xy()]
    descriptors = []
    for C in fibres.known_classes:
        D = C.covering.model
        marking = C.covering.marking
        if marking.translation is not None or marking.target != E0:
            raise ValueError('the reference cocycle backend requires the assembled canonical marking')
        f = R([_sage_rational(sage,c) for c in D.f])
        xi = R([_sage_rational(sage,c) for c in D.xi])
        alpha = f.roots(ring=QQbar,multiplicities=False)[0]
        A = f.derivative()(alpha)
        s = (f.derivative(2)(alpha)/2-f[2])/3
        u,shift,_,_ = tuple(_sage_rational(sage,c) for c in marking.change.urst)
        j = tuple(u*u*e+shift for e in roots)
        beta = tuple(alpha+A/(v-s) for v in j)
        delta = tuple((s-v)*xi(b)/xi(alpha) for v,b in zip(j,beta))
        square_roots = [delta[0].sqrt(),delta[1].sqrt()]
        signed_product = -f[4]*A*_sage_rational(sage,D.norm_witness)/xi(alpha)**2
        square_roots.append(signed_product/(square_roots[0]*square_roots[1]))
        if any(v*v!=d for v,d in zip(square_roots,delta)):
            raise ArithmeticError('the signed algebraic square roots are inconsistent')
        start = len(numbers)
        numbers.extend((alpha,*square_roots))
        descriptors.append((C.key,D,marking.change.urst,start))

    K,values,embedding = number_field_elements_from_algebraics(
        numbers, minimal=True, same_field=True, name='a')
    if any(embedding(v)!=a for v,a in zip(values,numbers)):
        raise ArithmeticError('common number-field conversion changed an algebraic value')
    if cap is not None and K.degree()>cap:
        raise RuntimeError(f'common field degree {K.degree()} exceeds the requested cap {cap}')
    automorphisms = tuple(K.automorphisms())
    if len(automorphisms)==K.degree():
        L = K
    else:
        L,inclusion = K.galois_closure('g',map=True)
        values = [inclusion(v) for v in values]
        automorphisms = tuple(L.automorphisms())
    n = int(L.degree())
    if cap is not None and n>cap:
        raise RuntimeError(f'normal closure degree {n} exceeds the requested cap {cap}')
    if len(automorphisms)!=n:
        raise ArithmeticError('normal closure did not supply every automorphism')
    def pack(v):
        coefficients = list(L(v).list())
        return tuple(sp.Rational(str(c)) for c in coefficients+[0]*(n-len(coefficients)))
    primitive = L.gen()
    automorphisms = tuple(sorted(automorphisms,key=lambda a:(a(primitive)!=primitive,pack(a(primitive)))))
    lookup = {a(primitive):i for i,a in enumerate(automorphisms)}
    generators,seen = [],{0}
    while len(seen)<n:
        generators.append(next(i for i in range(n) if i not in seen))
        seen,todo = {0},[0]
        while todo:
            i = todo.pop()
            for g in generators:
                j = lookup[automorphisms[i](automorphisms[g](primitive))]
                if j not in seen:
                    seen.add(j);todo.append(j)
    polynomial = tuple(sp.Rational(str(c)) for c in L.defining_polynomial().monic().list())
    field = ExactGaloisField(polynomial,tuple(pack(a(primitive)) for a in automorphisms),tuple(generators))
    T = FourTorsionCertificate(field,E0.ainvs,tuple(pack(v) for v in values[:3]),
            ((pack(values[3]),pack(values[4])),(pack(values[5]),pack(values[6]))))
    covers = tuple(CoverTrivialization(key,D.f,D.xi,D.norm_witness,urst,
                    pack(values[start]),tuple(pack(v) for v in values[start+1:start+4]))
                    for key,D,urst,start in descriptors)
    evidence = GaloisGroupEvidence(T,covers)
    if not evidence.verify():
        raise ArithmeticError('the independent exact Galois/covering certificate audit failed')
    return evidence
