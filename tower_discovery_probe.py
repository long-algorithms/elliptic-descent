"""Diagnostic exact field-tower construction; not an accepted production change."""
import argparse
import json
import sys
from pathlib import Path
from time import perf_counter


def discover(fibres, report):
    from sage.all import QQ, PolynomialRing, EllipticCurve
    import sympy as sp
    from elliptic_descent.galois_covering import (
        ExactGaloisField, FourTorsionCertificate, CoverTrivialization, GaloisGroupEvidence, evidence_record,
    )
    class Tower:
        def __init__(self):
            self.K=QQ
            self.values=[]
            self.extensions=0
        def append(self,value):
            self.values.append(self.K(value))
            return len(self.values)-1
        def root(self, coefficients, tag):
            R=PolynomialRing(self.K,'z')
            p=R(coefficients)
            if p.degree()<1:
                raise ValueError('root polynomial has no positive degree')
            fac=p.factor()
            q=min((q for q,e in fac),key=lambda q:q.degree()).monic()
            if q.degree()==1:
                value=-q[0]
            else:
                old=self.K
                self.extensions+=1
                report('adjoin', tag=tag, old_degree=int(old.degree()), relative_degree=int(q.degree()))
                L=old.extension(q,'b'+str(self.extensions))
                if old is QQ:
                    new=L
                    value=L.gen()
                    moved=[L(v) for v in self.values]
                    new_coefficients=[L(c) for c in p.list()]
                else:
                    new=L.absolute_field('a'+str(self.extensions))
                    from_new,to_new=new.structure()
                    value=to_new(L.gen())
                    moved=[to_new(L(v)) for v in self.values]
                    new_coefficients=[to_new(L(c)) for c in p.list()]
                    if any(from_new(to_new(L(v)))!=L(v) for v in (old.gen(),old(1))):
                        raise ArithmeticError('relative/absolute structure maps failed')
                self.K=new
                self.values=moved
                p=PolynomialRing(new,'z')(new_coefficients)
                report('extended',tag=tag,degree=int(new.degree()))
            if p(value)!=0:
                raise ArithmeticError('adjoined value does not satisfy its defining equation')
            return self.append(value)
    E0=fibres.selmer.curves.original
    rational=lambda v:QQ(str(v))
    ainvs=tuple(rational(c) for c in E0.ainvs)
    b2,b4,b6,_=E0.b_invariants
    hq=tuple(rational(c) for c in (b6/4,b4/2,b2/4,1))
    tower=Tower()
    root_indices=[]
    for i in range(3):
        R=PolynomialRing(tower.K,'z');z=R.gen();p=R(hq)
        for j in root_indices:
            p=p.quo_rem(z-tower.values[j])[0]
        root_indices.append(tower.root(p.list(),'cubic-root-'+str(i)))
    i_index=tower.root((1,0,1),'fourth-root-of-unity')
    for a,b in ((0,1),(0,2),(1,2)):
        d=tower.values[root_indices[a]]-tower.values[root_indices[b]]
        tower.root((-d,0,1),'torsion-difference-'+str(a)+str(b))
    K=tower.K
    roots=tuple(tower.values[i] for i in root_indices)
    E=EllipticCurve(K,ainvs)
    basis=[]
    for index in (0,1):
        e=roots[index]
        w=[K(0)]*3
        for j in range(3):
            if j!=index:
                d=e-roots[j]
                w[j]=d.sqrt()
                if w[j]*w[j]!=d:
                    raise ArithmeticError('torsion square root failed')
        pairs=w[0]*w[1]+w[0]*w[2]+w[1]*w[2]
        x=e+pairs;v=sum(w)*pairs
        point=E(x,v-(ainvs[0]*x+ainvs[2])/2)
        if 2*point!=E(e,-(ainvs[0]*e+ainvs[2])/2):
            raise ArithmeticError('torsion halving failed')
        basis.append(point)
    if len({a*basis[0]+b*basis[1] for a in range(4) for b in range(4)})!=16:
        raise ArithmeticError('torsion basis does not span E[4]')
    tower.values=list(roots)+[c for P in basis for c in P.xy()]
    descriptors=[]
    report('torsion-field',degree=int(K.degree()))
    for C in fibres.known_classes:
        D=C.covering.model; marking=C.covering.marking
        if marking.translation is not None or marking.target!=E0:
            raise ValueError('unexpected marking')
        fcoeff=tuple(rational(c) for c in D.f)
        xicoeff=tuple(rational(c) for c in D.xi)
        alpha_index=tower.root(fcoeff,'branch-'+str(C.key))
        def data():
            K=tower.K;R=PolynomialRing(K,'z');f=R(fcoeff);xi=R(xicoeff)
            alpha=tower.values[alpha_index]
            A=f.derivative()(alpha)
            s=(f.derivative(2)(alpha)/2-f[2])/3
            u,shift,_,_=tuple(K(rational(c)) for c in marking.change.urst)
            j=tuple(u*u*e+shift for e in tower.values[:3])
            beta=tuple(alpha+A/(v-s) for v in j)
            delta=tuple((s-v)*xi(b)/xi(alpha) for v,b in zip(j,beta))
            signed=-f[4]*A*K(rational(D.norm_witness))/xi(alpha)**2
            return delta,signed
        for j in (0,1):
            delta,_=data()
            tower.root((-delta[j],0,1),'cover-square-'+str(C.key)+'-'+str(j))
        delta,signed=data()
        r0,r1=tower.values[alpha_index+1:alpha_index+3]
        r2=signed/(r0*r1)
        tower.append(r2)
        if any(r*r!=d for r,d in zip((r0,r1,r2),delta)):
            raise ArithmeticError('signed cover square-root check failed')
        descriptors.append((C.key,D,marking.change.urst,alpha_index))
        report('cover-field',key=C.key,degree=int(tower.K.degree()))
    K=tower.K;values=tower.values
    report('automorphisms-start',degree=int(K.degree()))
    aut=tuple(K.automorphisms())
    if len(aut)!=K.degree():
        report('normal-closure-start',degree=int(K.degree()),automorphisms=len(aut))
        L,inclusion=K.galois_closure('normal',map=True)
        values=[inclusion(v) for v in values]
        K=L;aut=tuple(K.automorphisms())
    n=int(K.degree())
    if len(aut)!=n:
        raise ArithmeticError('incomplete Galois automorphisms')
    def pack(v):
        c=list(K(v).list())
        return tuple(sp.Rational(str(a)) for a in c+[0]*(n-len(c)))
    primitive=K.gen()
    aut=tuple(sorted(aut,key=lambda a:(a(primitive)!=primitive,pack(a(primitive)))))
    lookup={a(primitive):i for i,a in enumerate(aut)}
    generators=[];seen={0}
    while len(seen)<n:
        generators.append(next(i for i in range(n) if i not in seen))
        seen={0};todo=[0]
        while todo:
            i=todo.pop()
            for g in generators:
                j=lookup[aut[i](aut[g](primitive))]
                if j not in seen:
                    seen.add(j);todo.append(j)
    f=tuple(sp.Rational(str(c)) for c in K.defining_polynomial().monic().list())
    field=ExactGaloisField(f,tuple(pack(a(primitive)) for a in aut),tuple(generators))
    T=FourTorsionCertificate(field,E0.ainvs,tuple(pack(v) for v in values[:3]),
                             ((pack(values[3]),pack(values[4])),(pack(values[5]),pack(values[6]))))
    covers=tuple(CoverTrivialization(key,D.f,D.xi,D.norm_witness,urst,
                                    pack(values[start]),tuple(pack(v) for v in values[start+1:start+4]))
                 for key,D,urst,start in descriptors)
    evidence=GaloisGroupEvidence(T,covers)
    Path('tower-unverified-evidence.json').write_text(json.dumps(evidence_record(evidence)))
    report('portable-verification-start',degree=n)
    if not evidence.verify():
        raise ArithmeticError('portable exact evidence failed')
    report('portable-verification-passed',degree=n)
    return evidence


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('fixture')
    args=parser.parse_args()
    sys.path[:0]=[str(args.source/'src'),str(args.source/'benchmarks')]
    from elliptic_descent import two_selmer,cassels_tate_pairing,assemble_four_selmer_fibres
    from elliptic_descent.four_selmer_group import four_selmer_group
    from four_selmer_group import EXPECTED
    from four_selmer_fibres import FIXTURES
    start=perf_counter()
    def report(stage,**fields):
        record={'stage':stage,'seconds':perf_counter()-start,**fields}
        print(json.dumps(record),flush=True)
        with open('tower-progress.jsonl','a') as f:f.write(json.dumps(record)+'\n')
    S=two_selmer(FIXTURES[args.fixture],proof=True)
    F=assemble_four_selmer_fibres(S,cassels_tate_pairing(S),proof=True)
    assert F.certified and F.verify_evidence()
    report('parent-verified',classes=F.cardinality)
    evidence=discover(F,report)
    G=four_selmer_group(F,backend='galois',evidence=evidence)
    assert G.certified and G.verify_evidence() and G.invariant_factors==EXPECTED[args.fixture]
    for a in G.elements:
        assert G.from_coordinates(a.coordinates)==a and 2*a==G.include_two_selmer(G.projection(a))
        for b in G.elements:
            assert (a+b)-b==a
    report('group-verified',invariant_factors=G.invariant_factors,classes=G.cardinality)

if __name__=='__main__':main()
