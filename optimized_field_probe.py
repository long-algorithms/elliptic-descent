"""Diagnostic reduced primitive element; no proof/check is bypassed."""
import argparse
import json
import sys
from pathlib import Path
from time import perf_counter


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('fixture')
    args=parser.parse_args()
    sys.path[:0]=[str(args.source/'src'),str(args.source/'benchmarks')]
    from elliptic_descent import two_selmer,cassels_tate_pairing,assemble_four_selmer_fibres
    from elliptic_descent.four_selmer_group import four_selmer_group
    from elliptic_descent.galois_covering import evidence_record
    import elliptic_descent.backends.sage_four_group as backend
    from four_selmer_group import EXPECTED
    from four_selmer_fibres import FIXTURES
    start=perf_counter()
    def report(stage,**fields):
        record={'stage':stage,'seconds':perf_counter()-start,**fields}
        print(json.dumps(record),flush=True)
        with open('reduced-progress.jsonl','a') as f:f.write(json.dumps(record)+'\n')
    def reduce_field(K,values):
        oldbits=max(int(abs(c)).bit_length() for c in K.defining_polynomial().list())
        report('reduction-start',degree=int(K.degree()),polynomial_bits=oldbits)
        new,from_new,to_new=K.optimized_representation(name='small')
        if new.degree()!=K.degree() or from_new(to_new(K.gen()))!=K.gen() or to_new(from_new(new.gen()))!=new.gen():
            raise ArithmeticError('primitive-element isomorphisms are not mutual inverses')
        moved=[to_new(v) for v in values]
        if any(from_new(w)!=v for w,v in zip(moved,values)):
            raise ArithmeticError('primitive-element change did not preserve all retained data')
        bits=max(int(abs(c)).bit_length() for c in new.defining_polynomial().list())
        report('reduction-complete',degree=int(new.degree()),polynomial_bits=bits)
        return new,moved
    def save_evidence(E):
        Path('reduced-unverified-evidence.json').write_text(json.dumps(evidence_record(E)))
        report('portable-check-start',degree=E.torsion.field.degree)
    source=Path(backend.__file__).read_text()
    anchor='    automorphisms = tuple(K.automorphisms())\n'
    assert source.count(anchor)==1
    source=source.replace(anchor,'    K, values = _reduce_field(K, values)\n'+anchor)
    anchor='    if not evidence.verify():\n'
    assert source.count(anchor)==1
    source=source.replace(anchor,'    _save_evidence(evidence)\n'+anchor)
    namespace=dict(backend.__dict__,_reduce_field=reduce_field,_save_evidence=save_evidence)
    exec(compile(source,'<reduced-primitive-element-diagnostic>','exec'),namespace)
    S=two_selmer(FIXTURES[args.fixture],proof=True)
    F=assemble_four_selmer_fibres(S,cassels_tate_pairing(S),proof=True)
    assert F.certified and F.verify_evidence()
    report('parent-verified',classes=F.cardinality)
    evidence=namespace['discover_group_evidence'](F)
    report('portable-check-passed')
    G=four_selmer_group(F,backend='galois',evidence=evidence)
    assert G.certified and G.verify_evidence() and G.invariant_factors==EXPECTED[args.fixture]
    for a in G.elements:
        assert G.from_coordinates(a.coordinates)==a and 2*a==G.include_two_selmer(G.projection(a))
        for b in G.elements:
            assert (a+b)-b==a
    report('group-verified',invariant_factors=G.invariant_factors,classes=G.cardinality)

if __name__=='__main__':main()
