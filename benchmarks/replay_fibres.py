"""Replay exact global-equivalence diagrams and kernel labels without Sage.

This does not rerun Selmer or relative local enumeration. The recorded parent
classes are trusted inputs; all recorded geometric identifications are checked.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import sympy as sp

from elliptic_descent.cubic_local import CubicCover
from elliptic_descent.cubic_arithmetic import kummer_algebra, rational_two_torsion, rational_halves_of_two_torsion
from elliptic_descent.conic_parametrization import ConicParametrization
from elliptic_descent.selmer_quartic import QuarticPresentation
from elliptic_descent.geometry import construct_four_cover
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange
from elliptic_descent.equivalence import ScalarSquareEquivalence
from elliptic_descent.fibre_actions import quartic_deck_action, FibreEquivalence, zero_fibre_kummer
from elliptic_descent.four_selmer_fibres import action_orbits
from elliptic_descent.quotient import row_basis, normal_form
from elliptic_descent.two_selmer import F2Labels


def rationals(v):
    return tuple(sp.Rational(c) for c in v)


def rational_point(P):
    return None if P is None else rationals(P)


def require(value, message):
    if not value:
        raise ValueError(message)


def replay_record(record):
    r=record
    E=WeierstrassModel(rationals(r['original_ainvs']))
    W=WeierstrassModel(rationals(r['working_ainvs']))
    change=WeierstrassChange(W,rationals(r['working_to_original']))
    require(change.verify() and change.target==E,'wrong original-model marking')
    dimension=r['selmer_dimension']
    require(type(dimension) is int and dimension>=0,'invalid Selmer dimension')
    B=r['pairing_matrix']
    require(len(B)==dimension and all(len(row)==dimension for row in B),'matrix shape')
    require(all(type(B[i][j]) is int and B[i][j] in (0,1) and B[i][j]==B[j][i]
                and (i!=j or B[i][j]==0) for i in range(dimension) for j in range(dimension)),'nonalternating pairing')
    evaluate=lambda a,b:sum(((a>>i)&1)*B[i][j]*((b>>j)&1) for i in range(dimension) for j in range(dimension))%2
    radical=tuple(c for c in range(1<<dimension) if all(evaluate(c,1<<j)==0 for j in range(dimension)))
    require(tuple(r['image_classes'])==radical,'wrong pairing image')
    entries=r['torsion']; points=tuple(rational_point(e['point']) for e in entries)
    require(points==rational_two_torsion(E),'missing rational 2-torsion points')
    for e,P in zip(entries,points):
        require(tuple(rational_point(Q) for Q in e['halves'])==rational_halves_of_two_torsion(E,P),
                'wrong rational halves')
    image=row_basis(e['coordinates'] for e in entries)
    require(tuple(r['torsion_image_basis'])==image,'wrong Kummer image span')
    require(all(type(e['coordinates']) is int and 0<=e['coordinates']<1<<dimension for e in entries),'invalid torsion labels')
    require(all((e['coordinates']==0)==bool(e['halves']) for e in entries),'wrong torsion stabilizers')
    require(tuple(f['coordinates'] for f in r['fibres'])==radical,'missing or duplicate image fibre')
    require(tuple(e['coordinates'] for e in r['exclusions'])==tuple(c for c in range(1<<dimension) if c not in radical),
            'missing excluded class')
    require(all(evaluate(e['coordinates'],e['test_vector'])==1 for e in r['exclusions']),'invalid pairing exclusion')
    total,witnesses,zero_models,zero_cover=0,0,None,None
    for f in r['fibres']:
        auxiliary=CubicCover(W,rationals(f['delta']))
        P=QuarticPresentation(auxiliary,ConicParametrization(tuple(rationals(row) for row in f['conic_matrix'])))
        require(P.verify() and P.quartic.f==rationals(f['quartic']),'invalid intermediate presentation')
        marking=P.marking.retarget(change)
        require(marking.change.urst==rationals(f['jacobian_to_original']),'wrong retained target map')
        models=tuple(construct_four_cover(P.quartic.f,rationals(m['xi']),norm_witness=sp.Rational(m['norm_witness']))
                     for m in f['models'])
        require(len(models)==1<<dimension and all(D.check_construction() for D in models),'relative model list')
        require(tuple(rational_point(a['point']) for a in f['actions'])==points,'incomplete action group')
        permutations=[]
        for a,T in zip(f['actions'],points):
            g=quartic_deck_action(P,change.inverse().apply(T))
            require(g.matrix==tuple(rationals(row) for row in a['matrix']) and g.y_scale==sp.Rational(a['y_scale'])
                    and g.multiplier==rationals(a['cubic_multiplier']) and g.point==rational_point(a['working_point']),
                    'wrong intermediate torsion action')
            permutation=tuple(a['permutation']);permutations.append(permutation)
            require(len(a['witnesses'])==len(models) and set(permutation)==set(range(len(models))),'action row')
            for i,(j,w) in enumerate(zip(permutation,a['witnesses'])):
                t=g.transport(models[i])
                require(t.raw.xi==rationals(w['raw']['xi']) and t.raw.norm_witness==sp.Rational(w['raw']['norm_witness']),
                        'wrong transported norm pair')
                evidence=FibreEquivalence(t,models[j],ScalarSquareEquivalence(sp.Rational(w['scalar']),rationals(w['multiplier'])))
                require(evidence.verify() and evidence.coordinate_matrix==sp.Matrix([rationals(row) for row in w['coordinate_matrix']]),
                        'failed explicit global equivalence')
                witnesses+=1
        orbits=action_orbits(permutations,len(models))
        require(orbits==tuple(tuple(o) for o in f['orbits']),'wrong global orbit partition')
        require(all(len(o)==1<<len(image) for o in orbits),'wrong stabilizer/orbit size')
        for i,P0 in enumerate(points):
            identity=tuple(range(len(models)))
            require((permutations[i]==identity)==(entries[i]['coordinates']==0),'wrong action kernel')
            for j,Q in enumerate(points):
                k=points.index(E.add(P0,Q))
                require(tuple(permutations[i][permutations[j][v]] for v in range(len(models)))==permutations[k],
                        'action violates the group law')
        total+=len(orbits)
        if f['coordinates']==0:
            zero_models=models
            zero_cover=SimpleNamespace(ambient_label=0,quartic=P.quartic,marking=marking,
                        selmer=SimpleNamespace(curves=SimpleNamespace(working=W,to_working=change.inverse())))
            zero_actions=permutations
            zero_orbits=orbits
    require(total==r['global_cardinality']==r['predicted_cardinality']==len(radical)*(1<<(dimension-len(image))),
            'global exact-sequence cardinality mismatch')
    require(zero_models is not None,'missing zero fibre')
    A=kummer_algebra(W); arithmetic=SimpleNamespace(algebra=A)
    labels=tuple(r['kernel_labels'])
    require(len(labels)==len(zero_models) and set(labels)==set(range(1<<dimension)),'kernel labels not bijective')
    require(len(r['kernel_classes'])==len(labels),'missing kernel class witnesses')
    for label,D,e in zip(labels,zero_models,r['kernel_classes']):
        delta=zero_fibre_kummer(zero_cover,D,arithmetic)
        root=rationals(e['square_root']); representative=rationals(e['representative'])
        require(e['coordinates']==label and delta==rationals(e['delta'])
                and A.multiply(representative,A.multiply(root,root))==delta,'wrong zero-fibre Kummer identity')
    for permutation,e in zip(zero_actions,entries):
        require(all(labels[j]==labels[i]^e['coordinates'] for i,j in enumerate(permutation)),
                'wrong actual torsion action on the kernel labels')
    for i in range(len(labels)):
        for j in range(len(labels)):
            require(any(i in o and j in o for o in zero_orbits)==(normal_form(labels[i]^labels[j],image)==0),
                    'kernel quotient disagrees with actual global equivalence')
    require(tuple(r['zero_class_key'])==(0,next(i for i,o in enumerate(zero_orbits) if labels.index(0) in o)),
            'wrong global zero class')
    keys=tuple((f['coordinates'],i) for f in r['fibres'] for i in range(len(f['orbits'])))
    require(tuple(tuple(v) for v in r['global_keys'])==keys,'global keys do not enumerate the quotient')
    return {'fixture':r['fixture'],'global_classes':total,'fibres':len(r['fibres']),'equivalence_witnesses':witnesses,
            'verified':True,'scope':'exact global-equivalence and kernel geometry, not new Selmer completeness'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',type=Path)
    args=parser.parse_args()
    reports=[replay_record(json.loads(line)) for line in args.path.read_text().splitlines() if line.strip()]
    require(len({r['fixture'] for r in reports})==len(reports),'duplicate curve records')
    print(json.dumps({'passed':True,'reports':reports,'global_classes':sum(r['global_classes'] for r in reports),
                      'equivalence_witnesses':sum(r['equivalence_witnesses'] for r in reports)},sort_keys=True))


if __name__=='__main__':
    main()
