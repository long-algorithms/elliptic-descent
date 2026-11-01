"""Replay group operations tied to recorded marked coverings, without Sage.

The parent fibre replay checks geometric orbits and inclusion diagrams. It does
not independently recompute Selmer arithmetic or all local conditions. General
addition is checked through the exact normal field and geometric cocycles.
"""
import argparse
import json
from pathlib import Path

import sympy as sp

from replay_fibres import replay_record as replay_fibres, require, rationals
from elliptic_descent.four_selmer_group import FiniteClassLaw, forced_class_law
from elliptic_descent.galois_covering import evidence_from_record


def sequence_data(parent):
    """Recover the ACTUAL inclusion from the recorded zero-fibre labels/orbits."""
    keys=tuple(tuple(k) for k in parent['global_keys'])
    zero=tuple(parent['zero_class_key'])
    fibre=next(f for f in parent['fibres'] if f['coordinates']==0)
    labels=parent['kernel_labels']
    if len(labels)!=sum(len(o) for o in fibre['orbits']):
        raise ValueError('kernel labels and relative diagrams disagree')
    inclusion=[]
    for c in range(1<<parent['selmer_dimension']):
        indices=[i for i,label in enumerate(labels) if label==c]
        if len(indices)!=1:
            raise ValueError('missing or duplicate kernel-inclusion label')
        orbits=[j for j,o in enumerate(fibre['orbits']) if indices[0] in o]
        if len(orbits)!=1:
            raise ValueError('kernel diagram is not in exactly one orbit')
        inclusion.append((0,orbits[0]))
    return keys,zero,tuple(k[0] for k in keys),tuple(inclusion)


def replay_group_record(record, *, check_parent=True):
    parent=record['fibres']
    parent_report=replay_fibres(parent) if check_parent else None
    keys,zero,projections,inclusion=sequence_data(parent)
    law=FiniteClassLaw(tuple(tuple(k) for k in record['keys']),
                      tuple(tuple(row) for row in record['addition_table']),record['zero_index'])
    require(law.verify() and law.keys==keys and law.keys[law.zero]==zero,'invalid group law or keys')
    method=record['method']
    field_degree=0
    if method=='geometric_galois_cocycles':
        E=evidence_from_record(record['galois_evidence'])
        require(E.verify(),'invalid field/torsion/covering cocycles')
        require(E.torsion.ainvs==rationals(parent['original_ainvs']),'wrong original elliptic curve')
        require(tuple(C.key for C in E.coverings)==keys,'cocycles are not attached to all actual classes')
        for C in E.coverings:
            f=next(f for f in parent['fibres'] if f['coordinates']==C.key[0])
            m=f['models'][f['orbits'][C.key[1]][0]]
            require(C.f==rationals(f['quartic']) and C.xi==rationals(m['xi'])
                    and C.norm_witness==sp.Rational(m['norm_witness'])
                    and C.jacobian_to_original==rationals(f['jacobian_to_original']),
                    'cocycle certificate belongs to a different marked covering')
        T=E.class_table
        require(law.zero==T.zero_index and all(law.table[i][j]==T.add_indices(i,j)
                for i in range(len(keys)) for j in range(len(keys))), 'addition does not match geometric cocycles')
        field_degree=E.torsion.field.degree
    else:
        forced=forced_class_law(keys,zero,projections,inclusion)
        require(forced is not None and forced==(method,law) and record['galois_evidence'] is None,
                'a nonunique extension was assigned an arbitrary group law')
    index={k:i for i,k in enumerate(keys)}
    require(all(law.keys[law.table[i][i]]==inclusion[p] for i,p in enumerate(projections)),
            'doubling is not inclusion after projection')
    require(all(projections[law.table[i][j]]==a^b for i,a in enumerate(projections)
                for j,b in enumerate(projections)),'projection is not a homomorphism')
    require(all(law.keys[law.table[index[a]][index[b]]]==inclusion[i^j]
                for i,a in enumerate(inclusion) for j,b in enumerate(inclusion)),
            'kernel inclusion is not a homomorphism')
    factors,gens,coords=law.decomposition
    require(tuple(record['invariant_factors'])==factors,'incorrect invariant factors')
    require(tuple(tuple(k) for k in record['generator_keys'])==tuple(keys[i] for i in gens),
            'generators are not the actual recorded covering classes')
    require(len(coords)==len(keys),'invariant-factor coordinates omit classes')
    require(record['complete'] is True and record['verified'] is True,'incomplete or unverified record')
    return {'fixture':record['fixture'],'verified':True,'method':method,
            'classes':len(keys),'sums_checked':len(keys)**2,'field_degree':field_degree,
            'invariant_factors':factors,'parent_replay':parent_report,
            'scope':'exact recorded group/cover evidence; not a new class-group computation'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records',type=Path)
    args=parser.parse_args()
    rows=[json.loads(line) for line in args.records.read_text().splitlines() if line.strip()]
    seen=set()
    for r in rows:
        key=(r['fixture'],r['repetition'])
        require(key not in seen,'duplicate group benchmark record')
        seen.add(key)
        print(json.dumps(replay_group_record(r),sort_keys=True),flush=True)


if __name__=='__main__':
    main()
