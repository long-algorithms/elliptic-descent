"""Replay exact Fisher triples/local evidence from benchmark JSONL without Sage.

This does not re-establish the parent Selmer basis or its class-group proof.
It checks every formula, support prime and local contribution, and reconstructs
matrix/radical accounting from the verified pair values.
"""
import argparse
import json
from pathlib import Path

import sympy as sp

from elliptic_descent.arithmetic import affine_f2
from elliptic_descent.cassels_tate_forms import FisherTriple, change_parameters, substitute_generator
from elliptic_descent.cassels_tate_local import FisherEvaluation, PairingLocalEvidence
from elliptic_descent.weierstrass import WeierstrassModel
from elliptic_descent.cassels_tate import working_invariants
from elliptic_descent.quotient import row_basis


def replay(record):
    n = record['selmer_dimension']
    basis = tuple(record['basis_labels'])
    if (type(n) is not int or n < 0 or len(basis) != n
            or not all(record.get(k) is True for k in ('complete','certified','verified'))
            or record.get('four_selmer_complete') is not False or record.get('assumptions') != []
            or record.get('proof') != 'certified' or record.get('max_local_trials') is not None
            or any(type(b) is not int or b < 0 for b in basis) or len(row_basis(basis)) != n):
        raise ValueError('wrong or incomplete benchmark scope')
    rows = [[0]*n for _ in range(n)]
    expected = {(1 << i, 1 << j) for i in range(n) for j in range(i+1, n)}
    seen = set()
    E = WeierstrassModel(tuple(map(sp.Rational,record['working_ainvs'])))
    for e in record['entries']:
        key = e['left_coordinates'], e['right_coordinates']
        if key not in expected or key in seen:
            raise ValueError('missing, duplicate or mislabelled entry')
        seen.add(key)
        if e.get('reason') != 'Fisher_Theorem_3_1':
            raise ValueError('an off-diagonal basis entry needs formula evidence')
        F = FisherTriple(tuple(tuple(map(sp.Rational,q)) for q in e['normalized_quartics']),
                         tuple(map(sp.Rational,e['product_square_root'])))
        local = tuple(PairingLocalEvidence(r['place'], None if r['point'] is None else
                      tuple(map(sp.Rational,r['point'])),r['bit'],r['trials'],r['kind']) for r in e['local'])
        evidence = FisherEvaluation(F,tuple(e['support']),local)
        if (F.quartics[0].invariants != working_invariants(E)
                or not evidence.verify() or evidence.bit != e['bit']
                or F.gamma != tuple(map(sp.Rational,e['gamma']))):
            raise ValueError('local/form evidence mismatch')
        labels=[]
        if len(e['normalizations']) != 3:
            raise ValueError('expected three normalization records')
        for N,Q,coordinates in zip(e['normalizations'],F.quartics,(key[0],key[1],key[0]^key[1])):
            transformed=change_parameters(tuple(map(sp.Rational,N['original_quartic'])),
                tuple(tuple(map(sp.Rational,row)) for row in N['parameter_matrix']))
            if transformed != Q.f or substitute_generator(Q.z,-3,-E.b_invariants[0]/4) != tuple(map(sp.Rational,N['theta_invariant'])):
                raise ValueError('normalization identity mismatch')
            c=N['coordinates']; label=0
            if type(c) is not int or c != coordinates or not 0 <= c < 1 << n:
                raise ValueError('bad class coordinates')
            for i,b in enumerate(basis):
                if c >> i & 1: label ^= b
            if label != N['ambient_label']: raise ValueError('basis label mismatch')
            labels.append(label)
        if len(labels)!=3 or labels[0]^labels[1]!=labels[2]:
            raise ValueError('triple labels do not sum to zero')
        if labels[:2]!=[e['left_label'],e['right_label']]: raise ValueError('entry label mismatch')
        i,j=key[0].bit_length()-1,key[1].bit_length()-1
        rows[i][j]=rows[j][i]=evidence.bit
    if seen!=expected or rows != record['matrix']:
        raise ValueError('matrix assembly mismatch')
    columns=tuple(sum(rows[i][j]<<i for i in range(n)) for j in range(n))
    radical=affine_f2(columns,0).basis
    ambient=[]
    for c in radical:
        label=0
        for i,b in enumerate(basis):
            if c >> i & 1: label ^= b
        ambient.append(label)
    if (tuple(record['radical_basis'])!=radical or record['pairing_rank']!=n-len(radical)
            or record.get('radical_ambient_basis') != ambient
            or record.get('pari_pairing_rank') != n-len(radical)):
        raise ValueError('radical accounting mismatch')
    return {'fixture':record['fixture'],'entries_verified':len(seen),'pairing_rank':n-len(radical)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',type=Path)
    args=parser.parse_args()
    reports=[replay(json.loads(line)) for line in args.path.read_text().splitlines() if line.strip()]
    print(json.dumps({'passed':True,'records':reports,
                     'scope':'formula/local evidence replay, not an independent Selmer computation'},sort_keys=True))


if __name__=='__main__': main()
