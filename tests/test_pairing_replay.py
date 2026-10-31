"""Replay rejects modified formula, local evidence and matrix/label metadata."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest
import sympy as sp

from elliptic_descent.cassels_tate_forms import FisherTriple, change_parameters, substitute_generator
from elliptic_descent.cassels_tate_local import evaluate_fisher
from elliptic_descent.weierstrass import WeierstrassModel

spec=importlib.util.spec_from_file_location('replay_pairing',Path(__file__).resolve().parents[1]/'benchmarks/replay_pairing.py')
replay_module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay_module)


def record():
    # Normalize Fisher's published example onto the working 571a1 invariants.
    paper=((-64,-164,-52,68,-11),(-3,-52,-232,-60,-4),(-53,102,32,-78,-31))
    quartics=tuple(change_parameters(f,((sp.Rational(1,2),0),(0,1))) for f in paper)
    m=tuple(c/512 for c in substitute_generator(tuple(sp.Rational(c,9) for c in (936032,-8656,20)),4,0))
    F=FisherTriple(quartics,m)
    evaluation=evaluate_fisher(F)
    assert evaluation.verify() and evaluation.bit==1
    E=WeierstrassModel((0,-1,1,-929,-10595))
    entry={'left_coordinates':1,'right_coordinates':2,'left_label':1,'right_label':2,
           'bit':1,'reason':'Fisher_Theorem_3_1',
           'normalized_quartics':[[str(c) for c in q.f] for q in F.quartics],
           'product_square_root':[str(c) for c in m],'gamma':[str(c) for c in F.gamma],
           'support':list(evaluation.support),
           'local':[{'place':r.place,'point':None if r.point is None else [str(c) for c in r.point],
                     'kind':r.kind,'trials':r.trials,'bit':r.bit} for r in evaluation.local],
           'normalizations':[{'coordinates':i,'ambient_label':i,
                 'original_quartic':[str(c) for c in q.f],
                 'parameter_matrix':[['1','0'],['0','1']],
                 'theta_invariant':[str(c) for c in substitute_generator(q.z,-3,-E.b_invariants[0]/4)]}
                 for i,q in enumerate(F.quartics,1)]}
    # The coordinate metadata is a fixture labelling, not a Sage class-group proof.
    return {'fixture':'published_571a1','selmer_dimension':2,'basis_labels':[1,2],
            'complete':True,'certified':True,'verified':True,'proof':'certified','assumptions':[],
            'max_local_trials':None,'four_selmer_complete':False,
            'working_ainvs':[str(c) for c in E.ainvs], 'entries':[entry],
            'matrix':[[0,1],[1,0]],'radical_basis':[],'radical_ambient_basis':[],
            'pairing_rank':2,'pari_pairing_rank':2}


def test_formula_local_and_matrix_replay():
    result=replay_module.replay(record())
    assert result['entries_verified']==1 and result['pairing_rank']==2


@pytest.mark.parametrize('mutation', ['bit','support','gamma','square_root','point','matrix',
    'radical','label','coordinates','normalization','duplicate','missing','assumption','scope'])
def test_corrupt_replay_record_fails(mutation):
    r=deepcopy(record()); e=r['entries'][0]
    if mutation=='bit':e['bit']=0
    elif mutation=='support':e['support']=e['support'][:-1]
    elif mutation=='gamma':e['gamma'][0]='1'
    elif mutation=='square_root':e['product_square_root'][0]='1'
    elif mutation=='point':e['local'][0]['point']=['0','0']
    elif mutation=='matrix':r['matrix']=[[0,0],[0,0]]
    elif mutation=='radical':r['radical_basis']=[1]
    elif mutation=='label':e['normalizations'][0]['ambient_label']=3
    elif mutation=='coordinates':e['normalizations'][0]['coordinates']=2
    elif mutation=='normalization':e['normalizations'][0]['parameter_matrix']=[['2','0'],['0','2']]
    elif mutation=='duplicate':r['entries'].append(deepcopy(e))
    elif mutation=='missing':r['entries']=[]
    elif mutation=='assumption':r['assumptions']=['GRH']
    else:r['four_selmer_complete']=True
    with pytest.raises(ValueError):replay_module.replay(r)
