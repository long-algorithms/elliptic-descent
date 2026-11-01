"""Real covering-record replay and API tests, explicitly conditional on old arithmetic.

The lightweight parent adapter below reconstructs actual maps from PR5 records.
It is NOT a Sage substitute: it deliberately remains uncertified and records
that class-group and Selmer exhaustiveness have not been rerun by this replay.
"""
from copy import deepcopy
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
# Benchmark modules deliberately also work as standalone scripts.
sys.path.insert(0,str(ROOT/'benchmarks'))
from replay_four_group import replay_group_record, sequence_data
from replay_fibres import replay_record as replay_parent, rationals

from elliptic_descent import FourSelmerGroup, four_selmer_group
from elliptic_descent.galois_covering import evidence_from_record
from elliptic_descent.geometry import construct_four_cover
from elliptic_descent.marked_cover import MarkedQuarticCover
from elliptic_descent.two_cover import QuarticTwoCover
from elliptic_descent.weierstrass import WeierstrassModel, WeierstrassChange


@lru_cache(maxsize=2)
def load(name):
    return json.loads((Path(__file__).parent/'data'/name).read_text())


class RecordedFibres:
    """Geometric parent replay with an explicit, non-removable arithmetic caveat."""
    def __init__(self,parent):
        self.record=parent
        keys,zero,projections,inclusion=sequence_data(parent)
        self.selmer=SimpleNamespace(cardinality=1<<parent['selmer_dimension'],
            curves=SimpleNamespace(original=WeierstrassModel(rationals(parent['original_ainvs']))))
        classes=[]
        for key in keys:
            f=next(f for f in parent['fibres'] if f['coordinates']==key[0])
            source=QuarticTwoCover(rationals(f['quartic']))
            marking=MarkedQuarticCover(source,WeierstrassChange(source.ainvs,rationals(f['jacobian_to_original'])))
            m=f['models'][f['orbits'][key[1]][0]]
            D=construct_four_cover(source.f,rationals(m['xi']),norm_witness=rationals((m['norm_witness'],))[0])
            classes.append(SimpleNamespace(key=key,covering=marking.compose(D)))
        self.known_classes=tuple(classes)
        lookup={C.key:C for C in classes}
        self.zero_class=lookup[zero]
        self.kernel_map=SimpleNamespace(apply=lambda c:lookup[inclusion[c]])
        self.cardinality=len(keys)
        self.complete=True  # Exhaustive list asserted by the retained PR5 record.
        self.certified=False  # This replay does NOT rerun its arithmetic backend.
        self.assumptions=('parent Selmer/class-group exhaustiveness imported from PR5 record',)
        self._report=None

    def projection(self,key):
        return key[0]

    def verify_evidence(self,**kwargs):
        if self._report is None:
            self._report=replay_parent(self.record)
        return self._report['verified']


@pytest.mark.parametrize('name,expected', [('group_split_full.json',(2,2)),('group_cyclic_four.json',(4,))])
def test_geometric_replay_and_public_api_on_actual_global_classes(name,expected):
    r=load(name)
    report=replay_group_record(r)
    assert report['verified'] and tuple(report['invariant_factors'])==expected
    parent=RecordedFibres(r['fibres'])
    E=None if r['galois_evidence'] is None else evidence_from_record(r['galois_evidence'])
    G=four_selmer_group(parent,evidence=E)
    assert G.complete and G.verify_evidence() and not G.certified and G.assumptions
    assert G.invariant_factors==expected
    for a in G.elements:
        assert G.from_coordinates(a.coordinates)==a
        assert a-a==G.zero and 4*a==G.zero
        assert a.covering.target==parent.selmer.curves.original
        assert 2*a==G.include_two_selmer(G.projection(a))
        for b in G.elements:
            assert (a+b)-b==a
            assert G.projection(a+b)==G.projection(a)^G.projection(b)
    for g,n in zip(G.generators,G.invariant_factors):
        assert g.order()==n and n*g==G.zero
    assert len(set(G.elements))==4
    with pytest.raises(ValueError):G.element((99,0))
    with pytest.raises(ValueError):G.from_coordinates((1,2,3))
    with pytest.raises(TypeError):G.multiply(0.5,G.zero)
    H=four_selmer_group(parent,evidence=E)
    with pytest.raises(ValueError):G.zero+H.zero


@pytest.mark.parametrize('mutation',['table','generator','model','field','missing'])
def test_portable_records_do_not_silently_relabel_coverings(mutation):
    r=deepcopy(load('group_split_full.json'))
    if mutation=='table':r['addition_table'][1][2]=0
    elif mutation=='generator':r['generator_keys'][0]=r['keys'][0]
    elif mutation=='model':r['galois_evidence']['coverings'][1]['xi'][0]='999'
    elif mutation=='field':r['galois_evidence']['field']['automorphism_images'].pop()
    else:r['galois_evidence']['coverings'].pop()
    with pytest.raises(ValueError):replay_group_record(r,check_parent=False)
