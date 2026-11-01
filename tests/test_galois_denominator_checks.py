"""The verifier checks equivalent identities without dense-field inversions.

Its denominators A=f'(alpha) and xi(alpha) are separately checked nonzero.
Clearing them preserves the signed identities and the rejection of tampering.
"""
from dataclasses import replace
import json
from pathlib import Path

from elliptic_descent.galois_covering import evidence_from_record


def certificate():
    record = json.loads((Path(__file__).parent/'data/group_split_full.json').read_text())
    return evidence_from_record(record['galois_evidence'])


def test_cover_checks_do_not_invert_nonrational_field_denominators(monkeypatch):
    E = certificate()
    assert E.torsion.verify()
    value_type = type(E.torsion.field.one)
    original = value_type.__truediv__
    attempted = []

    def rational_division_only(self, other):
        if hasattr(other, 'to_list') and len(other.to_list()) > 1:
            attempted.append(other)
            raise ArithmeticError('verification attempted a dense field inverse')
        return original(self, other)

    monkeypatch.setattr(value_type, '__truediv__', rational_division_only)
    for C in E.coverings:
        assert C.verify(E.torsion)
        assert not replace(C, norm_witness=-C.norm_witness).verify(E.torsion)
    assert not attempted


def test_cleared_checks_agree_with_the_original_quotient_identities():
    E = certificate()
    assert E.torsion.verify()
    F = E.torsion.field
    from elliptic_descent.galois_covering import polynomial_value
    for C in E.coverings:
        f, xi, alpha, r, _, A, _, beta, xia, _, roots = C._data(E.torsion)
        assert A and xia
        product = roots[0]*roots[1]*roots[2]
        rhs = -f[4]*A*r
        assert product == rhs/(xia*xia)
        assert product*xia*xia == rhs
        P = C.origin(E.torsion)
        for t in (alpha, *beta):
            value = polynomial_value(xi,t,F.zero)*polynomial_value(P,t,F.zero)**2
            rhs = xia*(alpha-t)
            assert value == rhs/A
            assert value*A == rhs
