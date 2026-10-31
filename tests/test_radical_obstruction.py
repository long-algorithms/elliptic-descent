"""Independent small-prime enumeration and tamper checks for mod-p^2 certificates."""
from dataclasses import replace
from itertools import product

import pytest
import sympy as sp

from elliptic_descent.geometry import U
from elliptic_descent.cubic_local import CubicCover, cubic_local_solubility
from elliptic_descent.local import QuadraticSystem
from elliptic_descent.local_obstructions import radical_plane_obstruction
from elliptic_descent.local_complete import padic_solubility_complete
from elliptic_descent.status import LocalStatus


def direct_primitive_residues(system, p, depth):
    """All projective chart representatives, without radicals or polynomial gcds."""
    m = p**depth
    for pivot in range(4):
        free = [i for i in range(4) if i != pivot]
        ranges = [range(0, m, p) if i < pivot else range(m) for i in free]
        for values in product(*ranges):
            point = [0]*4
            point[pivot] = 1
            for i, v in zip(free, values):
                point[i] = v
            if all(v % m == 0 for v in system.values(point)):
                yield tuple(point)


@pytest.mark.parametrize('p', [3, 5, 7])
def test_obstructions_against_exhaustive_projective_residues(p):
    x,y,z,w = U
    q = (x*x+p*z*z+p*p*y*w, y*y+p*w*w+p*p*x*z)
    system = QuadraticSystem.from_quadrics(q)
    evidence = radical_plane_obstruction(system, p)
    assert evidence is not None and evidence.status == LocalStatus.INSOLUBLE
    assert evidence.verify(system)
    assert next(direct_primitive_residues(system,p,1),None) is not None
    assert next(direct_primitive_residues(system,p,2),None) is None
    # Conjugate away the coordinate-aligned radical; no special basis is assumed.
    N = sp.Matrix([[1,2,1,0],[0,1,1,1],[0,0,1,2],[0,0,0,1]])
    sub = dict(zip(U,N*sp.Matrix(U)))
    changed = tuple(sp.expand(v.subs(sub,simultaneous=True)) for v in q)
    other = radical_plane_obstruction(changed,p)
    assert other is not None and other.verify(changed)
    assert not other.verify(q)
    bad = replace(other, certificate=replace(other.certificate,basis=((0,0,0,0),)*4))
    assert not bad.verify(changed)
    assert not replace(other,status=LocalStatus.SOLUBLE).verify(changed)


def test_sufficient_test_does_not_misclassify_common_binary_zeros():
    x,y,z,w = U
    for p in (3,5,7):
        # Restricted forms have a common point at infinity / finite point.
        for q in ((x*x+p*z*z,y*y+p*z*z),
                  (x*x+p*w*w,y*y+p*w*w),
                  (x*x+p*z*z,x*y+p*w*w)):
            assert radical_plane_obstruction(q,p) is None
    assert radical_plane_obstruction((x*x+2*z*z,y*y+2*w*w),2) is None


def test_571a1_obstructions_have_small_independently_checkable_certificates():
    E = (0,-1,1,-929,-10595)
    for delta in ((-1953439,-58168,3184),(-8764919,-260912,14284),
                  (25541581,760580,-41632),(135527537,4034412,-220868)):
        D = CubicCover(E,delta)
        local = cubic_local_solubility(D)
        assert local.status == LocalStatus.INSOLUBLE and local.verify(D)
        report = local.finite[-1]
        assert report.place == 571 and report.nodes == 1 and report.depth == 2
        inner = report.certificate.inner
        assert inner.certificate.kind == 'radical_plane_mod_p2'
        assert inner.verify(report.certificate.change.target,max_nodes=1)
        assert report.verify(D,max_nodes=1)
        # It must not evade an explicitly zero search budget.
        assert padic_solubility_complete(D,571,max_nodes=0).status == LocalStatus.UNKNOWN
