import pytest

from elliptic_descent.geometry import U
from elliptic_descent.local import QuadraticSystem
from elliptic_descent.local_complete import _Budget, _initial_solutions


def test_degenerate_residue_fibres_respect_the_work_budget():
    from elliptic_descent.local_complete import _Limit
    u,v,w,z=U
    system=QuadraticSystem.from_quadrics((u*u,v*v))
    # A whole plane modulo p must not be materialized before checking the limit.
    with pytest.raises(_Limit):
        list(_initial_solutions(system,1_000_000_007,_Budget(3)))
