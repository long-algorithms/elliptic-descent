"""Release audit: a same-length list with duplicate/missing classes is not complete."""
from dataclasses import replace

import pytest

from elliptic_descent.accounting import valid_accounting
from elliptic_descent.arithmetic import AffineF2


@pytest.mark.parametrize('labels', [(), (0,), (2, 0), (0, 2, 4, 6)])
def test_correct_partial_and_complete_quotient_ledgers(labels):
    assert valid_accounting(AffineF2(3, 0, (2, 4)), labels, 4, (1,), 3)


@pytest.mark.parametrize('labels', [(0, 0, 4, 6), (0, 2, 4, 5), (8,), (-1,), (True,), ('0',)])
def test_duplicates_and_invalid_labels_fail_closed(labels):
    assert not valid_accounting(AffineF2(3, 0, (2, 4)), labels, 4, (1,), 3)


@pytest.mark.parametrize('space,total,scalars,dimension', [
    (AffineF2(3, 0, (2, 2)), 4, (1,), 3),
    (AffineF2(3, 0, (2, 4)), 8, (1,), 3),
    (AffineF2(3, 0, (2, 4)), 4, (1, 1), 3),
    (AffineF2(3, 1, (2, 4)), 4, (1,), 3),
    (AffineF2(3, 0, (1, 4)), 4, (1,), 3),
    (AffineF2(3, 0, (2, 4)), 4, (1,), 2),
    (AffineF2(3, None, (2,)), 0, (1,), 0),
    (AffineF2(3, None, ()), 1, (1,), 0),
    (AffineF2(-1, 0, ()), 1, (), 0),
])
def test_corrupted_affine_space_metadata(space, total, scalars, dimension):
    assert not valid_accounting(space, (), total, scalars, dimension)


def test_empty_and_nonzero_affine_fibres():
    empty = AffineF2(4, None, ())
    assert valid_accounting(empty, (), 0, (1, 2), 0)
    assert not valid_accounting(empty, (0,), 0, (1, 2), 0)
    assert valid_accounting(AffineF2(4, 8, (2, 4)), (8, 10, 12, 14), 4, (1,), 3)
    assert not valid_accounting(AffineF2(4, 8, (2, 4)), (0,), 4, (1,), 3)
