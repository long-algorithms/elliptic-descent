"""Check a descent's finite enumeration ledger without redoing class groups.

Membership, independence and duplicate checks are exact F2 operations. They
certify internal accounting only: the backend must still establish that the
recorded affine space is the correct exhaustive arithmetic search space.
"""
from .arithmetic import AffineF2
from .quotient import normal_form, row_basis


def valid_accounting(space, labels, total_count, scalar_basis, norm_dimension):
    """Validate a possibly partial, nonredundant traversal of an affine quotient.

    An empty fibre has particular=None and no basis. For a nonempty fibre,
    the independent scalar image and quotient dimensions add to norm_dimension.
    All vectors are native integer bitsets, as produced by the solver.
    """
    if not isinstance(space, AffineF2):
        return False
    if type(space.nvariables) is not int or space.nvariables < 0:
        return False
    if type(total_count) is not int or total_count < 0:
        return False
    if type(norm_dimension) is not int or not 0 <= norm_dimension <= space.nvariables:
        return False
    bound = 1 << space.nvariables

    def bitset(value):
        return type(value) is int and 0 <= value < bound

    labels, scalars, basis = tuple(labels), tuple(scalar_basis), tuple(space.basis)
    if not all(bitset(v) for v in (*labels, *scalars, *basis)):
        return False
    if len(set(labels)) != len(labels) or len(labels) > total_count:
        return False
    scalar_rows, rows = row_basis(scalars), row_basis(basis)
    if len(scalar_rows) != len(scalars) or len(rows) != len(basis):
        return False
    if space.particular is None:
        return not labels and not basis and total_count == 0 and norm_dimension == 0
    if not bitset(space.particular) or total_count != 1 << len(basis):
        return False
    if norm_dimension != len(basis)+len(scalar_rows):
        return False
    # Representatives must be in the canonical complement of the scalar image.
    if any(normal_form(v, scalar_rows) != v for v in (space.particular, *basis)):
        return False
    return all(normal_form(v ^ space.particular, rows) == 0 for v in labels)
