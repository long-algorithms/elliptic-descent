"""Exact F2 quotient by rational square classes; no geometric assertions here."""
from .arithmetic import AffineF2, nonnegative_integer


def row_basis(vectors):
    rows = {}
    for value in vectors:
        value = nonnegative_integer(value, "vector")
        for pivot in sorted(rows, reverse=True):
            if value >> pivot & 1:
                value ^= rows[pivot]
        if value:
            pivot = value.bit_length()-1
            for k in rows:
                if rows[k] >> pivot & 1:
                    rows[k] ^= value
            rows[pivot] = value
    return tuple(rows[p] for p in sorted(rows, reverse=True))


def normal_form(value, rows):
    value = nonnegative_integer(value, "vector")
    for row in rows:
        if value >> (row.bit_length()-1) & 1:
            value ^= row
    return value


def quotient_fibre(fibre, scalar_vectors):
    """A canonical transversal of an affine fibre modulo a kernel subspace."""
    scalars = row_basis(scalar_vectors)
    if fibre.particular is None:
        return AffineF2(fibre.nvariables, None, ()), scalars
    kernel = row_basis(fibre.basis)
    if any(normal_form(v, kernel) for v in scalars):
        raise ValueError("rational-scalar space is not contained in the norm kernel")
    particular = normal_form(fibre.particular, scalars)
    basis = row_basis(normal_form(v, scalars) for v in fibre.basis)
    return AffineF2(fibre.nvariables, particular, basis), scalars
