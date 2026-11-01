"""Explicit resource-limit exceptions, distinct from arithmetic failures."""
from .arithmetic import nonnegative_integer


class FieldDegreeLimit(RuntimeError):
    """An exact discovered field exceeds a requested acceptance bound.

    This is not an insolubility certificate. The bound does not limit the
    time/memory spent discovering the field. RuntimeError compatibility with
    the original low-level group interface is retained.
    """

    def __init__(self, degree, limit):
        self.degree = nonnegative_integer(degree, 'field degree')
        self.limit = nonnegative_integer(limit, 'field degree limit')
        super().__init__(f'field degree {self.degree} exceeds the requested cap {self.limit}')
