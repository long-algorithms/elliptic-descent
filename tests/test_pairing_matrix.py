"""Regression for repeated matrix construction, with no simulated Selmer arithmetic."""
from itertools import product


def test_matrix_operations_materialize_once_and_preserve_exact_linear_algebra():
    from elliptic_descent.cassels_tate import _MatrixOperations
    from elliptic_descent.two_selmer import F2Labels

    # Isolate the linear-algebra interface, not a simulated Selmer computation.
    class CoordinateSpace:
        labels = F2Labels(4, (1, 2, 4, 8))
        dimension = 4
        ambient_label = labels.ambient

    class MatrixProbe(_MatrixOperations):
        selmer = CoordinateSpace()
        calls = 0

        @property
        def matrix(self):
            self.calls += 1
            return ((0,1,0,0), (1,0,0,0), (0,0,0,0), (0,0,0,0))

    probe = MatrixProbe()
    for a,b in product(range(16), repeat=2):
        before = probe.calls
        assert probe.evaluate(a,b) == ((a & 1)*((b >> 1) & 1) + ((a >> 1) & 1)*(b & 1)) % 2
        assert probe.calls == before+1
    before = probe.calls
    assert probe.radical_basis == (4,8) and probe.calls == before+1
    for a in range(16):
        before = probe.calls
        assert probe.in_radical(a) == (a % 4 == 0)
        assert probe.calls == before+1
