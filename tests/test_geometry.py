import pytest
import sympy as sp

from elliptic_descent import construct_four_cover


def test_known_quadrics_and_point():
    D = construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0))
    u0, u1, u2, u3 = sp.symbols("u0:4")
    assert sp.expand(D.quadrics[0] - (-2*u0*u1 - 2*u1*u3 - u2**2 + 2*u2*u3)) == 0
    assert sp.expand(D.quadrics[1] - (-2*u0*u2 - u1**2 - 2*u2*u3 + u3**2)) == 0
    eq, image = D.evaluate((1, 0, 0, 0))
    assert eq == (0, 0)
    X, Y, Z = image
    assert sp.expand(Y**2 - sum(D.f[i] * X**i * Z**(4-i) for i in range(5))) == 0
    assert D.check_algebra()
    assert D.check_map_identity()


@pytest.mark.parametrize("f,xi", [
    ((1, -1, 0, 0, 1), (1, 0, 0, 0)),
    ((1, -1, 0, 0, 4), (1, 0, 0, 0)),
    ((1, 1, 0, 0, -4), (0, 1, 0, 0)),
    ((4, 0, -5, 0, 1), (1, 0, 0, 0)),
])
def test_universal_norm_identity(f, xi):
    D = construct_four_cover(f, xi)
    assert D.check_algebra()
    assert D.check_map_identity()


def test_rational_square_change_of_xi():
    f = (1, -1, 0, 0, 1)
    D = construct_four_cover(f, (0, -1, 0, 0))
    E = construct_four_cover(f, (0, -9, 0, 0))
    assert E.norm_witness == 81 * D.norm_witness
    assert all(sp.expand(y - 9*x) == 0 for x,y in zip(D.q,E.q))
    assert sp.expand(E.map_xyz[1] - 81 * D.map_xyz[1]) == 0


def test_invalid_input():
    with pytest.raises(ValueError, match="squarefree"):
        construct_four_cover((1, 0, 2, 0, 1), (1, 0, 0, 0))
    with pytest.raises(ValueError, match="leading coefficient"):
        construct_four_cover((1, 0, 0, 1, 0), (1, 0, 0, 0))
    with pytest.raises(ValueError, match="norm witness"):
        construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0), norm_witness=2)
    with pytest.raises(TypeError, match="floating-point"):
        construct_four_cover((1.0, -1, 0, 0, 1), (1, 0, 0, 0))
    with pytest.raises(ValueError, match="invertible"):
        construct_four_cover((1, -1, 0, 0, 1), (0, 0, 0, 0))


def test_projective_scaling_and_zero_rejection():
    D = construct_four_cover((1, -1, 0, 0, 1), (0, -1, 0, 0))
    eq, image = D.evaluate((1, 0, 0, 0))
    eq2, image2 = D.evaluate((3, 0, 0, 0))
    assert eq2 == tuple(9 * v for v in eq)
    assert image2 == (9 * image[0], 81 * image[1], 9 * image[2])
    with pytest.raises(ValueError, match="projective point"):
        D.evaluate((0, 0, 0, 0))
