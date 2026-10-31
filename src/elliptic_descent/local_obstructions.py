"""Small exact mod-p^2 obstructions for two quadrics, without residue enumeration.

At an odd prime, both reductions may have a common two-dimensional radical.
If their quotient binary forms have no common geometric projective zero, all
residue points lie in that radical. Divide the restricted forms by p: another
empty binary intersection proves that none of those residues lifts to mod p^2.
See docs/two-selmer.md. This optional certificate never replaces the fallback.
"""
from dataclasses import dataclass

import sympy as sp

from .arithmetic import prime
from .local import LocalResult, QuadraticSystem
from .status import LocalStatus


def _no_binary_zero(matrices, p):
    """A sufficient, exact test over Fp-bar, including the point at infinity."""
    t = sp.Symbol('t')
    forms = []
    for H in matrices:
        # Integral Hessians have even diagonal; p is odd when a division is used.
        if any(c.is_Integer is not True for c in H) or H[0, 0] % 2 or H[1, 1] % 2:
            return False
        forms.append(sp.Poly(H[0, 0]//2*t*t+H[0, 1]*t+H[1, 1]//2, t, modulus=p))
    if all(f.nth(2) == 0 for f in forms):
        return False  # Both forms vanish at (1:0).
    return sp.gcd(*forms).degree() == 0


def _radical_basis(matrices, p):
    """Column basis [complement | common radical], or None unless its rank is 2."""
    rows = [[int(v) % p for v in row] for H in matrices for row in H.tolist()]
    pivots = []
    for j in range(4):
        k = next((k for k in range(len(pivots), len(rows)) if rows[k][j]), None)
        if k is None:
            continue
        i = len(pivots)
        rows[i], rows[k] = rows[k], rows[i]
        inverse = pow(rows[i][j], -1, p)
        rows[i] = [v*inverse % p for v in rows[i]]
        for k in range(len(rows)):
            if k != i:
                c = rows[k][j]
                rows[k] = [(v-c*w) % p for v, w in zip(rows[k], rows[i])]
        pivots.append(j)
    if len(pivots) != 2:
        return None
    columns = [sp.eye(4)[:, j] for j in pivots]
    for j in range(4):
        if j in pivots:
            continue
        v = [0]*4
        v[j] = 1
        for i, k in enumerate(pivots):
            c = -rows[i][j] % p
            v[k] = c-p if 2*c > p else c
        columns.append(sp.Matrix(v))
    return sp.Matrix.hstack(*columns)


@dataclass(frozen=True)
class RadicalPlaneCertificate:
    system_key: str
    place: int
    basis: tuple
    kind: str = 'radical_plane_mod_p2'

    def verify(self, quadrics):
        system = QuadraticSystem.from_quadrics(quadrics)
        if system.key != self.system_key:
            return False
        p = prime(self.place)
        if p == 2:
            return False
        try:
            N = sp.Matrix(self.basis)
        except (TypeError, ValueError):
            return False
        if (N.shape != (4, 4) or any(v.is_Integer is not True for v in N)
                or N.det() % p == 0):
            return False
        transformed = [N.T*H*N for H in system.matrices()]
        # Check the actual certificate, not the detector's nullspace calculation.
        if any(H[i, j] % p for H in transformed for i in range(4) for j in (2, 3)):
            return False
        quotient = [H[:2, :2] for H in transformed]
        restricted = [H[2:, 2:]/p for H in transformed]
        return _no_binary_zero(quotient, p) and _no_binary_zero(restricted, p)


class RadicalPlaneResult(LocalResult):
    def verify(self, quadrics, *, max_nodes=1_000_000):
        c = self.certificate
        return bool(isinstance(c, RadicalPlaneCertificate)
                    and self.status == LocalStatus.INSOLUBLE
                    and self.place == c.place and c.verify(quadrics))


def radical_plane_obstruction(quadrics, p):
    """Return a proved obstruction, or None when this sufficient test does not apply."""
    system, p = QuadraticSystem.from_quadrics(quadrics), prime(p)
    if p == 2:
        return None
    N = _radical_basis(system.matrices(), p)
    if N is None:
        return None
    certificate = RadicalPlaneCertificate(system.key, p, tuple(map(tuple, N.tolist())))
    if not certificate.verify(system):
        return None
    return RadicalPlaneResult(p, LocalStatus.INSOLUBLE, certificate, nodes=1, depth=2)
