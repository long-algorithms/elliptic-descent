"""Elementary integral model improvements with exact change-of-model evidence.

The safe operations saturate the two equation rows at p and move common or
square linear factors of the reductions into a p-scaled coordinate.
Strict improvements decrease the determinant-pencil discriminant valuation.
At most four neutral square-factor steps are permitted between strict improvements.
This is not a full minimization/reduction implementation.
"""
from dataclasses import dataclass
from math import gcd

import sympy as sp

from .geometry import U, _rat
from .arithmetic import prime
from .local import LocalCertificate, LocalResult, MONOMIALS, QuadraticSystem, _point
from .status import LocalStatus


def _coefficients(H):
    if any(c.is_Integer is not True for c in H) or any(H[i,i] % 2 for i in range(4)):
        raise ArithmeticError("expected an integral quadratic-form Hessian")
    return tuple(int(H[i,j]/2) if i == j else int(H[i,j]) for i,j in MONOMIALS)


def _primitive(row):
    divisor = gcd(*row)
    if not divisor:
        raise ValueError("dependent equations do not define a smooth intersection")
    if next(c for c in row if c) < 0:
        divisor = -divisor
    return tuple(c//divisor for c in row), divisor


def _square_factor(row, p):
    """L when a nonzero row reduces to a scalar times L^2, else None."""
    if p == 2:
        if any(c % p for c,(i,j) in zip(row,MONOMIALS) if i != j):
            return None
        linear = tuple(c % p for c,(i,j) in zip(row,MONOMIALS) if i == j)
        return linear if any(linear) else None
    H = QuadraticSystem((row,row)).matrices()[0]
    pivot = next((i for i in range(4) if H[i,i] % p), None)
    if pivot is None:
        return None
    inverse = pow(int(H[pivot,pivot]) % p,-1,p)
    linear = tuple(int(H[pivot,j])*inverse % p for j in range(4))
    if all((H[i,j]-H[pivot,pivot]*linear[i]*linear[j]) % p == 0
           for i in range(4) for j in range(4)):
        return linear
    return None


def _vp(n,p):
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


@dataclass(frozen=True)
class IntegralModelChange:
    source_key: str
    target: QuadraticSystem
    coordinate_matrix: tuple
    equation_matrix: tuple
    steps: int

    def verify(self, source):
        source = QuadraticSystem.from_quadrics(source)
        if source.key != self.source_key:
            return False
        try:
            N = sp.Matrix([[ _rat(c) for c in row] for row in self.coordinate_matrix])
            M = sp.Matrix([[ _rat(c) for c in row] for row in self.equation_matrix])
            if N.shape != (4,4) or M.shape != (2,2) or not N.det() or not M.det():
                return False
            H = [N.T*A*N for A in source.matrices()]
            expected = [sum((M[i,j]*H[j] for j in range(2)), sp.zeros(4)) for i in range(2)]
            # Compare matrices over Q; do not truncate fractional forged data.
            return all(A == B for A,B in zip(expected,self.target.matrices()))
        except (TypeError, ValueError):
            return False


def improve_integral_model(source, p, *, _smooth_checked=False):
    """Return an exactly equivalent Q_p model; never discard a residue class.

    In the common-factor step the coordinate matrix has determinant p, while
    both new equation rows are divisible by p. Thus the discriminant valuation
    drops by at least 12. Row saturation also strictly decreases it. A single
    row that is a scalar square permits a non-increasing step; at most four
    such neutral steps are allowed between strict decreases. Smooth input
    therefore permits only finitely many steps. The transformations are
    invertible over Q_p (indeed Q), so no branch is lost.
    """
    source = QuadraticSystem.from_quadrics(source)
    p = prime(p)
    if not _smooth_checked:
        A,B = source.matrices()
        t = sp.Symbol("t")
        regular = False
        for shift in range(5):
            M = A+shift*B
            if M.det():
                regular = bool(sp.Poly((t*M+B).det(),t).discriminant())
                break
        if not regular:
            raise NotImplementedError("integral model improvement requires a smooth intersection")
    rows = list(source.coefficients)
    N, M, steps = sp.eye(4), sp.eye(2), 0
    neutral_steps = 0
    while True:
        pivot = next(i for i,c in enumerate(rows[0]) if c % p)
        multiplier = rows[1][pivot]*pow(rows[0][pivot] % p,-1,p) % p
        difference = tuple(b-multiplier*a for a,b in zip(rows[0],rows[1]))
        if not any(c % p for c in difference):
            if not any(difference):
                break  # A bounded caller may intentionally supply singular data.
            rows[1], divisor = _primitive(difference)
            M[1,:] = (M[1,:]-multiplier*M[0,:])/divisor
            steps += 1
            neutral_steps = 0
            continue
        expressions = [sum(c*U[i]*U[j] for c,(i,j) in zip(row,MONOMIALS)) for row in rows]
        common = sp.gcd(*(sp.Poly(q,*U,modulus=p) for q in expressions))
        common_factor = common.total_degree() == 1
        if common_factor:
            linear = [int(common.coeff_monomial(u)) % p for u in U]
        else:
            if neutral_steps >= 4:
                break
            squares = [_square_factor(row,p) for row in rows]
            index = next((i for i,linear in enumerate(squares) if linear is not None),None)
            if index is None:
                break
            linear = list(squares[index])
            if index == 1:
                rows.reverse()
                M.row_swap(0,1)
        pivot = next(i for i,c in enumerate(linear) if c)
        inverse = pow(linear[pivot],-1,p)
        linear = [(c*inverse) % p for c in linear]
        linear = [c-p if c > p//2 else c for c in linear]
        B = sp.eye(4)
        for i,c in enumerate(linear):
            B[pivot,i] = p if i == pivot else -c
        matrices = QuadraticSystem(tuple(rows)).matrices()
        new_rows, divided_valuation = [], 0
        for i,A in enumerate(matrices):
            transformed = _coefficients(B.T*A*B)
            if (common_factor or i == 0) and any(c % p for c in transformed):
                raise ArithmeticError("linear-factor model transformation failed")
            row, divisor = _primitive(transformed)
            divided_valuation += _vp(abs(divisor),p)
            new_rows.append(row)
            M[i,:] = M[i,:]/divisor
        if divided_valuation < 1:
            raise ArithmeticError("model improvement increased the discriminant valuation")
        neutral_steps = neutral_steps+1 if divided_valuation == 1 else 0
        rows, N = new_rows, N*B
        steps += 1
    return IntegralModelChange(source.key, QuadraticSystem(tuple(rows)),
                               tuple(tuple(c for c in row) for row in N.tolist()),
                               tuple(tuple(c for c in row) for row in M.tolist()), steps)


@dataclass(frozen=True)
class ChangedModelCertificate:
    system_key: str
    place: int
    change: IntegralModelChange
    inner: LocalResult
    kind: str = "model_change"


class ChangedModelResult(LocalResult):
    """A local decision on an exactly equivalent model, with checked maps."""

    def verify(self, quadrics, *, max_nodes=1_000_000):
        c = self.certificate
        return bool(isinstance(c,ChangedModelCertificate)
                    and c.system_key == QuadraticSystem.from_quadrics(quadrics).key
                    and c.place == self.place == c.inner.place
                    and c.inner.status == self.status != LocalStatus.UNKNOWN
                    and c.change.verify(quadrics)
                    and c.inner.verify(c.change.target,max_nodes=max_nodes))


def transport_local_result(source, change, result):
    source = QuadraticSystem.from_quadrics(source)
    if not change.verify(source):
        raise ArithmeticError("local model-change certificate failed")
    if result.status == LocalStatus.UNKNOWN:
        return LocalResult(result.place,result.status,nodes=result.nodes,depth=result.depth,
                           reason=f"after {change.steps} integral model changes: {result.reason}")
    if result.certificate.kind == "rational_point":
        point = _point(tuple(sp.Matrix(change.coordinate_matrix)*sp.Matrix(result.certificate.point)))
        if source.values(point) != (0,0):
            raise ArithmeticError("rational point transport failed")
        certificate = LocalCertificate(source.key,result.place,"rational_point",point)
        return LocalResult(result.place,result.status,certificate,result.nodes,result.depth)
    certificate = ChangedModelCertificate(source.key,result.place,change,result)
    return ChangedModelResult(result.place,result.status,certificate,result.nodes,result.depth,
                              reason=f"after {change.steps} integral model changes")
