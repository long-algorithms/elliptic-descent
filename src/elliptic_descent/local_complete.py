"""Terminating reference p-adic decisions for smooth quadric intersections.

Complete mode has no depth/node cutoff. Finite-field elimination replaces a
p^3 scan; higher digits solve exact affine linear systems. Compactness plus
smooth Hensel lifting proves termination (docs/relative-descent.md).
"""
from dataclasses import dataclass
from itertools import count, product

import sympy as sp

from .arithmetic import nonnegative_integer, prime
from .geometry import U
from .local import (LocalCertificate, LocalResult, MONOMIALS, QuadraticSystem,
                    _digits, _hensel, bad_primes)
from .real import real_squareclass_solubility
from .local_models import improve_integral_model, transport_local_result
from .status import LocalStatus


class _Limit(Exception):
    pass


class _Budget:
    def __init__(self, limit):
        self.limit = limit
        self.used = 0

    def tick(self):
        if self.limit is not None and self.used >= self.limit:
            raise _Limit
        self.used += 1


def _roots(poly, variable, p):
    poly = sp.Poly(poly, variable, modulus=p)
    if poly.is_zero:
        return range(p)
    return sorted(int(root) % p for root in poly.ground_roots())


def _bivariate_points(first, second, y, z, p, budget):
    """All common F_p zeros of two quadratics; degenerate fibres included."""
    q, r = (sp.Poly(expr, z) for expr in (first, second))
    A, B, C = (sp.Poly(q.nth(i), y, modulus=p).as_expr() for i in (2, 1, 0))
    D, E, F = (sp.Poly(r.nth(i), y, modulus=p).as_expr() for i in (2, 1, 0))
    # Total degree two makes A,D constant, including after specialization.
    A, D = int(A) % p, int(D) % p
    if not A and D:
        A, B, C, D, E, F = D, E, F, A, B, C
    if A:
        b, c = sp.expand(A*E-D*B), sp.expand(A*F-D*C)
        eliminant = sp.expand(A*c*c-B*b*c+C*b*b)
    else:
        eliminant = sp.expand(B*F-E*C)
    for value in _roots(eliminant, y, p):
        budget.tick()
        bb, cc, ee, ff = (int(expr.subs(y, value)) % p for expr in (B, C, E, F))
        if A:
            linear, constant = (A*ee-D*bb) % p, (A*ff-D*cc) % p
            if linear:
                zs = [(-constant*pow(linear, -1, p)) % p]
            elif constant:
                continue
            else:
                zs = _roots(A*z*z+bb*z+cc, z, p)
        elif bb:
            zs = [(-cc*pow(bb, -1, p)) % p]
        elif ee:
            zs = [(-ff*pow(ee, -1, p)) % p]
        elif cc or ff:
            continue
        else:
            zs = range(p)
        for value_z in zs:
            budget.tick()
            # Check every reconstructed root; eliminants can have extra roots.
            if ((A*value_z*value_z+bb*value_z+cc) % p == 0
                    and (D*value_z*value_z+ee*value_z+ff) % p == 0):
                yield value, value_z


def _affine_points(equations, variables, p, budget):
    n = len(variables)
    if n == 3:
        x, y, z = variables
        for value in range(p):
            budget.tick()
            specialized = tuple(sp.expand(q.subs(x, value)) for q in equations)
            for pair in _bivariate_points(*specialized, y, z, p, budget):
                yield (value, *pair)
    elif n == 2:
        yield from _bivariate_points(*equations, *variables, p, budget)
    elif n == 1:
        x = variables[0]
        g = sp.gcd(*(sp.Poly(q, x, modulus=p) for q in equations))
        for value in _roots(g.as_expr(), x, p):
            budget.tick()
            yield (value,)
    else:
        budget.tick()
        if all(int(q) % p == 0 for q in equations):
            yield ()


def _initial_solutions(system, p, budget):
    expressions = tuple(sum(c*U[i]*U[j] for c, (i, j) in zip(row, MONOMIALS))
                        for row in system.coefficients)
    for pivot in range(4):
        substitutions = {U[i]: 0 for i in range(pivot)}
        substitutions[U[pivot]] = 1
        equations = tuple(sp.expand(q.subs(substitutions)) for q in expressions)
        for suffix in _affine_points(equations, U[pivot+1:], p, budget):
            yield pivot, (0,)*pivot+(1,)+suffix


def _linear_solutions(matrix, target, p):
    """All solutions of a 2x3 system, with no probabilistic rank decisions."""
    rows = [[int(x) % p for x in row]+[int(rhs) % p]
            for row, rhs in zip(matrix, target)]
    pivots = []
    for j in range(3):
        k = next((k for k in range(len(pivots), 2) if rows[k][j]), None)
        if k is None:
            continue
        i = len(pivots)
        rows[i], rows[k] = rows[k], rows[i]
        inverse = pow(rows[i][j], -1, p)
        rows[i] = [(x*inverse) % p for x in rows[i]]
        for k in range(2):
            if k != i:
                multiple = rows[k][j]
                rows[k] = [(a-multiple*b) % p for a, b in zip(rows[k], rows[i])]
        pivots.append(j)
    if any(not any(row[:3]) and row[3] for row in rows):
        return
    free = tuple(j for j in range(3) if j not in pivots)
    for values in _digits(p, len(free)):
        point = [0]*3
        for j, value in zip(free, values):
            point[j] = value
        for i, j in enumerate(pivots):
            point[j] = (rows[i][3]-sum(rows[i][k]*point[k] for k in free)) % p
        yield tuple(point)


def _next_solutions(system, frontier, modulus, p, budget):
    for pivot, point in frontier:
        free = tuple(i for i in range(4) if i != pivot)
        jacobian = system.jacobian(point)
        matrix = [[row[i] for i in free] for row in jacobian]
        target = [-v//modulus for v in system.values(point)]
        for digits in _linear_solutions(matrix, target, p):
            budget.tick()
            lifted = list(point)
            for i, digit in zip(free, digits):
                lifted[i] += modulus*digit
            yield pivot, tuple(lifted)


def padic_solubility_complete(quadrics, p, *, max_depth=None, max_nodes=None,
                              _smooth_checked=False, _model_improved=False):
    """Decide Q_p-solubility of a smooth intersection, or return bounded unknown.

    Node counts are algebraic search work units, not all primitive residues.
    None means no artificial cutoff. The supported smooth domain is checked
    before an unbounded search, unless the caller has already checked it.
    """
    system, p = QuadraticSystem.from_quadrics(quadrics), prime(p)
    max_depth = None if max_depth is None else nonnegative_integer(max_depth, "max_depth")
    max_nodes = None if max_nodes is None else nonnegative_integer(max_nodes, "max_nodes")
    if not _smooth_checked and max_depth is None and max_nodes is None:
        bad_primes(system)  # Nonzero determinant discriminant implies smoothness.
        _smooth_checked = True
    if not _model_improved and max_depth != 0 and max_nodes != 0:
        try:
            change = improve_integral_model(system,p,_smooth_checked=_smooth_checked)
        except NotImplementedError:
            change = None  # Bounded searches can also examine singular test models.
        if change is not None and change.steps:
            report = padic_solubility_complete(change.target,p,max_depth=max_depth,
                                               max_nodes=max_nodes,_smooth_checked=True,
                                               _model_improved=True)
            return transport_local_result(system,change,report)
    budget = _Budget(max_nodes)
    frontier, modulus, depth = [], 1, 0
    try:
        for depth in count(1):
            if max_depth is not None and depth > max_depth:
                return LocalResult(p, LocalStatus.UNKNOWN, nodes=budget.used,
                                   depth=max_depth, reason="depth budget exhausted")
            solutions = (_initial_solutions(system, p, budget) if depth == 1
                         else _next_solutions(system, frontier, modulus, p, budget))
            new_frontier = []
            for pivot, point in solutions:
                if any(value % (modulus*p) for value in system.values(point)):
                    raise ArithmeticError("invalid lifted residue")
                certificate = _hensel(system, point, p, pivot)
                if certificate:
                    return LocalResult(p, LocalStatus.SOLUBLE, certificate, budget.used, depth)
                new_frontier.append((pivot, point))
            if not new_frontier:
                certificate = LocalCertificate(system.key, p, "empty_modulus", depth=depth)
                return LocalResult(p, LocalStatus.INSOLUBLE, certificate, budget.used, depth)
            frontier, modulus = new_frontier, modulus*p
    except _Limit:
        return LocalResult(p, LocalStatus.UNKNOWN, nodes=budget.used, depth=depth,
                           reason="node budget exhausted")


@dataclass(frozen=True)
class CompleteLocalScreen:
    real: object
    finite: tuple
    good_reduction: object

    @property
    def status(self):
        if self.real.status == LocalStatus.INSOLUBLE or any(r.status == LocalStatus.INSOLUBLE for r in self.finite):
            return LocalStatus.INSOLUBLE
        if any(r.certificate and r.certificate.kind == "rational_point" for r in self.finite):
            return LocalStatus.SOLUBLE
        required = set(self.good_reduction.primes) if self.good_reduction is not None else None
        if (self.real.status == LocalStatus.SOLUBLE and required is not None
                and {r.place for r in self.finite if r.status == LocalStatus.SOLUBLE} >= required):
            return LocalStatus.SOLUBLE
        return LocalStatus.UNKNOWN

    def verify(self, model, *, max_nodes=1_000_000):
        if self.status == LocalStatus.UNKNOWN or not self.real.verify(model):
            return False
        if self.real.status == LocalStatus.INSOLUBLE:
            return self.status == LocalStatus.INSOLUBLE
        if self.good_reduction != bad_primes(model):
            return False
        return all(r.verify(model, max_nodes=max_nodes) for r in self.finite if r.status != LocalStatus.UNKNOWN)


def complete_local_screen(model, *, max_depth=None, max_nodes=None):
    real = real_squareclass_solubility(model)
    if real.status == LocalStatus.INSOLUBLE:
        return CompleteLocalScreen(real, (), None)
    bound = bad_primes(model)
    finite = []
    for p in bound.primes:
        report = padic_solubility_complete(model, p, max_depth=max_depth, max_nodes=max_nodes,
                                           _smooth_checked=True)
        finite.append(report)
        if report.status == LocalStatus.INSOLUBLE or (report.certificate and report.certificate.kind == "rational_point"):
            break
    return CompleteLocalScreen(real, tuple(finite), bound)
