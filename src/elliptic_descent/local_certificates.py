"""Independent checker for empty projective residue trees.

Unlike the solver, enumerate every initial projective residue directly; no
resultants or finite-field polynomial factorization are used. Later levels
use explicit rank-0/rank-1 compatibility checks, not the solver's RREF code.
"""
from itertools import product


def verify_empty_residue_tree(system, p, depth, max_nodes):
    """Return True only after all residue branches have been excluded.

    False includes insufficient verification resources. A unit Jacobian minor
    already proves solubility and therefore immediately refutes the certificate.
    """
    if depth < 1 or p**3+p**2+p+1 > max_nodes:
        return False
    frontier, used = [], 0
    for pivot in range(4):
        for suffix in product(range(p),repeat=3-pivot):
            if used >= max_nodes:
                return False
            used += 1
            point = (0,)*pivot+(1,)+suffix
            if all(v % p == 0 for v in system.values(point)):
                frontier.append((pivot,point))
    modulus = p
    for level in range(1,depth):
        if not frontier:
            return True
        following = []
        for pivot,point in frontier:
            if used >= max_nodes:
                return False
            used += 1
            indices = tuple(i for i in range(4) if i != pivot)
            jacobian = system.jacobian(point)
            A = [[row[i] % p for i in indices] for row in jacobian]
            rhs = [(-v//modulus) % p for v in system.values(point)]
            if any((A[0][i]*A[1][j]-A[0][j]*A[1][i]) % p
                   for i in range(3) for j in range(i+1,3)):
                return False  # Smooth F_p point, hence a Q_p point.
            row_index = next((i for i in (0,1) if any(A[i])),None)
            if row_index is None:
                if any(rhs):
                    continue
                free = (0,1,2)
                pivot_digit = None
            else:
                row = A[row_index]
                pivot_digit = next(i for i,c in enumerate(row) if c)
                inverse = pow(row[pivot_digit],-1,p)
                other = 1-row_index
                multiplier = A[other][pivot_digit]*inverse % p
                if rhs[other] != multiplier*rhs[row_index] % p:
                    continue
                free = tuple(i for i in range(3) if i != pivot_digit)
            for values in product(range(p),repeat=len(free)):
                if used >= max_nodes:
                    return False
                used += 1
                digits = [0]*3
                for i,value in zip(free,values):
                    digits[i] = value
                if pivot_digit is not None:
                    digits[pivot_digit] = ((rhs[row_index]-sum(row[i]*digits[i] for i in free))*inverse) % p
                lifted = list(point)
                for i,digit in zip(indices,digits):
                    lifted[i] += modulus*digit
                if any(v % (p*modulus) for v in system.values(lifted)):
                    return False
                following.append((pivot,tuple(lifted)))
        frontier,modulus = following,modulus*p
    return not frontier
