"""Exact quadratic square tests and Hilbert symbols over Q_p and R.

Inputs are nonzero rationals. No floating approximation or p-adic precision
is used: valuations and unit residues determine the quadratic square class.
Values are bits: 0 means +1 (zero in Q/Z), 1 means -1 (one half in Q/Z).
"""
from .arithmetic import prime
from .geometry import _rat


def valuation_unit(value, p, modulus=None):
    value = _rat(value)
    if not value:
        raise ValueError('zero has no multiplicative square class')
    p = prime(p)
    n, d, valuation = int(value.p), int(value.q), 0
    while n % p == 0:
        n //= p
        valuation += 1
    while d % p == 0:
        d //= p
        valuation -= 1
    modulus = p if modulus is None else modulus
    return valuation, (n * pow(d, -1, modulus)) % modulus


def local_square(value, place):
    value = _rat(value)
    if not value:
        raise ValueError('use nonzero values for open square-class tests')
    if place == 'real':
        return value > 0
    p = prime(place)
    v, unit = valuation_unit(value, p, 8 if p == 2 else p)
    return v % 2 == 0 and (unit == 1 if p == 2 else pow(unit, (p-1)//2, p) == 1)


def hilbert_bit(left, right, place):
    a, b = _rat(left), _rat(right)
    if not a or not b:
        raise ValueError('Hilbert symbols require nonzero arguments')
    if place == 'real':
        return int(bool(a < 0 and b < 0))
    p = prime(place)
    alpha, u = valuation_unit(a, p, 8 if p == 2 else p)
    beta, v = valuation_unit(b, p, 8 if p == 2 else p)
    if p == 2:
        return (((u-1)//2)*((v-1)//2) + alpha*((v*v-1)//8) + beta*((u*u-1)//8)) % 2
    lu = int(pow(u, (p-1)//2, p) == p-1)
    lv = int(pow(v, (p-1)//2, p) == p-1)
    return (alpha*beta*((p-1)//2) + beta*lu + alpha*lv) % 2
