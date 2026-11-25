"""H(lam): slanted primitive units on the upper chain of Delta(lam) = conv S(lam).

Theorem C of the paper: dim E_lam - U(lam) >= ceil(H(lam) / 2). Hence the extended Courant
property holds at lam whenever H(lam) >= 2 r(lam) - 3. The C programs in ../c sweep this
criterion up to 10^9; this module is the slow reference implementation.
"""
from math import gcd, isqrt

from .lattice import mult


def upper_chain(lam):
    """Vertices of the upper chain of conv S(lam), left to right."""
    pts = []
    x = 0
    while True:
        rem = lam - (x + 1) ** 2
        if rem < 1:
            break
        c = isqrt(rem) - 1
        if c < 0:
            break
        while len(pts) >= 2:
            (x1, y1), (x2, y2) = pts[-2], pts[-1]
            if (x2 - x1) * (c - y1) - (y2 - y1) * (x - x1) >= 0:
                pts.pop()
            else:
                break
        pts.append((x, c))
        x += 1
    return pts


def H(lam):
    ch = upper_chain(lam)
    return sum(gcd(x2 - x1, y1 - y2) for (x1, y1), (x2, y2) in zip(ch, ch[1:])
               if x2 != x1 and y2 != y1)


def slack(lam):
    """ceil(H/2) - (r - 1); ECP holds at lam when this is >= 0."""
    return (H(lam) + 1) // 2 - (mult(lam) - 1)
