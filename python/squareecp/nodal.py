"""Numerical nodal-domain counts for explicit polynomials P(X, Y) on (-1, 1)^2.

These are evidence for the lower-bound constructions, not proofs: a cell of an n x n grid is
a wall cell if its four corner values do not share a strict sign, and the domains are the
4-connected components of the remaining positive and negative cells. Positive and negative
cells never share an edge, so domains cannot leak into each other at crossings; a domain
thinner than about one cell can be lost, so counts are compared at several resolutions.
"""
import numpy as np
from scipy import ndimage


def count(P, n):
    """Nodal domains of the polynomial function P (vectorised callable) on an n-cell grid."""
    g = np.linspace(-1, 1, n + 1)
    X, Y = np.meshgrid(g, g, indexing="ij")
    V = P(X, Y)
    tol = 1e-12 * np.abs(V).max()
    pos, neg = V > tol, V < -tol
    cpos = pos[:-1, :-1] & pos[1:, :-1] & pos[:-1, 1:] & pos[1:, 1:]
    cneg = neg[:-1, :-1] & neg[1:, :-1] & neg[:-1, 1:] & neg[1:, 1:]
    return ndimage.label(cpos)[1] + ndimage.label(cneg)[1]


def stable_count(P, resolutions=(1200, 2399, 3001)):
    """The common count at several resolutions, or None if they disagree."""
    c = {count(P, n) for n in resolutions}
    return c.pop() if len(c) == 1 else None


# Explicit extremal examples (Section "Exact values and sharpness"): (lam, P(X, Y), nodal count).
# At lam = 5, 8, 10, 13, 20, 25, 34 the count equals kappa(lam); at 17 and 18 it equals U(lam).
EXAMPLES = [
    (5,  "X", 2),
    (8,  "X*Y", 4),
    (10, "X**2 + Y**2 - 3/2", 5),
    (13, "X*Y*(X + Y - 1/2)", 7),
    (17, "(X**2 + Y**2 - 3/2)*(Y - X)", 8),
    (18, "(X**2 - 9/25)*(Y**2 - 9/25)", 9),
    (20, "(X - 17/20)*(Y - 17/20)*(X**2 + Y**2 - 3/2)", 12),
    (25, "(X**2 - 9/100)*(Y**2 - 9/100)*(Y - X - 1/10)", 14),
    (34, "(X**2 - 289/400)*(Y**2 - 289/400)*(X**2 + Y**2 - 3/2)", 21),
]
