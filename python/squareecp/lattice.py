"""Spectral data of the square and exact lattice-polygon primitives.

A combination of Dirichlet eigenfunctions sin(mx) sin(ky) of (0, pi)^2 with
m^2 + k^2 <= lam equals sin x sin y P(cos x, cos y) with supp P in

    S(lam) = {(a, b) in Z_{>=0}^2 : (a + 1)^2 + (b + 1)^2 <= lam}.

Polygons are tuples of integer vertices in counterclockwise order (a point is a 1-tuple,
a segment a 2-tuple). All arithmetic is exact integer arithmetic.
"""
from math import gcd, isqrt


# --- spectrum -----------------------------------------------------------------------------
def support(lam):
    """S(lam), the exponent set of P."""
    n = isqrt(lam) + 1
    return [(a, b) for a in range(n) for b in range(n) if (a + 1) ** 2 + (b + 1) ** 2 <= lam]


def dim(lam):
    """dim E_lam = #S(lam)."""
    return len(support(lam))


def mult(lam):
    """Multiplicity r(lam) = #{(m, k) in Z_{>=1}^2 : m^2 + k^2 = lam}."""
    c = 0
    for m in range(1, isqrt(lam) + 1):
        k2 = lam - m * m
        k = isqrt(k2)
        if k >= 1 and k * k == k2:
            c += 1
    return c


def kappa(lam):
    """kappa(lam) = 1 + #{eigenvalues < lam}, counted with multiplicity."""
    return dim(lam) - mult(lam) + 1


def eigenvalues(upto):
    """Distinct Dirichlet eigenvalues m^2 + k^2 <= upto, in increasing order."""
    return sorted({m * m + k * k for m in range(1, isqrt(upto) + 1)
                   for k in range(1, isqrt(upto) + 1) if m * m + k * k <= upto})


# --- polygons ----------------------------------------------------------------------------
def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull(points):
    """Convex hull (monotone chain), vertices counterclockwise, collinear points dropped."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return tuple(pts)
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return tuple(lower[:-1] + upper[:-1])


def area2(h):
    """Twice the area."""
    if len(h) < 3:
        return 0
    return abs(sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2) in zip(h, h[1:] + h[:1])))


def boundary_points(h):
    """Number of lattice points on the boundary (for a segment: on the segment)."""
    if len(h) == 1:
        return 1
    if len(h) == 2:
        (x1, y1), (x2, y2) = h
        return gcd(abs(x2 - x1), abs(y2 - y1)) + 1
    return sum(gcd(abs(x2 - x1), abs(y2 - y1)) for (x1, y1), (x2, y2) in zip(h, h[1:] + h[:1]))


def interior(h):
    """I(h), the number of interior lattice points (Pick); 0 for points and segments."""
    if len(h) < 3:
        return 0
    return (area2(h) - boundary_points(h) + 2) // 2


def msum(h1, h2):
    """Minkowski sum."""
    return hull([(a + c, b + d) for a, b in h1 for c, d in h2])


def mv(h1, h2):
    """Mixed area, normalised so that two unit segments in different directions give 1."""
    return (area2(msum(h1, h2)) - area2(h1) - area2(h2)) // 2


def inside(h, big):
    """True if every vertex of h lies in the convex polygon big (counterclockwise)."""
    for p in h:
        for (x1, y1), (x2, y2) in zip(big, big[1:] + big[:1]):
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) < 0:
                return False
    return True


def _on_segment(p, h):
    (x1, y1), (x2, y2) = h
    return ((x2 - x1) * (p[1] - y1) == (y2 - y1) * (p[0] - x1)
            and min(x1, x2) <= p[0] <= max(x1, x2) and min(y1, y2) <= p[1] <= max(y1, y2))


def lattice_points(h):
    """All lattice points of the polygon (or segment, or point) h in the first quadrant."""
    if len(h) == 1:
        return [h[0]]
    W = max(p[0] for p in h)
    Hh = max(p[1] for p in h)
    box = [(a, b) for a in range(W + 1) for b in range(Hh + 1)]
    if len(h) == 2:
        return [p for p in box if _on_segment(p, h)]
    return [p for p in box if inside((p,), h)]
