"""The bound U(lam) of Theorem A, computed exactly.

U(lam) is the maximum over multisets D_1, ..., D_k of factor polygons (Definition "The bound" of the
paper) with D_1 + ... + D_k inside Delta(lam) = conv S(lam) of

    1 + sum_i (I(D_i) + 1) + sum_{i<j} MV(D_i, D_j)
      + max_T [ min(floor(cap(S(lam)) / 2), sum_{i in T} eta(D_i)) - |T| ].

The search is exhaustive (Appendix A of the paper).
"""
import functools
import itertools

from .lattice import area2, hull, interior, msum, mv, inside, lattice_points, support


def cap(A):
    """Corner-sign number cap(A) of a finite exponent set A (Definition "corner-sign number")."""
    dx = max(a for a, b in A)
    dy = max(b for a, b in A)
    top_even = not any(a % 2 for a, b in A if b == dy)
    right_even = not any(b % 2 for a, b in A if a == dx)

    def beta(d, u, v, s):
        # roots in (-1, 1) of a degree-<=d polynomial with sign u at -1, v at +1, leading sign s
        if d == 0:
            return 0
        return d - (v != s) - (u != (s if d % 2 == 0 else -s))
    best = 0
    for pp, pm, mp, mm in itertools.product((-1, 1), repeat=4):   # signs at (1,1),(1,-1),(-1,1),(-1,-1)
        for R, L, T, B in itertools.product((-1, 1), repeat=4):
            if top_even and L != R:
                continue
            if right_even and B != T:
                continue
            best = max(best, beta(dy, pm, pp, R) + beta(dy, mm, mp, L)
                       + beta(dx, mp, pp, T) + beta(dx, mm, pm, B))
    return best


def factor_polygons(big):
    """All down-closed lattice polygons in big, other than points and axis segments of
    lattice length >= 2 (the set D_lam). A down-closed set is a non-increasing profile."""
    W = max(p[0] for p in big)
    Hh = max(p[1] for p in big)
    out = set()

    def rec(x, prev, cols):
        if x > W or prev < 0:
            if cols:
                h = hull([(i, j) for i, c in enumerate(cols) for j in range(c + 1)])
                if len(h) >= 2 and inside(h, big):
                    out.add(h)
            return
        rec(W + 1, -1, cols)                # stop here
        for c in range(prev + 1):
            rec(x + 1, c, cols + [c])
    rec(0, Hh, [])
    good = []
    for h in out:
        if len(h) == 2:
            (a, b), (c, d) = h
            if abs(c - a) + abs(d - b) != 1:
                continue
        good.append(h)
    return sorted(good, key=lambda h: (-area2(h), h))


UNIT_TRIANGLE = hull([(0, 0), (1, 0), (0, 1)])


def half_ends(h):
    """h(D), the half-ends bound in the definition of U."""
    wx = max(p[0] for p in h)
    wy = max(p[1] for p in h)
    if h == UNIT_TRIANGLE:
        return 1
    if h == hull([(0, 0), (1, 0), (0, wy)]) or h == hull([(0, 0), (wx, 0), (0, 1)]):
        return max(wx, wy)
    return wx + wy


def eta(h):
    """eta(D) = min(h(D), floor(cap(D)/2))."""
    return min(half_ends(h), cap(lattice_points(h)) // 2)


def U(lam, return_maximiser=False):
    """The bound U(lam) (exact). Optionally also a maximising multiset of polygons."""
    S = support(lam)
    big = hull(S)
    polys = factor_polygons(big)
    data = [(eta(h), interior(h)) for h in polys]
    half_cap = cap(S) // 2

    def gain(etas):
        g, run = 0, 0
        for j, x in enumerate(sorted(etas, reverse=True), 1):
            run += x
            g = max(g, min(half_cap, run) - j)
        return g

    @functools.lru_cache(maxsize=None)
    def rec(start, cur, etas):
        val, arg = gain(etas), ()
        for i in range(start, len(polys)):
            h = polys[i]
            nxt = msum(cur, h) if cur else h
            if not inside(nxt, big):
                continue
            e, I = data[i]
            add = I + 1 + (mv(cur, h) if cur else 0)
            v, a = rec(i, nxt, tuple(sorted(etas + (e,))))
            if add + v > val:
                val, arg = add + v, (i,) + a
        return val, arg
    val, arg = rec(0, (), ())
    if return_maximiser:
        return 1 + val, [polys[i] for i in arg]
    return 1 + val


def describe(h):
    named = {hull([(0, 0), (1, 0)]): "vline", hull([(0, 0), (0, 1)]): "hline",
             UNIT_TRIANGLE: "line", hull([(0, 0), (1, 0), (0, 1), (1, 1)]): "hyperbola",
             hull([(0, 0), (2, 0), (0, 2)]): "conic"}
    return named.get(h, "poly" + str(list(h)))
