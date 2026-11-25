"""LB(lam), an exact lower bound for the deficiency dim E_lam - U(lam) (Proposition "LB" in Section 7).

    LB(lam) = min over down-closed lattice polygons D in Delta(lam) of
              miss(D) + ceil((n(D) + phi(o(D))) / 2),

where miss(D) = #(S(lam) minus D), n(D) is the number of slanted primitive units of the
upper chain of D, o(D) the number of parity classes occurring an odd number of times among
them, and phi(0, 1, 2, 3) = (0, 1, 2, 1). Computed by a convex-chain dynamic programme.
"""
import math

from .lattice import support

PHI = (0, 1, 2, 1)


def column_heights(lam):
    S = support(lam)
    W = max(a for a, b in S)
    return [max(b for a, b in S if a == x) for x in range(W + 1)]


def lower_bound(lam):
    Hc = column_heights(lam)
    W = len(Hc) - 1
    Hy = Hc[0]
    INF = 10 ** 18

    def inS(x, y):
        return 0 <= x <= W and 0 <= y <= Hc[x]
    tail = [0] * (W + 2)
    for x in range(W, -1, -1):
        tail[x] = tail[x + 1] + Hc[x] + 1
    # state (x, y, parity bits); costs are doubled so that the half-units stay integral
    best = {}
    for x0 in range(W + 1):
        for y0 in range(Hc[x0] + 1):
            c = 2 * (tail[x0 + 1] + (Hc[x0] - y0))
            best[(x0, y0, 0)] = min(best.get((x0, y0, 0), INF), c)
    dirs = [(a, b) for a in range(1, W + 1) for b in range(0, Hy + 1) if math.gcd(a, b) == 1]
    dirs.sort(key=lambda d: -d[1] / d[0])            # decreasing slope keeps the chain convex
    for a, b in dirs:
        slanted = 1 if b >= 1 else 0
        cls = (1 if a % 2 else 0) | (2 if b % 2 else 0)
        for key in sorted(best, key=lambda k: -k[0]):
            cost = best[key]
            cx, cy, u = key
            while True:
                nx, ny = cx - a, cy + b
                if not inS(nx, ny):
                    break
                for col in range(nx, cx):             # columns under the new edge
                    cost += 2 * (Hc[col] - (cy + (b * (cx - col)) // a))
                cost += slanted
                if slanted:
                    u ^= 1 << (cls - 1)
                k = (nx, ny, u)
                if cost < best.get(k, INF):
                    best[k] = cost
                cx, cy = nx, ny
    total = min(best[k] + PHI[bin(k[2]).count("1")] for k in best if k[0] == 0)
    return -(-total // 2)
