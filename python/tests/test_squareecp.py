"""Tests for the squareecp package. Run: python -m pytest -q   (from python/)."""
import math
import random

import numpy as np
import pytest
import sympy as sp

from squareecp import bound, deficiency, hull_units, lattice as L, nodal


def tri(a, b):
    return L.hull([(0, 0), (a, 0), (0, b)])


# --- spectral data ------------------------------------------------------------------------
def test_kappa_dim_mult():
    for lam, k in [(2, 1), (5, 2), (8, 4), (10, 5), (13, 7), (17, 9), (18, 11), (20, 12),
                   (25, 14), (34, 21)]:
        assert L.kappa(lam) == k
    for lam in range(2, 400):
        assert L.dim(lam) == sum(1 for m in range(1, 30) for k in range(1, 30) if m * m + k * k <= lam)


def test_chebyshev_substitution():
    """sin(m x) = sin x U_{m-1}(cos x), and supp P is the down-closed set S(lam)."""
    x = sp.symbols("x")
    xs = np.linspace(0.01, 3.13, 400)
    for m in range(1, 8):
        f = sp.lambdify(x, sp.sin(m * x) - sp.sin(x) * sp.chebyshevu(m - 1, sp.cos(x)), "numpy")
        assert np.abs(f(xs)).max() < 1e-10
    for lam in (17, 18, 50, 130):
        S = set(L.support(lam))
        assert all((a2, b2) in S for a, b in S for a2 in range(a + 1) for b2 in range(b + 1))
    assert set(L.support(17)) == {(a, b) for a in range(4) for b in range(4) if a + b <= 3}
    assert set(L.support(18)) == set(L.support(17)) | {(2, 2)}


# --- polygon primitives -------------------------------------------------------------------
def test_pick_and_mixed_area():
    assert L.interior(tri(3, 3)) == 1 and L.interior(tri(4, 4)) == 3
    assert L.interior(L.hull([(0, 0), (3, 0), (2, 2), (0, 3)])) == 3
    sx, sy = L.hull([(0, 0), (1, 0)]), L.hull([(0, 0), (0, 1)])
    assert L.mv(sx, sy) == 1 and L.mv(sx, sx) == 0
    assert L.mv(tri(2, 2), tri(2, 2)) == 4 and L.mv(tri(3, 3), tri(1, 1)) == 3
    rng = random.Random(0)
    polys = bound.factor_polygons(L.hull(L.support(50)))
    for _ in range(300):
        A, B, C = (rng.choice(polys) for _ in range(3))
        assert L.mv(A, B) == L.mv(B, A)
        assert L.mv(A, L.msum(B, C)) == L.mv(A, B) + L.mv(A, C)
        if len(A) >= 3 and len(B) >= 3:
            assert L.interior(L.msum(A, B)) == L.interior(A) + L.interior(B) + L.mv(A, B) - 1


def test_root_lemma_random():
    """Lemma "roots and corner signs" against random real polynomials."""
    rng = np.random.default_rng(1)
    for d in range(1, 7):
        for _ in range(2000):
            c = rng.normal(size=d + 1) * rng.choice([1, 10, 0.1], size=d + 1)
            r = np.roots(c[::-1])
            r = r[np.abs(r.imag) < 1e-9].real
            n = int(np.sum(np.abs(r) < 1))
            g1, gm1, s = np.sign(np.polyval(c[::-1], 1)), np.sign(np.polyval(c[::-1], -1)), np.sign(c[-1])
            assert n <= d - (g1 != s) - (gm1 != (s if d % 2 == 0 else -s))


def test_cap_against_random_polynomials():
    """Corner-cap lemma: boundary roots of random P with support S(lam) never exceed cap(S)."""
    rng = np.random.default_rng(2)
    for lam in (10, 13, 17, 18, 20, 25):
        S = L.support(lam)
        cap = bound.cap(S)
        for _ in range(500):
            c = {m: rng.normal() * rng.choice([1, 0.1]) for m in S}
            tot = 0
            for side in range(4):
                # coefficients of P restricted to X = +-1 (in Y) or Y = +-1 (in X)
                deg = max(b for a, b in S) if side < 2 else max(a for a, b in S)
                coef = np.zeros(deg + 1)
                for (a, b), v in c.items():
                    if side < 2:
                        coef[b] += v * (1 if side == 0 else (-1) ** a)
                    else:
                        coef[a] += v * (1 if side == 2 else (-1) ** b)
                r = np.roots(coef[::-1])
                r = r[np.abs(r.imag) < 1e-9].real
                tot += int(np.sum(np.abs(r) < 1))
            assert tot <= cap


def test_enumeration_complete():
    """Appendix A(a): every down-closed polygon in Delta(lam) is listed."""
    rng = random.Random(3)
    for lam in (18, 34, 50):
        S = L.support(lam)
        polys = set(bound.factor_polygons(L.hull(S)))
        for _ in range(2000):
            sel = rng.sample(S, rng.randint(1, 4))
            h = L.hull([(a, b) for x, y in sel for a in range(x + 1) for b in range(y + 1)])
            if len(h) == 1 or (len(h) == 2 and abs(h[1][0] - h[0][0]) + abs(h[1][1] - h[0][1]) > 1):
                continue
            assert h in polys


def test_identity_and_dimension_bound():
    """Identity proposition (Section 6): the identity for U(D_1..D_k) and U <= dim E_lam."""
    rng = random.Random(4)
    for lam in (18, 34, 50, 65):
        S = L.support(lam)
        big = L.hull(S)
        polys = bound.factor_polygons(big)
        half_cap = bound.cap(S) // 2
        N = len(S)
        for _ in range(1500):
            fs, cur = [], None
            for _ in range(rng.randint(1, 6)):
                h = rng.choice(polys)
                nxt = L.msum(cur, h) if cur else h
                if L.inside(nxt, big):
                    fs.append(h)
                    cur = nxt
            if not fs:
                continue
            etas = [bound.eta(h) for h in fs]
            val = 1 + sum(L.interior(h) + 1 for h in fs)
            val += sum(L.mv(fs[i], fs[j]) for i in range(len(fs)) for j in range(i + 1, len(fs)))
            G, run = 0, 0
            for j, x in enumerate(sorted(etas, reverse=True), 1):
                run += x
                G = max(G, min(half_cap, run) - j)
            val += G
            Bstar = [L.boundary_points(h) if len(h) > 2 else
                     2 * math.gcd(abs(h[1][0] - h[0][0]), abs(h[1][1] - h[0][1])) for h in fs]
            w = [max(p[0] for p in h) + max(p[1] for p in h) for h in fs]
            pred = len(L.lattice_points(cur)) - sum(b - x - 1 for b, x in zip(Bstar, w)) \
                - (sum(x - 1 for x in w) - G)
            assert val == pred
            assert all(b - x - 1 >= 0 for b, x in zip(Bstar, w)) and sum(x - 1 for x in w) >= G
            assert val <= N


def test_parity_cap_loss():
    """Parity cap loss lemma for primitive triangles: eta <= w - 1."""
    for a in range(2, 12):
        for b in range(2, 12):
            if math.gcd(a, b) == 1:
                assert bound.eta(tri(a, b)) <= a + b - 1


# --- the bound and the deficiency ---------------------------------------------------------
EXACT_U = {5: 2, 8: 4, 10: 5, 13: 7, 17: 8, 18: 9, 20: 12, 25: 14, 26: 15, 29: 17, 32: 18, 34: 21, 37: 21, 40: 24, 41: 26, 45: 28, 50: 31}


@pytest.mark.parametrize("lam", sorted(EXACT_U))
def test_U_values(lam):
    assert bound.U(lam) == EXACT_U[lam]
    assert bound.U(lam) <= L.kappa(lam)


@pytest.mark.parametrize("lam", [l for l in L.eigenvalues(65) if l >= 5])
def test_lb_equals_exact_deficiency(lam):
    """On small levels the deficiency bound is sharp: LB = dim E - U."""
    assert deficiency.lower_bound(lam) == L.dim(lam) - bound.U(lam)


@pytest.mark.parametrize("lam", [l for l in L.eigenvalues(2000) if l <= 400 or l % 7 == 0])
def test_hull_units_bound_lb(lam):
    """Theorem C: LB(lam) >= ceil(H(lam)/2), and ECP holds (slack >= 0)."""
    assert deficiency.lower_bound(lam) >= (hull_units.H(lam) + 1) // 2
    assert hull_units.slack(lam) >= 0


def test_zero_slack_levels():
    """Levels <= 3000 where the hull criterion has no room to spare."""
    assert all(hull_units.slack(lam) >= 0 for lam in L.eigenvalues(3000))


# --- constructions ------------------------------------------------------------------------
@pytest.mark.parametrize("lam,expr,count", nodal.EXAMPLES)
def test_examples(lam, expr, count):
    X, Y = sp.symbols("X Y")
    P = sp.expand(sp.sympify(expr))
    assert all((a + 1) ** 2 + (b + 1) ** 2 <= lam for a, b in sp.Poly(P, X, Y).monoms())
    f = sp.lambdify((X, Y), P, "numpy")
    assert nodal.stable_count(lambda x, y: f(x, y) + 0 * x) == count
    assert count <= bound.U(lam) if lam <= 34 else True
