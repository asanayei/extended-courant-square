"""Figures of the paper. Run from this directory:  ../../.venv/bin/python make_figures.py

Writes fig_*.pdf into the parent directory (paper/). Everything is computed from the
package in ../../python and from data/ ; nothing is drawn by hand.
"""
import json
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Polygon, Rectangle
from scipy import ndimage

sys.path.insert(0, "../../python")
from squareecp import lattice as L, hull_units as Hu, nodal  # noqa: E402

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9.5,
                     "axes.edgecolor": "#52514e", "axes.linewidth": 0.8,
                     "xtick.color": "#52514e", "ytick.color": "#52514e",
                     "pdf.fonttype": 42})
POS, NEG, INK, MUTED, GRID = "#a9c9f0", "#f5bba3", "#0b0b0b", "#8a8984", "#e6e5e1"
BLUE, ORANGE = "#2a78d6", "#eb6834"
cmap = ListedColormap([NEG, POS])
OUT = "../"


def count_domains(F):
    """Grid count of the nodal domains of P = F in the box (sanity check of the pictures)."""
    c = nodal.stable_count(F)
    assert c is not None
    return c


def sign_panel(ax, F, coords="square", n=1200, lim=None, lw=1.1):
    if coords == "square":
        t = np.linspace(0, np.pi, n + 2)[1:-1]
        x, y = np.meshgrid(t, t)
        V = np.sin(x) * np.sin(y) * F(np.cos(x), np.cos(y)); ext = (0, np.pi)
    else:
        a = lim or 1.0
        t = np.linspace(-a, a, n + 2)[1:-1]
        x, y = np.meshgrid(t, t); V = F(x, y); ext = (-a, a)
    ax.imshow((V > 0).astype(float), origin="lower", extent=ext * 2, cmap=cmap, vmin=0, vmax=1,
              interpolation="nearest", rasterized=True)
    ax.contour(x, y, V, levels=[0], colors=INK, linewidths=lw)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(INK); s.set_linewidth(1.0)
    return V


def tidy(ax, grid="y"):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.6); ax.set_axisbelow(True)


EX = [(10, lambda X, Y: X**2 + Y**2 - 1.5, 5),
      (13, lambda X, Y: X*Y*(X + Y - .5), 7),
      (17, lambda X, Y: (X**2 + Y**2 - 1.5)*(Y - X), 8),
      (18, lambda X, Y: (X**2 - .36)*(Y**2 - .36), 9),
      (20, lambda X, Y: (X - .85)*(Y - .85)*(X**2 + Y**2 - 1.5), 12),
      (34, lambda X, Y: (X**2 - .7225)*(Y**2 - .7225)*(X**2 + Y**2 - 1.5), 21)]


# ---- Figure 1: gallery of extremal combinations ---------------------------------------------
def fig_gallery():
    fig, axs = plt.subplots(2, 3, figsize=(6.3, 4.55))
    for ax, (lam, F, nu) in zip(axs.flat, EX):
        sign_panel(ax, F)
        assert count_domains(F) == nu, (lam, count_domains(F))
        k = L.kappa(lam)
        tag = rf"$\nu={nu}=\kappa$" if nu == k else rf"$\nu={nu}$,  $\kappa={k}$"
        ax.set_title(rf"$\lambda={lam}$:   {tag}", fontsize=10, pad=5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01, wspace=0.06, hspace=0.16)
    fig.savefig(OUT + "fig_gallery.pdf", dpi=300)


# ---- Figure 2: the reduction ---------------------------------------------------------------
def fig_reduction():
    lam = 20
    F = EX[4][1]
    fig, axs = plt.subplots(1, 3, figsize=(6.5, 2.55), gridspec_kw={"width_ratios": [1, 1, 1.02]})
    sign_panel(axs[0], F)
    axs[0].set_title(r"$f=\sin x\,\sin y\;P(\cos x,\cos y)$", fontsize=9.5)
    axs[0].set_xlabel(r"on $(0,\pi)^2$", fontsize=9)
    # the algebraic curve, seen beyond the box
    ax = axs[1]; a = 1.7
    t = np.linspace(-a, a, 900); X, Y = np.meshgrid(t, t); V = F(X, Y)
    inside = (np.abs(X) < 1) & (np.abs(Y) < 1)
    img = np.where(inside, (V > 0).astype(float), np.nan)
    ax.imshow(img, origin="lower", extent=(-a, a, -a, a), cmap=cmap, vmin=0, vmax=1,
              interpolation="nearest", rasterized=True)
    ax.contour(X, Y, V, levels=[0], colors=MUTED, linewidths=0.9)
    ax.contour(X, Y, np.where(inside, V, np.nan), levels=[0], colors=INK, linewidths=1.2)
    ax.add_patch(Rectangle((-1, -1), 2, 2, fill=False, ec=INK, lw=1.0))
    ax.set_xlim(a, -a); ax.set_ylim(a, -a); ax.set_aspect("equal")  # cos is decreasing
    ax.set_xticks([-1, 1]); ax.set_yticks([-1, 1]); ax.tick_params(labelsize=7.5, length=2)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(r"the curve $P(X,Y)=0$", fontsize=9.5)
    ax.set_xlabel(r"and the box $B=(-1,1)^2$", fontsize=9)
    # exponents
    ax = axs[2]
    S = L.support(lam)
    used = {(3, 1), (1, 3), (1, 1), (3, 0), (1, 2), (1, 0), (2, 1), (0, 3), (0, 1), (2, 0), (0, 2), (0, 0)}
    assert used <= set(S)
    th = np.linspace(0, np.pi/2, 200)
    ax.plot(math.sqrt(lam)*np.cos(th) - 1, math.sqrt(lam)*np.sin(th) - 1, color=ORANGE, lw=1.3)
    ax.add_patch(Polygon(L.hull(S), closed=True, fc="#eef3fb", ec=BLUE, lw=1.0))
    xs, ys = zip(*S)
    ax.scatter(xs, ys, s=26, facecolors="white", edgecolors=INK, linewidths=0.9, zorder=3)
    xs, ys = zip(*used)
    ax.scatter(xs, ys, s=26, color=INK, zorder=4)
    ax.set_xlim(-0.6, 4.4); ax.set_ylim(-0.6, 4.4); ax.set_aspect("equal")
    ax.set_xticks(range(5)); ax.set_yticks(range(5)); ax.tick_params(labelsize=7.5, length=2)
    ax.set_xlabel(r"exponent $a$ of $X$", fontsize=9); ax.set_ylabel(r"exponent $b$ of $Y$", fontsize=9, labelpad=2)
    ax.set_title(r"$\mathrm{supp}\,P\subseteq S(20)$", fontsize=9.5)
    tidy(ax, None)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.17, wspace=0.16)
    fig.savefig(OUT + "fig_reduction.pdf", dpi=300)


# ---- Figure 3: the Euler count, and why maximisers are degenerate --------------------------
def fig_euler():
    F0 = lambda X, Y: (X**2 - .36)*(Y**2 - .36)
    fig, axs = plt.subplots(1, 3, figsize=(6.5, 2.75))
    titles = [r"$P=(X^2-\frac{9}{25})(Y^2-\frac{9}{25})$:  $\nu=9$", r"$P+\delta$:  $\nu=5$", r"$P-\delta$:  $\nu=6$"]
    for ax, e, want, ttl in zip(axs, (0, .012, -.012), (9, 5, 6), titles):
        G = lambda X, Y, e=e: F0(X, Y) + e
        sign_panel(ax, G, "box")
        assert count_domains(G) == want
        ax.set_title(ttl, fontsize=9.5)
    ends = [(a*.6, b) for a in (1, -1) for b in (1, -1)] + [(b, a*.6) for a in (1, -1) for b in (1, -1)]
    cr = [(a*.6, b*.6) for a in (1, -1) for b in (1, -1)]
    ax = axs[0]
    ax.scatter(*zip(*ends), s=34, color=ORANGE, edgecolors="white", linewidths=0.8, zorder=4, clip_on=False,
               label=r"ends: $\varepsilon=8$")
    ax.scatter(*zip(*cr), s=34, marker="s", color=BLUE, edgecolors="white", linewidths=0.8, zorder=4,
               label=r"nodes: $h_p=4$")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=2, frameon=False, fontsize=8.5,
              handletextpad=0.1, columnspacing=0.9)
    fig.subplots_adjust(left=0.015, right=0.985, top=0.9, bottom=0.12, wspace=0.08)
    fig.savefig(OUT + "fig_euler.pdf", dpi=300)


# ---- Figure 4: the corner-sign argument -----------------------------------------------------
def fig_corner():
    fig, axs = plt.subplots(1, 3, figsize=(6.5, 2.5), gridspec_kw={"width_ratios": [1, 1, 0.95]})
    t = np.linspace(-1.75, 1.75, 500)
    cubics = [(lambda t: (t + .7)*t*(t - .65), [-.7, 0, .65], "signs agree: three roots\nfit in $(-1,1)$"),
              (lambda t: (t + .6)*(t - .15)*(t - 1.35), [-.6, .15, 1.35],
               "$g(1)<0<g(+\\infty)$: one root\nis pushed outside")]
    for ax, (g, roots, note) in zip(axs[:2], cubics):
        ax.axvspan(-1, 1, color="#eef3fb", lw=0)
        ax.axhline(0, color=MUTED, lw=0.7)
        ax.plot(t, g(t), color=BLUE, lw=1.8)
        inside = [r for r in roots if abs(r) < 1]; outside = [r for r in roots if abs(r) > 1]
        ax.scatter(inside, [0]*len(inside), s=30, color=INK, zorder=4)
        if outside:
            ax.scatter(outside, [0]*len(outside), s=34, facecolors="white", edgecolors=ORANGE, linewidths=1.5, zorder=4)
        for e in (-1, 1):
            ax.scatter([e], [g(e)], s=30, color=ORANGE, edgecolors="white", linewidths=0.8, zorder=5)
            ax.annotate(r"$+$" if g(e) > 0 else r"$-$", (e, g(e)), xytext=(7*e, 6 if g(e) > 0 else -9),
                        textcoords="offset points", ha="center", va="center", fontsize=11, color="#b84a1c")
        ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.5, 1.5)
        ax.set_xticks([-1, 1]); ax.set_yticks([]); ax.tick_params(labelsize=8, length=2)
        ax.set_xlabel(note, fontsize=8.5, labelpad=3)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    axs[0].set_title(r"a cubic side restriction $g$", fontsize=9.5)
    axs[1].set_title(r"a sign mismatch at $t=1$", fontsize=9.5)
    # the square with the corner signs wanted by each side (odd degree, shared leading signs)
    ax = axs[2]
    ax.add_patch(Rectangle((-1, -1), 2, 2, fc="#f7f6f3", ec=INK, lw=1.1))
    want = {(1, 1): ("+", "+"), (1, -1): ("-", "+"), (-1, 1): ("+", "-"), (-1, -1): ("-", "-")}
    for (cx, cy), (v, h) in want.items():
        bad = v != h
        ax.text(cx*1.19, cy*0.80, rf"${v}$", ha="center", va="center", fontsize=12, color=BLUE)
        ax.text(cx*0.80, cy*1.19, rf"${h}$", ha="center", va="center", fontsize=12, color="#b84a1c")
        ax.scatter([cx], [cy], s=60 if bad else 26, zorder=4, clip_on=False,
                   facecolors="white" if bad else INK, edgecolors=INK, linewidths=1.3)
    ax.text(1.19, 0, "sides\n$X=\\pm1$", ha="left", va="center", fontsize=8, color=BLUE)
    ax.text(0, 1.22, "sides $Y=\\pm1$", ha="center", va="bottom", fontsize=8, color="#b84a1c")
    ax.text(0, 0, "two corners\ncannot satisfy\nboth sides", ha="center", va="center", fontsize=8)
    ax.set_xlim(-1.45, 1.95); ax.set_ylim(-1.45, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(r"odd degree: $12\to10$ ends", fontsize=9.5)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.2, wspace=0.12)
    fig.savefig(OUT + "fig_corner.pdf")


# ---- Figure 5: from the curve to lattice polygons -------------------------------------------
def lattice_axes(ax, lam, size):
    S = L.support(lam)
    ax.add_patch(Polygon(L.hull(S), closed=True, fill=False, ec=MUTED, lw=1.0, ls=(0, (3, 2))))
    ax.set_xlim(-0.5, size + 0.9); ax.set_ylim(-0.5, size + 0.9); ax.set_aspect("equal")
    ax.set_xticks(range(size + 1)); ax.set_yticks(range(size + 1)); ax.tick_params(labelsize=7, length=2)
    tidy(ax, None)
    return S


def fig_polygons():
    rows = [
        (17, EX[2][1], [([(0, 0), (2, 0), (0, 2)], BLUE, r"$D_1$ (conic)", (1.05, 1.05)),
                        ([(0, 0), (1, 0), (0, 1)], ORANGE, r"$D_2$ (line)", (0.25, 0.25))],
         [(0, 0), (3, 0), (0, 3)], r"missed $0$ + chain $1$ + ends $1$ $=2$"),
        (18, EX[3][1], [([(0, 0), (1, 0)], BLUE, r"$D_1=D_2$ (lines $X=\pm\frac{3}{5}$)", (0.5, -0.02)),
                        ([(0, 0), (0, 1)], ORANGE, r"$D_3=D_4$ (lines $Y=\pm\frac{3}{5}$)", (0.05, 0.75))],
         [(0, 0), (2, 0), (2, 2), (0, 2)], r"missed $2$ + chain $0$ + ends $0$ $=2$"),
    ]
    fig, axs = plt.subplots(2, 3, figsize=(6.5, 4.5))
    for (lam, F, facs, total, acct), row in zip(rows, axs):
        sign_panel(row[0], F, "box")
        row[0].set_ylabel(rf"$\lambda={lam}$", fontsize=10.5, rotation=0, labelpad=22, va="center")
        ax = row[1]; S = lattice_axes(ax, lam, 3)
        xs, ys = zip(*S); ax.scatter(xs, ys, s=12, color=MUTED, zorder=3)
        for P, col, lab, (tx, ty) in facs:
            if len(P) == 2:
                ax.plot(*zip(*P), color=col, lw=3.2, solid_capstyle="round", zorder=4)
                ax.plot([], [], color=col, lw=3.2, label=lab)
            else:
                ax.add_patch(Polygon(P, closed=True, fc=col, alpha=0.22, ec="none", zorder=2))
                ax.add_patch(Polygon(P, closed=True, fill=False, ec=col, lw=1.6, zorder=4))
                ax.plot([], [], color=col, lw=3.2, label=lab)
        ax.legend(loc="upper right", frameon=False, fontsize=7.5, handlelength=1.2, borderaxespad=0.1)
        ax = row[2]; S = lattice_axes(ax, lam, 3)
        ax.add_patch(Polygon(total, closed=True, fc=POS, ec=BLUE, lw=1.5, zorder=2))
        Dpts = set(L.lattice_points(L.hull(total))) if hasattr(L, "lattice_points") else None
        if Dpts is None:
            xm = max(p[0] for p in total); ym = max(p[1] for p in total)
            if len(total) == 3:
                Dpts = {(a, b) for a in range(xm + 1) for b in range(ym + 1) if a + b <= xm}
            else:
                Dpts = {(a, b) for a in range(xm + 1) for b in range(ym + 1)}
        hit = [p for p in S if p in Dpts]; miss = [p for p in S if p not in Dpts]
        ax.scatter(*zip(*hit), s=14, color=INK, zorder=4)
        if miss:
            ax.scatter(*zip(*miss), s=40, facecolors="white", edgecolors=ORANGE, linewidths=1.6, zorder=4)
        ax.set_xlabel(acct, fontsize=8)
    axs[0, 0].set_title("nodal set in $B$", fontsize=9.5)
    axs[0, 1].set_title("polygons of the factors", fontsize=9.5)
    axs[0, 2].set_title(r"their sum $D\subseteq\Delta_\lambda$", fontsize=9.5)
    fig.subplots_adjust(left=0.09, right=0.97, top=0.94, bottom=0.09, wspace=0.25, hspace=0.3)
    fig.savefig(OUT + "fig_polygons.pdf", dpi=300)


# ---- Figure 6: the hull and the circle ------------------------------------------------------
def hull_panel(ax, lam, ms):
    R = math.sqrt(lam)
    Sp = [(a + 1, b + 1) for a, b in L.support(lam)]
    ch = [(x + 1, y + 1) for x, y in Hu.upper_chain(lam)]
    th = np.linspace(0, np.pi/2, 400)
    ax.plot(R*np.cos(th), R*np.sin(th), color=MUTED, lw=0.8)
    xs, ys = zip(*Sp); ax.scatter(xs, ys, s=ms, color=INK, linewidths=0, zorder=2)
    nH = 0
    for (x1, y1), (x2, y2) in zip(ch, ch[1:]):
        slanted = x1 != x2 and y1 != y2
        ax.plot([x1, x2], [y1, y2], color=BLUE if slanted else MUTED, lw=1.8 if slanted else 1.2, zorder=3,
                solid_capstyle="round")
        if slanted:
            g = math.gcd(x2 - x1, y1 - y2); nH += g
            px = [x1 + (x2 - x1)*i/g for i in range(g + 1)]; py = [y1 + (y2 - y1)*i/g for i in range(g + 1)]
            ax.scatter(px, py, s=ms*2.6, color=BLUE, zorder=4, linewidths=0)
    arc = [(m, k) for m, k in Sp if m*m + k*k == lam]
    ax.scatter(*zip(*arc), s=ms*9, facecolors="none", edgecolors=ORANGE, linewidths=1.4, zorder=5)
    assert nH == Hu.H(lam) and len(arc) == L.mult(lam)
    ax.set_aspect("equal"); ax.set_xlim(0, R*1.04); ax.set_ylim(0, R*1.04)
    ax.set_xlabel(r"$m$", labelpad=1); ax.set_ylabel(r"$k$", rotation=0, labelpad=6)
    ax.tick_params(labelsize=7.5, length=2); tidy(ax, None)
    ax.set_title(rf"$\lambda={lam}$:  $r={len(arc)}$,  $H={nH}$", fontsize=9.5)


def fig_hull():
    fig = plt.figure(figsize=(6.5, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.25, 1], hspace=0.42, wspace=0.22,
                          left=0.085, right=0.985, top=0.95, bottom=0.085)
    hull_panel(fig.add_subplot(gs[0, 0]), 130, 9)
    hull_panel(fig.add_subplot(gs[0, 1]), 1105, 2.2)
    ax = fig.add_subplot(gs[1, :])
    lam = 32045; R = math.sqrt(lam)
    ch = [(x + 1, y + 1) for x, y in Hu.upper_chain(lam)]
    pts = []
    for (x1, y1), (x2, y2) in zip(ch, ch[1:]):
        g = math.gcd(x2 - x1, abs(y1 - y2))
        pts += [(x1 + (x2 - x1)*i/g, y1 + (y2 - y1)*i/g) for i in range(g)]
    pts.append(ch[-1])
    ang = [math.degrees(math.atan2(b, a)) for a, b in pts]; dep = [R - math.hypot(a, b) for a, b in pts]
    ax.plot(ang, dep, color=BLUE, lw=1.0)
    ax.scatter(ang, dep, s=5, color=BLUE, zorder=3, linewidths=0)
    arc = [(a, b) for a, b in ch if a*a + b*b == lam]
    ax.scatter([math.degrees(math.atan2(b, a)) for a, b in arc], [0]*len(arc), s=30, facecolors="white",
               edgecolors=ORANGE, linewidths=1.4, zorder=4, clip_on=False)
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.set_xlim(0, 90); ax.set_ylim(max(dep)*1.08, -0.004)
    ax.set_xticks(range(0, 91, 15)); ax.tick_params(labelsize=7.5, length=2)
    ax.set_xlabel("angle (degrees)", labelpad=1); ax.set_ylabel("depth below circle", fontsize=8.5)
    ax.set_title(rf"$\lambda={lam}$, unrolled:  $r={len(arc)}$,  $H={Hu.H(lam)}$", fontsize=9.5)
    tidy(ax, None)
    fig.savefig(OUT + "fig_hull.pdf")


# ---- Figure 7: supply and demand ------------------------------------------------------------
def fig_supply():
    w = np.loadtxt("../../data/H_window_minima_1e9.txt")
    mid = np.sqrt(w[:, 0]*w[:, 1]); loH = w[:, 2]
    rec = json.load(open("../../data/rmax_records.json"))
    X = np.logspace(0.3, 13, 900)
    rmax = np.array([max(d for n, d in rec if n <= x) for x in X])
    fig, ax = plt.subplots(figsize=(6.3, 3.5))
    ax.axvspan(1, 1e9, color="#eef3fb", lw=0)
    ax.plot(X, np.maximum(2*rmax - 3, 0.5), color=ORANGE, lw=1.8, drawstyle="steps-post")
    ax.plot(mid, loH, "o", ms=3.4, color=BLUE, mec="white", mew=0.4)
    bb = 0.33/4*X**(1/3) - 4
    ax.plot(X[bb > 0.8], bb[bb > 0.8], "--", color=BLUE, lw=1.2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(2, 1e13); ax.set_ylim(0.8, 2e4)
    ax.text(2.5e5, 600, "smallest $H(\\lambda)$ in each\nhalf-octave of $\\lambda$", color=INK, fontsize=8.5, ha="right")
    ax.text(2e9, 4000, r"$2\,\max_{\mu\leq\lambda}r(\mu)-3$", color=INK, fontsize=8.5, ha="left")
    ax.text(2.2e6, 2.6, "$\\frac{0.33}{4}\\lambda^{1/3}-4$  (Balog–Bárány, for large $\\lambda$)", color=INK, fontsize=8.5, ha="left", va="center")
    ax.text(4, 6e3, "checked at every eigenvalue", color="#1c5cab", fontsize=8.5, ha="left")
    ax.set_xlabel(r"eigenvalue $\lambda$"); ax.set_ylabel("number of segments")
    tidy(ax, "both")
    fig.tight_layout(pad=0.4); fig.savefig(OUT + "fig_supply.pdf")


if __name__ == "__main__":
    which = sys.argv[1:] or ["gallery", "reduction", "euler", "corner", "polygons", "hull", "supply"]
    for name in which:
        globals()["fig_" + name](); print("wrote fig_%s.pdf" % name)
