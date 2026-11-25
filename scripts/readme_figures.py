"""Figures in the README. Usage: python scripts/readme_figures.py  (writes docs/figures/*.png)"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage

repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
out = os.path.join(repo, "docs", "figures")
os.makedirs(out, exist_ok=True)
sys.path.insert(0, os.path.join(repo, "python"))
from squareecp import nodal, hull_units, lattice as L

THEMES = {
    "light": dict(bg="#ffffff", ink="#1f2328", muted="#59636e", grid="#d1d9e0",
                  pos="#2a78d6", neg="#eb6834", pos_soft="#c8dcf4", neg_soft="#f9d3c3", accent="#e34948"),
    "dark":  dict(bg="#0d1117", ink="#e6edf3", muted="#9198a1", grid="#3d444d",
                  pos="#3987e5", neg="#d95926", pos_soft="#1d3557", neg_soft="#5a2a16", accent="#e66767"),
}

def style(t):
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 11, "text.color": t["ink"],
                         "axes.labelcolor": t["ink"], "xtick.color": t["muted"], "ytick.color": t["muted"],
                         "axes.edgecolor": t["grid"], "svg.fonttype": "none"})

# ---------- 1. nodal domains of extremal combinations ----------
PICKS = [(18, 9, "U(18) = 9"), (20, 12, "κ(20) = 12"), (34, 21, "κ(34) = 21")]
def nodal_fig(t, name):
    style(t)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.7), facecolor=t["bg"])
    n = 900
    g = np.linspace(0, np.pi, n)
    X, Y = np.meshgrid(g, g, indexing="xy")
    for ax, (lam, count, lab) in zip(axes, PICKS):
        expr = dict((l, e) for l, e, _ in nodal.EXAMPLES)[lam]
        P = eval("lambda X, Y: " + expr)
        F = np.sin(X) * np.sin(Y) * P(np.cos(X), np.cos(Y))
        # check the count on the polynomial grid
        assert nodal.stable_count(eval("lambda X, Y: " + expr)) == count, lam
        img = np.zeros(F.shape + (3,))
        from matplotlib.colors import to_rgb
        img[F > 0] = to_rgb(t["pos_soft"]); img[F <= 0] = to_rgb(t["neg_soft"])
        ax.imshow(img, origin="lower", extent=(0, np.pi, 0, np.pi), interpolation="bilinear")
        ax.contour(X, Y, F, levels=[0], colors=[t["ink"]], linewidths=1.3)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_color(t["muted"]); s.set_linewidth(1)
        ax.set_title(f"λ = {lam}:  {count} nodal domains", fontsize=12, pad=8, color=t["ink"])
        ax.set_facecolor(t["bg"])
    fig.tight_layout(w_pad=2.2)
    fig.savefig(os.path.join(out, name[:-4]+".png"), facecolor=t["bg"], dpi=200)
    plt.close(fig)

# ---------- 2. the hull and H(λ) ----------
def hull_fig(t, name, lam=130):
    style(t)
    fig, ax = plt.subplots(figsize=(5.2, 5.6), facecolor=t["bg"])
    ax.set_facecolor(t["bg"])
    R = np.sqrt(lam)
    pts = [(m, k) for m in range(1, int(R) + 1) for k in range(1, int(R) + 1) if m * m + k * k <= lam]
    on = [(m, k) for m, k in pts if m * m + k * k == lam]
    th = np.linspace(0, np.pi / 2, 300)
    ax.plot(R * np.cos(th), R * np.sin(th), color=t["muted"], lw=1, ls=(0, (4, 3)))
    ax.scatter(*zip(*pts), s=14, color=t["muted"], zorder=3, linewidths=0)
    # hull chain in (m, k) coordinates = exponent chain shifted by (1, 1)
    ch = [(a + 1, b + 1) for a, b in hull_units.upper_chain(lam)]
    ax.plot(*zip(*ch), color=t["ink"], lw=1.4, zorder=2)
    H = 0
    for (x1, y1), (x2, y2) in zip(ch, ch[1:]):
        if x1 != x2 and y1 != y2:
            from math import gcd
            d = gcd(abs(x2 - x1), abs(y2 - y1))
            H += d
            ax.plot([x1, x2], [y1, y2], color=t["pos"], lw=3.2, zorder=4, solid_capstyle="round")
            for j in range(d + 1):
                ax.scatter(x1 + (x2 - x1) * j / d, y1 + (y2 - y1) * j / d, s=34, color=t["pos"],
                           zorder=5, edgecolors=t["bg"], linewidths=1.2)
    ax.scatter(*zip(*on), s=90, facecolors="none", edgecolors=t["neg"], linewidths=1.8, zorder=6)
    r = L.mult(lam)
    assert H == hull_units.H(lam)
    ax.set_xlim(0, R + 0.8); ax.set_ylim(0, R + 0.8); ax.set_aspect("equal")
    ax.set_xlabel("m"); ax.set_ylabel("k")
    ax.set_xticks(range(0, int(R) + 1, 2)); ax.set_yticks(range(0, int(R) + 1, 2))
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(f"λ = {lam}:  H = {H} slanted units,  r = {r}", fontsize=12, color=t["ink"])
    ax.text(0.5, -0.17, "blue: slanted hull units (H)     orange rings: m² + k² = λ (r)",
            transform=ax.transAxes, ha="center", color=t["muted"], fontsize=9.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, name[:-4]+".png"), facecolor=t["bg"], dpi=200)
    plt.close(fig)

# ---------- 3. slack ⌈H/2⌉ − (r − 1) ----------
def slack_data(X=200000):
    lams = L.eigenvalues(X)
    return np.array(lams), np.array([hull_units.slack(l) for l in lams])

def slack_fig(t, name, lams, sl):
    style(t)
    fig, ax = plt.subplots(figsize=(8.2, 3.9), facecolor=t["bg"])
    ax.set_facecolor(t["bg"])
    ax.scatter(lams, sl, s=2, color=t["pos"], alpha=0.35, linewidths=0, rasterized=True)
    z = sl == 0
    ax.scatter(lams[z], sl[z], s=26, color=t["neg"], zorder=4, edgecolors=t["bg"], linewidths=0.8)
    ax.axhline(0, color=t["ink"], lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("eigenvalue λ")
    ax.set_ylabel(r"$\lceil H/2\rceil - (r-1)$")
    ax.grid(axis="y", color=t["grid"], lw=0.6)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(f"The criterion of Theorem C at all {len(lams):,} eigenvalues λ ≤ {lams[-1]:,}",
                 fontsize=12, color=t["ink"], loc="left")
    ax.text(0.02, 0.86, f"never below 0; equal to 0 at {z.sum()} levels (orange), the last at λ = {lams[z].max()}",
            transform=ax.transAxes, ha="left", color=t["muted"], fontsize=9.5)
    assert (sl >= 0).all()
    fig.tight_layout()
    fig.savefig(os.path.join(out, name[:-4]+".png"), facecolor=t["bg"], dpi=200)
    plt.close(fig)

lams, sl = slack_data()
for mode, t in THEMES.items():
    nodal_fig(t, f"nodal-{mode}.svg")
    hull_fig(t, f"hull-{mode}.svg")
    slack_fig(t, f"slack-{mode}.svg", lams, sl)

