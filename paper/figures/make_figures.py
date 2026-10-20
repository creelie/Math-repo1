#!/usr/bin/env python3
"""Figures of the paper (PNG, white background, 256 colours, at most 1600 px wide).

usage: python3 paper/figures/make_figures.py [RESULTS_DIR]

fig_depth.png   improper lifts by binary level for the rows of Table 2
fig_shape.png   the function H on a two-parameter slice, against the level 900
fig_shift.png   the shifts forbidden by the exceptional speeds (shift lemma)
fig_budget.png  log p for the primes 17..523, verified primes marked
"""
import glob, json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "results")
plt.rcParams.update({"font.family": "serif", "font.size": 11, "mathtext.fontset": "cm"})


def save(fig, name):
    path = os.path.join(HERE, name)
    fig.savefig(path, dpi=260, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    im = Image.open(path).convert("RGB")
    gray = np.asarray(im.convert("L"))
    ys, xs = np.nonzero(gray < 250)
    m = 24
    im = im.crop((max(xs.min() - m, 0), max(ys.min() - m, 0),
                  min(xs.max() + m, im.width), min(ys.max() + m, im.height)))
    if im.width > 1600:
        im = im.resize((1600, round(im.height * 1600 / im.width)), Image.LANCZOS)
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(path, optimize=True)
    print(name, im.size)


# ---------------------------------------------------------------- fig_depth
rows = [
    (r"$p=19$: $1^7,2,\dots,9$", [3984537312, 11836500480, 21372422400, 0]),
    (r"$p=23$: $1^5,2,\dots,11$", [21976320, 593280, 0, 0]),
    (r"$p=23$: $1^3,2^2,3^2,4,\dots,11$", [25459728, 3335616, 50688, 0]),
    (r"$p=37$: $1^7,3,4,7,\dots,16$", [29756160, 82615680, 0, 0]),
    (r"$p=41$: $1^6,3,4,7,\dots,18$", [7687200, 11865600, 0, 0]),
][::-1]
levels = [16, 32, 64, 128]
fig = plt.figure(figsize=(8.0, 5.6))
ax = fig.add_subplot(projection="3d")
colors = cm.viridis(np.linspace(0.12, 0.85, len(rows)))
for r, (label, counts) in enumerate(rows):
    for j, c in enumerate(counts):
        h = math.log10(c) if c else 0.0
        if c:
            ax.bar3d(j - 0.3, r - 0.3, 0, 0.6, 0.6, h, color=colors[r], alpha=0.95, shade=True,
                     edgecolor="k", linewidth=0.3)
        else:
            ax.bar3d(j - 0.3, r - 0.3, 0, 0.6, 0.6, 0.12, color="#d9534f", alpha=0.9, shade=True,
                     edgecolor="k", linewidth=0.3)
ax.set_xticks(range(len(levels)))
ax.set_xticklabels([str(l) for l in levels])
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([lab for lab, _ in rows], fontsize=8.5, ha="left")
ax.set_xlabel("binary level", labelpad=6)
ax.set_zlabel(r"$\log_{10}$ (improper lifts)", labelpad=4)
ax.set_zlim(0, 11)
ax.tick_params(axis="y", pad=-2)
ax.view_init(elev=24, azim=-58)
ax.text2D(0.02, 0.93, "red: no improper lift left", transform=ax.transAxes, fontsize=9, color="#b52b27")
save(fig, "fig_depth.png")

# ---------------------------------------------------------------- fig_shape
def q(u): return 1 + 2 * u - u * u
def a(u): return q(u) / u ** (2 / 15)
def w(u): return u * u / q(u)
s = np.linspace(0.035, 0.15, 90)
t = np.linspace(0.035, 0.15, 90)
S, T = np.meshgrid(s, t)
Hs = (0.5 + 7 * w(S) + 7 * w(T)) * 2 * a(S) ** 7 * a(T) ** 7
fig = plt.figure(figsize=(7.6, 5.8))
ax = fig.add_subplot(projection="3d")
ax.plot_surface(S, T, Hs, cmap=cm.coolwarm, linewidth=0, antialiased=True, alpha=0.95, rstride=1, cstride=1)
ax.plot_surface(S, T, np.full_like(S, 900.0), color="#607d8b", alpha=0.35, linewidth=0)
t0 = 0.07245
H0 = (1 + 2 * t0 + 27 * t0 * t0) * q(t0) ** 13 / t0 ** (28 / 15)
ax.plot([t0, t0], [t0, t0], [900, H0], color="k", linewidth=1.2)
ax.scatter([t0], [t0], [H0], color="k", s=22, depthshade=False)
ax.text2D(0.03, 0.90, r"one coordinate $1$, seven equal to $s$, seven equal to $t$", transform=ax.transAxes, fontsize=10)
ax.text2D(0.03, 0.84, r"minimum at $s=t\approx0.0725$, $H\approx%.0f$; grey plane $H=900$" % H0,
          transform=ax.transAxes, fontsize=10)
ax.set_xlabel(r"$s$", labelpad=6)
ax.set_ylabel(r"$t$", labelpad=6)
ax.set_zlabel(r"$H$", labelpad=4)
ax.set_zlim(880, 1700)
ax.view_init(elev=30, azim=-55)
save(fig, "fig_shape.png")

# ---------------------------------------------------------------- fig_shift
def dist(x):
    x -= math.floor(x)
    return min(x, 1 - x)
d, t0s = 16, 0.013
speeds = [1, 3, 2, 4, 6]
fig = plt.figure(figsize=(9.0, 6.6))
ax = fig.add_subplot(projection="3d")
theta = np.linspace(0, 2 * np.pi, 200)
good = [qq for qq in range(d) if all(dist((t0s + qq / d) * u) >= 1 / 16 for u in speeds)]
for level, u in enumerate(speeds):
    z = level * 1.0
    ax.plot(np.cos(theta), np.sin(theta), z, color="#888888", linewidth=0.8)
    g = math.gcd(d, u)
    bad = [qq for qq in range(d) if dist((t0s + qq / d) * u) < 1 / 16]
    for qq in range(d):
        ang = 2 * np.pi * qq / d
        col = "#d9534f" if qq in bad else "#2e7d32"
        ax.scatter(np.cos(ang), np.sin(ang), z, color=col, s=30 if qq in bad else 14, depthshade=False)
    ax.text(1.35, -0.6, z, r"$u\equiv%d$ (mod 16): forbids %d $\leq %d$" % (
        u, len(bad), g * math.ceil(2 * (d // g) / 16)), fontsize=10)
qq = good[0]
ang = 2 * np.pi * qq / d
ax.plot([np.cos(ang)] * 2, [np.sin(ang)] * 2, [0, len(speeds) - 1], color="#1565c0", linewidth=2.2)
ax.text(np.cos(ang) - 0.2, np.sin(ang), len(speeds) - 0.6, r"shift $q=%d$ is good for every speed" % qq,
        fontsize=10, color="#1565c0")
ax.set_xlim(-1.1, 2.6)
ax.set_ylim(-1.1, 1.1)
ax.set_axis_off()
ax.view_init(elev=16, azim=-80)
ax.set_box_aspect((3.4, 2.0, 3.0))
save(fig, "fig_shift.png")

# ---------------------------------------------------------------- fig_budget
def primes(lo, hi):
    return [n for n in range(lo, hi + 1) if n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))]
closed = set()
for f in glob.glob(os.path.join(RESULTS, "summary_*.json")):
    sm = json.load(open(f))
    if sm.get("status") == "GATE_CLOSED":
        closed.add(int(sm["gate"]["p"]))
ps = primes(17, 523)
fig = plt.figure(figsize=(8.4, 4.8))
ax = fig.add_subplot(projection="3d")
for i, p in enumerate(ps):
    done = p in closed
    ax.bar3d(i, 0, 0, 0.8, 0.8, math.log(p), color="#1565c0" if done else "#cfd8dc",
             edgecolor="k" if done else "#90a4ae", linewidth=0.2, alpha=1.0 if done else 0.85, shade=True)
ticks = [i for i, p in enumerate(ps) if p in (17, 101, 179, 239, 307, 401, 523)]
ax.set_xticks(ticks)
ax.set_xticklabels([str(ps[i]) for i in ticks], fontsize=9)
ax.set_yticks([])
ax.set_zticks([0, 2, 4, 6])
ax.set_zlabel(r"$\log p$", labelpad=0)
ax.set_xlabel("prime $p$", labelpad=16)
ax.view_init(elev=26, azim=-72)
ax.set_box_aspect((3.6, 0.7, 1.1))
tot = sum(math.log(p) for p in closed if p >= 17)
ax.text2D(0.02, 0.80, r"blue: verified primes, $\sum\log p=%.2f$; needed $483.54$" % tot,
          transform=ax.transAxes, fontsize=10)
save(fig, "fig_budget.png")
