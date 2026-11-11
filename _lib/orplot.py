"""Shared plotting helpers so every figure on the site looks the same.

Usage inside a .qmd code cell:

    from _lib.orplot import new_axes, shade_polygon, COLORS
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

# Colorblind-safe palette, used consistently across the site:
#   disjuncts / pieces -> blue, orange
#   relaxations        -> red (weak), green (strong / hull)
COLORS = {
    "p1": "#4C72B0",
    "p2": "#DD8452",
    "weak": "#C44E52",
    "hull": "#55A868",
    "grid": "#BBBBBB",
    "text": "#333333",
}

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 160,
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 12,
    "legend.frameon": False,
})


def new_axes(ncols: int = 1, xlim=(0, 10), ylim=(0, 10), size: float = 4.2):
    """Square 2-D axes with a light grid, labelled x1 / x2."""
    fig, axes = plt.subplots(1, ncols, figsize=(size * ncols, size), squeeze=False)
    axes = axes[0]
    for ax in axes:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.grid(True, color=COLORS["grid"], lw=0.5, alpha=0.5)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
    return fig, (axes if ncols > 1 else axes[0])


def shade_polygon(ax, vertices, color, label=None, alpha=0.35, hatch=None, lw=1.6):
    """Fill a polygon given its vertices in order."""
    v = np.asarray(vertices, dtype=float)
    ax.fill(v[:, 0], v[:, 1], color=color, alpha=alpha, label=label, hatch=hatch, lw=0)
    closed = np.vstack([v, v[:1]])
    ax.plot(closed[:, 0], closed[:, 1], color=color, lw=lw)


def convex_hull(points):
    """Vertices of the 2-D convex hull, in counter-clockwise order (monotone chain)."""
    pts = sorted(set(map(tuple, np.asarray(points, dtype=float))))
    if len(pts) <= 2:
        return np.array(pts)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def mark_point(ax, xy, color, label=None, marker="o", size=60):
    ax.scatter([xy[0]], [xy[1]], color=color, s=size, marker=marker,
               zorder=5, edgecolor="white", linewidth=1.2, label=label, clip_on=False)
