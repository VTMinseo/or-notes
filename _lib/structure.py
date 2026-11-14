"""Schematic drawings of the block structure of a constraint matrix.

    draw_structure(ax, "benders")  -> block-diagonal + linking VARIABLES (a column)
    draw_structure(ax, "dw")       -> block-diagonal + linking CONSTRAINTS (a row)

Coloring follows the site palette: linking part orange (COLORS["p2"]),
diagonal blocks blue (COLORS["p1"]), objective row grey.
"""
from __future__ import annotations

import matplotlib.patches as mpatches

from _lib.orplot import COLORS

GREY = "#9a9a9a"


def _cell(ax, x, y, w, h, color, label, alpha=0.75, fontsize=12, text_color="white"):
    ax.add_patch(mpatches.FancyBboxPatch((x + 0.04, y + 0.04), w - 0.08, h - 0.08,
                                         boxstyle="round,pad=0,rounding_size=0.06",
                                         fc=color, ec=color, alpha=alpha, lw=1.2))
    if label:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=fontsize, color=text_color, fontweight="bold")


def draw_structure(ax, kind: str = "benders", K: int = 3, show_rhs: bool = True):
    """Draw the constraint-matrix structure exploited by Benders or Dantzig–Wolfe."""
    blue, orange = COLORS["p1"], COLORS["p2"]
    bw, bh = 1.6, 1.15          # size of one diagonal block
    obj_h = 0.55
    sub = "₁₂₃₄₅₆₇₈₉"

    if kind == "benders":
        lw = 1.2                                 # width of the linking column x
        width = lw + K * bw
        height = obj_h + K * bh
        top = height
        # objective row
        _cell(ax, 0, top - obj_h, lw, obj_h, GREY, "$c^\\top$", fontsize=11)
        for k in range(K):
            _cell(ax, lw + k * bw, top - obj_h, bw, obj_h, GREY, f"$d_{k + 1}^\\top$", fontsize=11)
        # rows
        for k in range(K):
            y = top - obj_h - (k + 1) * bh
            _cell(ax, 0, y, lw, bh, orange, f"$A_{k + 1}$")
            _cell(ax, lw + k * bw, y, bw, bh, blue, f"$B_{k + 1}$")
            if show_rhs:
                ax.text(width + 0.25, y + bh / 2, f"$\\geq\\; b_{k + 1}$", va="center", fontsize=11)
        # column headers
        ax.text(lw / 2, top + 0.22, "$x$", ha="center", fontsize=13, color=orange, fontweight="bold")
        for k in range(K):
            ax.text(lw + k * bw + bw / 2, top + 0.22, f"$y_{k + 1}$", ha="center", fontsize=13)
        # bracket for the linking column
        ax.annotate("", xy=(0.05, -0.25), xytext=(lw - 0.05, -0.25),
                    arrowprops=dict(arrowstyle="-", color=orange, lw=2))
        ax.text(lw / 2, -0.5, "linking\nvariables", ha="center", va="top",
                fontsize=10, color=orange)
        ax.text(lw + K * bw / 2, -0.5, "fix $x$ $\\Rightarrow$ the $y_k$ separate into\n"
                f"{K} independent LPs (subproblems)", ha="center", va="top", fontsize=10)
    elif kind == "dw":
        width = K * bw
        link_h = 0.8
        height = obj_h + link_h + K * bh
        top = height
        for k in range(K):
            _cell(ax, k * bw, top - obj_h, bw, obj_h, GREY, f"$c_{k + 1}^\\top$", fontsize=11)
            _cell(ax, k * bw, top - obj_h - link_h, bw, link_h, orange, f"$D_{k + 1}$")
        if show_rhs:
            ax.text(width + 0.25, top - obj_h - link_h / 2, "$\\leq\\; d$", va="center", fontsize=11)
        for k in range(K):
            y = top - obj_h - link_h - (k + 1) * bh
            _cell(ax, k * bw, y, bw, bh, blue, f"$B_{k + 1}$")
            if show_rhs:
                ax.text(width + 0.25, y + bh / 2, f"$\\leq\\; b_{k + 1}$", va="center", fontsize=11)
        for k in range(K):
            ax.text(k * bw + bw / 2, top + 0.22, f"$x_{k + 1}$", ha="center", fontsize=13)
        # bracket for the linking rows (left side)
        y0, y1 = top - obj_h - link_h + 0.05, top - obj_h - 0.05
        ax.annotate("", xy=(-0.25, y0), xytext=(-0.25, y1),
                    arrowprops=dict(arrowstyle="-", color=orange, lw=2))
        ax.text(-0.4, (y0 + y1) / 2, "linking\nconstraints", ha="right", va="center",
                fontsize=10, color=orange)
        ax.text(width / 2, -0.25, "price out the linking rows $\\Rightarrow$ the $x_k$ separate\n"
                f"into {K} independent subproblems (pricing)", ha="center", va="top", fontsize=10)
    else:
        raise ValueError(kind)

    # outline of the whole matrix; "0" in two empty off-diagonal blocks
    ax.add_patch(mpatches.Rectangle((0, 0), width, height, fill=False, ec="#666666", lw=1.0))
    x0 = lw if kind == "benders" else 0.0
    for (col, row) in [(0, K - 1), (K - 1, 0)]:          # bottom-left, top-right
        ax.text(x0 + col * bw + bw / 2, (K - 1 - row) * bh + bh / 2, "0", ha="center",
                va="center", fontsize=15, color="#aaaaaa")
    ax.set_xlim(-2.0 if kind == "dw" else -0.4, width + 1.3)
    ax.set_ylim(-1.3, height + 0.6)
    ax.set_aspect("equal")
    ax.axis("off")
