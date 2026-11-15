"""Shared toy network for the Network Optimization pages, so every page shows the same numbers.

Six nodes s, a, b, c, d, t and nine directed arcs. Each arc has a length/cost
(for shortest paths and min-cost flow) and a capacity (for max flow).

    shortest s -> t path:  s - a - d - c - t, length 8
    maximum s -> t flow:   18, minimum cut S = {s, a, b} (arcs a->c, a->d, b->d)

Also: a small bipartite graph for the matching pages, and plotting helpers.
"""
from __future__ import annotations

from collections import deque
import heapq

import numpy as np
from scipy.optimize import linprog

NODES = ["s", "a", "b", "c", "d", "t"]
POS = {"s": (0.0, 1.0), "a": (1.0, 2.0), "b": (1.0, 0.0),
       "c": (2.2, 2.0), "d": (2.2, 0.0), "t": (3.2, 1.0)}
# (tail, head, cost/length, capacity)
ARCS = [("s", "a", 2, 12), ("s", "b", 4, 10), ("a", "b", 1, 3), ("a", "c", 7, 6),
        ("a", "d", 3, 4), ("b", "d", 3, 8), ("d", "c", 2, 5), ("c", "t", 1, 11),
        ("d", "t", 5, 9)]
COST = {(u, v): w for u, v, w, _ in ARCS}
CAP = {(u, v): k for u, v, _, k in ARCS}


def incidence_matrix(nodes=NODES, arcs=ARCS):
    """Node-arc incidence matrix: +1 where an arc leaves a node, -1 where it enters."""
    N = np.zeros((len(nodes), len(arcs)), dtype=int)
    idx = {v: i for i, v in enumerate(nodes)}
    for j, (u, v, *_) in enumerate(arcs):
        N[idx[u], j] = 1
        N[idx[v], j] = -1
    return N


# ---------------------------------------------------------------- shortest paths
def dijkstra_steps(source="s", nodes=NODES, arcs=ARCS):
    """Dijkstra's algorithm; returns one snapshot per node made permanent."""
    out = {v: [] for v in nodes}
    for u, v, w, _ in arcs:
        out[u].append((v, w))
    dist = {v: np.inf for v in nodes}
    pred = {v: None for v in nodes}
    dist[source] = 0
    done, heap, steps = [], [(0, source)], []
    while heap:
        d, u = heapq.heappop(heap)
        if u in done or d > dist[u]:
            continue
        done.append(u)
        for v, w in out[u]:
            if dist[u] + w < dist[v]:
                dist[v], pred[v] = dist[u] + w, u
                heapq.heappush(heap, (dist[v], v))
        steps.append(dict(permanent=list(done), dist=dict(dist), pred=dict(pred), current=u))
    return steps


def path_to(pred, target):
    path = [target]
    while pred[path[-1]] is not None:
        path.append(pred[path[-1]])
    return path[::-1]


def shortest_path_lp(source="s", target="t", nodes=NODES, arcs=ARCS):
    """Shortest path as an LP (send one unit); returns arc flows and node potentials."""
    N = incidence_matrix(nodes, arcs)
    rhs = np.zeros(len(nodes))
    rhs[nodes.index(source)], rhs[nodes.index(target)] = 1, -1
    c = np.array([w for _, _, w, _ in arcs], float)
    res = linprog(c, A_eq=N, b_eq=rhs, bounds=[(0, None)] * len(arcs), method="highs")
    pot = -res.eqlin.marginals
    pot = pot - pot[nodes.index(source)]
    return dict(x=res.x, value=res.fun, potential=dict(zip(nodes, pot)))


# ---------------------------------------------------------------- max flow
def max_flow_steps(source="s", sink="t", nodes=NODES, arcs=ARCS):
    """Edmonds-Karp (shortest augmenting paths). One snapshot per augmentation.

    Each snapshot: flow per arc, the augmenting path as a list of (u, v, direction)
    with direction +1 (forward on arc u->v) or -1 (cancelling flow on arc v->u),
    the bottleneck, and the flow value. The last snapshot has path=None and the
    source side S of a minimum cut.
    """
    cap = {(u, v): k for u, v, _, k in arcs}
    flow = {e: 0 for e in cap}

    def residual_neighbors(u):
        for (x, y), k in cap.items():
            if x == u and flow[(x, y)] < k:
                yield y, (x, y, +1), k - flow[(x, y)]
            if y == u and flow[(x, y)] > 0:
                yield x, (x, y, -1), flow[(x, y)]

    steps, value = [], 0
    while True:
        prev = {source: None}
        q = deque([source])
        while q and sink not in prev:
            u = q.popleft()
            for v, arc, r in residual_neighbors(u):
                if v not in prev:
                    prev[v] = (u, arc, r)
                    q.append(v)
        if sink not in prev:
            steps.append(dict(flow=dict(flow), path=None, bottleneck=0, value=value,
                              cut=sorted(prev, key=nodes.index)))
            return steps
        path, v = [], sink
        while prev[v] is not None:
            u, arc, r = prev[v]
            path.append((arc, r))
            v = u
        path = path[::-1]
        delta = min(r for _, r in path)
        steps.append(dict(flow=dict(flow), path=[a for a, _ in path], bottleneck=delta,
                          value=value, cut=None))
        for (x, y, sgn), _ in path:
            flow[(x, y)] += sgn * delta
        value += delta


def min_cost_flow(amount, source="s", sink="t", nodes=NODES, arcs=ARCS):
    """Send `amount` units from source to sink at minimum cost (LP). Returns flows,
    cost and node potentials (duals of flow conservation)."""
    N = incidence_matrix(nodes, arcs)
    rhs = np.zeros(len(nodes))
    rhs[nodes.index(source)], rhs[nodes.index(sink)] = amount, -amount
    c = np.array([w for _, _, w, _ in arcs], float)
    res = linprog(c, A_eq=N, b_eq=rhs, bounds=[(0, k) for *_, k in arcs], method="highs")
    if res.status != 0:
        return None
    pot = res.eqlin.marginals
    pot = pot - pot[nodes.index(source)]
    return dict(x=res.x, cost=res.fun, potential=dict(zip(nodes, -pot)))


# ---------------------------------------------------------------- bipartite graph
WORKERS = ["Ana", "Ben", "Cai", "Dee"]
JOBS = ["J1", "J2", "J3", "J4"]
# which worker can do which job, with a value (for the assignment problem).
# Ana, Ben and Dee can only do J1 or J2, so at most 3 workers can be matched.
BIP_EDGES = [("Ana", "J1", 7), ("Ana", "J2", 5), ("Ben", "J1", 6), ("Cai", "J1", 3),
             ("Cai", "J2", 8), ("Cai", "J3", 4), ("Cai", "J4", 6), ("Dee", "J2", 2)]
BIP_POS = {**{w: (0.0, 3.0 - i) for i, w in enumerate(WORKERS)},
           **{j: (2.0, 3.0 - i) for i, j in enumerate(JOBS)}}


def bipartite_matching_lp(weighted=False, edges=None):
    """Max-cardinality (or max-weight) bipartite matching by LP; returns edges, value,
    and dual node prices. Integral because the incidence matrix is totally unimodular."""
    edges = edges or BIP_EDGES
    nodes = WORKERS + JOBS
    A = np.zeros((len(nodes), len(edges)))
    for k, (u, v, _) in enumerate(edges):
        A[nodes.index(u), k] = A[nodes.index(v), k] = 1
    w = np.array([wt if weighted else 1.0 for *_, wt in edges])
    res = linprog(-w, A_ub=A, b_ub=np.ones(len(nodes)), bounds=[(0, None)] * len(edges),
                  method="highs")
    chosen = [(u, v) for (u, v, _), x in zip(edges, res.x) if x > 0.5]
    return dict(x=res.x, matching=chosen, value=-res.fun,
                price=dict(zip(nodes, -res.ineqlin.marginals)))


def min_vertex_cover(edges=None):
    """Smallest set of nodes touching every edge (brute force; the graph is tiny)."""
    import itertools as it
    edges = edges or BIP_EDGES
    nodes = WORKERS + JOBS
    for k in range(len(nodes) + 1):
        for C in it.combinations(nodes, k):
            if all(u in C or v in C for u, v, _ in edges):
                return list(C)


# ---------------------------------------------------------------- drawing
def draw_network(ax, label=None, flow=None, highlight=(), dim=(), color_hl="#c44e52",
                 node_color=None, pos=POS, arcs=ARCS, curved=()):
    """Draw the toy network. label: 'cost', 'cap', 'flow', or a dict arc->text."""
    import matplotlib.patches as mpatches
    node_color = node_color or {}
    for (u, v, w, k) in arcs:
        p, q = np.array(pos[u]), np.array(pos[v])
        e = (u, v)
        col = color_hl if e in highlight else ("#bbbbbb" if e in dim else "#555555")
        lw = 2.8 if e in highlight else 1.4
        rad = 0.25 if e in curved else 0.0
        arr = mpatches.FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=14,
                                       color=col, lw=lw, shrinkA=13, shrinkB=13,
                                       connectionstyle=f"arc3,rad={rad}", zorder=2)
        ax.add_patch(arr)
        if label is not None:
            if isinstance(label, dict):
                txt = label.get(e, "")
            elif label == "cost":
                txt = f"{w}"
            elif label == "cap":
                txt = f"{k}"
            elif label == "flow":
                txt = f"{flow[e]:g}/{k}"
            else:
                txt = ""
            m = (p + q) / 2
            d = q - p
            nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d)
            off = nrm * (0.13 + rad * 0.5)
            ax.text(*(m + off), txt, ha="center", va="center", fontsize=9,
                    color=col if e in highlight else "black",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85),
                    zorder=4)
    for v, (x, y) in pos.items():
        fc = node_color.get(v, "white")
        ax.add_patch(mpatches.Circle((x, y), 0.15, fc=fc, ec="black", lw=1.3, zorder=5))
        ax.text(x, y, v, ha="center", va="center", fontsize=11, zorder=6)
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 0.35, max(xs) + 0.35)
    ax.set_ylim(min(ys) - 0.4, max(ys) + 0.4)
    ax.set_aspect("equal")
    ax.axis("off")
