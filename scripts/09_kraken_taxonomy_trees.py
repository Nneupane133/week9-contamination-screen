#!/usr/bin/env python3
"""
Draw a Kraken2 --report-minimizer-data report as a taxonomy tree, two ways:

  reports/figures/kraken_tree_<sample>.png      rectangular cladogram + side panel
  reports/figures/kraken_circular_<sample>.png  radial cladogram + rings

The tree is the NCBI taxonomy restricted to taxa that received reads, NOT a
phylogeny: branch lengths carry no meaning, every step is one retained node.
Unary chains (a taxon with exactly one child, e.g. Metazoa -> ... -> Homo) are
collapsed into a single edge so the figure shows only the places where the
classified reads actually split.

Per tip the side panel / rings show
  * clade reads            reads assigned to the taxon or anything below it
  * distinct minimizers    Kraken2's estimate of distinct minimizers hitting it
  * minimizers per read    distinct minimizers / clade reads. A real genome is
                           hit across many different minimizers; a spurious
                           call is the same few minimizers over and over.
Dots on internal nodes are sized by the reads Kraken2 parked AT that node
(could not resolve further down), which the tip panel would otherwise hide.

usage: 09_kraken_taxonomy_trees.py REPORT [REPORT ...] [--outdir DIR]
"""
import argparse, math, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
from matplotlib import cm

WK = "/home/aubnxn001/week9"

INK, INK2, MUTED, GUIDE = "#1a1a1a", "#444444", "#8a8a8a", "#e6e6e6"
# tab10 first; beyond 10 groups fall back to darker, still-legible extras
PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b",
           "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#393b79", "#637939",
           "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d", "#31a354",
           "#756bb1", "#636363", "#ad494a", "#5254a3"]
FIXED = {"Chordata": "#1f77b4", "Ascomycota": "#d62728",
         "Basidiomycota": "#843c39", "Apicomplexa": "#ff7f0e",
         "Euglenozoa": "#9467bd", "Pseudomonadota": "#e377c2",
         "Bacteroidota": "#17becf", "Bacillota": "#8c6d31",
         "Actinomycetota": "#2ca02c", "Cyanobacteriota": "#637939",
         "Uroviricota": "#8c564b", "Artverviricota": "#7b4173",
         "Kitrinoviricota": "#bcbd22"}
RANK_NAME = {"R": "root", "D": "domain", "K": "kingdom", "P": "phylum",
             "C": "class", "O": "order", "F": "family", "G": "genus",
             "S": "species"}


class Node:
    def __init__(self, name, rank, taxid, clade, direct, mins, distinct):
        self.name, self.rank, self.taxid = name, rank, taxid
        self.clade, self.direct = clade, direct
        self.mins, self.distinct = mins, distinct
        self.children, self.parent = [], None
        self.collapsed = []          # nodes folded into the edge above this one
        self.group = None

    @property
    def is_leaf(self):
        return not self.children


def parse(path):
    root, unclassified, stack = None, 0, []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            raw = f[7]
            name = raw.strip()
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            n = Node(name, f[5], f[6], int(f[1]), int(f[2]),
                     int(f[3]), int(f[4]))
            if f[5] == "U":
                unclassified = n.clade
                continue
            while stack and stack[-1][0] >= depth:
                stack.pop()
            if stack:
                n.parent = stack[-1][1]
                n.parent.children.append(n)
            else:
                root = n
            stack.append((depth, n))
    return root, unclassified


def collapse_unary(n):
    """Fold single-child chains into the deepest node; returns the new subtree."""
    chain = []
    while len(n.children) == 1:
        chain.append(n)
        n = n.children[0]
    n.collapsed = chain
    n.children = [collapse_unary(c) for c in n.children]
    for c in n.children:
        c.parent = n
    return n


def walk(n):
    yield n
    for c in n.children:
        yield from walk(c)


def assign_groups(root):
    """Colour group = phylum if there is one on the path, else kingdom, else domain."""
    def full_path(n):
        names = []
        while n is not None:
            names = [(x.rank, x.name) for x in n.collapsed + [n]] + names
            n = n.parent
        return names
    for n in walk(root):
        ranks = {}
        for r, nm in full_path(n):
            if r in ("P", "K", "R2"):
                ranks.setdefault(r[0], nm)
        n.group = ranks.get("P") or ranks.get("K") or ranks.get("R") or n.name
    groups = []
    for n in walk(root):
        if n.is_leaf and n.group not in groups:
            groups.append(n.group)
    # common groups keep one colour across samples so figures can be compared
    colours = {g: FIXED[g] for g in groups if g in FIXED}
    spare = [c for c in PALETTE if c not in FIXED.values()]
    for g in groups:
        if g not in colours:
            colours[g] = spare[len(colours) % len(spare)] if spare else MUTED
    return colours


def layout(root):
    """x = depth (one step per retained node), y = tip order."""
    leaves = []

    def place(n, d):
        n.x = d
        if n.is_leaf:
            n.y = len(leaves)
            leaves.append(n)
        else:
            # big clades first, like a ladderized tree
            n.children.sort(key=lambda c: -c.clade)
            for c in n.children:
                place(c, d + 1)
            n.y = (n.children[0].y + n.children[-1].y) / 2
    place(root, 0)
    return leaves


def tip_label(n):
    rank = RANK_NAME.get(n.rank[0], n.rank)
    if n.rank[1:]:
        rank = f"{rank}+{n.rank[1:]}"
    return rank


def short(name, width=34):
    return name if len(name) <= width else name[:width - 1] + "…"


def mpr(n):
    return n.distinct / n.clade if n.clade else 0.0


def direct_size(reads):
    return 6 + 9 * math.log10(reads + 1)


# ---------------------------------------------------------------- rectangular

def draw_rect(root, leaves, colours, sample, unclassified, out):
    nl = len(leaves)
    maxd = max(l.x for l in leaves)
    fig_h = max(6.0, 0.24 * nl + 2.4)
    fig = plt.figure(figsize=(15, fig_h), dpi=200)
    top = 1 - 1.15 / fig_h
    bot = 0.9 / fig_h
    ax = fig.add_axes([0.02, bot, 0.56, top - bot])
    panels = [fig.add_axes([0.60 + i * 0.13, bot, 0.115, top - bot], sharey=ax)
              for i in range(3)]

    def Y(n):
        return nl - 1 - n.y

    for n in walk(root):
        for c in n.children:
            ax.plot([n.x, n.x], [Y(n), Y(c)], color=INK, lw=0.9,
                    solid_capstyle="butt")
            ax.plot([n.x, c.x], [Y(c), Y(c)], color=INK, lw=0.9)
            # folded single-child taxa that still hold reads sit on the edge
            k = len(c.collapsed)
            for i, f in enumerate(c.collapsed):
                if f.direct:
                    ax.scatter(n.x + (c.x - n.x) * (i + 1) / (k + 1), Y(c),
                               s=direct_size(f.direct), color=MUTED, zorder=3,
                               linewidths=0)
        if not n.is_leaf:
            if n.direct:
                ax.scatter(n.x, Y(n), s=direct_size(n.direct), color=INK,
                           zorder=3, linewidths=0)
            if n is not root:
                ax.text(n.x - 0.06, Y(n) + 0.12, n.name, fontsize=6.4,
                        color=INK2, ha="right", va="bottom")
                ax.text(n.x - 0.06, Y(n) - 0.12, f"{n.clade:,}", fontsize=6,
                        color=MUTED, ha="right", va="top")

    for l in leaves:
        y = Y(l)
        ax.plot([l.x, maxd + 5.2], [y, y], color=GUIDE, lw=0.6, zorder=0)
        if l.direct:
            ax.scatter(l.x, y, s=direct_size(l.direct), color=colours[l.group],
                       zorder=3, linewidths=0)
        style = "italic" if l.rank.startswith(("S", "G")) else "normal"
        t = ax.text(l.x + 0.18, y, short(l.name, 40), fontsize=8.2, color=colours[l.group],
                    va="center", style=style)
        t.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.6))
        ax.annotate(f"{tip_label(l)}  {l.taxid}", xy=(1, 0.5),
                    xycoords=t, xytext=(4, 0), textcoords="offset points",
                    fontsize=6, color=MUTED, va="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.6))

    ax.set_xlim(-0.3, maxd + 5.2)
    ax.set_ylim(-0.8, nl - 0.2)
    ax.axis("off")

    cols = [("clade reads", [l.clade for l in leaves], "{:,}"),
            ("distinct minimizers", [l.distinct for l in leaves], "{:,}"),
            ("distinct minimizers / read", [mpr(l) for l in leaves], "{:.2g}")]
    for p, (title, vals, fmt) in zip(panels, cols):
        vmax = max(vals) or 1
        lo = min(v for v in vals if v > 0) if any(v > 0 for v in vals) else 1
        lo = min(lo, 1)
        for l, v in zip(leaves, vals):
            y = Y(l)
            w = (math.log10(v) - math.log10(lo) + 0.15) / \
                (math.log10(vmax) - math.log10(lo) + 0.15) if v > 0 else 0
            p.barh(y, w, height=0.66, color=colours[l.group], alpha=0.85,
                   linewidth=0)
            p.text(w + 0.02, y, fmt.format(v), fontsize=6.2, color=INK2,
                   va="center")
        p.set_xlim(0, 1.45)
        p.set_title(title, fontsize=8.5, color=INK, pad=4)
        p.set_yticks([])
        for s in ("top", "right", "left"):
            p.spines[s].set_visible(False)
        p.set_xticks([0, 1])
        p.set_xticklabels([f"{lo:.2g}" if lo < 1 else "1", f"{vmax:,.3g}"],
                          fontsize=6.5, color=MUTED)
        p.set_xlabel("log scale", fontsize=6.5, color=MUTED, labelpad=1)
        p.tick_params(length=2, colors=MUTED)

    total = root.clade + unclassified
    fig.text(0.5, 1 - 0.32 / fig_h,
             f"{sample}: Kraken2 classified reads on the NCBI taxonomy",
             ha="center", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.5, 1 - 0.62 / fig_h,
             f"{root.clade:,} of {total:,} reads classified "
             f"({100 * root.clade / total:.2f}%). Cladogram, not a phylogeny: "
             "one step per branching taxon, single-child chains folded into the edge (grey dots)",
             ha="center", fontsize=8.5, color=INK2)

    handles = [Line2D([], [], marker="o", ls="", color=INK,
                      markersize=math.sqrt(direct_size(r)), label=f"{r:,}")
               for r in (1, 100, 10000)]
    handles.append(Line2D([], [], marker="o", ls="", color=MUTED, markersize=4,
                          label="folded taxon"))
    leg = ax.legend(handles=handles, title="reads Kraken2 assigned at that node (could not resolve lower)",
                    loc="lower left", fontsize=6.5, title_fontsize=7,
                    frameon=False, ncol=4, handletextpad=0.2,
                    bbox_to_anchor=(0, -0.075))
    leg._legend_box.align = "left"
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)


# ------------------------------------------------------------------- circular

def draw_circular(root, leaves, colours, sample, unclassified, out):
    nl = len(leaves)
    maxd = max(l.x for l in leaves)
    gap = 0.22 * 2 * math.pi / 8          # opening for the radial scale
    span = 2 * math.pi - gap
    theta0 = math.pi / 2 - gap / 2

    def ang(y):
        return theta0 - span * (y + 0.5) / nl

    r_in, r_tree = 0.42, 1.0              # tree occupies [r_in, r_tree]

    def R(x):
        return r_in + (r_tree - r_in) * x / maxd

    k = max(1.0, (nl / 45) ** 0.5)        # grow the canvas with the tip count
    fig = plt.figure(figsize=(13.5 * k, 11 * k), dpi=200)
    ax = fig.add_axes([0.0, 0.03, 0.80, 0.92])
    ax.set_aspect("equal")
    ax.axis("off")

    def pol(r, a):
        return r * math.cos(a), r * math.sin(a)

    for n in walk(root):
        if n.is_leaf:
            continue
        a1, a2 = ang(n.children[0].y), ang(n.children[-1].y)
        steps = max(2, int(abs(a2 - a1) / 0.02))
        arc = [pol(R(n.x), a1 + (a2 - a1) * i / steps) for i in range(steps + 1)]
        ax.plot(*zip(*arc), color=INK, lw=0.8)
        for c in n.children:
            a = ang(c.y)
            ax.plot(*zip(pol(R(n.x), a), pol(R(c.x), a)), color=INK, lw=0.8)
        if n.direct:
            ax.scatter(*pol(R(n.x), ang(n.y)), s=direct_size(n.direct),
                       color=INK, zorder=3, linewidths=0)

    # rings, outward from the tree
    ring_w, ring_gap = 0.13, 0.025
    r_grp = r_tree + 0.06
    rings = [("clade reads", [l.clade for l in leaves], "Blues"),
             ("distinct minimizers", [l.distinct for l in leaves], "Greens"),
             ("distinct minimizers / read", [mpr(l) for l in leaves], "Oranges")]
    dth = math.degrees(span / nl)
    for l in leaves:
        a = ang(l.y)
        ax.plot(*zip(pol(R(l.x), a), pol(r_grp, a)), color=GUIDE, lw=0.6,
                zorder=0)
        if l.direct:
            ax.scatter(*pol(R(l.x), a), s=direct_size(l.direct),
                       color=colours[l.group], zorder=3, linewidths=0)
        ax.add_patch(Wedge((0, 0), r_grp + 0.05, math.degrees(a) - dth / 2 + 0.3,
                           math.degrees(a) + dth / 2 - 0.3, width=0.05,
                           color=colours[l.group], linewidth=0))

    r = r_grp + 0.05 + ring_gap
    cbar_specs = []
    for title, vals, cmap_name in rings:
        pos = [v for v in vals if v > 0]
        norm = LogNorm(vmin=min(pos), vmax=max(pos)) if pos else None
        cmap = matplotlib.colormaps[cmap_name]
        ax.add_patch(Wedge((0, 0), r + ring_w, math.degrees(theta0 - span),
                           math.degrees(theta0), width=ring_w,
                           color="#f4f4f4", linewidth=0, zorder=0))
        for l, v in zip(leaves, vals):
            if v <= 0:
                continue
            frac = 0.15 + 0.85 * norm(v) if norm.vmax > norm.vmin else 1.0
            a = math.degrees(ang(l.y))
            ax.add_patch(Wedge((0, 0), r + ring_w * frac, a - dth / 2 + 0.3,
                               a + dth / 2 - 0.3, width=ring_w * frac,
                               color=cmap(0.35 + 0.6 * norm(v)), linewidth=0))
        ax.text(*pol(r + ring_w / 2, math.pi / 2), str(len(cbar_specs) + 1),
                fontsize=8, fontweight="bold", color=INK2, ha="center",
                va="center")
        cbar_specs.append((title, norm, cmap))
        r += ring_w + ring_gap

    for l in leaves:
        a = ang(l.y)
        deg = math.degrees(a) % 360
        flip = 90 < deg < 270
        style = "italic" if l.rank.startswith(("S", "G")) else "normal"
        ax.text(*pol(r + 0.02, a), short(l.name), fontsize=7.2,
                color=colours[l.group], style=style,
                rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                ha="right" if flip else "left", va="center")

    # clade labels for the deepest-branching internal nodes
    for n in walk(root):
        if not n.is_leaf and n.x <= 2 and n is not root:
            ax.text(*pol(R(n.x) + 0.06, ang(n.y)), n.name, fontsize=6.5,
                    color=INK2, ha="center", va="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.5,
                              alpha=0.85))

    total = root.clade + unclassified
    ax.text(0, 0.06, sample, ha="center", va="center", fontsize=11,
            fontweight="bold", color=INK)
    ax.text(0, -0.05, f"{root.clade:,} classified\nof {total:,} reads\n"
            f"({100 * root.clade / total:.2f}%)", ha="center", va="center",
            fontsize=7.5, color=INK2, linespacing=1.3)
    lim = r + 0.75 / k ** 0.5
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    lax = fig.add_axes([0.80, 0.08, 0.19, 0.84])
    lax.axis("off")
    used = [g for g in colours if any(l.group == g for l in leaves)]
    leg1 = lax.legend(handles=[Patch(color=colours[g], label=g) for g in used],
                      title="tip colour: phylum (else kingdom)",
                      loc="upper left", fontsize=7.5, title_fontsize=8,
                      frameon=False)
    leg1._legend_box.align = "left"
    lax.add_artist(leg1)
    handles = [Line2D([], [], marker="o", ls="", color=INK,
                      markersize=math.sqrt(direct_size(v)), label=f"{v:,}")
               for v in (1, 100, 10000)]
    leg2 = lax.legend(handles=handles, title="reads assigned at that node",
                      loc="center left", fontsize=7.5, title_fontsize=8,
                      frameon=False, bbox_to_anchor=(0, 0.42))
    leg2._legend_box.align = "left"
    lax.add_artist(leg2)
    for i, (title, norm, cmap) in enumerate(cbar_specs):
        cax = fig.add_axes([0.81, 0.24 - i * 0.07, 0.14, 0.012])
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                          orientation="horizontal")
        cb.ax.tick_params(labelsize=6, colors=INK2, length=2)
        cb.set_label(f"{title} (ring {i + 1}, log)", fontsize=7, color=INK2,
                     labelpad=2)
        cb.outline.set_visible(False)

    fig.text(0.5, 0.012,
             "Cladogram of the NCBI taxonomy restricted to taxa with reads "
             "(not a phylogeny; radius = rank steps, single-child chains "
             "collapsed). Rings, inside out: tip colour, clade reads, distinct "
             "minimizers, distinct minimizers per read.",
             ha="center", fontsize=7, color=INK2)
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--outdir", default=os.path.join(WK, "reports", "figures"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    for path in a.reports:
        sample = os.path.basename(path).split(".")[0]
        root, unclassified = parse(path)
        root = collapse_unary(root)
        colours = assign_groups(root)
        leaves = layout(root)
        for kind, fn in (("tree", draw_rect), ("circular", draw_circular)):
            out = os.path.join(a.outdir, f"kraken_{kind}_{sample}.png")
            fn(root, leaves, colours, sample, unclassified, out)
            print(out)


if __name__ == "__main__":
    main()
