#!/usr/bin/env python3
"""
STEP 7: the two class-facing deliverables.

  reports/detection_table.<level>.{tsv,md}      concise detection table
  reports/figures/read_provenance_<sample>.png  read-provenance tree

WHAT THIS ADDS ON TOP OF 05_score_detections.py
-----------------------------------------------
05 measures. This script adjudicates, and it does three things 05 cannot:

1. PANEL CORRECTIONS. Two candidates are mis-measured by construction because
   of how the competitive reference was built, not because of anything in the
   data. Both are corrected here, explicitly and with a stated reason, rather
   than quietly:

     * PhiX174 appears TWICE in the reference -- once as the decoy and once
       inside UniVec_Core, which carries the whole phi-X174 genome as
       `gnl|uv|J02482.1:1-5386-49` (5,375/5,375 identical bases). Every PhiX
       read therefore has two perfect placements, gets MAPQ 0, and is discarded
       by the MAPQ>=20 filter. Measured at MAPQ 0 with the two copies summed,
       PhiX is a complete genome at ~14x.

     * UniVec_Core is a COLLECTION of ~3,000 unrelated vector and adapter
       entries, not an organism. Breadth over the whole 688 kb collection is
       meaningless -- an adapter that is present at a million-fold depth still
       covers 0.2% of the collection. It is measured on the single entry that
       carries its reads instead.

   Both corrections are the same underlying mistake: treating a reference whose
   internal structure matters as if it were one genome.

2. OBSERVED VS EXPECTED BREADTH. The 0.5% breadth floor in 05 is
   reference-length dependent and therefore unreachable for a large reference
   at low abundance: 18,000 reads cannot cover 0.5% of the 3.1 Gb human genome
   no matter how real they are. The length-free version of the same question is

       evenness = observed breadth / (1 - exp(-depth))

   the denominator being the breadth a randomly-placed library of that depth
   would give (Lander-Waterman). ~1 means the reads are spread over the genome
   the way a real population's reads are; <<1 means they are piled on a few
   conserved or repetitive loci, which is what a cross-mapping artifact looks
   like. This separates the two cases that a fixed breadth threshold conflates.

3. ATTRIBUTION. The final bin and the source category come from
   metadata/attribution.tsv, where each is tied to named evidence. Where the
   final bin disagrees with 05's automatic bin, both are printed.
"""
import argparse, csv, math, os, re, textwrap
from collections import defaultdict

WK = "/home/aubnxn001/week9"

SAMPLES = ("BYJ10_RNA", "SRR39197041_DNA")
HOSTS = {"BYJ10_RNA": "Gallus_gallus_GRCg7b",
         "SRR39197041_DNA": "Aedes_aegypti_AaegL5"}
HOST_LABEL = {"BYJ10_RNA": "Gallus gallus (chicken)",
              "SRR39197041_DNA": "Aedes aegypti (mosquito)"}
SUBTITLE = {"BYJ10_RNA": "jejunum RNA-seq, 2x150 bp, flowcell 23C55KLT4 lane 5",
            "SRR39197041_DNA": "whole-genome DNA, 2x151 bp, flowcell H5F2GDSX3 lane 4"}

# The PhiX174 genome, verbatim, inside UniVec_Core. Verified with
#   minimap2 -x asm20 univec.fa phix.fa  ->  5375/5375 identical bases.
PHIX_IN_UNIVEC = "UniVec_Core::gnl|uv|J02482.1:1-5386-49"

# Panel corrections: (sample, source) -> how to measure it instead.
#   mode "merge_at_mq0"  sum this source with `extra` contigs and read at MAPQ 0
#   mode "single_contig" measure on one named contig of a collection reference
CORRECTIONS = {
    ("SRR39197041_DNA", "PhiX174"): dict(
        mode="merge_at_mq0", extra=[PHIX_IN_UNIVEC], ref_len=5386,
        label="PhiX174 (+ its duplicate inside UniVec_Core)",
        note="identical copy in UniVec_Core forces MAPQ 0; summed and read at MAPQ 0"),
    ("SRR39197041_DNA", "UniVec_Core"): dict(
        mode="single_contig", contig="UniVec_Core::gnl|uv|NGB00360.1:1-58",
        label="Illumina PCR primer (UniVec NGB00360)",
        note="collection reference; measured on the one 58 bp entry that carries the reads"),
    ("BYJ10_RNA", "UniVec_Core"): dict(
        mode="single_contig", contig="UniVec_Core::gnl|uv|NGB00846.1:1-65",
        label="NEBNext adaptor (UniVec NGB00846)",
        note="collection reference; measured on the one 65 bp entry that carries the reads"),
}

BINS = ("noise", "trace", "consequential", "biological")

# Ordinal ramp, light -> dark, one hue: the four bins are an ORDERED scale, not
# four categories. Steps 250/400/550/700 of the blue ramp; the lightest clears
# 2:1 on the light surface, which is the ordinal floor.
BIN_COLOR = {"noise": "#86b6ef", "trace": "#3987e5",
             "consequential": "#1c5cab", "biological": "#0d366b"}
BIN_INK = {"noise": "#0b0b0b", "trace": "#ffffff",
           "consequential": "#ffffff", "biological": "#ffffff"}

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
RULE = "#e1e0d9"
HOSTC = "#c3c2b7"


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def load_coverage(path):
    """samtools coverage -> {contig: dict}, and {source: aggregated dict}."""
    per, agg = {}, defaultdict(lambda: dict(length=0, covbases=0, reads=0, depthsum=0.0))
    if not os.path.exists(path):
        return per, agg
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            name, length = f[0], int(f[2])
            d = dict(length=length, reads=int(f[3]), covbases=int(f[4]),
                     depth=float(f[6]), mapq=float(f[8]))
            per[name] = d
            a = agg[name.split("::", 1)[0]]
            a["length"] += length
            a["covbases"] += d["covbases"]
            a["reads"] += d["reads"]
            a["depthsum"] += d["depth"] * length
    return per, agg


def load_tsv(path, key):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get(key):
                out[r[key]] = r
    return out


def load_attribution(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            out[(f[0], f[1])] = dict(display=f[2], bin=f[3], source=f[4], why=f[5])
    return out


def load_host_competition(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] in ("source", ""):
                continue
            try:
                out[f[0]] = float(f[3])
            except ValueError:
                pass
    return out


def total_reads(flagstat):
    """Primary alignments = reads that went in. See 05 for why not 'in total'."""
    if not os.path.exists(flagstat):
        return None
    for line in open(flagstat):
        m = re.match(r"(\d+) \+ \d+ primary$", line.strip())
        if m:
            return int(m.group(1))
    return None


# --------------------------------------------------------------------------
# assembly of the per-sample row set
# --------------------------------------------------------------------------
def evenness(breadth, depth):
    """observed breadth / Lander-Waterman expectation, 1 - exp(-depth)."""
    if depth is None or depth <= 0:
        return None
    exp = 1.0 - math.exp(-depth)
    return breadth / exp if exp > 0 else None


def build_rows(sample, level, attrib):
    per0, agg0 = load_coverage(f"{WK}/results/coverage/{level}/{sample}.coverage.mq0.tsv")
    per20, agg20 = load_coverage(f"{WK}/results/coverage/{level}/{sample}.coverage.mq20.tsv")
    det = load_tsv(f"{WK}/results/coverage/{level}/{sample}.detections.tsv", "source")
    hc = load_host_competition(f"{WK}/results/host_competition/{level}/{sample}.host_competition.tsv")
    total = total_reads(f"{WK}/results/coverage/{level}/{sample}.flagstat.txt")
    host = HOSTS[sample]

    rows = []
    for src, d in det.items():
        if src == host:
            continue
        corr = CORRECTIONS.get((sample, src))
        label = attrib.get((sample, src), {}).get("display", src.replace("_", " "))
        note = ""
        a20, a0 = agg20.get(src, {}), agg0.get(src, {})

        if corr and corr["mode"] == "merge_at_mq0":
            L = corr["ref_len"]
            reads = a0.get("reads", 0) + sum(per0.get(c, {}).get("reads", 0) for c in corr["extra"])
            cov = max([a0.get("covbases", 0)] +
                      [per0.get(c, {}).get("covbases", 0) for c in corr["extra"]])
            cov = min(cov, L)
            dsum = a0.get("depthsum", 0.0) + sum(
                per0.get(c, {}).get("depth", 0.0) * per0.get(c, {}).get("length", 0)
                for c in corr["extra"])
            depth, breadth = dsum / L, cov / L
            reads20 = reads                       # MAPQ is meaningless here, say so
            measured, note = "MAPQ 0 (merged)", corr["note"]
            label = corr["label"]
        elif corr and corr["mode"] == "single_contig":
            c = corr["contig"]
            c20, c0 = per20.get(c, {}), per0.get(c, {})
            L = c20.get("length") or c0.get("length") or 1
            reads20, reads = c20.get("reads", 0), c0.get("reads", 0)
            depth, breadth = c20.get("depth", 0.0), c20.get("covbases", 0) / L
            measured, note = "MAPQ>=20 (single entry)", corr["note"]
            label = corr["label"]
            # Measuring one entry of a collection leaves the rest of the
            # collection's reads unaccounted for. Emit them as their own row so
            # the provenance tree still partitions the library exactly, instead
            # of letting them vanish into the "ambiguous" residual.
            claimed = c20.get("reads", 0)
            if (sample, "PhiX174") in CORRECTIONS:
                claimed += per20.get(PHIX_IN_UNIVEC, {}).get("reads", 0)
            rest = a20.get("reads", 0) - claimed
            if rest > 0:
                rows.append(dict(
                    sample=sample, source=f"{src}__rest",
                    display=f"{src}, all other entries (vectors, other adapters)",
                    ref_len=a20.get("length", 0), reads=rest, reads_all=rest,
                    read_fraction=(rest / total if total else float("nan")),
                    depth=float("nan"), breadth=float("nan"), evenness=None,
                    mapq_survival=None, host_competing=None,
                    bin_auto="noise", bin_called="noise", source_called="reagent",
                    why=("the remainder of the UniVec_Core collection once the one "
                         "entry that carries the bulk of the reads is counted "
                         "separately; spread thinly over ~3,000 unrelated vector "
                         "and adapter sequences"),
                    measured_at="MAPQ>=20", correction=""))
        else:
            L = a20.get("length") or a0.get("length") or 1
            reads20, reads = a20.get("reads", 0), a0.get("reads", 0)
            # 05 switches a host-integrated element to MAPQ 0, since the host is
            # legitimately competing for its reads. Honour the same choice here.
            if d.get("measured_at") == "mq0":
                depth = a0.get("depthsum", 0.0) / L
                breadth = a0.get("covbases", 0) / L
                reads20, measured = reads, "MAPQ 0 (host-integrated)"
            else:
                depth = a20.get("depthsum", 0.0) / L
                breadth = a20.get("covbases", 0) / L
                measured = "MAPQ>=20"

        att = attrib.get((sample, src), {})
        rows.append(dict(
            sample=sample, source=src, display=label,
            ref_len=L, reads=reads20, reads_all=reads,
            read_fraction=(reads20 / total if total else float("nan")),
            depth=depth, breadth=breadth,
            evenness=evenness(breadth, depth),
            # Meaningless for a row read at MAPQ 0 -- it would always be 1.0.
            mapq_survival=(None if measured.startswith("MAPQ 0")
                           else (reads20 / reads if reads else None)),
            host_competing=hc.get(src),
            bin_auto=d.get("bin", "-"),
            bin_called=att.get("bin", d.get("bin", "-")),
            source_called=att.get("source", "-"),
            why=att.get("why", ""),
            measured_at=measured, correction=note,
        ))

    rows.sort(key=lambda r: (-r["reads"], r["display"]))

    host_row = dict(
        sample=sample, source=host, display=HOST_LABEL[sample],
        reads=agg20.get(host, {}).get("reads", 0),
        reads_all=agg0.get(host, {}).get("reads", 0),
        ref_len=agg20.get(host, {}).get("length", 0),
        bin_called="host", source_called="host", why="reference organism for this library",
    )
    # Everything samtools saw, at MAPQ 0, is a true partition: one primary
    # alignment per read. What is left over never aligned to the panel at all.
    placed0 = sum(a["reads"] for a in agg0.values())
    return dict(sample=sample, total=total, host=host_row, rows=rows,
                placed_mq0=placed0, unplaced=(total - placed0) if total else 0)


# --------------------------------------------------------------------------
# table
# --------------------------------------------------------------------------
def pct(x, nd=4):
    return "-" if x is None or x != x else f"{x * 100:.{nd}f}%"


def fmt_pct(frac):
    """Percent with enough significant figures to stay non-zero down to 1e-8."""
    if frac is None or frac != frac:
        return "-"
    if not frac:
        return "0%"
    p = frac * 100
    if p >= 99.999:
        return "100%"
    nd = 1 if p >= 10 else max(1, min(6, 2 - int(math.floor(math.log10(p)))))
    return f"{p:.{nd}f}%"


def num(x, nd=3):
    if x is None or x != x:
        return "-"
    return f"{x:,.0f}" if abs(x) >= 1000 else f"{x:.{nd}g}"


def write_table(data, level, out):
    cols = ["sample", "organism", "bin_called", "bin_auto", "likely_source",
            "reads", "read_fraction", "depth", "breadth", "evenness",
            "mapq_survival", "host_competing", "ref_len", "measured_at",
            "panel_correction", "evidence"]
    tsv = f"{out}/detection_table.{level}.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for d in data:
            for r in d["rows"]:
                w.writerow([
                    r["sample"], r["display"], r["bin_called"], r["bin_auto"],
                    r["source_called"], r["reads"],
                    f"{r['read_fraction']:.3e}", f"{r['depth']:.6g}",
                    f"{r['breadth']:.6g}",
                    f"{r['evenness']:.3f}" if r["evenness"] is not None else "-",
                    f"{r['mapq_survival']:.3f}" if r["mapq_survival"] is not None else "-",
                    f"{r['host_competing']:.3f}" if r["host_competing"] is not None else "-",
                    r["ref_len"], r["measured_at"], r["correction"], r["why"]])
    print(f"wrote {tsv}")

    order = {b: i for i, b in enumerate(BINS)}
    md = [f"# Detection table ({level})", "",
          "`bin` is the final call; where the automatic rule in "
          "`05_score_detections.py` disagreed, its verdict is shown as "
          "`auto: <bin>` and the reason is in the evidence column.",
          "",
          "`evenness` is observed breadth divided by the breadth a randomly "
          "placed library of that depth would give (`1 - e^-depth`). Near 1 "
          "means the reads are spread across the reference like a real "
          "population's; near 0 means they are piled on a few loci.",
          "",
          "† measured at MAPQ 0 rather than MAPQ>=20, because the sequence is "
          "integrated in the host genome or present twice in the reference "
          "panel — in both cases a MAPQ filter throws away its real reads. "
          "MAPQ survival is left blank for those rows, where it would be 1.0 "
          "by construction.", ""]
    for d in data:
        s = d["sample"]
        md += [f"## {s} — {HOST_LABEL[s]}", "",
               f"*{SUBTITLE[s]}; {d['total']:,} reads; host claims "
               f"{d['host']['reads']:,} ({d['host']['reads'] / d['total']:.1%}) "
               f"at MAPQ>=20.*", "",
               "| organism | bin | likely source | reads | read fraction | depth | breadth | evenness | MAPQ surv. | why |",
               "|---|---|---|---:|---:|---:|---:|---:|---:|---|"]
        for r in sorted(d["rows"], key=lambda r: (-order.get(r["bin_called"], 0), -r["reads"])):
            b = f"**{r['bin_called']}**"
            if r["bin_auto"] != r["bin_called"] and r["bin_auto"] != "-":
                b += f"<br>`auto: {r['bin_auto']}`"
            md.append("| {} | {} | {} | {:,} | {} | {} | {} | {} | {} | {} |".format(
                r["display"] + (" †" if r["measured_at"].startswith("MAPQ 0") else ""),
                b, r["source_called"], r["reads"],
                fmt_pct(r["read_fraction"]), num(r["depth"]), fmt_pct(r["breadth"]),
                num(r["evenness"], 2), num(r["mapq_survival"], 2), r["why"]))
        md.append("")
    mdp = f"{out}/detection_table.{level}.md"
    open(mdp, "w").write("\n".join(md) + "\n")
    print(f"wrote {mdp}")
    return mdp


# --------------------------------------------------------------------------
# read-provenance tree
# --------------------------------------------------------------------------
def draw_tree(d, level, outdir, min_reads=1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    s, total = d["sample"], d["total"]
    shown = [r for r in d["rows"] if r["reads"] >= min_reads]
    hidden = [r for r in d["rows"] if r["reads"] < min_reads]
    nonhost_unique = sum(r["reads"] for r in d["rows"])
    ambiguous = d["placed_mq0"] - d["host"]["reads"] - nonhost_unique

    # (depth, label, reads, bin, source) -- depth 0 root, 1 branch, 2 leaf
    items = [(1, f"Host — {HOST_LABEL[s]}", d["host"]["reads"], "host", "")]
    items.append((1, "Non-host, uniquely placed (MAPQ>=20)", nonhost_unique, "", ""))
    # A row whose reads had to be read at MAPQ 0 -- because the sequence is
    # integrated in the host, or because the panel carries it twice -- is not a
    # MAPQ>=20 count and is marked so the branch heading cannot mislead.
    for r in shown:
        mark = " †" if r["measured_at"].startswith("MAPQ 0") else ""
        items.append((2, r["display"] + mark, r["reads"],
                      r["bin_called"], r["source_called"]))
    if hidden:
        items.append((2, f"{len(hidden)} further references, 0 reads each",
                      0, "noise", "absent"))
    items.append((1, "Mapped, but ambiguous between references (MAPQ<20)",
                  ambiguous, "", ""))
    items.append((1, "Unplaced — no alignment to any reference in the panel",
                  d["unplaced"], "", ""))

    n = len(items)
    fig_h = 2.35 + 0.295 * n
    fig, ax = plt.subplots(figsize=(14.2, fig_h), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    top = 1 - (1.35 / fig_h)
    row_h = 0.295 / fig_h
    x_name = {0: 0.012, 1: 0.038, 2: 0.078}
    # CHIP_W is fixed so the "likely source" column never collides with the
    # widest chip ("consequential").
    X_READS, X_PCT, X_CHIP = 0.545, 0.645, 0.665
    CHIP_W = 0.112
    X_SRC = X_CHIP + CHIP_W + 0.022

    ax.text(0.012, 1 - (0.32 / fig_h), s.replace("_", "  "),
            fontsize=15, fontweight="bold", color=INK, va="center")
    ax.text(0.012, 1 - (0.64 / fig_h),
            f"{SUBTITLE[s]}   ·   where every read went",
            fontsize=9.5, color=INK2, va="center")

    hdr = top + 0.62 * row_h
    for x, t, ha in ((X_READS, "reads", "right"), (X_PCT, "% of library", "right"),
                     (X_CHIP, "significance", "left"), (X_SRC, "likely source", "left")):
        ax.text(x, hdr, t, fontsize=8, color=MUTED, va="center", ha=ha)
    ax.plot([0.012, 0.988], [hdr - 0.45 * row_h] * 2, color=RULE, lw=0.8)

    # root
    y0 = top - 0.1 * row_h
    ax.text(x_name[0], y0, "LIBRARY", fontsize=10.5, fontweight="bold",
            color=INK, va="center")
    ax.text(X_READS, y0, f"{total:,}", fontsize=10, color=INK,
            va="center", ha="right")
    ax.text(X_PCT, y0, "100%", fontsize=10, color=INK, va="center", ha="right")

    ys, prev_branch_y = [], y0
    branch_xs = {0: x_name[0] + 0.009, 1: x_name[1] + 0.009}
    spans = defaultdict(list)          # parent depth -> child y positions

    for i, (depth, label, reads, bin_, src) in enumerate(items):
        y = y0 - (i + 1) * row_h
        ys.append(y)
        spans[depth - 1].append(y)
        frac = reads / total if total else 0

        # elbow connector from the parent spine
        xs_ = branch_xs[depth - 1]
        ax.plot([xs_, xs_], [y + row_h * 0.5, y], color=RULE, lw=0.9,
                solid_capstyle="butt", zorder=1)
        ax.plot([xs_, x_name[depth] - 0.006], [y, y], color=RULE, lw=0.9, zorder=1)

        bold = depth == 1
        ax.text(x_name[depth], y, label,
                fontsize=10 if bold else 9.3,
                fontweight="bold" if bold else "normal",
                color=INK if bold else INK2, va="center")
        ax.text(X_READS, y, f"{reads:,}", fontsize=9.5 if bold else 9,
                color=INK if bold else INK2, va="center", ha="right",
                fontfamily="DejaVu Sans")
        ax.text(X_PCT, y, fmt_pct(frac),
                fontsize=9.5 if bold else 9, color=INK if bold else INK2,
                va="center", ha="right")

        if bin_:
            face = HOSTC if bin_ == "host" else BIN_COLOR[bin_]
            ink = INK if bin_ == "host" else BIN_INK[bin_]
            ax.add_patch(FancyBboxPatch(
                (X_CHIP, y - 0.33 * row_h), CHIP_W, 0.66 * row_h,
                boxstyle="round,pad=0,rounding_size=0.004",
                linewidth=0, facecolor=face, zorder=2))
            ax.text(X_CHIP + CHIP_W / 2, y, bin_, fontsize=8, color=ink,
                    va="center", ha="center", zorder=3)
        if src:
            ax.text(X_SRC, y, src, fontsize=8.6, color=INK2, va="center")

    # vertical spines
    for pdepth, child_ys in spans.items():
        if pdepth not in branch_xs or not child_ys:
            continue
        xs_ = branch_xs[pdepth]
        ax.plot([xs_, xs_], [max(child_ys), min(child_ys)], color=RULE, lw=0.9, zorder=1)

    foot = min(ys) - 1.05 * row_h
    ax.plot([0.012, 0.988], [foot + 0.45 * row_h] * 2, color=RULE, lw=0.8)
    notes = [r for r in shown if r["correction"]]
    txt = ("Counts partition the library: each read has one primary alignment. "
           "“Uniquely placed” = that alignment survives MAPQ>=20, i.e. no other "
           "reference explains the read as well.")
    if any(r["measured_at"].startswith("MAPQ 0") for r in shown):
        txt += ("  † counted at MAPQ 0, not MAPQ>=20: the reference is either "
                "integrated in the host genome or present twice in the panel, so "
                "a MAPQ filter would discard its real reads.")
    if notes:
        txt += "  " + "  ".join(f"[{r['display'].split(' (')[0]}] {r['correction']}."
                                for r in notes)
    ax.text(0.012, foot, "\n".join(textwrap.wrap(txt, 168)),
            fontsize=7.6, color=MUTED, va="top", linespacing=1.5)

    os.makedirs(outdir, exist_ok=True)
    p = f"{outdir}/read_provenance_{s}.png"
    fig.savefig(p, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    print(f"wrote {p}")
    return p


BEGIN = "<!-- BEGIN generated: detection table -->"
END = "<!-- END generated: detection table -->"


def inject_readme(mdp, readme):
    """Drop the generated table into README between the two markers.

    The README is the thing people read; keeping the table generated rather
    than pasted means it cannot drift from the numbers in results/.
    """
    if not os.path.exists(readme):
        return
    body = open(mdp).read()
    body = "\n".join(l for l in body.split("\n") if not l.startswith("# "))
    text = open(readme).read()
    if BEGIN not in text or END not in text:
        print(f"note: markers not found in {readme}; table not injected")
        return
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    open(readme, "w").write(f"{head}{BEGIN}\n{body.strip()}\n{END}{tail}")
    print(f"injected table into {readme}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", default="full")
    ap.add_argument("--out", default=f"{WK}/reports")
    ap.add_argument("--readme", default=f"{WK}/README.md",
                    help="inject the table between the generated-table markers")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    attrib = load_attribution(f"{WK}/metadata/attribution.tsv")
    data = [build_rows(s, a.level, attrib) for s in SAMPLES]
    mdp = write_table(data, a.level, a.out)
    for d in data:
        draw_tree(d, a.level, f"{a.out}/figures")
    if a.readme and a.level == "full":
        inject_readme(mdp, a.readme)


if __name__ == "__main__":
    main()
