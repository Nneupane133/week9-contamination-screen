#!/usr/bin/env python3
"""
STEP 4: turn alignments + classifications into the detection table.

Computes, per candidate organism:
    read fraction  = reads uniquely assigned to it / total reads in the library
    depth          = mean per-base depth across its ENTIRE reference
    breadth        = fraction of its reference bases covered at >= 1x

and three auxiliary statistics that exist specifically to catch false positives:

    mapq_survival  = reads at MAPQ>=20 / reads at MAPQ>=0
        A real organism's reads are unique to it, so they keep a high MAPQ in a
        competitive reference. A hit that is really cross-mapping from the host
        or from a conserved region collapses when ambiguous reads are removed.
        Low survival means "these reads were never really yours".

    gini           = Gini coefficient of per-base depth (0 = perfectly even)
        Breadth and mean depth cannot distinguish an even genome-wide 0.5x from
        a 25x spike over 2 % of the genome, yet those mean opposite things.
        A genuine organism gives low-to-moderate Gini; a conserved-locus or
        database-error artifact gives Gini near 1.

    top1pct_share  = share of all aligned bases falling in the most-covered 1 %
                     of the reference
        A blunt, easily-explained version of the same idea, for the write-up.

The four bins and why the thresholds sit where they do are documented in
assign_bin() below.
"""
import argparse, gzip, os, re, sys
from collections import defaultdict

HOSTS = ("Gallus_gallus_GRCg7b", "Aedes_aegypti_AaegL5")


def read_coverage(path):
    """samtools coverage -> {source: dict of aggregated stats}"""
    agg = defaultdict(lambda: dict(length=0, covbases=0, numreads=0, depth_sum=0.0,
                                   mapq_w=0.0, contigs=0))
    if not os.path.exists(path):
        return agg
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            rname = f[0]
            src = rname.split("::", 1)[0] if "::" in rname else rname
            length = int(f[2])              # endpos == sequence length
            covb = int(f[4])
            nreads = int(f[3])
            meandepth = float(f[6])
            meanmapq = float(f[8])
            a = agg[src]
            a["length"] += length
            a["covbases"] += covb
            a["numreads"] += nreads
            a["depth_sum"] += meandepth * length
            a["mapq_w"] += meanmapq * nreads
            a["contigs"] += 1
    return agg


def gini_and_spike(depth_path, sources_wanted):
    """Per-base depth -> Gini coefficient and top-1 % share, per source."""
    vals = defaultdict(list)
    if not os.path.exists(depth_path):
        return {}, {}
    op = gzip.open if depth_path.endswith(".gz") else open
    with op(depth_path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3:
                continue
            src = f[0].split("::", 1)[0]
            if src not in sources_wanted:
                continue
            vals[src].append(int(f[2]))
    gini, spike = {}, {}
    for src, v in vals.items():
        if not v:
            continue
        v.sort()
        n = len(v)
        total = sum(v)
        if total == 0:
            gini[src], spike[src] = float("nan"), 0.0
            continue
        # Gini over the per-base depth vector, zeros included.
        cum = 0
        weighted = 0
        for i, x in enumerate(v, start=1):
            weighted += i * x
        gini[src] = (2 * weighted) / (n * total) - (n + 1) / n
        k = max(1, n // 100)
        spike[src] = sum(v[-k:]) / total
    return gini, spike


def read_total_reads(flagstat):
    if not os.path.exists(flagstat):
        return None
    with open(flagstat) as fh:
        for line in fh:
            m = re.match(r"(\d+) \+ \d+ in total", line)
            if m:
                return int(m.group(1))
    return None


def read_bracken(path):
    """Bracken species table -> {name: (new_est_reads, fraction)}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_name = hdr.index("name"); i_new = hdr.index("new_est_reads")
            i_frac = hdr.index("fraction_total_reads")
        except ValueError:
            return out
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_name, i_new, i_frac):
                continue
            out[f[i_name]] = (int(f[i_new]), float(f[i_frac]))
    return out


def read_minimizer_report(path):
    """kraken2 --report-minimizer-data -> {taxname: (taxon_reads, distinct_minimizers)}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            try:
                taxon_reads = int(f[2]); distinct = int(f[3])
            except ValueError:
                continue
            name = f[7].strip()
            out[name] = (taxon_reads, distinct)
    return out


def assign_bin(read_frac, depth, breadth, mapq_survival, gini, mapped_reads,
               checkv_completeness=None):
    """
    noise / trace / consequential / biological.

    The ordering matters: the credibility filters run FIRST, so that a hit can
    never be promoted to a higher bin on the strength of numbers that the
    credibility checks have already shown to be artefactual.

    noise -- not credible.
        Fewer than 50 unique reads, or breadth under 0.5 %, or most of the
        reads evaporate when ambiguous alignments are dropped
        (mapq_survival < 0.25), or essentially all coverage sits in one spike
        (gini > 0.95 while breadth < 5 %). Any of these means the evidence is
        consistent with cross-mapping or a database artifact.

    trace -- real but harmless.
        Survives the credibility filters, but read fraction < 0.1 %, depth
        < 1x and breadth < 10 %. Something really is there; there is not
        enough of it to change any downstream result.

    consequential -- would change an analysis.
        read fraction >= 0.1 %, or depth >= 1x, or breadth >= 10 %.
        At this level the organism can shift expression estimates, pull
        mis-assembled contigs into an assembly, or seed false variant calls.

    biological -- worth reporting in its own right.
        breadth >= 50 % at depth >= 1x (a near-complete genome is present), or
        a viral contig reconstructed to >= 90 % CheckV completeness. At that
        point the detection is a finding about the sample's biology, not a
        quality-control footnote.
    """
    if mapped_reads < 50:
        return "noise", "fewer than 50 uniquely-mapped reads"
    if breadth < 0.005:
        return "noise", f"breadth {breadth:.3%} below 0.5% of the reference"
    if mapq_survival is not None and mapq_survival < 0.25:
        return "noise", (f"only {mapq_survival:.0%} of reads survive MAPQ>=20; "
                         "alignments are ambiguous with host/decoy")
    if gini is not None and gini == gini and gini > 0.95 and breadth < 0.05:
        return "noise", (f"coverage concentrated in a spike (Gini {gini:.2f}, "
                         f"breadth {breadth:.2%})")

    if checkv_completeness is not None and checkv_completeness >= 90:
        return "biological", f"viral genome reconstructed to {checkv_completeness:.0f}% completeness"
    if breadth >= 0.50 and depth >= 1.0:
        return "biological", f"near-complete genome present (breadth {breadth:.1%}, depth {depth:.2f}x)"

    if read_frac >= 0.001 or depth >= 1.0 or breadth >= 0.10:
        return "consequential", (f"read fraction {read_frac:.3%}, depth {depth:.2f}x, "
                                 f"breadth {breadth:.1%}")
    return "trace", (f"credible but low: read fraction {read_frac:.4%}, "
                     f"depth {depth:.3f}x, breadth {breadth:.2%}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--level", default="test")
    ap.add_argument("--wk", default="/home/aubnxn001/week9")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    cov = f"{a.wk}/results/coverage/{a.level}/{a.sample}"
    k2 = f"{a.wk}/results/kraken2/{a.level}/{a.sample}"

    mq20 = read_coverage(f"{cov}.coverage.mq20.tsv")
    mq0 = read_coverage(f"{cov}.coverage.mq0.tsv")
    total_reads = read_total_reads(f"{cov}.flagstat.txt")
    bracken = read_bracken(f"{k2}.bracken.species.tsv")
    minim = read_minimizer_report(f"{k2}.minimizer.report")

    sources = set(mq20) | set(mq0)
    nonhost = {s for s in sources if not s.startswith(HOSTS)}
    gini, spike = gini_and_spike(f"{cov}.depth.nonhost.tsv.gz", nonhost)

    rows = []
    for src in sorted(sources):
        a20 = mq20.get(src, dict(length=0, covbases=0, numreads=0, depth_sum=0.0, mapq_w=0.0, contigs=0))
        a0 = mq0.get(src, a20)
        L = a20["length"] or a0["length"]
        if L == 0:
            continue
        breadth = a20["covbases"] / L
        depth = a20["depth_sum"] / L
        reads20, reads0 = a20["numreads"], a0["numreads"]
        surv = (reads20 / reads0) if reads0 else None
        frac = (reads20 / total_reads) if total_reads else float("nan")
        g = gini.get(src)
        is_host = src.startswith(HOSTS)
        if is_host:
            bin_, why = "host", "reference organism for this library"
        else:
            bin_, why = assign_bin(frac, depth, breadth, surv, g, reads20)
        rows.append(dict(
            sample=a.sample, source=src, role="host" if is_host else "nonhost",
            ref_len=L, contigs=a20["contigs"],
            reads_mq20=reads20, reads_mq0=reads0,
            read_fraction=frac, depth=depth, breadth=breadth,
            mapq_survival=surv if surv is not None else float("nan"),
            mean_mapq=(a20["mapq_w"] / reads20) if reads20 else float("nan"),
            gini=g if g is not None else float("nan"),
            top1pct_share=spike.get(src, float("nan")),
            bin=bin_, rationale=why,
        ))

    rows.sort(key=lambda r: (-r["read_fraction"] if r["read_fraction"] == r["read_fraction"] else 0))
    out = a.out or f"{a.wk}/results/coverage/{a.level}/{a.sample}.detections.tsv"
    cols = ["sample", "source", "role", "ref_len", "contigs", "reads_mq20", "reads_mq0",
            "read_fraction", "depth", "breadth", "mapq_survival", "mean_mapq",
            "gini", "top1pct_share", "bin", "rationale"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"total reads in library (flagstat): {total_reads}")
    print(f"wrote {out}  ({len(rows)} references)")
    for r in rows:
        if r["role"] == "nonhost" and r["bin"] != "noise":
            print(f"  {r['bin']:>14}  {r['source']:<42} frac={r['read_fraction']:.5%} "
                  f"depth={r['depth']:.3f}x breadth={r['breadth']:.2%}")


if __name__ == "__main__":
    main()
