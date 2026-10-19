#!/usr/bin/env python3
"""
STEP 6: merge every line of evidence into one detection table.

Columns come from five independent places, which is the point: a detection is
only as good as the agreement between methods that fail differently.

    sendsketch (RefSeq + nt)  presence / ANI, database-free
    Kraken2                   direct read count + DISTINCT MINIMIZERS
    Bracken                   species-level re-estimate (and see the caveat it exposes)
    competitive mapping       read fraction, depth, breadth, Gini, MAPQ survival
    host-competition test     is the sequence integrated in the host?
    control samples           does it track the lane, the batch, or the sample?
    assembly / CheckV         is a near-complete genome actually reconstructable?
"""
import argparse, csv, glob, gzip, os, re, sys
from collections import defaultdict

WK = "/home/aubnxn001/week9"
CONTROLS = {
    "BYJ10_RNA": [("AYJ9_RNA", "same lane"), ("SYJ7_RNA", "diff flowcell"),
                  ("SYJ10_RNA", "diff flowcell")],
    "SRR39197041_DNA": [("SRR39197043_DNA", "same lane")],
}


def load_manifest(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rows[f[0]] = dict(name=f[0], role=f[1], taxid=f[2], source=f[3],
                              id=f[4], restrict=f[5])
    return rows


def kraken_by_taxid(path):
    """-> {taxid: (clade_reads, direct_reads, distinct_minimizers, rank, name)}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            try:
                out[f[6].strip()] = (int(f[1]), int(f[2]), int(f[3]), f[5], f[7].strip())
            except ValueError:
                continue
    return out


def kraken_total(path):
    """Total reads processed, from the kraken2 stderr log."""
    log = path.replace(".minimizer.report", ".kraken2.log")
    if os.path.exists(log):
        m = re.search(r"(\d+) sequences \(", open(log).read())
        if m:
            return int(m.group(1))
    return None


def sendsketch_by_taxid(path):
    """-> {taxid: (ANI, matches, unique, complt)}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 33 or f[0] in ("WKID", ""):
                continue
            try:
                out[f[20].strip()] = (f[2], int(f[14]), int(f[15]), f[4])
            except (ValueError, IndexError):
                continue
    return out


def load_detections(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["source"]] = r
    return out


def checkv_best(sample, level):
    """Best (most complete) CheckV record per contig, keyed by contig id."""
    p = f"{WK}/results/assembly/{level}/{sample}/checkv/quality_summary.tsv"
    rows = []
    if not os.path.exists(p):
        return rows
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                comp = float(r.get("completeness") or "nan")
            except ValueError:
                comp = float("nan")
            r["_completeness"] = comp
            rows.append(r)
    rows.sort(key=lambda r: (r["_completeness"] if r["_completeness"] == r["_completeness"] else -1),
              reverse=True)
    return rows


def fmt(x, pct=False, nd=4):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "-"
    if v != v:
        return "-"
    return f"{v*100:.{nd}f}%" if pct else f"{v:.{nd}g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", default="full")
    ap.add_argument("--out", default=f"{WK}/reports")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    man = load_manifest(f"{WK}/metadata/candidates.tsv")
    all_rows = []

    for sample in ("BYJ10_RNA", "SRR39197041_DNA"):
        det = load_detections(f"{WK}/results/coverage/{a.level}/{sample}.detections.tsv")
        k2p = f"{WK}/results/kraken2/{a.level}/{sample}.minimizer.report"
        k2 = kraken_by_taxid(k2p)
        ktot = kraken_total(k2p)
        ss_rs = sendsketch_by_taxid(f"{WK}/results/sendsketch/{a.level}/{sample}.refseq.tsv")
        ss_nt = sendsketch_by_taxid(f"{WK}/results/sendsketch/{a.level}/{sample}.nt.tsv")
        ctrl = {}
        for cname, rel in CONTROLS.get(sample, []):
            cp = f"{WK}/results/kraken2/controls/{cname}.minimizer.report"
            ctrl[cname] = (kraken_by_taxid(cp), kraken_total(cp), rel)

        for name, m in man.items():
            if m["restrict"] not in ("both", sample):
                continue
            if m["role"] == "host":
                continue
            d = det.get(name, {})
            tid = m["taxid"]
            kc, kd, kmin, krank, kname = k2.get(tid, (0, 0, 0, "", ""))
            row = dict(
                sample=sample, organism=name, role=m["role"], taxid=tid,
                # --- screen A
                sendsketch_refseq_ani=ss_rs.get(tid, ("-",))[0],
                sendsketch_refseq_matches=ss_rs.get(tid, ("-", "-"))[1] if tid in ss_rs else "-",
                sendsketch_nt_ani=ss_nt.get(tid, ("-",))[0],
                # --- screen B
                kraken_clade_reads=kc, kraken_direct_reads=kd,
                kraken_distinct_minimizers=kmin,
                kraken_reads_per_minimizer=(f"{kd/kmin:.4g}" if kmin else "-"),
                kraken_frac=(f"{kd/ktot:.3e}" if ktot else "-"),
                # --- measurement
                ref_len=d.get("ref_len", "-"), measured_at=d.get("measured_at", "-"),
                reads_mq20=d.get("reads_mq20", "-"), reads_mq0=d.get("reads_mq0", "-"),
                read_fraction=d.get("read_fraction", ""), depth=d.get("depth", ""),
                breadth=d.get("breadth", ""), mapq_survival=d.get("mapq_survival", ""),
                host_competing=d.get("host_competing", ""), gini=d.get("gini", ""),
                top1pct_share=d.get("top1pct_share", ""),
                bin=d.get("bin", "not-measured"), rationale=d.get("rationale", ""),
            )
            for cname, (ck, ctotal, rel) in ctrl.items():
                cd = ck.get(tid, (0, 0, 0, "", ""))[1]
                row[f"ctrl_{cname}_reads"] = cd
                row[f"ctrl_{cname}_frac"] = f"{cd/ctotal:.3e}" if ctotal else "-"
            all_rows.append(row)

    cols = []
    for r in all_rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    tsv = f"{a.out}/detection_table.{a.level}.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="-")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"wrote {tsv}  ({len(all_rows)} rows)")

    # Compact markdown view, non-noise first.
    order = {"biological": 0, "consequential": 1, "trace": 2, "noise": 3, "not-measured": 4}
    all_rows.sort(key=lambda r: (r["sample"], order.get(r["bin"], 9),
                                 -float(r["read_fraction"] or 0)))
    md = [f"# Detection table ({a.level})", ""]
    for sample in ("BYJ10_RNA", "SRR39197041_DNA"):
        md += [f"## {sample}", "",
               "| organism | bin | read fraction | depth | breadth | MAPQ surv. | host-competing | Gini | kraken2 reads | distinct minimizers |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for r in all_rows:
            if r["sample"] != sample:
                continue
            md.append("| {} | **{}** | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                r["organism"].replace("_", " "), r["bin"],
                fmt(r["read_fraction"], pct=True), fmt(r["depth"]),
                fmt(r["breadth"], pct=True, nd=2), fmt(r["mapq_survival"], nd=2),
                fmt(r["host_competing"], nd=2), fmt(r["gini"], nd=2),
                r["kraken_direct_reads"], r["kraken_distinct_minimizers"]))
        md.append("")
    mdp = f"{a.out}/detection_table.{a.level}.md"
    open(mdp, "w").write("\n".join(md) + "\n")
    print(f"wrote {mdp}")

    for sample in ("BYJ10_RNA", "SRR39197041_DNA"):
        cv = checkv_best(sample, a.level)
        if cv:
            print(f"\n{sample}: top CheckV contigs")
            for r in cv[:5]:
                print(f"  {r['contig_id'][:40]:<40} len={r['contig_length']:>7} "
                      f"quality={r['checkv_quality']:<18} completeness={r.get('completeness','NA')}")


if __name__ == "__main__":
    main()
