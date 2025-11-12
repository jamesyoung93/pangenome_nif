#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO enrichment + cellular map for diazotrophy predictors.

Inputs (current directory):
  - llm_context_features.csv  (from step 6; preferred)  [gene_family, direction/effect_size, GO, ...]
  - narrative_top_features.csv (from step 8; optional)  [subset with GO carried over]
Optional:
  - go-basic.obo (for GO ID -> name/namespace). If missing, labels use GO:IDs.

Outputs:
  - enriched_go_positive.csv / enriched_go_negative.csv
  - go_cell_map.svg
"""

import argparse, math, re, os, sys, random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------- small helpers ----------

def read_go_obo(path="go-basic.obo"):
    """Return dicts: id->name, id->namespace if OBO is present; else empty."""
    if not os.path.exists(path):
        return {}, {}
    name, ns = {}, {}
    cur_id, cur_name, cur_ns = None, None, None
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("[Term]"):
                if cur_id:
                    if cur_name: name[cur_id] = cur_name
                    if cur_ns: ns[cur_id] = cur_ns
                cur_id = cur_name = cur_ns = None
            elif line.startswith("id: GO:"):
                cur_id = line.strip().split("id: ")[1]
            elif line.startswith("name: "):
                cur_name = line.strip().split("name: ",1)[1]
            elif line.startswith("namespace: "):
                cur_ns = line.strip().split("namespace: ",1)[1]
        if cur_id:
            if cur_name: name[cur_id] = cur_name
            if cur_ns: ns[cur_id] = cur_ns
    return name, ns

def parse_go_field(val):
    """Split a GO field like 'GO:0005524; GO:0055114' into a set of GO IDs."""
    if pd.isna(val) or not str(val).strip():
        return set()
    tokens = re.split(r"[;, \t]+", str(val))
    return {t for t in tokens if t.startswith("GO:")}

def bh_fdr(pvals):
    """Benjamini-Hochberg; returns array of q-values."""
    p = np.array(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, int); ranks[order] = np.arange(1, n+1)
    q = p * n / np.maximum(ranks, 1)
    # enforce monotonicity
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty(n, float); out[order] = np.clip(q_sorted, 0, 1)
    return out

def fisher_right_tail(a, b, c, d):
    """
    One-sided Fisher exact p-value (enrichment in group A): a/(a+b) vs c/(c+d).
    Implemented via hypergeometric tail; no SciPy required.
    Table:
        group = Positive     Not-Positive
    GO       a               c
    not GO   b               d
    """
    from math import comb
    N = a + b + c + d
    K = a + c           # GO in universe
    n = a + b           # size of Positive group
    # P[X >= a] where X ~ Hypergeom(N, K, n)
    def hypergeom_pmf(x):
        return comb(K, x) * comb(N - K, n - x) / comb(N, n)
    p = 0.0
    for x in range(a, min(K, n) + 1):
        p += hypergeom_pmf(x)
    return min(1.0, p)

# Compartment classifier from simple keywords (robust without OBO)
COMPARTMENT_RULES = [
    ("membrane", ["membrane","transmembrane","periplasm","secretion","transport","porin","flagell"]),
    ("dna_region", ["dna","replication","recombination","repair","transcription"]),
    ("ribosome", ["ribosom","translation","rrna","trna","elongation"]),
    ("redox", ["oxidoreduct","electron","respirat","dehydrogenase","cytochrome","ferredoxin","flavodoxin","quinone"]),
    ("metabolism", ["metabolic","biosynth","catabolic","cofactor","sulfur","iron","molybden","heme"]),
    ("stress", ["stress","oxygen","reactive","peroxid","superoxide","catalase"]),
    ("other", [""])
]

def term_compartment(go_id, go_name):
    name = (go_name or "").lower()
    for comp, keys in COMPARTMENT_RULES:
        if any(k in name for k in keys):
            return comp
    return "other"

def load_inputs(prefer_narrative=None):
    # Prefer the richer LLM merge
    if prefer_narrative and os.path.exists(prefer_narrative):
        df = pd.read_csv(prefer_narrative)
    elif os.path.exists("llm_context_features.csv"):
        df = pd.read_csv("llm_context_features.csv")
    else:
        raise SystemExit("No llm_context_features.csv or narrative_top_features.csv found.")
    # Expected columns from Steps 5/6/8
    # gene_family | direction | effect_size | GO | (importance...)  (Step 6/8 produce these)  # 
    for c in ("gene_family","GO"):
        if c not in df.columns:
            raise SystemExit(f"Missing required column {c}")
    # Normalize direction from effect_size if needed
    if "direction" not in df.columns and "effect_size" in df.columns:
        df["direction"] = np.where(df["effect_size"]>0, "Positive",
                            np.where(df["effect_size"]<0, "Negative", "Neutral"))
    return df

def make_gf_go_table(df):
    # One row per GF with set of GO IDs
    gf_go = (df[["gene_family","GO"]].copy())
    gf_go["GOset"] = gf_go["GO"].map(parse_go_field)
    gf_go = (gf_go.groupby("gene_family", as_index=False)
                  .agg(GOset=("GOset", lambda s: set().union(*s))))
    return gf_go

def compute_enrichment(df, gf_go, direction_label="Positive", min_count=2, q_cut=0.1):
    pos_gf = set(df.loc[df["direction"]==direction_label, "gene_family"])
    all_gf = set(gf_go["gene_family"])
    # Universe: GFs with any GO
    gf_has_go = {r.gene_family for r in gf_go.itertuples() if len(r.GOset)>0}
    pos_gf = pos_gf & gf_has_go
    U = gf_go[gf_go["gene_family"].isin(gf_has_go)].copy()

    # Build GO term -> set of GFs carrying it
    go_to_gf = {}
    for r in U.itertuples():
        for go in r.GOset:
            go_to_gf.setdefault(go, set()).add(r.gene_family)

    rows = []
    m = len(pos_gf)                 # size of positive set
    M = len(gf_has_go)              # universe size
    for go, gfs in go_to_gf.items():
        a = len(gfs & pos_gf)       # GO n positive
        b = m - a                   # not GO in positive
        c = len(gfs) - a            # GO in not-positive
        d = (M - m) - c             # not GO in not-positive
        if a < min_count:
            continue
        p = fisher_right_tail(a,b,c,d)
        rows.append((go, a, len(gfs), p))
    if not rows:
        return pd.DataFrame(columns=["GO","pos_hits","term_size","p","q","direction"])
    out = pd.DataFrame(rows, columns=["GO","pos_hits","term_size","p"])
    out["q"] = bh_fdr(out["p"].values)
    out["direction"] = direction_label
    out = out.sort_values(["q","p","pos_hits"], ascending=[True, True, False]).reset_index(drop=True)
    return out[(out["q"] <= q_cut) | (out["pos_hits"]>=3)].copy()

def draw_cell_map(pos_df, neg_df, go_name_map, outfile="go_cell_map.svg"):
    # Build plotting table
    pos_df = pos_df.copy(); neg_df = neg_df.copy()
    pos_df["name"] = [go_name_map.get(g, g) for g in pos_df["GO"]]
    neg_df["name"] = [go_name_map.get(g, g) for g in neg_df["GO"]]
    pos_df["comp"] = [term_compartment(g, n) for g,n in zip(pos_df["GO"], pos_df["name"])]
    neg_df["comp"] = [term_compartment(g, n) for g,n in zip(neg_df["GO"], neg_df["name"])]
    pos_df["size"] = -np.log10(np.clip(pos_df["q"].values, 1e-300, 1.0)) * 150
    neg_df["size"] = -np.log10(np.clip(neg_df["q"].values, 1e-300, 1.0)) * 150

    rng = np.random.default_rng(42)
    # Canvas
    fig = plt.figure(figsize=(10, 7))
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw a simple Gram-negative-like outline (outer and inner envelope)
    # Outer ellipse
    outer = plt.Circle((0,0), 4.5, fill=False, linewidth=2)
    ax.add_artist(outer)
    # Inner membrane (ring)
    inner = plt.Circle((0,0), 3.8, fill=False, linewidth=1, linestyle="--")
    ax.add_artist(inner)
    # "Nucleoid" region
    nuc = plt.Circle((0,0), 2.0, fill=False, linewidth=1)
    ax.add_artist(nuc)

    # Regions and bounding boxes for scattering
    boxes = {
        "membrane":  (-4.5, -4.5, 9.0, 1.4),   # top arc band
        "redox":     (-3.5,  1.2,  7.0, 1.4),
        "metabolism":(-3.2, -1.0,  6.4, 1.6),
        "dna_region":(-2.0, -0.6,  4.0, 1.2),
        "ribosome":  (-3.4, -2.8,  6.8, 1.4),
        "stress":    (-3.8,  2.8,  7.6, 1.2),
        "other":     (-1.5, -4.0,  3.0, 1.0),
    }

    def scatter_terms(df, marker):
        coords = []
        for comp in boxes:
            sub = df[df["comp"]==comp]
            if not len(sub): continue
            x0,y0,w,h = boxes[comp]
            xs = rng.uniform(x0, x0+w, size=len(sub))
            ys = rng.uniform(y0, y0+h, size=len(sub))
            coords.append((sub, xs, ys))
        for sub, xs, ys in coords:
            ax.scatter(xs, ys, s=sub["size"].values, marker=marker, alpha=0.7)
            # labels
            for x,y,name in zip(xs, ys, sub["name"].values):
                ax.text(x, y, name, fontsize=8, ha="center", va="center")

    # Positive as circles, Negative as squares (no explicit colors)
    if len(pos_df):
        scatter_terms(pos_df, marker="o")
    if len(neg_df):
        scatter_terms(neg_df, marker="s")

    # Legend proxies (default colors)
    ax.scatter([], [], s=80, marker="o", label="Positive-enriched")
    ax.scatter([], [], s=80, marker="s", label="Negative-enriched")
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    fig.suptitle("Enriched GO terms mapped onto a stylized cell", fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--use", default="llm_context_features.csv",
                    help="Which table to use: llm_context_features.csv or narrative_top_features.csv")
    ap.add_argument("--q", type=float, default=0.10, help="FDR cutoff")
    ap.add_argument("--min-count", type=int, default=2, help="Min positive hits for a GO term")
    args = ap.parse_args()

    df = load_inputs(prefer_narrative=args.use if os.path.exists(args.use) else None)
    gf_go = make_gf_go_table(df)

    pos = compute_enrichment(df, gf_go, "Positive", min_count=args.min_count, q_cut=args.q)
    neg = compute_enrichment(df, gf_go, "Negative", min_count=args.min_count, q_cut=args.q)

    # Optional GO names if go-basic.obo is available
    go_name_map, go_ns_map = read_go_obo("go-basic.obo")
    pos["name"] = [go_name_map.get(g, g) for g in pos["GO"]]
    neg["name"] = [go_name_map.get(g, g) for g in neg["GO"]]

    pos.to_csv("enriched_go_positive.csv", index=False)
    neg.to_csv("enriched_go_negative.csv", index=False)

    # Draw cell map
    draw_cell_map(pos, neg, go_name_map, outfile="go_cell_map.svg")

    # Brief console summary
    def headlines(tag, dfx):
        print(f"\nTop enriched GO ({tag}):")
        for i, r in dfx.head(12).iterrows():
            print(f"  {r['name']} [{r['GO']}]  hits={r['pos_hits']}/{r['term_size']}  p={r['p']:.2e}  q={r['q']:.2e}")
    if len(pos): headlines("Positive", pos)
    if len(neg): headlines("Negative", neg)

if __name__ == "__main__":
    main()
