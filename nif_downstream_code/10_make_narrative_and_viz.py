#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make a module-based narrative + figure pack from top_features_narrative.tsv

Fixes vs prior version:
- Uses explicit color mapping + Line2D proxy legend (no empty scatter calls).
- Avoids matplotlib 'c' parsing edge cases on some builds.

Usage:
  python3 10_make_narrative_and_viz.py --table results/baseline/top_features_narrative.tsv
  python3 10_make_narrative_and_viz.py --table results/baseline/top_features_narrative.tsv --top-k 60 --outdir results/baseline/narrative_viz
"""

import argparse
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


REQUIRED_COLS = ["gene_family", "consensus_score", "models_supporting", "product_mode"]
OPTIONAL_COLS = [
    "effect_size", "n_genomes_with_feature", "module_guess",
    "diaz_rate_present", "diaz_rate_absent",
    "uniprot_ids", "representative_protein"
]


def _ensure_cols(df: pd.DataFrame):
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise RuntimeError(f"Input TSV missing required columns: {missing}")
    for c in OPTIONAL_COLS:
        if c not in df.columns:
            df[c] = np.nan
    if "module_guess" not in df.columns or df["module_guess"].isna().all():
        df["module_guess"] = "Unclassified/other"
    return df


def _read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = _ensure_cols(df)
    df["consensus_score"] = pd.to_numeric(df["consensus_score"], errors="coerce").fillna(0.0)
    df["models_supporting"] = pd.to_numeric(df["models_supporting"], errors="coerce").fillna(0).astype(int)
    df["effect_size"] = pd.to_numeric(df["effect_size"], errors="coerce")
    df["n_genomes_with_feature"] = pd.to_numeric(df["n_genomes_with_feature"], errors="coerce")
    df["product_mode"] = df["product_mode"].fillna("").astype(str)
    df["gene_family"] = df["gene_family"].astype(str)
    df["module_guess"] = df["module_guess"].fillna("Unclassified/other").astype(str)
    return df


def _topk(df: pd.DataFrame, k: int) -> pd.DataFrame:
    df = df.sort_values(["consensus_score", "models_supporting"], ascending=[False, False]).copy()
    return df.head(k).copy()


def _savefig(fig, out_png: Path, out_pdf: Path):
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def panelA_module_composition(df: pd.DataFrame, outdir: Path):
    counts = df.groupby("module_guess")["gene_family"].count().sort_values(ascending=False)
    weights = df.groupby("module_guess")["consensus_score"].sum().sort_values(ascending=False)

    modules = list(counts.index)
    counts = counts.reindex(modules).fillna(0)
    weights = weights.reindex(modules).fillna(0)

    fig = plt.figure(figsize=(10, 5.5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.barh(range(len(modules)), counts.values)
    ax1.set_yticks(range(len(modules)))
    ax1.set_yticklabels(modules)
    ax1.invert_yaxis()
    ax1.set_xlabel("Number of gene families")
    ax1.set_title("Panel A1: Module composition (counts)")

    ax2.barh(range(len(modules)), weights.values)
    ax2.set_yticks(range(len(modules)))
    ax2.set_yticklabels([])
    ax2.invert_yaxis()
    ax2.set_xlabel("Sum of consensus importance")
    ax2.set_title("Panel A2: Module composition (importance-weighted)")

    fig.tight_layout()
    _savefig(fig, outdir / "panelA_module_composition.png", outdir / "panelA_module_composition.pdf")


def panelB_effect_vs_consensus(df: pd.DataFrame, outdir: Path, annotate_n: int = 12):
    # x: effect size (if missing, set to 0)
    x = df["effect_size"].copy()
    if x.isna().all():
        x = pd.Series(np.zeros(len(df)), index=df.index)
    else:
        x = x.fillna(0.0)

    # y: consensus score
    y = df["consensus_score"].copy().fillna(0.0)

    # bubble size: prevalence (if missing, use 1)
    size = df["n_genomes_with_feature"].copy()
    if size.isna().all():
        size = pd.Series(np.ones(len(df)), index=df.index)
    size = size.fillna(1.0)

    # gentle size scaling (log), with clamping for readability
    size2 = (np.log10(size.values + 1.0) + 0.2) * 180.0
    size2 = np.clip(size2, 35.0, 900.0)

    # categorical colors (explicit mapping avoids Matplotlib c-parsing issues)
    modules = df["module_guess"].astype("category")
    codes = modules.cat.codes.values
    cmap = plt.get_cmap("tab20")
    point_colors = [cmap(int(c) % 20) for c in codes]

    fig = plt.figure(figsize=(9, 6.5))
    ax = fig.add_subplot(111)

    ax.scatter(
        x.values,
        y.values,
        s=size2,
        color=point_colors,
        alpha=0.8,
        edgecolors="none"
    )

    ax.axvline(0.0, linestyle="--", linewidth=1)
    ax.set_xlabel("Effect size: P(diaz|present) - P(diaz|absent)")
    ax.set_ylabel("Consensus importance (mean normalized)")
    ax.set_title("Panel B: Directionality vs importance (bubble size ~ prevalence)")

    # annotate top points by consensus score
    top_annot = df.sort_values("consensus_score", ascending=False).head(annotate_n)
    for _, r in top_annot.iterrows():
        ax.text(
            float(r["effect_size"]) if pd.notna(r["effect_size"]) else 0.0,
            float(r["consensus_score"]),
            r["gene_family"],
            fontsize=8,
            ha="left",
            va="bottom"
        )

    # legend using Line2D proxies (robust on all Matplotlib builds)
    handles = []
    labels = []
    for code, name in enumerate(modules.cat.categories):
        handles.append(Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            markersize=7,
            markerfacecolor=cmap(int(code) % 20),
            markeredgecolor="none",
            alpha=0.8
        ))
        labels.append(name)

    ax.legend(handles, labels, title="Module", loc="best", fontsize=8, title_fontsize=9, frameon=True)

    fig.tight_layout()
    _savefig(fig, outdir / "panelB_effect_vs_consensus.png", outdir / "panelB_effect_vs_consensus.pdf")


def panelC_module_network(df: pd.DataFrame, outdir: Path, max_per_module: int = 12):
    df2 = df.sort_values(["module_guess", "consensus_score"], ascending=[True, False]).copy()
    df2 = df2.groupby("module_guess").head(max_per_module).copy()

    modules = list(df2["module_guess"].unique())
    families = list(df2["gene_family"].unique())

    fig = plt.figure(figsize=(11, max(6, 0.35 * (len(families) + len(modules)))))
    ax = fig.add_subplot(111)
    ax.set_axis_off()

    mod_y = np.linspace(0.95, 0.05, len(modules)) if len(modules) > 1 else np.array([0.5])
    fam_y = np.linspace(0.98, 0.02, len(families)) if len(families) > 1 else np.array([0.5])

    mod_pos = {m: (0.08, mod_y[i]) for i, m in enumerate(modules)}
    fam_pos = {g: (0.92, fam_y[i]) for i, g in enumerate(families)}

    for m, (x, y) in mod_pos.items():
        ax.text(x, y, m, ha="left", va="center", fontsize=10, fontweight="bold")
    for g, (x, y) in fam_pos.items():
        ax.text(x, y, g, ha="right", va="center", fontsize=8)

    s = df2["consensus_score"].values
    if np.nanmax(s) <= 0:
        s_scaled = np.ones_like(s)
    else:
        s_scaled = 0.5 + 3.0 * (s / (np.nanmax(s) + 1e-12))

    for (i, row), lw in zip(df2.iterrows(), s_scaled):
        m = row["module_guess"]
        g = row["gene_family"]
        x0, y0 = mod_pos[m]
        x1, y1 = fam_pos[g]
        ax.plot([x0 + 0.18, x1 - 0.18], [y0, y1], linewidth=float(lw), alpha=0.35)

    ax.set_title("Panel C: Module ? gene-family network (top families per module)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    _savefig(fig, outdir / "panelC_module_network.png", outdir / "panelC_module_network.pdf")


def panelD_top_by_module(df: pd.DataFrame, outdir: Path, top_per_module: int = 6):
    df2 = df.sort_values(["module_guess", "consensus_score"], ascending=[True, False]).copy()
    df2 = df2.groupby("module_guess").head(top_per_module).copy()

    modules = list(df2["module_guess"].unique())

    fig = plt.figure(figsize=(12, max(5, 0.55 * len(modules))))
    ax = fig.add_subplot(111)
    ax.set_axis_off()
    ax.set_title("Panel D: Top families per module (product-mode + directionality)", fontsize=12, fontweight="bold")

    y = 0.95
    line_h = 0.055
    for m in modules:
        sub = df2[df2["module_guess"] == m].copy()
        ax.text(0.02, y, m, fontsize=11, fontweight="bold", ha="left", va="top")
        y -= line_h

        for _, r in sub.iterrows():
            prod = (r["product_mode"] or "").strip() or "NA"
            eff = r["effect_size"]
            eff_s = f"{eff:.2f}" if pd.notna(eff) else "NA"
            it = f"{r['gene_family']}: {prod} (?={eff_s})"
            wrapped = textwrap.fill(it, width=120)
            ax.text(0.04, y, wrapped, fontsize=9, ha="left", va="top")
            y -= line_h * (1 + wrapped.count("\n"))

        y -= line_h * 0.35

    fig.tight_layout()
    _savefig(fig, outdir / "panelD_top_by_module.png", outdir / "panelD_top_by_module.pdf")


def write_auto_narrative(df: pd.DataFrame, out_md: Path):
    n = len(df)
    n_mod = df["module_guess"].nunique()
    pos_frac = float((df["effect_size"] > 0).mean()) if df["effect_size"].notna().any() else np.nan

    lines = []
    lines.append("# Auto-generated narrative (from top_features_narrative.tsv)\n")
    lines.append(f"- Top-K families analyzed: **{n}**")
    lines.append(f"- Modules represented: **{n_mod}**")
    if not np.isnan(pos_frac):
        lines.append(f"- Fraction of families with positive directionality (effect_size > 0): **{pos_frac:.2f}**")
    lines.append("")

    for m, sub in df.groupby("module_guess", sort=False):
        sub = sub.sort_values("consensus_score", ascending=False)
        lines.append(f"## {m}\n")
        lines.append(f"Families: **{len(sub)}**")
        if sub["effect_size"].notna().any():
            lines.append(f"Median effect size: **{sub['effect_size'].median():.2f}**")
        lines.append("")
        lines.append("Top examples:")
        for _, r in sub.head(6).iterrows():
            prod = (r["product_mode"] or "").strip() or "NA"
            eff = r["effect_size"]
            eff_s = f"{eff:.2f}" if pd.notna(eff) else "NA"
            sup = int(r["models_supporting"]) if pd.notna(r["models_supporting"]) else 0
            lines.append(f"- **{r['gene_family']}**: {prod} (?={eff_s}, models={sup})")
        lines.append("")

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="Path to top_features_narrative.tsv")
    ap.add_argument("--top-k", type=int, default=60, help="Subset the table to top K by consensus_score")
    ap.add_argument("--outdir", default="", help="Output directory (default: <table_dir>/narrative_viz)")
    args = ap.parse_args()

    table_path = Path(args.table)
    if not table_path.exists():
        raise FileNotFoundError(f"Missing input table: {table_path}")

    df = _read_table(table_path)
    df = _topk(df, args.top_k)

    outdir = Path(args.outdir) if args.outdir else (table_path.parent / "narrative_viz")
    outdir.mkdir(parents=True, exist_ok=True)

    panelA_module_composition(df, outdir)
    panelB_effect_vs_consensus(df, outdir)
    panelC_module_network(df, outdir)
    panelD_top_by_module(df, outdir)
    write_auto_narrative(df, outdir / "auto_narrative.md")

    print("\nDONE")
    print(f"Outputs written to: {outdir.resolve()}")
    print("Panels:")
    print(" - panelA_module_composition.(png|pdf)")
    print(" - panelB_effect_vs_consensus.(png|pdf)")
    print(" - panelC_module_network.(png|pdf)")
    print(" - panelD_top_by_module.(png|pdf)")
    print(" - auto_narrative.md\n")


if __name__ == "__main__":
    main()
