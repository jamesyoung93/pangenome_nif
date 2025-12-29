#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Top-feature consensus + functional annotation table for biological narrative

What it does
------------
Given a pipeline run directory (e.g., results/baseline), this script:
  1) Loads all per-model feature importance files: feature_importance_*.csv
  2) Normalizes importances within each model so scales are comparable
  3) Builds a consensus ranking across models (mean normalized importance + rank aggregation)
  4) Merges directionality stats if present (feature_directionality_full.csv)
  5) Attaches functional summaries per gene family from protein_family_cds.tsv:
       - modal "product"
       - top 3 product strings with counts
       - top GO IDs with counts
       - example proteins / locus tags
  6) Optionally maps representative RefSeq protein accessions (WP_*) to UniProt via refseq_to_uniprot.tsv
  7) Writes:
       - <run_dir>/top_features_narrative.tsv
       - <run_dir>/module_summary.tsv
       - <run_dir>/narrative_outline.md  (module-binned outline to help writing)

Expected inputs (typical in your repo)
-------------------------------------
Run directory (e.g., results/baseline/) containing:
  - feature_importance_Random_Forest.csv
  - feature_importance_Gradient_Boosting.csv
  - feature_importance_Logistic_Regression.csv
  - feature_importance_XGBoost.csv (optional)
  - feature_directionality_full.csv (optional)

And mapping tables (provide via args or defaults):
  - protein_family_cds.tsv
  - refseq_to_uniprot.tsv (optional but recommended)
  - gene_family_info.csv (optional but recommended)

Usage examples
--------------
  python3 09_top_features_narrative.py --run-dir results/baseline --top-k 60
  python3 09_top_features_narrative.py --run-dir results/genus_purity_0p80 --top-k 60

If your mapping files are in /mnt/data (as in this environment), defaults will work.
On HPC, pass explicit paths to your TSVs as needed.

Outputs
-------
  results/<run>/top_features_narrative.tsv
  results/<run>/module_summary.tsv
  results/<run>/narrative_outline.md
"""

import argparse
from pathlib import Path
import glob
import re
from collections import Counter

import numpy as np
import pandas as pd


# ----------------------------
# I/O helpers
# ----------------------------

def _first_existing(paths):
    for p in paths:
        if p and Path(p).exists():
            return str(Path(p))
    return ""


def _safe_read_csv(path, sep=None):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    if sep is None:
        return pd.read_csv(path)
    return pd.read_csv(path, sep=sep, low_memory=False)


# ----------------------------
# Core: importances consensus
# ----------------------------

def load_feature_importances(run_dir: Path):
    imp_files = sorted(glob.glob(str(run_dir / "feature_importance_*.csv")))
    if not imp_files:
        raise FileNotFoundError(f"No feature_importance_*.csv found in: {run_dir}")

    frames = []
    for fn in imp_files:
        model = Path(fn).stem.replace("feature_importance_", "")
        df = _safe_read_csv(fn)
        if not {"gene_family", "importance"}.issubset(df.columns):
            continue
        df = df[["gene_family", "importance"]].copy()
        df["model"] = model
        df["importance"] = pd.to_numeric(df["importance"], errors="coerce").fillna(0.0)
        frames.append(df)

    if not frames:
        raise RuntimeError(f"Found feature_importance files, but none had expected columns in: {run_dir}")

    long = pd.concat(frames, ignore_index=True)

    # Normalize within each model (scale-free)
    long["importance_norm"] = long.groupby("model")["importance"].transform(
        lambda s: s / (s.sum() + 1e-12)
    )

    # Rank within each model: 1 = most important
    long["rank_in_model"] = long.groupby("model")["importance_norm"].rank(
        method="average", ascending=False
    )

    # Consensus: mean normalized importance, rank-sum (Borda-ish), support count
    agg = (long.groupby("gene_family")
           .agg(
               consensus_score=("importance_norm", "mean"),
               consensus_rank_sum=("rank_in_model", "sum"),
               models_supporting=("model", "nunique"),
           )
           .reset_index())

    # Wide view: model columns of normalized importance
    wide = (long.pivot_table(index="gene_family", columns="model", values="importance_norm", aggfunc="mean")
            .reset_index())

    consensus = agg.merge(wide, on="gene_family", how="left")
    return consensus, long


def add_models_where_topk(consensus: pd.DataFrame, long: pd.DataFrame, topm_per_model: int):
    tmp = long.copy()
    tmp["is_topm"] = tmp.groupby("model")["importance_norm"].rank(method="first", ascending=False) <= topm_per_model
    sup = (tmp.groupby("gene_family")["is_topm"].sum()
           .reset_index()
           .rename(columns={"is_topm": f"models_where_top{topm_per_model}"}))
    return consensus.merge(sup, on="gene_family", how="left")


# ----------------------------
# Directionality merge
# ----------------------------

def load_directionality(run_dir: Path):
    fn = run_dir / "feature_directionality_full.csv"
    if not fn.exists():
        return pd.DataFrame()
    df = _safe_read_csv(fn)
    keep = [c for c in [
        "gene_family", "direction", "effect_size",
        "diaz_rate_present", "diaz_rate_absent",
        "n_genomes_with_feature", "p_value"
    ] if c in df.columns]
    return df[keep].copy()


# ----------------------------
# Annotation summarization
# ----------------------------

def _mode_or_blank(vals):
    vals = [v for v in vals if isinstance(v, str) and v.strip()]
    if not vals:
        return ""
    return Counter(vals).most_common(1)[0][0]


def _top_k_with_counts(vals, k=3):
    vals = [v for v in vals if isinstance(v, str) and v.strip()]
    if not vals:
        return ""
    return "; ".join([f"{p} ({n})" for p, n in Counter(vals).most_common(k)])


def _collect_go(vals, k=10):
    gos = []
    for x in vals:
        if not isinstance(x, str):
            continue
        parts = re.split(r"[;,|\s]+", x.strip())
        gos.extend([p for p in parts if p.startswith("GO:")])
    if not gos:
        return ""
    return "; ".join([f"{g} ({n})" for g, n in Counter(gos).most_common(k)])


def summarize_family_annotations(protein_family_cds_tsv: Path, gene_families):
    cds = _safe_read_csv(protein_family_cds_tsv, sep="\t")

    if "gene_family" not in cds.columns:
        raise RuntimeError(f"'gene_family' column not found in {protein_family_cds_tsv}")

    cds = cds[cds["gene_family"].isin(set(gene_families))].copy()

    # Ensure expected columns exist
    for col in ["product", "go_ids", "protein_accession", "locus_tag"]:
        if col not in cds.columns:
            cds[col] = np.nan

    cds["product"] = cds["product"].fillna("").astype(str)
    cds["go_ids"] = cds["go_ids"].fillna("").astype(str)

    ann = (cds.groupby("gene_family")
           .agg(
               product_mode=("product", _mode_or_blank),
               product_top3=("product", lambda s: _top_k_with_counts(list(s), k=3)),
               go_top=("go_ids", lambda s: _collect_go(list(s), k=10)),
               example_proteins=("protein_accession", lambda s: "; ".join(list(pd.Series(s.dropna().astype(str)).unique())[:6])),
               example_locus_tags=("locus_tag", lambda s: "; ".join(list(pd.Series(s.dropna().astype(str)).unique())[:6])),
           )
           .reset_index())
    return ann


def attach_representatives(consensus: pd.DataFrame, gene_family_info_csv: Path):
    if not gene_family_info_csv.exists():
        return consensus

    info = _safe_read_csv(gene_family_info_csv)
    if not {"gene_family", "representative"}.issubset(info.columns):
        return consensus

    info = info.copy()
    # Representative headers are typically like "GCF_xxx|WP_yyyy"
    info["representative_protein"] = info["representative"].astype(str).str.split("|").str[-1]
    keep_cols = ["gene_family", "representative", "representative_protein"]
    if "num_genomes" in info.columns:
        keep_cols.append("num_genomes")

    return consensus.merge(info[keep_cols], on="gene_family", how="left")


def attach_uniprot(consensus: pd.DataFrame, refseq_to_uniprot_tsv: Path):
    if not refseq_to_uniprot_tsv.exists():
        return consensus

    if "representative_protein" not in consensus.columns:
        return consensus

    m = _safe_read_csv(refseq_to_uniprot_tsv, sep="\t")
    # tolerate unknown headers
    if m.shape[1] >= 2:
        m = m.iloc[:, :2].copy()
        m.columns = ["refseq", "uniprot"]
    else:
        return consensus

    tmp = consensus[["gene_family", "representative_protein"]].merge(
        m, left_on="representative_protein", right_on="refseq", how="left"
    )

    u = (tmp.groupby("gene_family")["uniprot"]
         .apply(lambda s: "; ".join(list(pd.Series(s.dropna().astype(str)).unique())[:6]))
         .reset_index()
         .rename(columns={"uniprot": "uniprot_ids"}))

    return consensus.merge(u, on="gene_family", how="left")


# ----------------------------
# Module guessing (for narrative scaffolding)
# ----------------------------

def guess_module(product_mode: str) -> str:
    s = (product_mode or "").lower()

    # Oxygen/ROS protection
    ros_kw = ["flavodiiron", "peroxiredoxin", "catalase", "superoxide dismutase", "sod", "rubredoxin", "alkyl hydroperoxide"]
    if any(k in s for k in ros_kw):
        return "O2/ROS protection"

    # Redox / electron transfer
    redox_kw = ["ferredoxin", "flavodoxin", "thioredoxin", "glutaredoxin", "oxidoreductase", "dehydrogenase", "redox"]
    if any(k in s for k in redox_kw):
        return "Redox/electron supply"

    # Cofactor/Fe-S/metallocluster biogenesis
    cofac_kw = ["fe-s", "iron-sulfur", "isc", "suf", "nifu", "nifs", "chaperone", "assembly protein", "metallocluster"]
    if any(k in s for k in cofac_kw):
        return "Cofactor/Fe-S biogenesis"

    # Membrane/envelope remodeling
    env_kw = ["glycolipid", "lipid", "polysaccharide", "capsule", "envelope", "membrane", "peptidoglycan"]
    if any(k in s for k in env_kw):
        return "Envelope remodeling"

    # Transport
    tr_kw = ["transporter", "permease", "abc transporter", "symporter", "antiporter", "porin", "channel"]
    if any(k in s for k in tr_kw):
        return "Transport"

    # Regulation
    reg_kw = ["transcriptional regulator", "response regulator", "sensor kinase", "two-component", "sigma factor", "repressor", "activator"]
    if any(k in s for k in reg_kw):
        return "Regulation"

    return "Unclassified/other"


def write_narrative_outline_md(out_md: Path, top_df: pd.DataFrame):
    lines = []
    lines.append("# Narrative outline (auto-binned)\n")
    lines.append("This is an automatically generated scaffold. Refine module assignments and wording.\n")

    for module, sub in top_df.groupby("module_guess", sort=False):
        lines.append(f"## {module}\n")
        for _, r in sub.iterrows():
            gf = r.get("gene_family", "")
            prod = r.get("product_mode", "")
            eff = r.get("effect_size", np.nan)
            score = r.get("consensus_score", np.nan)
            support = r.get("models_supporting", np.nan)
            lines.append(f"- {gf}: {prod} | effect_size={eff:.3f} | consensus={score:.4g} | models={int(support) if pd.notna(support) else ''}")
        lines.append("")

    out_md.write_text("\n".join(lines), encoding="utf-8")


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Build a consensus top-feature table with functional summaries for narrative writing."
    )
    ap.add_argument("--run-dir", required=True, help="Pipeline run dir, e.g., results/baseline")
    ap.add_argument("--top-k", type=int, default=60, help="How many top families to output by consensus score")
    ap.add_argument("--support-topm-per-model", type=int, default=60,
                    help="Count in how many models a family lands in top M (per-model)")

    ap.add_argument("--protein-family-cds", default="",
                    help="Path to protein_family_cds.tsv (tab-separated)")
    ap.add_argument("--refseq-to-uniprot", default="",
                    help="Path to refseq_to_uniprot.tsv (tab-separated, 2 columns)")
    ap.add_argument("--gene-family-info", default="",
                    help="Path to gene_family_info.csv (to get representative proteins)")

    ap.add_argument("--out-tsv", default="",
                    help="Output TSV path (default: <run-dir>/top_features_narrative.tsv)")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")

    # Defaults that work in this environment; override on HPC if needed
    protein_family_cds_default = _first_existing([
        args.protein_family_cds,
        "protein_family_cds.tsv",
        str(Path("/mnt/data") / "protein_family_cds.tsv"),
    ])
    refseq_to_uniprot_default = _first_existing([
        args.refseq_to_uniprot,
        "refseq_to_uniprot.tsv",
        str(Path("/mnt/data") / "refseq_to_uniprot.tsv"),
    ])
    gene_family_info_default = _first_existing([
        args.gene_family_info,
        "gene_family_info.csv",
        str(run_dir.parent / "gene_family_info.csv"),
        str(Path.cwd() / "gene_family_info.csv"),
    ])

    if not protein_family_cds_default:
        raise FileNotFoundError(
            "Could not locate protein_family_cds.tsv. Provide --protein-family-cds PATH."
        )

    # Load and aggregate importances
    consensus, long = load_feature_importances(run_dir)
    consensus = add_models_where_topk(consensus, long, topm_per_model=args.support_topm_per_model)

    # Attach representative protein + UniProt mapping (optional)
    if gene_family_info_default:
        consensus = attach_representatives(consensus, Path(gene_family_info_default))
    if refseq_to_uniprot_default:
        consensus = attach_uniprot(consensus, Path(refseq_to_uniprot_default))

    # Merge directionality if available
    direction = load_directionality(run_dir)
    if not direction.empty:
        consensus = consensus.merge(direction, on="gene_family", how="left")

    # Select top K by consensus score; tie-break by rank-sum then support
    consensus = consensus.sort_values(
        ["consensus_score", "models_supporting", "consensus_rank_sum"],
        ascending=[False, False, True]
    )
    top = consensus.head(args.top_k).copy()

    # Summarize function annotations
    ann = summarize_family_annotations(Path(protein_family_cds_default), top["gene_family"].tolist())
    top = top.merge(ann, on="gene_family", how="left")

    # Add module guess for narrative scaffolding
    top["module_guess"] = top["product_mode"].apply(guess_module)

    # Output paths
    out_tsv = Path(args.out_tsv) if args.out_tsv else (run_dir / "top_features_narrative.tsv")
    out_summary = run_dir / "module_summary.tsv"
    out_md = run_dir / "narrative_outline.md"

    # Write outputs
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    top.to_csv(out_tsv, sep="\t", index=False)

    module_summary = (top.groupby("module_guess")
                      .agg(
                          n_families=("gene_family", "count"),
                          mean_consensus=("consensus_score", "mean"),
                          median_effect_size=("effect_size", "median"),
                      )
                      .reset_index()
                      .sort_values(["n_families", "mean_consensus"], ascending=[False, False]))
    module_summary.to_csv(out_summary, sep="\t", index=False)

    write_narrative_outline_md(out_md, top.sort_values(["module_guess", "consensus_score"], ascending=[True, False]))

    # Console hints
    print("\nDONE")
    print(f"Top-feature table:     {out_tsv}")
    print(f"Module summary:        {out_summary}")
    print(f"Narrative outline:     {out_md}\n")

    # Print quick preview
    preview_cols = [c for c in [
        "gene_family", "consensus_score", "models_supporting",
        f"models_where_top{args.support_topm_per_model}",
        "effect_size", "product_mode", "module_guess"
    ] if c in top.columns]
    print("Preview:")
    print(top[preview_cols].head(min(12, len(top))).to_string(index=False))


if __name__ == "__main__":
    main()
