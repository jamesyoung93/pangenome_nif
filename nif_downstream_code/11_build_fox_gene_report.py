#!/usr/bin/env python3
import argparse
import glob
import os
import re
from collections import Counter

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


GF_RE = re.compile(r"^GF_\d+$")


def read_table(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    sep = "\t" if ext in [".tsv", ".txt"] else ","
    return pd.read_csv(path, sep=sep, dtype=str)


def coerce_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def pick_gene_family_col(df: pd.DataFrame) -> str:
    best_col, best_score = df.columns[0], -1.0
    for col in df.columns:
        vals = df[col].astype(str).head(500)
        score = vals.str.match(GF_RE).mean()
        if score > best_score:
            best_col, best_score = col, score
    return best_col


def pick_numeric_col(df: pd.DataFrame, preferred_names) -> str | None:
    lower_map = {c.lower(): c for c in df.columns}
    for p in preferred_names:
        if p.lower() in lower_map:
            return lower_map[p.lower()]
    for col in df.columns[1:]:
        s = coerce_num(df[col])
        if s.notna().mean() >= 0.50:
            return col
    return None


def find_run_dirs(root: str) -> list[str]:
    """
    Find directories under any results* folder that contain directionality or importance outputs.
    """
    runs = []
    for entry in os.listdir(root):
        base = os.path.join(root, entry)
        if not os.path.isdir(base):
            continue
        if not entry.startswith("results"):
            continue
        for dirpath, _, filenames in os.walk(base):
            has_imp = any(f.startswith("feature_importance_") and f.endswith(".csv") for f in filenames)
            has_dir = any(f.startswith("feature_directionality") and f.endswith(".csv") for f in filenames)
            has_pur = ("gene_family_purity_stats.csv" in filenames)
            if has_imp or has_dir or has_pur:
                runs.append(dirpath)

    # de-dup preserving order
    seen = set()
    out = []
    for r in runs:
        if r not in seen:
            out.append(r)
            seen.add(r)
    return out


def load_directionality(run_dir: str) -> pd.DataFrame:
    cand = [
        os.path.join(run_dir, "feature_directionality_full.csv"),
        os.path.join(run_dir, "feature_directionality_summary.csv"),
    ]
    path = next((p for p in cand if os.path.exists(p)), None)
    if not path:
        return pd.DataFrame(columns=["gene_family", "effect_size"])

    df = read_table(path)
    gf_col = pick_gene_family_col(df)

    eff_col = pick_numeric_col(df, preferred_names=["effect_size", "delta", "difference"])
    if eff_col is None:
        for c in df.columns:
            cl = c.lower()
            if "effect" in cl or "delta" in cl:
                eff_col = c
                break

    out = pd.DataFrame({
        "gene_family": df[gf_col].astype(str),
        "effect_size": coerce_num(df[eff_col]) if eff_col else np.nan,
    })
    return out.dropna(subset=["gene_family"])


def load_purity(run_dir: str, purity_threshold: float) -> pd.DataFrame:
    path = os.path.join(run_dir, "gene_family_purity_stats.csv")
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_family", "diazotroph_pct", "is_pure"])

    df = read_table(path)
    gf_col = pick_gene_family_col(df)

    pct_col = pick_numeric_col(df, preferred_names=["diazotroph_pct", "pct_positive", "pos_frac", "positive_frac"])
    if pct_col is None:
        for c in df.columns:
            cl = c.lower()
            if "pct" in cl and ("diaz" in cl or "pos" in cl):
                pct_col = c
                break

    is_pure_col = None
    for c in df.columns:
        if c.lower() == "is_pure":
            is_pure_col = c
            break

    out = pd.DataFrame({
        "gene_family": df[gf_col].astype(str),
        "diazotroph_pct": coerce_num(df[pct_col]) if pct_col else np.nan,
    })

    if is_pure_col:
        # accept True/False strings
        out["is_pure"] = df[is_pure_col].astype(str).str.lower().isin(["true", "1", "yes"])
    else:
        # derive from threshold if explicit boolean not present
        pt = purity_threshold
        out["is_pure"] = (out["diazotroph_pct"] >= pt) | (out["diazotroph_pct"] <= (1 - pt))

    # carry counts if present
    for c in df.columns:
        cl = c.lower()
        if cl in ("total_genomes", "diazotroph_count", "present_count", "n_present"):
            out[c] = coerce_num(df[c])

    return out.dropna(subset=["gene_family"])


def load_used_set(run_dir: str) -> set[str]:
    path = os.path.join(run_dir, "final_gene_families_used.csv")
    if not os.path.exists(path):
        return set()
    df = read_table(path)
    gf_col = pick_gene_family_col(df)
    return set(df[gf_col].dropna().astype(str).tolist())


def load_importances_rankpct(run_dir: str, top_per_model: int | None = None) -> pd.DataFrame:
    """
    Robust consensus importance:
    - compute abs importance
    - compute within-file percentile rank (highest -> 1.0)
    - aggregate later across run-model hits
    """
    paths = sorted(glob.glob(os.path.join(run_dir, "feature_importance_*.csv")))
    paths = [p for p in paths if os.path.basename(p) != "feature_importance_top.csv"]

    out = []
    for p in paths:
        model = os.path.basename(p).replace("feature_importance_", "").replace(".csv", "")
        df = read_table(p)
        gf_col = pick_gene_family_col(df)

        imp_col = pick_numeric_col(df, preferred_names=["importance", "gain", "weight", "abs_coef", "coef", "coefficient"])
        if imp_col is None:
            if len(df.columns) >= 2:
                imp_col = df.columns[1]
            else:
                continue

        tmp = pd.DataFrame({
            "gene_family": df[gf_col].astype(str),
            "importance_raw": coerce_num(df[imp_col]),
        }).dropna(subset=["gene_family", "importance_raw"])

        tmp["importance_abs"] = tmp["importance_raw"].abs()

        # optional truncation for speed (note: percentile ranks then reflect truncated file)
        if top_per_model is not None and top_per_model > 0:
            tmp = tmp.sort_values("importance_abs", ascending=False).head(top_per_model).copy()

        # percentile rank (highest importance_abs -> 1.0)
        tmp["importance_rank_pct"] = tmp["importance_abs"].rank(pct=True, ascending=True)

        # best rank (1 is best)
        tmp = tmp.sort_values("importance_abs", ascending=False)
        tmp["rank_in_model"] = np.arange(1, len(tmp) + 1)

        tmp["model"] = model
        out.append(tmp[["gene_family", "model", "importance_abs", "importance_rank_pct", "rank_in_model"]])

    if not out:
        return pd.DataFrame(columns=["gene_family", "model", "importance_abs", "importance_rank_pct", "rank_in_model"])
    return pd.concat(out, ignore_index=True)


def load_module_map(run_dir: str) -> pd.DataFrame:
    """
    If module_summary.tsv exists, ingest it (gene_family -> module/category).
    """
    path = os.path.join(run_dir, "module_summary.tsv")
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_family", "module"])

    df = read_table(path)
    gf_col = pick_gene_family_col(df)

    mod_col = None
    for c in df.columns:
        if c.lower() in ("module", "functional_module", "category", "bin"):
            mod_col = c
            break
    if mod_col is None:
        # heuristic: first non-family column
        cand = [c for c in df.columns if c != gf_col]
        if cand:
            mod_col = cand[0]

    out = pd.DataFrame({
        "gene_family": df[gf_col].astype(str),
        "module": df[mod_col].astype(str) if mod_col else "",
    }).dropna(subset=["gene_family"])

    return out


def load_family_level_annotations(root: str, families: set[str]) -> pd.DataFrame:
    """
    Prefer protein_family_cds.tsv for family -> gene/product/GO fields.
    """
    path = os.path.join(root, "protein_family_cds.tsv")
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_family"])

    df = read_table(path)
    gf_col = pick_gene_family_col(df)
    df = df.rename(columns={gf_col: "gene_family"})
    df = df[df["gene_family"].isin(families)].copy()

    keep = ["gene_family"]
    for c in df.columns:
        cl = c.lower()
        if any(k in cl for k in ["locus", "gene", "product", "desc", "annotation", "uniprot", "swiss", "go", "kegg", "pfam", "cog", "ec"]):
            if c not in keep:
                keep.append(c)
    return df[keep].drop_duplicates("gene_family")


def load_membership_examples(root: str, families: set[str], max_members: int) -> pd.DataFrame:
    """
    Use genome_protein_family_map.tsv to provide example proteins/genomes per family.
    """
    path = os.path.join(root, "genome_protein_family_map.tsv")
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene_family"])

    df = read_table(path)
    gf_col = pick_gene_family_col(df)

    prot_col = None
    gen_col = None
    for c in df.columns:
        cl = c.lower()
        if prot_col is None and ("protein" in cl or "prot" in cl or "accession" in cl):
            prot_col = c
        if gen_col is None and ("genome" in cl or "assembly" in cl or "gcf_" in cl):
            gen_col = c

    df = df.rename(columns={gf_col: "gene_family"})
    df = df[df["gene_family"].isin(families)].copy()

    if prot_col is None:
        cand = [c for c in df.columns if c != "gene_family"]
        prot_col = cand[0] if cand else None

    rows = []
    for gf, sub in df.groupby("gene_family", sort=False):
        prots = sub[prot_col].dropna().astype(str).unique().tolist() if prot_col and prot_col in sub.columns else []
        gens = sub[gen_col].dropna().astype(str).unique().tolist() if gen_col and gen_col in sub.columns else []
        rows.append({
            "gene_family": gf,
            "n_member_proteins": len(prots),
            "n_member_genomes": len(gens),
            "example_proteins": ";".join(prots[:max_members]),
            "example_genomes": ";".join(gens[:max_members]),
        })
    return pd.DataFrame(rows)


def mode_or_empty(values: list[str]) -> str:
    vals = [v for v in values if v and v.lower() != "nan"]
    if not vals:
        return ""
    return Counter(vals).most_common(1)[0][0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="Project root (default: .)")
    ap.add_argument("--outdir", default="fox_report", help="Output dir (default: fox_report)")
    ap.add_argument("--purity_threshold", type=float, default=0.90, help="Purity threshold (default: 0.90)")
    ap.add_argument("--top_per_model", type=int, default=0,
                    help="Optional: keep only top N per model importance file (0 = keep all; recommended)")

    ap.add_argument("--max_members", type=int, default=10, help="Max example proteins/genomes per family (default: 10)")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    run_dirs = find_run_dirs(root)
    if not run_dirs:
        raise SystemExit(f"No run directories found under {root} (expected results* folders).")

    all_dir = []
    all_purity = []
    all_imp = []
    all_mod = []
    used_any = set()

    for rd in run_dirs:
        run_name = os.path.relpath(rd, root)

        ddf = load_directionality(rd)
        if not ddf.empty:
            ddf["run"] = run_name
            all_dir.append(ddf)

        pdf = load_purity(rd, purity_threshold=args.purity_threshold)
        if not pdf.empty:
            pdf["run"] = run_name
            all_purity.append(pdf)

        topn = None if args.top_per_model == 0 else args.top_per_model
        idf = load_importances_rankpct(rd, top_per_model=topn)
        if not idf.empty:
            idf["run"] = run_name
            all_imp.append(idf)

        mdf = load_module_map(rd)
        if not mdf.empty:
            mdf["run"] = run_name
            all_mod.append(mdf)

        used_any |= load_used_set(rd)

    dir_long = pd.concat(all_dir, ignore_index=True) if all_dir else pd.DataFrame(columns=["gene_family", "effect_size", "run"])
    pur_long = pd.concat(all_purity, ignore_index=True) if all_purity else pd.DataFrame(columns=["gene_family", "diazotroph_pct", "is_pure", "run"])
    imp_long = pd.concat(all_imp, ignore_index=True) if all_imp else pd.DataFrame(columns=["gene_family", "model", "importance_rank_pct", "rank_in_model", "run"])
    mod_long = pd.concat(all_mod, ignore_index=True) if all_mod else pd.DataFrame(columns=["gene_family", "module", "run"])

    # Directionality summary (only where computed)
    dir_sum = (dir_long.groupby("gene_family", as_index=False)
               .agg(effect_size_mean=("effect_size", "mean"),
                    effect_size_median=("effect_size", "median"),
                    effect_size_max=("effect_size", "max"),
                    n_runs_with_directionality=("run", "nunique")))
    dir_sum["is_positive_assoc"] = dir_sum["effect_size_mean"] > 0

    # Purity summary (works even when directionality is missing)
    pur_sum = (pur_long.groupby("gene_family", as_index=False)
               .agg(diazotroph_pct_mean=("diazotroph_pct", "mean"),
                    diazotroph_pct_max=("diazotroph_pct", "max"),
                    is_pure_any=("is_pure", "max"),
                    n_runs_with_purity=("run", "nunique")))

    pt = float(args.purity_threshold)
    pur_sum["is_pure_pos_any"] = pur_sum["diazotroph_pct_max"] >= pt
    pur_sum["is_pure_neg_any"] = pur_sum["diazotroph_pct_max"] <= (1 - pt)

    # Importance summary (robust: rank percentile)
    if not imp_long.empty:
        imp_sum = (imp_long.groupby("gene_family", as_index=False)
                   .agg(n_models_with_importance=("model", "nunique"),
                        n_run_model_hits=("run", "count"),
                        consensus_rank_pct_mean=("importance_rank_pct", "mean"),
                        consensus_rank_pct_median=("importance_rank_pct", "median"),
                        best_rank=("rank_in_model", "min")))
    else:
        imp_sum = pd.DataFrame(columns=[
            "gene_family", "n_models_with_importance", "n_run_model_hits",
            "consensus_rank_pct_mean", "consensus_rank_pct_median", "best_rank"
        ])

    # Module consensus (mode across runs)
    if not mod_long.empty:
        mod_sum = (mod_long.groupby("gene_family", as_index=False)
                   .agg(module_mode=("module", lambda s: mode_or_empty(list(s))),
                        n_runs_with_module=("run", "nunique")))
    else:
        mod_sum = pd.DataFrame(columns=["gene_family", "module_mode", "n_runs_with_module"])

    # Master table
    fams = set(dir_sum["gene_family"].tolist()) | set(pur_sum["gene_family"].tolist()) | set(imp_sum["gene_family"].tolist()) | used_any
    master = pd.DataFrame({"gene_family": sorted(fams)})

    master = master.merge(dir_sum, on="gene_family", how="left")
    master = master.merge(pur_sum, on="gene_family", how="left")
    master = master.merge(imp_sum, on="gene_family", how="left")
    master = master.merge(mod_sum, on="gene_family", how="left")

    master["used_in_any_final_set"] = master["gene_family"].isin(used_any)
    master["is_positive_assoc"] = master["is_positive_assoc"].fillna(False)

    master["n_models_with_importance"] = master["n_models_with_importance"].fillna(0).astype(int)
    master["n_run_model_hits"] = master["n_run_model_hits"].fillna(0).astype(int)

    master["is_pure_pos_any"] = master["is_pure_pos_any"].fillna(False)
    master["is_pure_neg_any"] = master["is_pure_neg_any"].fillna(False)

    # Tier definitions (fixed)
    # Tier 1: positive association (effect_size_mean > 0) AND (appears in any model importance OR used in any final set)
    tier1 = master[
        master["is_positive_assoc"] &
        ((master["n_models_with_importance"] > 0) | master["used_in_any_final_set"])
    ].copy()

    # Tier 2: purity-held-out positives (diazotroph_pct >= threshold), not used anywhere
    # Note: directionality is typically not computed for these, so Tier 2 does NOT require effect_size.
    tier2 = master[
        master["is_pure_pos_any"] &
        (~master["used_in_any_final_set"])
    ].copy()

    tier1 = tier1.sort_values(["consensus_rank_pct_mean", "effect_size_mean"], ascending=False, na_position="last")
    tier2 = tier2.sort_values(["diazotroph_pct_max", "diazotroph_pct_mean"], ascending=False, na_position="last")

    families_of_interest = set(tier1["gene_family"].tolist()) | set(tier2["gene_family"].tolist())

    # Attach family-level annotations (gene/product/GO) if available
    fam_annot = load_family_level_annotations(root, families_of_interest)
    if not fam_annot.empty:
        master = master.merge(fam_annot, on="gene_family", how="left")
        tier1 = tier1.merge(fam_annot, on="gene_family", how="left")
        tier2 = tier2.merge(fam_annot, on="gene_family", how="left")

    # Attach membership examples
    mem = load_membership_examples(root, families_of_interest, max_members=args.max_members)
    if not mem.empty:
        master = master.merge(mem, on="gene_family", how="left")
        tier1 = tier1.merge(mem, on="gene_family", how="left")
        tier2 = tier2.merge(mem, on="gene_family", how="left")

    # Write tables
    master.to_csv(os.path.join(outdir, "master_family_table.tsv"), sep="\t", index=False)
    tier1.to_csv(os.path.join(outdir, "tier1_positive_model_selected.tsv"), sep="\t", index=False)
    tier2.to_csv(os.path.join(outdir, "tier2_pure_positive_heldout.tsv"), sep="\t", index=False)

    # Plot (fixed): rank-percentile consensus avoids outliers; encode support
    if len(tier1) > 0 and "effect_size_mean" in tier1.columns and "consensus_rank_pct_mean" in tier1.columns:
        x = tier1["effect_size_mean"].astype(float)
        y = tier1["consensus_rank_pct_mean"].astype(float)

        # size encodes run-model hits (bounded for readability)
        hits = tier1["n_run_model_hits"].astype(float).fillna(0.0)
        s = 15 + 2.0 * np.sqrt(hits.clip(lower=0.0))

        # color encodes # models contributing
        c = tier1["n_models_with_importance"].astype(float).fillna(0.0)

        plt.figure()
        sc = plt.scatter(x, y, s=s, c=c)  # default colormap
        plt.xlabel("Effect size (mean across runs; directionality)")
        plt.ylabel("Consensus importance (mean percentile rank across run-model hits)")
        plt.title("Tier 1 families: directionality vs consensus importance (rank-percentile)")
        cb = plt.colorbar(sc)
        cb.set_label("# models with importance")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "tier1_effect_vs_importance_rankpct.png"), dpi=300)
        plt.close()

    # Manuscript scaffold
    with open(os.path.join(outdir, "manuscript_bullets.md"), "w", encoding="utf-8") as f:
        f.write("# FOX gene family report (auto-generated)\n\n")
        f.write(f"- Runs scanned: {len(run_dirs)}\n")
        f.write(f"- Tier 1 (positive + model-selected/used): {len(tier1)} families\n")
        f.write(f"- Tier 2 (purity-held-out positives = {pt:.2f}, not used): {len(tier2)} families\n\n")

        f.write("## Top Tier 1 families (first 30)\n\n")
        show = tier1.head(30).copy()
        cols = [c for c in [
            "gene_family", "effect_size_mean",
            "consensus_rank_pct_mean", "n_models_with_importance", "n_run_model_hits",
            "used_in_any_final_set", "module_mode",
            "diazotroph_pct_max", "diazotroph_pct_mean",
            "gene", "product", "go_ids",
            "n_member_genomes", "n_member_proteins"
        ] if c in show.columns]
        f.write("\t".join(cols) + "\n")
        for _, r in show[cols].iterrows():
            row = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, float) and np.isnan(v):
                    v = ""
                row.append(str(v))
            f.write("\t".join(row) + "\n")

        f.write("\n## Top Tier 2 purity-held-out families (first 30)\n\n")
        show2 = tier2.head(30).copy()
        cols2 = [c for c in [
            "gene_family", "diazotroph_pct_max", "diazotroph_pct_mean",
            "module_mode", "gene", "product", "go_ids",
            "n_member_genomes", "n_member_proteins"
        ] if c in show2.columns]
        f.write("\t".join(cols2) + "\n")
        for _, r in show2[cols2].iterrows():
            row = []
            for c in cols2:
                v = r.get(c, "")
                if isinstance(v, float) and np.isnan(v):
                    v = ""
                row.append(str(v))
            f.write("\t".join(row) + "\n")

    print(f"[OK] Wrote report to: {outdir}")
    print("Key outputs:")
    print(" - master_family_table.tsv")
    print(" - tier1_positive_model_selected.tsv")
    print(" - tier2_pure_positive_heldout.tsv")
    print(" - manuscript_bullets.md")
    print(" - tier1_effect_vs_importance_rankpct.png")


if __name__ == "__main__":
    main()
