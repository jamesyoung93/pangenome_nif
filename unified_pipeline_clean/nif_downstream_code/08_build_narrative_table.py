#!/usr/bin/env python3
# 08_build_narrative_table.py  (FIXED: aligned counts + guardrails)

import os
import pandas as pd
import numpy as np

TOP_N = int(os.environ.get("TOP_N", "100"))

IMP_RF       = "feature_importance_Random_Forest.csv"
DIR_FULL     = "feature_directionality_full.csv"
MEMBERS_META = "gf_members_with_metadata.csv"
LLM_MERGED   = "llm_context_features.csv"
MATRIX_CSV   = "gene_family_matrix.csv"

OUT_CSV      = "narrative_top_features.csv"

def load_rf_importance(path):
    df = pd.read_csv(path)
    need = {"gene_family", "importance"}
    if not need.issubset(df.columns):
        raise RuntimeError(f"{path} missing columns {need}")
    df = df.sort_values("importance", ascending=False).reset_index(drop=True)
    df["rf_rank"] = np.arange(1, len(df) + 1)
    return df[["gene_family","importance","rf_rank"]]

def load_directionality(path):
    df = pd.read_csv(path)
    need = {"gene_family","effect_size","p_value","direction"}
    if not need.issubset(df.columns):
        raise RuntimeError(f"{path} missing columns {need}")
    keep = ["gene_family","effect_size","p_value","direction",
            "diaz_rate_present","diaz_rate_absent","n_genomes_with_feature"]
    for c in ["diaz_rate_present","diaz_rate_absent","n_genomes_with_feature"]:
        if c not in df.columns:
            df[c] = np.nan
    # carry 'importance' across models if present
    if "importance" in df.columns:
        df = df.rename(columns={"importance":"importance_across_models"})
    else:
        df["importance_across_models"] = np.nan
    return df[keep + ["importance_across_models"]]

def load_members(path):
    mm = pd.read_csv(path)
    need = {"gene_family","assembly_accession","organism_name","genus","is_diazotroph","member_id"}
    missing = [c for c in need if c not in mm.columns]
    if missing:
        raise RuntimeError(f"{path} missing columns {missing}")
    mm["is_diazotroph"] = (
        mm["is_diazotroph"].astype(str).map({"True":1,"False":0,"1":1,"0":0}).fillna(0).astype(int)
    )
    return mm

def attach_annotations(df):
    if os.path.exists(LLM_MERGED):
        ann = pd.read_csv(LLM_MERGED)
        keep = [c for c in ["gene_family","GeneName","ProteinName","EC","GO","target","pident","evalue"] if c in ann.columns]
        if keep:
            return df.merge(ann[keep].drop_duplicates("gene_family"), on="gene_family", how="left")
    return df

def attach_annotations(top_gf):
    import os
    import pandas as pd
    import numpy as np

    LLM_MERGED = "llm_context_features.csv"
    SP_HITS    = "gf_reps_vs_sprot.tsv"
    CLUSTERS_TSV = "gene_families_clusters.tsv"

    def first_nonempty(series):
        for x in series:
            if pd.notna(x) and str(x).strip():
                return str(x).strip()
        return ""

    def join_unique(series):
        vals = [str(x).strip() for x in series if pd.notna(x) and str(x).strip()]
        return ";".join(sorted(set(vals)))

    # Prefer the merged LLM annotations
    if os.path.exists(LLM_MERGED):
        ann = pd.read_csv(LLM_MERGED)
        # Guard: if there are no expected columns, bail out gracefully
        expected_any = {"GeneName","ProteinName","EC","GO","target","pident","evalue"} & set(ann.columns)
        if "gene_family" in ann.columns and expected_any:
            agg = ann.groupby("gene_family", as_index=False).agg({
                # choose first non-empty strings
                **({ "GeneName": first_nonempty } if "GeneName"   in ann.columns else {}),
                **({ "ProteinName": first_nonempty } if "ProteinName" in ann.columns else {}),
                # consolidate GO/EC across rows
                **({ "GO": join_unique } if "GO" in ann.columns else {}),
                **({ "EC": join_unique } if "EC" in ann.columns else {}),
                # prefer best hit stats where present
                **({ "target": first_nonempty } if "target" in ann.columns else {}),
                **({ "pident": "max" } if "pident" in ann.columns else {}),
                **({ "evalue": "min" } if "evalue" in ann.columns else {}),
            })
            return top_gf.merge(agg, on="gene_family", how="left")
        return top_gf

    # Fallback: best Swiss-Prot IDs if available
    if os.path.exists(SP_HITS) and os.path.exists(CLUSTERS_TSV):
        hits = pd.read_csv(SP_HITS, sep="\t", header=None,
                           names=["query", "target", "evalue", "pident", "alnlen"])
        reps = pd.read_csv(CLUSTERS_TSV, sep="\t", header=None, names=["rep", "member"])
        rep_order = reps["rep"].drop_duplicates().reset_index(drop=True)
        gf_map = pd.DataFrame({"rep": rep_order, "gene_family": [f"GF_{i:05d}" for i in range(len(rep_order))]})
        best = (hits.sort_values(["query", "evalue", "pident"], ascending=[True, True, False])
                    .groupby("query", as_index=False).first())
        best = best.merge(gf_map, left_on="query", right_on="rep", how="left")
        ann_small = best[["gene_family", "target", "pident", "evalue"]].drop_duplicates("gene_family")
        return top_gf.merge(ann_small, on="gene_family", how="left")

    return top_gf


def aggregate_members(mm, gf_list):
    mm = mm[mm["gene_family"].isin(gf_list)].copy()

    # genome/genera counts per GF
    agg = (mm.groupby("gene_family", as_index=False)
             .agg(n_members=("member_id","count"),
                  n_genomes=("assembly_accession", pd.Series.nunique),
                  n_diazotroph=("is_diazotroph","sum"),
                  n_non_diazo=("is_diazotroph", lambda s: int((1-s).sum()))))

    denom = (agg["n_diazotroph"] + agg["n_non_diazo"]).replace(0, np.nan)
    agg["frac_diazotroph"] = (agg["n_diazotroph"] / denom).fillna(0.0)

    # genus breakdown (count a genus as diazo if any diazo member exists)
    g = (mm.dropna(subset=["genus"])
           .groupby(["gene_family","genus"])["is_diazotroph"].max())
    g = g.reset_index()
    n_gen = g.groupby("gene_family")["genus"].nunique().rename("n_genera")
    n_gen_d = g[g["is_diazotroph"] > 0].groupby("gene_family")["genus"].nunique().rename("n_genera_diazo")
    n_gen_n = g[g["is_diazotroph"] == 0].groupby("gene_family")["genus"].nunique().rename("n_genera_non")
    gen_join = (n_gen.to_frame()
                  .join(n_gen_d, how="left").join(n_gen_n, how="left")
                  .fillna(0).astype(int).reset_index())

    # compact organism roster: one line per genome
    def org_list(df):
        seen = set()
        vals = []
        for _, r in df.drop_duplicates(["assembly_accession"]).iterrows():
            acc = str(r["assembly_accession"])
            if acc in seen: 
                continue
            seen.add(acc)
            label = "diazo" if int(r["is_diazotroph"]) == 1 else "non"
            name = str(r.get("organism_name") or r.get("genus") or acc)
            if name.startswith("{") and r.get("genus"):
                name = str(r["genus"])
            vals.append(f"{name}|{label}")
        s = "; ".join(vals)
        return s if len(s) <= 6000 else (s[:6000] + " ...")

    orgs = (mm.groupby("gene_family", group_keys=False)[["assembly_accession","organism_name","genus","is_diazotroph"]].apply(org_list).rename("organisms").reset_index())

    out = (agg.merge(gen_join, on="gene_family", how="left")
              .merge(orgs, on="gene_family", how="left"))
    return out

def main():
    rf   = load_rf_importance(IMP_RF)
    dirf = load_directionality(DIR_FULL)
    mem  = load_members(MEMBERS_META)

    # Top-N by RF
    top = rf.sort_values("rf_rank").head(TOP_N).copy()
    top = top.merge(dirf, on="gene_family", how="left")

    # Aggregate membership (counts from members table)
    agg = aggregate_members(mem, set(top["gene_family"]))
    top = top.merge(agg, on="gene_family", how="left")

    # Add counts from matrix for alignment checks
    mat = pd.read_csv(MATRIX_CSV, index_col=0)
    gf_counts = (mat.sum(axis=0).rename("n_from_matrix")
                    .reset_index().rename(columns={"index":"gene_family"}))
    top = top.merge(gf_counts, on="gene_family", how="left")
    top["n_from_members"] = top["n_genomes"].fillna(0).astype(float)

    # Attach annotations (UniProt etc.) if present
    top = attach_annotations(top)

    # Signed importance (RF * sign(effect))
    top["signed_importance"] = top["importance"] * np.sign(top["effect_size"].fillna(0.0))

    # Order by RF rank then signed importance
    top = top.sort_values(["rf_rank","signed_importance","importance"], ascending=[True, False, False])

    cols = [
        "gene_family","rf_rank","importance","signed_importance",
        "direction","effect_size","p_value",
        "n_genomes","n_diazotroph","n_non_diazo","frac_diazotroph",
        "n_from_matrix","n_from_members",
        "n_genera","n_genera_diazo","n_genera_non",
        "organisms",
        "GeneName","ProteinName","EC","GO","target","pident","evalue",
        "importance_across_models"
    ]
    for c in cols:
        if c not in top.columns:
            top[c] = ""
    out = top[cols].fillna("")
    out.to_csv(OUT_CSV, index=False, encoding="utf-8")
    print(f"Wrote {OUT_CSV} with {len(out)} rows (TOP_N={TOP_N}).")

    # Quick console check for any residual misalignment
    if "n_from_matrix" in out and "n_from_members" in out:
        chk = out.copy()
        chk["diff"] = pd.to_numeric(chk["n_from_members"], errors="coerce").fillna(0) - \
                      pd.to_numeric(chk["n_from_matrix"], errors="coerce").fillna(0)
        bad = chk[chk["diff"] != 0].sort_values("diff", key=abs, ascending=False)
        if len(bad):
            nbad = len(bad)
            print(f"NOTE: {nbad} gene families still differ between members and matrix. Top 10:")
            print(bad[["gene_family","n_from_matrix","n_from_members","diff"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
