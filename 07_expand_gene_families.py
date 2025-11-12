#!/usr/bin/env python3
# 07_expand_gene_families.py  (FIXED: canonical GF mapping from Step 3)
import os
import pandas as pd
import numpy as np

# Inputs
CLUSTERS_TSV   = "gene_families_clusters.tsv"          # rep \t member
INFO_CSV       = "gene_family_info.csv"                 # canonical GF IDs from Step 3
METADATA_CSV   = "complete_genomes_with_proteins.csv"   # from Step 3
PA_MATRIX_CSV  = "gene_family_matrix.csv"               # from Step 3
SUMMARY_CSV    = "feature_directionality_summary.csv"   # optional, from Step 5

# Outputs
OUT_MEMBERS             = "gf_members.csv"
OUT_MEMBERS_METADATA    = "gf_members_with_metadata.csv"
OUT_MEMBER_SUMMARY      = "gf_member_summary.csv"
OUT_PRES_LONG           = "gf_presence_matrix_long.csv"
OUT_TOP_MEMBERS         = "top_features_members.csv"

TOP_N = 100

def load_canonical_rep_map(info_csv: str) -> pd.DataFrame:
    info = pd.read_csv(info_csv)
    need = {"gene_family", "representative"}
    if not need.issubset(info.columns):
        raise RuntimeError(f"{info_csv} missing columns {need}")
    rep_map = (info[["gene_family", "representative"]]
               .drop_duplicates()
               .rename(columns={"representative": "rep_id"}))
    return rep_map

def split_member(member_id: str):
    if "|" in member_id:
        acc, prot = member_id.split("|", 1)
        return acc, prot
    return member_id, ""

def main():
    # Canonical mapping from Step 3
    rep_map = load_canonical_rep_map(INFO_CSV)

    # Read all cluster edges, then **inner-join** to canonical reps
    reps = pd.read_csv(CLUSTERS_TSV, sep="\t", header=None, names=["rep_id", "member_id"])
    reps = reps.merge(rep_map, on="rep_id", how="inner")

    # Parse assembly_accession and protein_id
    accs, prots = zip(*[split_member(mid) for mid in reps["member_id"]])
    reps["assembly_accession"] = accs
    reps["protein_id"] = prots

    members = reps[["gene_family", "rep_id", "member_id", "assembly_accession", "protein_id"]].copy()
    members.to_csv(OUT_MEMBERS, index=False)
    print(f"Wrote {OUT_MEMBERS} ({len(members)} rows)")

    # Attach genome metadata
    meta = pd.read_csv(METADATA_CSV)
    if meta.duplicated("assembly_accession").any():
        meta = meta.drop_duplicates("assembly_accession", keep="first")
    members_meta = members.merge(meta, on="assembly_accession", how="left")
    members_meta.to_csv(OUT_MEMBERS_METADATA, index=False)
    print(f"Wrote {OUT_MEMBERS_METADATA} ({len(members_meta)} rows)")

    # Summaries (genomes, genera, diazo breakdown)
    mm = members_meta.copy()
    mm["is_diazotroph"] = mm.get("is_diazotroph", 0).astype(int)

    summary = (
        mm.groupby("gene_family", as_index=False)
          .agg(n_members=("member_id", "count"),
               n_genomes=("assembly_accession", pd.Series.nunique),
               n_diazotroph=("is_diazotroph", "sum"),
               n_non_diazo=("is_diazotroph", lambda s: int((1 - s).sum())),
               n_genera=("genus", lambda s: s.dropna().nunique()))
    )
    denom = (summary["n_diazotroph"] + summary["n_non_diazo"]).replace(0, np.nan)
    summary["frac_diazotroph"] = (summary["n_diazotroph"] / denom).fillna(0.0)
    summary.to_csv(OUT_MEMBER_SUMMARY, index=False)
    print(f"Wrote {OUT_MEMBER_SUMMARY} ({len(summary)} rows)")

    # Optional presence matrix (long)
    if os.path.exists(PA_MATRIX_CSV):
        pa = pd.read_csv(PA_MATRIX_CSV, index_col=0)
        pa_long = pa.stack().rename("present").reset_index()
        pa_long.columns = ["assembly_accession", "gene_family", "present"]
        pa_long.to_csv(OUT_PRES_LONG, index=False)
        print(f"Wrote {OUT_PRES_LONG} ({len(pa_long)} rows)")

    # Optional top-N members for quick browsing
    if os.path.exists(SUMMARY_CSV):
        summ = pd.read_csv(SUMMARY_CSV)
        if {"gene_family", "importance"}.issubset(summ.columns):
            top_gf = (summ.sort_values("importance", ascending=False)
                         .head(TOP_N)["gene_family"].tolist())
            top_members = members_meta[members_meta["gene_family"].isin(top_gf)].copy()
            top_members.to_csv(OUT_TOP_MEMBERS, index=False)
            print(f"Wrote {OUT_TOP_MEMBERS} for top {len(top_gf)} features ({len(top_members)} rows)")

    # Sanity: every GF here should be a column in the matrix
    if os.path.exists(PA_MATRIX_CSV):
        mat_cols = set(pd.read_csv(PA_MATRIX_CSV, index_col=0).columns)
        extra = set(members["gene_family"].unique()) - mat_cols
        if extra:
            print(f"WARNING: {len(extra)} GFs in members not in matrix (should be 0). Examples: {sorted(list(extra))[:5]}")

if __name__ == "__main__":
    main()
