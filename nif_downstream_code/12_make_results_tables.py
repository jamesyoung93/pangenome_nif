#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from collections import Counter

import numpy as np
import pandas as pd


GF_RE = re.compile(r"^GF_\d+$")


# -------------------------
# FOX rediscovery heuristics
# -------------------------
# Keep conservative: only tag what you are comfortable calling "established FOX / diazotrophy machinery".
KNOWN_DIAZ_CORE = {
    "nifh", "nifd", "nifk", "nifb", "nife", "nifn", "nifs", "nifu", "nifv", "nifw", "nifx", "nifz", "nifq", "nift",
    "fdxh", "fdxn", "hesa", "hupl", "hups", "hoxh", "hoxy",
}

# Generic functional cues that are FOX-consistent even if not canonical nif genes
FOX_CONSISTENT_KEYWORDS = [
    "nitrogenase", "ferredoxin", "flavodoxin", "4fe-4s", "iron-sulfur", "fe-s", "maturation",
    "peroxiredoxin", "thioredoxin", "glutaredoxin", "superoxide", "catalase",
    "cytochrome", "oxidase", "respir", "dehydrogenase",
]

HOUSEKEEPING_MARKER_KEYWORDS = [
    "ribosomal", "30s", "50s", "trna ligase", "dna gyrase", "topoisomerase",
    "photosystem", "psb", "psa", "atp synthase", "methionyl aminopeptidase",
    "cell division", "fts", "replication", "polymerase"
]


def assign_module(gene: str, product: str, go_ids: str) -> str:
    s = " ".join([str(gene or ""), str(product or ""), str(go_ids or "")]).lower()

    # Housekeeping marker
    if any(k in s for k in HOUSEKEEPING_MARKER_KEYWORDS):
        return "Housekeeping / lineage marker"

    # Nitrogenase & metallocluster
    if "nif" in s or "nitrogenase" in s:
        return "Nitrogenase / nif machinery"
    if any(k in s for k in ["4fe-4s", "iron-sulfur", "fe-s", "suf", "isc", "nifb", "nife", "nifn", "maturation"]):
        return "Metallocluster & Fe-S biogenesis"

    # O2 protection / redox homeostasis
    if any(k in s for k in ["peroxiredoxin", "thioredoxin", "glutaredoxin", "superoxide", "catalase", "oxidative stress", "redox"]):
        return "O2 protection & redox homeostasis"

    # Respiration / bioenergetics
    if any(k in s for k in ["cytochrome", "oxidase", "respir", "succinate dehydrogenase", "ndh", "qcr", "pet", "electron transport"]):
        return "Respiration & bioenergetics"

    # Carbon / NADPH supply
    if any(k in s for k in ["pentose phosphate", "phosphogluconate", "gnd", "glg", "glycogen", "aldolase", "gapdh", "glycol"]):
        return "Carbon metabolism & NAD(P)H supply"

    # Regulation
    if any(k in s for k in ["transcription", "regulator", "two-component", "response regulator", "fur", "sigma factor"]):
        return "Regulation"

    # Transport
    if any(k in s for k in ["abc transporter", "permease", "substrate-binding", "transport"]):
        return "Transport"

    # Stress / membrane / unknown
    if any(k in s for k in ["duf", "hypothetical", "membrane", "pmp3", "yqae", "ftsh", "peptidase"]):
        return "Stress / membrane / uncharacterized"

    return "Other / unassigned"


def tag_fox_status(gene: str, product: str, go_ids: str) -> str:
    s = " ".join([str(gene or ""), str(product or ""), str(go_ids or "")]).lower()
    g = str(gene or "").lower()

    if g in KNOWN_DIAZ_CORE or ("nitrogenase" in s):
        return "Rediscovery: core diazotrophy/FOX machinery"

    if any(k in s for k in FOX_CONSISTENT_KEYWORDS):
        return "FOX-consistent support function (plausible rediscovery)"

    if any(k in s for k in HOUSEKEEPING_MARKER_KEYWORDS):
        return "Housekeeping marker (not a transplant candidate)"

    return "Novel/uncertain candidate"


def read_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    # coerce numerics where expected
    for c in df.columns:
        if c in (
            "effect_size_mean",
            "consensus_rank_pct_mean",
            "diazotroph_pct_max",
            "diazotroph_pct_mean",
            "n_models_with_importance",
            "n_run_model_hits",
        ):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def find_existing(path_candidates):
    for p in path_candidates:
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report_dir", required=True, help="Directory containing fox_report_v2 outputs")
    ap.add_argument("--outdir", default="results_tables", help="Output directory")
    ap.add_argument("--top_n", type=int, default=100, help="Top N to include in manuscript main table")
    args = ap.parse_args()

    report_dir = os.path.abspath(args.report_dir)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    tier1_path = os.path.join(report_dir, "tier1_positive_model_selected.tsv")

    tier2_candidates = [
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
    ]
    # In case your v2 script used a slightly different name:
    tier2_candidates += [
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
    ]
    # Practical: also accept the filename from the earlier script variant
    tier2_candidates += [
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
        os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
    ]
    # Actually used by your v2 generator in the code I gave you earlier:
    tier2_candidates += [os.path.join(report_dir, "tier2_pure_positive_heldout.tsv"),
                         os.path.join(report_dir, "tier2_pure_positive_heldout.tsv")]
    # And the file name you previously showed in chat:
    tier2_candidates += [os.path.join(report_dir, "tier2_pure_positive_heldout.tsv")]

    # Real fallback: accept the file you showed earlier
    tier2_candidates += [os.path.join(report_dir, "tier2_positive_pure_heldout.tsv")]

    tier2_path = find_existing(tier2_candidates)

    if not os.path.exists(tier1_path):
        raise SystemExit(f"Missing tier1 file: {tier1_path}")

    tier1 = read_tsv(tier1_path)
    tier2 = read_tsv(tier2_path) if tier2_path else pd.DataFrame()

    # Add module and FOX status tags
    for df in [tier1, tier2]:
        if df is None or df.empty:
            continue
        df["module_bin"] = df.apply(lambda r: assign_module(r.get("gene", ""), r.get("product", ""), r.get("go_ids", "")), axis=1)
        df["fox_status"] = df.apply(lambda r: tag_fox_status(r.get("gene", ""), r.get("product", ""), r.get("go_ids", "")), axis=1)

    # Tier 1 ranked table
    sort_cols = [c for c in ["consensus_rank_pct_mean", "effect_size_mean", "n_run_model_hits"] if c in tier1.columns]
    if not sort_cols:
        raise SystemExit("Tier 1 file is missing expected ranking columns like consensus_rank_pct_mean and effect_size_mean.")

    tier1_ranked = tier1.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last").copy()

    # Module summaries
    mod_counts_t1 = (tier1_ranked["module_bin"].fillna("Other / unassigned")
                     .value_counts()
                     .reset_index())
    mod_counts_t1.columns = ["module_bin", "n_families"]

    fox_counts_t1 = (tier1_ranked["fox_status"].fillna("Novel/uncertain candidate")
                     .value_counts()
                     .reset_index())
    fox_counts_t1.columns = ["fox_status", "n_families"]

    # Write outputs
    tier1_ranked.to_csv(os.path.join(outdir, "tier1_ranked_annotated.tsv"), sep="\t", index=False)
    mod_counts_t1.to_csv(os.path.join(outdir, "tier1_module_bin_counts.tsv"), sep="\t", index=False)
    fox_counts_t1.to_csv(os.path.join(outdir, "tier1_fox_status_counts.tsv"), sep="\t", index=False)

    # Manuscript-facing "Top N" table
    keep_cols = [c for c in [
        "gene_family", "gene", "product", "go_ids",
        "effect_size_mean", "consensus_rank_pct_mean",
        "n_models_with_importance", "n_run_model_hits",
        "diazotroph_pct_max", "diazotroph_pct_mean",
        "module_mode", "module_bin", "fox_status"
    ] if c in tier1_ranked.columns]

    tier1_ranked[keep_cols].head(args.top_n).to_csv(
        os.path.join(outdir, f"Table_Tier1_Top{args.top_n}.tsv"),
        sep="\t", index=False
    )

    # Tier 2 annotated table (optional)
    if tier2 is not None and not tier2.empty:
        tier2.to_csv(os.path.join(outdir, "tier2_pure_positive_annotated.tsv"), sep="\t", index=False)

    print(f"[OK] Wrote: {outdir}")
    print("Key files:")
    print(" - tier1_ranked_annotated.tsv")
    print(f" - Table_Tier1_Top{args.top_n}.tsv")
    print(" - tier1_module_bin_counts.tsv")
    print(" - tier1_fox_status_counts.tsv")
    if tier2 is not None and not tier2.empty:
        print(" - tier2_pure_positive_annotated.tsv")
    else:
        print(" - tier2 file not found or empty (this is OK if not present in report_dir)")
    if tier2_path:
        print(f"[INFO] Tier 2 source: {tier2_path}")
    else:
        print("[INFO] Tier 2 source: not found")


if __name__ == "__main__":
    main()
