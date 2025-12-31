#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidated Pangenome Analysis Pipeline + Multi-Arm Comparative Experiment Matrix

This script is ASCII-safe (no smart quotes) to avoid Non-UTF-8 parsing errors on HPC.

Core pipeline steps:
  1. Prepare data and label genomes
  2. Download protein sequences
  3. Create gene families using MMseqs2
  4. Classification with blocked CV (genus/family/order)
  5. Feature directionality
  6-8. Optional wrappers

New rigorous experiment mode (multi-arm, matched folds):
  --mode comparative_experiment_matrix

Experiment design (requested):
  - Baseline: purity filtering only (NO genus/family/order breadth filters)
  - Arm 1: purity + min_genera >= 4
  - Arm 2: purity + min_genera >= 4 + min_families >= 3
  - Arm 3: purity + min_genera >= 4 + min_families >= 3 + min_orders >= 2

All arms run under each CV scheme:
  - Order-blocked CV
  - Family-blocked CV
  - Genus-blocked CV

Rigor and reproducibility:
  - CV folds are deterministic (seeded) and exported with participant lists
  - Class purity filtering uses y but is computed within TRAIN fold only (no leakage)
  - Breadth filters are unsupervised and computed within TRAIN fold only (X + taxonomy only)
  - Feature lists per fold per arm are saved (reproducibility and auditing)

Outputs (experiment mode):
  <output_root>/comparative_experiment_matrix/
    A_order/paired_fold_metrics.csv
    B_family/paired_fold_metrics.csv
    C_genus/paired_fold_metrics.csv
    paired_fold_metrics_all.csv
    summary_by_experiment_arm.csv
    paired_deltas_long.csv
    figure1_stratified_metrics.png
    figure2_pooled_delta_<metric>_<arm>.png   (one per metric per non-baseline arm)
    pooled_deltas_<metric>_<arm>.csv
    experiment_manifest.json

"""

import argparse
import os
import sys
import time
import warnings
import random
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional
import subprocess
import json
import ast
from datetime import datetime

import pandas as pd
import numpy as np

# SciPy (paired tests, chi-square in directionality)
try:
    from scipy import stats
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

# ML
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve
)
from sklearn.preprocessing import StandardScaler

# Plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import seaborn as sns
    HAS_SEABORN = True
except Exception:
    HAS_SEABORN = False

warnings.filterwarnings("ignore")

# Optional XGBoost
try:
    import xgboost as xgb
    HAS_XGB = True
except Exception:
    HAS_XGB = False

# Optional NCBI taxonomy lookup (Biopython)
try:
    from Bio import Entrez
    HAS_ENTREZ = True
except Exception:
    HAS_ENTREZ = False


# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

def _fmt_float_tag(x: float) -> str:
    s = f"{x:.2f}"
    return s.replace(".", "p")


class PipelineConfig:
    def __init__(self):
        # Inputs
        self.input_csv = "nif_hdk_hits_enriched_with_quality_checkm.csv"
        self.protein_dir = "downloads/proteins"

        # Outputs
        self.output_root = "results"
        self.run_variants = "both"  # pipeline: baseline|genus|both

        # Step 3 clustering
        self.mmseqs_identity = 0.8
        self.mmseqs_coverage = 0.8
        self.min_genomes_per_family = 40
        self.threads = 8

        # CV settings
        self.cv_folds = 5
        self.cv_seed = 42
        self.cv_group_col = "genus"  # pipeline step4: genus|family|order

        # Filters (pipeline Step 4)
        self.purity_threshold = 0.9
        self.min_genera_per_family = 3
        self.min_families_per_family = 1
        self.min_orders_per_family = 1
        self.max_genus_purity = 0.80  # pipeline genus variant only

        # Step 8
        self.top_n_features = 100

        # NCBI download parameters
        self.email = "james.young@jacks.sdstate.edu"
        self.download_jobs = 8
        self.download_sleep = 0.25
        self.download_retries = 3

        # Taxonomy lookup (Entrez)
        self.taxonomy_cache_path = None
        self.taxonomy_sleep = 0.34

        # Runtime
        self.non_interactive = False

        # Comparative experiment mode
        self.experiment_model = "xgboost"  # xgboost|random_forest|gradient_boosting|logistic_regression
        self.experiment_outdir = None
        self.paired_test = "wilcoxon"      # wilcoxon|ttest
        self.experiment_metrics = ["accuracy", "roc_auc", "f1", "precision", "recall"]


config = PipelineConfig()


# ============================================================================
# STEP 1: PREPARE DATA AND LABEL GENOMES
# ============================================================================

def parse_organism_name(name_str):
    try:
        parsed = ast.literal_eval(name_str) if pd.notna(name_str) and str(name_str).startswith("{") else {}
        return parsed.get("organism_name", "")
    except Exception:
        return ""


def extract_genus(organism_name):
    return organism_name.split()[0] if organism_name else ""


def step1_prepare_data():
    print("\n" + "=" * 80)
    print("STEP 1: DATA PREPARATION AND LABELING")
    print("=" * 80)

    print(f"\nLoading data from: {config.input_csv}")
    if not os.path.exists(config.input_csv):
        print(f"ERROR: Input file not found: {config.input_csv}")
        return False

    df = pd.read_csv(config.input_csv, on_bad_lines="skip", engine="python")
    print(f"Total genomes loaded: {len(df)}")

    if "organism_name" in df.columns:
        df["organism_full"] = df["organism_name"].apply(parse_organism_name)
        df["genus"] = df["organism_full"].apply(extract_genus)
    else:
        df["genus"] = ""

    for col in ["nifH_best_evalue", "nifD_best_evalue", "nifK_best_evalue"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "assembly_level" not in df.columns:
        print("ERROR: Missing 'assembly_level' column.")
        return False

    complete = df[df["assembly_level"] == "Complete Genome"].copy()
    print(f"Complete genomes: {len(complete)}")

    needed = ["nifH_best_evalue", "nifD_best_evalue", "nifK_best_evalue", "assembly_accession"]
    for c in needed:
        if c not in complete.columns:
            print(f"ERROR: Missing required column: {c}")
            return False

    complete["is_diazotroph"] = (
        (complete["nifH_best_evalue"] < 1e-50) &
        (complete["nifD_best_evalue"] < 1e-50) &
        (complete["nifK_best_evalue"] < 1e-50)
    )

    print(f"Diazotrophs: {int(complete['is_diazotroph'].sum())}")
    print(f"Non-diazotrophs: {int((~complete['is_diazotroph']).sum())}")
    print(f"Unique genera: {int(complete['genus'].nunique())}")

    complete.to_csv("complete_genomes_labeled.csv", index=False)
    complete["assembly_accession"].to_csv("genome_accessions.txt", index=False, header=False)

    return True


# ============================================================================
# STEP 2: DOWNLOAD PROTEINS (wrapper)
# ============================================================================

def step2_download_proteins():
    print("\n" + "=" * 80)
    print("STEP 2: DOWNLOAD PROTEIN SEQUENCES")
    print("=" * 80)

    if not os.path.exists("genome_accessions.txt"):
        print("ERROR: genome_accessions.txt not found. Run Step 1 first.")
        return False

    protein_dir = Path(config.protein_dir)
    protein_dir.mkdir(parents=True, exist_ok=True)

    existing = list(protein_dir.glob("*.faa"))
    print(f"Found {len(existing)} existing protein files")

    download_script = None
    search_paths = []
    if "__file__" in globals():
        search_paths.append(Path(__file__).parent.resolve() / "02_download_proteins.py")
    search_paths.append(Path.cwd() / "02_download_proteins.py")
    search_paths.append(Path("02_download_proteins.py"))

    for p in search_paths:
        if p.exists():
            download_script = p
            break

    if download_script:
        cmd = [
            "python3", str(download_script),
            "--nif-csv", config.input_csv,
            "--outdir", "downloads",
            "--jobs", str(config.download_jobs),
            "--sleep", str(config.download_sleep),
            "--retries", str(config.download_retries),
            "--email", config.email
        ]
        print("Running:", " ".join(cmd))
        _ = subprocess.run(cmd)
    else:
        if len(existing) == 0:
            print("ERROR: 02_download_proteins.py not found and no .faa files exist.")
            return False
        print("WARNING: 02_download_proteins.py not found. Using existing .faa files only.")

    return True


# ============================================================================
# STEP 3: CREATE GENE FAMILIES
# ============================================================================

def concatenate_protein_sequences(protein_dir, output_file):
    genome_protein_counts = {}
    with open(output_file, "w") as outf:
        for faa_file in Path(protein_dir).glob("*.faa"):
            genome_id = faa_file.stem
            count = 0
            with open(faa_file, "r") as inf:
                for line in inf:
                    if line.startswith(">"):
                        protein_id = line[1:].split()[0]
                        outf.write(f">{genome_id}|{protein_id}\n")
                        count += 1
                    else:
                        outf.write(line)
            genome_protein_counts[genome_id] = count
    return genome_protein_counts


def run_mmseqs2_clustering(input_fasta, output_prefix, identity=0.8, coverage=0.8, threads=8):
    db_file = f"{output_prefix}_DB"
    subprocess.run(["mmseqs", "createdb", input_fasta, db_file], check=True, capture_output=True)

    cluster_db = f"{output_prefix}_cluster"
    tmp_dir = f"{output_prefix}_tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    cmd = [
        "mmseqs", "cluster",
        db_file, cluster_db, tmp_dir,
        "--min-seq-id", str(identity),
        "-c", str(coverage),
        "--threads", str(threads),
        "--cluster-mode", "0",
        "--cov-mode", "0"
    ]
    subprocess.run(cmd, check=True)

    tsv_file = f"{output_prefix}_clusters.tsv"
    cmd = ["mmseqs", "createtsv", db_file, db_file, cluster_db, tsv_file, "--threads", str(threads)]
    subprocess.run(cmd, check=True)

    subprocess.run(["rm", "-rf", tmp_dir], capture_output=True)
    return tsv_file


def parse_clusters(tsv_file):
    clusters = defaultdict(set)
    with open(tsv_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                rep = parts[0]
                member = parts[1]
                genome_id = member.split("|")[0]
                clusters[rep].add(genome_id)
    cluster_sizes = {rep: len(genomes) for rep, genomes in clusters.items()}
    return clusters, cluster_sizes


def create_presence_absence_matrix(clusters, complete_genomes_df, min_genomes=40):
    all_genomes = set(complete_genomes_df["assembly_accession"])
    filtered_clusters = {rep: genomes for rep, genomes in clusters.items() if len(genomes) >= min_genomes}

    genome_list = sorted(all_genomes)
    cluster_list = sorted(filtered_clusters.keys())

    matrix = np.zeros((len(genome_list), len(cluster_list)), dtype=np.int8)
    for j, rep in enumerate(cluster_list):
        genomes_in_cluster = filtered_clusters[rep]
        for i, genome in enumerate(genome_list):
            if genome in genomes_in_cluster:
                matrix[i, j] = 1

    df_matrix = pd.DataFrame(
        matrix,
        index=genome_list,
        columns=[f"GF_{i:05d}" for i in range(len(cluster_list))]
    )
    return df_matrix, filtered_clusters


def step3_create_gene_families():
    print("\n" + "=" * 80)
    print("STEP 3: GENE FAMILY CLUSTERING")
    print("=" * 80)

    if not os.path.exists(config.protein_dir):
        print(f"ERROR: {config.protein_dir} not found. Run Step 2 first.")
        return False
    if not os.path.exists("complete_genomes_labeled.csv"):
        print("ERROR: complete_genomes_labeled.csv not found. Run Step 1 first.")
        return False

    protein_files = list(Path(config.protein_dir).glob("*.faa"))
    if len(protein_files) == 0:
        print("ERROR: No .faa files found.")
        return False

    complete_genomes = pd.read_csv("complete_genomes_labeled.csv")
    available_accessions = [f.stem for f in protein_files]
    complete_genomes = complete_genomes[complete_genomes["assembly_accession"].isin(available_accessions)].copy()

    min_genomes_per_family = max(10, min(config.min_genomes_per_family, len(complete_genomes) // 10))

    all_proteins_file = "all_proteins.faa"
    if not os.path.exists(all_proteins_file):
        _ = concatenate_protein_sequences(config.protein_dir, all_proteins_file)

    cluster_tsv = "gene_families_clusters.tsv"
    if not os.path.exists(cluster_tsv):
        cluster_tsv = run_mmseqs2_clustering(
            all_proteins_file,
            "gene_families",
            identity=config.mmseqs_identity,
            coverage=config.mmseqs_coverage,
            threads=config.threads
        )

    clusters, _ = parse_clusters(cluster_tsv)
    pa_matrix, filtered_clusters = create_presence_absence_matrix(clusters, complete_genomes, min_genomes=min_genomes_per_family)

    pa_matrix.to_csv("gene_family_matrix.csv")

    cluster_info = pd.DataFrame([
        {"gene_family": f"GF_{i:05d}", "representative": rep, "num_genomes": len(genomes)}
        for i, (rep, genomes) in enumerate(sorted(filtered_clusters.items()))
    ])
    cluster_info.to_csv("gene_family_info.csv", index=False)

    complete_genomes.to_csv("complete_genomes_with_proteins.csv", index=False)
    return True


# ============================================================================
# TAXONOMY HELPERS (family/order)
# ============================================================================

def _detect_taxid_column(df: pd.DataFrame) -> str:
    for c in ["taxid", "tax_id", "TaxId", "organism_taxid", "species_taxid", "taxon_id"]:
        if c in df.columns:
            return c
    return ""


def _parse_taxid_from_organism_field(val):
    if pd.isna(val):
        return np.nan
    s = str(val)
    try:
        if s.startswith("{") and s.endswith("}"):
            d = ast.literal_eval(s)
            if isinstance(d, dict):
                for k in ["taxid", "tax_id", "taxon_id", "TaxId", "taxonId"]:
                    if k in d and d[k] is not None:
                        return int(d[k])
                org = d.get("organism", {})
                if isinstance(org, dict):
                    for k in ["taxid", "tax_id", "taxon_id"]:
                        if k in org and org[k] is not None:
                            return int(org[k])
    except Exception:
        return np.nan
    return np.nan


def build_taxid_map_from_input_csv(input_csv: str) -> Dict[str, int]:
    if not os.path.exists(input_csv):
        return {}

    df = pd.read_csv(input_csv, on_bad_lines="skip", engine="python")
    if "assembly_accession" not in df.columns:
        return {}

    tax_col = _detect_taxid_column(df)
    if tax_col:
        df["__taxid__"] = pd.to_numeric(df[tax_col], errors="coerce")
    elif "organism_name" in df.columns:
        df["__taxid__"] = df["organism_name"].apply(_parse_taxid_from_organism_field)
    else:
        return {}

    tmp = df[["assembly_accession", "__taxid__"]].dropna()
    if len(tmp) == 0:
        return {}

    tmp["__taxid__"] = tmp["__taxid__"].astype(int)
    tmp = tmp.drop_duplicates("assembly_accession", keep="first")
    return dict(zip(tmp["assembly_accession"].astype(str), tmp["__taxid__"].astype(int)))


def ensure_taxid_in_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    if "taxid" in metadata.columns and int(metadata["taxid"].notna().sum()) >= max(10, int(0.8 * len(metadata))):
        return metadata
    tax_map = build_taxid_map_from_input_csv(config.input_csv)
    metadata["taxid"] = metadata.index.to_series().map(tax_map)
    return metadata


def _load_taxonomy_cache(path: Path) -> Dict[str, Dict[str, str]]:
    if path.exists():
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def _save_taxonomy_cache(cache: Dict[str, Dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(cache, f, indent=2)
    tmp.replace(path)


def _fetch_taxonomy_ranks_ncbi(taxid: int, cache: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    key = str(int(taxid))
    if key in cache:
        return cache[key]

    if not HAS_ENTREZ:
        raise RuntimeError("Biopython Entrez not available. Install biopython or provide family/order columns.")

    Entrez.email = config.email
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key

    handle = Entrez.efetch(db="taxonomy", id=key, retmode="xml")
    recs = Entrez.read(handle)
    handle.close()

    if not recs:
        cache[key] = {}
        return cache[key]

    rec = recs[0]
    lineage_ex = rec.get("LineageEx", [])
    ranks = {}

    for node in lineage_ex:
        r = node.get("Rank", "")
        n = node.get("ScientificName", "")
        if r and n:
            ranks[r] = n

    r0 = rec.get("Rank", "")
    n0 = rec.get("ScientificName", "")
    if r0 and n0:
        ranks[r0] = n0

    cache[key] = ranks
    time.sleep(config.taxonomy_sleep)
    return ranks


def ensure_family_order_columns(metadata: pd.DataFrame, ranks=("family", "order")) -> pd.DataFrame:
    if config.taxonomy_cache_path is None:
        config.taxonomy_cache_path = str(Path(config.output_root) / "taxonomy_cache.json")

    cache_path = Path(config.taxonomy_cache_path)
    cache = _load_taxonomy_cache(cache_path)

    metadata = ensure_taxid_in_metadata(metadata)

    for r in ranks:
        if r not in metadata.columns:
            metadata[r] = np.nan

    taxids = metadata["taxid"].dropna().astype(int).unique().tolist()
    if len(taxids) == 0:
        for r in ranks:
            metadata[r] = "Unknown"
        return metadata

    taxid_to_rank = {}
    for t in taxids:
        try:
            taxid_to_rank[str(t)] = _fetch_taxonomy_ranks_ncbi(t, cache)
        except Exception as e:
            print(f"WARNING: taxonomy lookup failed for taxid={t}: {e}")
            taxid_to_rank[str(t)] = {}

    for r in ranks:
        def _rank_for_row(taxid_val):
            if pd.isna(taxid_val):
                return "Unknown"
            d = taxid_to_rank.get(str(int(taxid_val)), {})
            return d.get(r, "Unknown") if d else "Unknown"
        metadata[r] = metadata["taxid"].apply(_rank_for_row)

    _save_taxonomy_cache(cache, cache_path)
    return metadata


def ensure_cv_group_column(metadata: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if group_col == "genus":
        if "genus" not in metadata.columns:
            raise RuntimeError("metadata missing genus column")
        metadata["genus"] = metadata["genus"].fillna("Unknown")
        return metadata

    if group_col in {"family", "order"}:
        metadata = ensure_family_order_columns(metadata, ranks=(group_col,))
        metadata[group_col] = metadata[group_col].fillna("Unknown")
        if metadata[group_col].nunique(dropna=True) <= 1:
            raise RuntimeError(
                f"CV grouping column '{group_col}' could not be populated (only Unknown). "
                f"Need taxids in {config.input_csv} and biopython Entrez."
            )
        return metadata

    raise ValueError(f"Unsupported group_col: {group_col}")


def ensure_family_column(metadata: pd.DataFrame) -> pd.DataFrame:
    if "family" in metadata.columns and metadata["family"].notna().any():
        metadata["family"] = metadata["family"].fillna("Unknown")
        return metadata
    metadata = ensure_family_order_columns(metadata, ranks=("family",))
    metadata["family"] = metadata["family"].fillna("Unknown")
    if metadata["family"].nunique(dropna=True) <= 1:
        raise RuntimeError("family column could not be populated (only Unknown).")
    return metadata


def ensure_order_column(metadata: pd.DataFrame) -> pd.DataFrame:
    if "order" in metadata.columns and metadata["order"].notna().any():
        metadata["order"] = metadata["order"].fillna("Unknown")
        return metadata
    metadata = ensure_family_order_columns(metadata, ranks=("order",))
    metadata["order"] = metadata["order"].fillna("Unknown")
    if metadata["order"].nunique(dropna=True) <= 1:
        raise RuntimeError("order column could not be populated (only Unknown).")
    return metadata


# ============================================================================
# SHARED: METRICS, CV, MODELS
# ============================================================================

def evaluate_model(y_true, y_pred, y_prob) -> Dict[str, float]:
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "roc_auc": roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.0
    }


def _choose_n_folds(meta: pd.DataFrame, group_col: str) -> int:
    n_groups = int(meta[group_col].nunique())
    n_folds = min(config.cv_folds, max(3, n_groups // 3))
    n_folds = min(n_folds, n_groups)
    return max(2, n_folds)


def group_cross_validation_split(metadata_df: pd.DataFrame, group_col="genus", n_folds=5, seed=42) -> List[List[str]]:
    rng = np.random.default_rng(seed)
    meta = metadata_df.copy()

    unique_groups = meta[group_col].dropna().unique().tolist()
    rng.shuffle(unique_groups)

    folds_groups = [[] for _ in range(n_folds)]
    for i, g in enumerate(unique_groups):
        folds_groups[i % n_folds].append(g)

    folds = []
    for fold_groups in folds_groups:
        idx_labels = meta.index[meta[group_col].isin(fold_groups)].unique().tolist()
        folds.append(idx_labels)

    return folds


def export_cv_splits(metadata_df: pd.DataFrame, folds: List[List[str]], outdir: Path, group_col="genus"):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_labels = metadata_df.index.tolist()
    summary = []

    for fold_idx, test_labels in enumerate(folds, start=1):
        test_labels = [lab for lab in test_labels if lab in metadata_df.index]
        test_set = set(test_labels)
        train_labels = [lab for lab in all_labels if lab not in test_set]

        fold_dir = outdir / f"fold_{fold_idx:02d}"
        fold_dir.mkdir(exist_ok=True)

        pd.Series(train_labels, name="assembly_accession").to_csv(fold_dir / "train_genomes.txt", index=False, header=False)
        pd.Series(test_labels, name="assembly_accession").to_csv(fold_dir / "test_genomes.txt", index=False, header=False)

        cols = []
        for c in [group_col, "genus", "family", "order", "is_diazotroph"]:
            if c in metadata_df.columns and c not in cols:
                cols.append(c)

        if cols:
            train_meta = metadata_df.loc[train_labels, cols].copy()
            train_meta.insert(0, "assembly_accession", train_meta.index)
            train_meta.to_csv(fold_dir / "train_metadata.csv", index=False)

            test_meta = metadata_df.loc[test_labels, cols].copy()
            test_meta.insert(0, "assembly_accession", test_meta.index)
            test_meta.to_csv(fold_dir / "test_metadata.csv", index=False)

        summary.append({
            "fold": fold_idx,
            "train_n": len(train_labels),
            "test_n": len(test_labels),
            "train_groups": int(metadata_df.loc[train_labels, group_col].nunique()),
            "test_groups": int(metadata_df.loc[test_labels, group_col].nunique()),
            "train_pos_rate": float(metadata_df.loc[train_labels, "is_diazotroph"].mean()) if "is_diazotroph" in metadata_df.columns else None,
            "test_pos_rate": float(metadata_df.loc[test_labels, "is_diazotroph"].mean()) if "is_diazotroph" in metadata_df.columns else None,
        })

    pd.DataFrame(summary).to_csv(outdir / "fold_summary.csv", index=False)


def _build_model(model_name: str):
    name = (model_name or "").lower().strip()

    if name == "xgboost":
        if not HAS_XGB:
            raise RuntimeError("XGBoost requested but not available.")
        return xgb.XGBClassifier(
            n_estimators=250,
            max_depth=5,
            learning_rate=0.08,
            subsample=0.85,
            colsample_bytree=0.85,
            random_state=42,
            seed=42,
            n_jobs=-1,
            use_label_encoder=False,
            eval_metric="logloss"
        )

    if name == "random_forest":
        return RandomForestClassifier(
            n_estimators=300,
            max_depth=12,
            min_samples_leaf=3,
            random_state=42,
            n_jobs=-1
        )

    if name == "gradient_boosting":
        return GradientBoostingClassifier(
            n_estimators=200,
            max_depth=4,
            learning_rate=0.08,
            random_state=42
        )

    if name == "logistic_regression":
        return LogisticRegression(max_iter=2000, random_state=42)

    raise ValueError(f"Unknown model: {model_name}")


def _predict_proba_01(model, X):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    if hasattr(model, "decision_function"):
        raw = model.decision_function(X)
        return (raw - raw.min()) / (raw.max() - raw.min() + 1e-12)
    pred = model.predict(X)
    return pred.astype(float)


# ============================================================================
# FEATURE FILTERS (purity + breadth)
# ============================================================================

def filter_pure_gene_families_supervised(X: pd.DataFrame, y: pd.Series, threshold: float) -> List[str]:
    present = (X > 0)
    tot = present.sum(axis=0).astype(float).replace(0, np.nan)
    diaz_mask = (y.astype(int) == 1)
    diaz_count = present.loc[diaz_mask].sum(axis=0).astype(float)
    diaz_pct = (diaz_count / tot).fillna(0.5)
    is_pure = (diaz_pct >= threshold) | (diaz_pct <= (1 - threshold))
    keep_cols = diaz_pct.index[~is_pure].tolist()
    return keep_cols


def _n_groups_present_per_feature(X: pd.DataFrame, group_series: pd.Series) -> pd.Series:
    X01 = (X > 0).astype(np.int8)
    g = group_series.fillna("Unknown").astype(str)
    G = pd.get_dummies(g, dtype=np.int8)
    counts_by_group = G.T.dot(X01)
    n_groups = (counts_by_group > 0).sum(axis=0)
    n_groups.index = X.columns
    return n_groups


def _filter_min_genera_unsupervised(X: pd.DataFrame, meta: pd.DataFrame, k: int) -> List[str]:
    if k <= 1:
        return list(X.columns)
    n = _n_groups_present_per_feature(X, meta["genus"])
    return n.index[n >= k].tolist()


def _filter_min_families_unsupervised(X: pd.DataFrame, meta: pd.DataFrame, k: int) -> List[str]:
    if k <= 1:
        return list(X.columns)
    n = _n_groups_present_per_feature(X, meta["family"])
    return n.index[n >= k].tolist()


def _filter_min_orders_unsupervised(X: pd.DataFrame, meta: pd.DataFrame, k: int) -> List[str]:
    if k <= 1:
        return list(X.columns)
    n = _n_groups_present_per_feature(X, meta["order"])
    return n.index[n >= k].tolist()


# ============================================================================
# COMPARATIVE EXPERIMENT MATRIX (multi-arm per CV scheme)
# ============================================================================

def _paired_p_value(baseline_vals: np.ndarray, arm_vals: np.ndarray, test_name: str) -> float:
    if not HAS_SCIPY:
        return np.nan
    if len(baseline_vals) < 2:
        return np.nan
    if np.allclose(arm_vals - baseline_vals, 0):
        return 1.0
    if test_name == "ttest":
        return float(stats.ttest_rel(arm_vals, baseline_vals).pvalue)
    try:
        return float(stats.wilcoxon(arm_vals, baseline_vals).pvalue)
    except Exception:
        return np.nan


def _apply_breadth_filters_unsupervised(
    X_train_cp: pd.DataFrame,
    meta_train: pd.DataFrame,
    min_genera: int,
    min_families: int,
    min_orders: int
) -> List[str]:
    cols = list(X_train_cp.columns)

    if min_genera > 1:
        cols = _filter_min_genera_unsupervised(X_train_cp[cols], meta_train, min_genera)
    if min_families > 1:
        cols = _filter_min_families_unsupervised(X_train_cp[cols], meta_train, min_families)
    if min_orders > 1:
        cols = _filter_min_orders_unsupervised(X_train_cp[cols], meta_train, min_orders)

    return cols


def _run_multiarm_experiment(
    X_full: pd.DataFrame,
    y_full: pd.Series,
    meta: pd.DataFrame,
    group_col: str,
    experiment_name: str,
    arms: List[Dict],
    outdir: Path
) -> pd.DataFrame:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Ensure grouping column exists
    meta = ensure_cv_group_column(meta, group_col)

    # Ensure taxonomy columns needed by arms
    # All arms need genus for min_genera checks; family/order only if thresholds >1
    if "genus" not in meta.columns:
        raise RuntimeError("metadata missing genus column")
    meta["genus"] = meta["genus"].fillna("Unknown")

    need_family = any(a.get("min_families", 1) > 1 for a in arms) or (group_col == "family")
    need_order = any(a.get("min_orders", 1) > 1 for a in arms) or (group_col == "order")

    if need_family:
        meta = ensure_family_column(meta)
    if need_order:
        meta = ensure_order_column(meta)

    n_folds = _choose_n_folds(meta, group_col)
    folds = group_cross_validation_split(meta, group_col=group_col, n_folds=n_folds, seed=config.cv_seed)

    splits_dir = outdir / f"cv_splits_{group_col}_{n_folds}fold_seed_{config.cv_seed}"
    export_cv_splits(meta, folds, splits_dir, group_col=group_col)

    model_name = config.experiment_model
    rows = []
    all_labels = X_full.index.tolist()

    for fold_idx, test_labels in enumerate(folds, start=1):
        test_labels = [lab for lab in test_labels if lab in X_full.index]
        test_set = set(test_labels)
        train_labels = [lab for lab in all_labels if lab not in test_set]

        X_train0 = X_full.loc[train_labels]
        X_test0 = X_full.loc[test_labels]
        y_train = y_full.loc[train_labels].astype(int)
        y_test = y_full.loc[test_labels].astype(int)

        meta_train = meta.loc[train_labels]
        meta_test = meta.loc[test_labels]

        # Train-fold supervised purity filter (common to baseline and all arms)
        purity_cols = filter_pure_gene_families_supervised(X_train0, y_train, threshold=config.purity_threshold)
        X_train_cp = X_train0[purity_cols]
        X_test_cp = X_test0[purity_cols]

        fold_dir = splits_dir / f"fold_{fold_idx:02d}"
        fold_dir.mkdir(exist_ok=True)

        for arm in arms:
            arm_id = arm["arm_id"]
            min_g = int(arm.get("min_genera", 1))
            min_f = int(arm.get("min_families", 1))
            min_o = int(arm.get("min_orders", 1))

            kept_cols = _apply_breadth_filters_unsupervised(
                X_train_cp=X_train_cp,
                meta_train=meta_train,
                min_genera=min_g,
                min_families=min_f,
                min_orders=min_o
            )

            # Save feature list for reproducibility
            pd.Series(kept_cols, name="gene_family").to_csv(fold_dir / f"features_{arm_id}.txt", index=False, header=False)

            # Fit and eval
            model = _build_model(model_name)
            if model_name.lower().strip() == "logistic_regression":
                scaler = StandardScaler()
                X_train_use = scaler.fit_transform(X_train_cp[kept_cols])
                X_test_use = scaler.transform(X_test_cp[kept_cols])
            else:
                X_train_use = X_train_cp[kept_cols]
                X_test_use = X_test_cp[kept_cols]

            model.fit(X_train_use, y_train)
            y_pred = model.predict(X_test_use)
            y_prob = _predict_proba_01(model, X_test_use)

            m = evaluate_model(y_test, y_pred, y_prob)

            rows.append({
                "experiment": experiment_name,
                "group_col": group_col,
                "fold": fold_idx,
                "arm": arm_id,
                "n_folds": n_folds,
                "n_train": len(train_labels),
                "n_test": len(test_labels),
                "train_pos_rate": float(y_train.mean()) if len(y_train) else np.nan,
                "test_pos_rate": float(y_test.mean()) if len(y_test) else np.nan,
                "n_features": int(len(kept_cols)),
                "purity_threshold": float(config.purity_threshold),
                "min_genera": min_g,
                "min_families": min_f,
                "min_orders": min_o,
                **m
            })

        # quick fold log
        fold_sub = pd.DataFrame([r for r in rows if r["experiment"] == experiment_name and r["fold"] == fold_idx])
        if not fold_sub.empty:
            auc_line = ", ".join([f"{a}:{v:.3f}" for a, v in zip(fold_sub["arm"].tolist(), fold_sub["roc_auc"].tolist())])
            print(f"[{experiment_name}] fold {fold_idx}/{n_folds}  {auc_line}")

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "paired_fold_metrics.csv", index=False)

    # Summary means per arm
    metrics = ["accuracy", "roc_auc", "f1", "precision", "recall", "n_features"]
    summ = (df.groupby(["experiment", "arm"])[metrics].agg(["mean", "std", "count"]))
    summ.to_csv(outdir / "paired_summary_by_arm.csv")

    return df


def plot_comparative_inference_multiarm(
    paired_df: pd.DataFrame,
    metrics: List[str],
    outdir: Path,
    paired_test: str,
    arms_order: List[str],
    baseline_arm: str = "baseline"
):
    """
    Figure 1: Stratified multi-metric plot
      - columns: Order CV, Family CV, Genus CV
      - rows: metrics
      - x-axis within each panel: baseline + arms
      - dashed horizontal line: pooled baseline mean for that metric (across all experiments and folds)
      - p-values: per arm vs baseline, paired within fold

    Figure 2: Pooled delta plots per metric per arm (arm - baseline), pooled across experiments and folds.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    exp_order = ["Order CV", "Family CV", "Genus CV"]

    df = paired_df.copy()
    if "experiment" not in df.columns:
        raise ValueError("paired_df missing 'experiment' column")

    # Ensure categorical ordering
    df["experiment"] = pd.Categorical(df["experiment"], categories=exp_order, ordered=True)
    df["arm"] = pd.Categorical(df["arm"], categories=arms_order, ordered=True)

    # Baseline means per metric across all experiments/folds
    baseline_means = {}
    for m in metrics:
        if m in df.columns:
            baseline_means[m] = float(df.loc[df["arm"] == baseline_arm, m].mean())

    # Figure 1 grid
    n_rows = len(metrics)
    n_cols = len(exp_order)

    if HAS_SEABORN:
        sns.set_context("paper", font_scale=1.1)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.6 * n_cols, 2.8 * n_rows), sharey="row")
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = np.array([axes])
    elif n_cols == 1:
        axes = np.array([[ax] for ax in axes])

    for i, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        for j, exp in enumerate(exp_order):
            ax = axes[i, j]
            sub = df[df["experiment"] == exp].copy()
            if sub.empty:
                ax.set_title(exp)
                continue

            # Plot
            if HAS_SEABORN:
                sns.boxplot(data=sub, x="arm", y=metric, ax=ax, showfliers=False)
                sns.stripplot(data=sub, x="arm", y=metric, ax=ax, jitter=0.12, alpha=0.8)
            else:
                # simple matplotlib fallback
                vals = [sub[sub["arm"] == a][metric].values for a in arms_order]
                ax.boxplot(vals, labels=arms_order, showfliers=False)
                for k, a in enumerate(arms_order, start=1):
                    yv = sub[sub["arm"] == a][metric].values
                    ax.scatter(np.zeros_like(yv) + k, yv, alpha=0.7)

            # Baseline mean line
            if metric in baseline_means and not np.isnan(baseline_means[metric]):
                ax.axhline(baseline_means[metric], linestyle="--", linewidth=1.2)

            # Titles / labels
            if i == 0:
                ax.set_title(exp, fontweight="bold")
            ax.set_xlabel("")
            if j == 0:
                ax.set_ylabel(metric)
            else:
                ax.set_ylabel("")

            # Paired p-values per arm vs baseline
            if HAS_SCIPY:
                piv = sub.pivot_table(index="fold", columns="arm", values=metric, aggfunc="mean").dropna()
                p_lines = []
                if baseline_arm in piv.columns:
                    b = piv[baseline_arm].values
                    for arm in arms_order:
                        if arm == baseline_arm:
                            continue
                        if arm in piv.columns:
                            p = _paired_p_value(b, piv[arm].values, paired_test)
                            if not np.isnan(p):
                                p_lines.append(f"{arm} p={p:.3g}")
                if p_lines:
                    ax.text(
                        0.02, 0.98,
                        "\n".join(p_lines),
                        transform=ax.transAxes,
                        ha="left",
                        va="top",
                        fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.2", alpha=0.12)
                    )

            ax.tick_params(axis="x", labelrotation=20)

    plt.tight_layout()
    fig1_path = outdir / "figure1_stratified_metrics.png"
    plt.savefig(fig1_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Wrote Figure 1: {fig1_path}")

    # Figure 2 pooled deltas per metric per arm
    delta_rows = []

    for exp in exp_order:
        sub_exp = df[df["experiment"] == exp].copy()
        for metric in metrics:
            if metric not in sub_exp.columns:
                continue
            piv = sub_exp.pivot_table(index="fold", columns="arm", values=metric, aggfunc="mean").dropna()
            if baseline_arm not in piv.columns:
                continue
            b = piv[baseline_arm]
            for arm in arms_order:
                if arm == baseline_arm:
                    continue
                if arm in piv.columns:
                    d = (piv[arm] - b)
                    for fold_idx, val in d.items():
                        delta_rows.append({
                            "experiment": exp,
                            "fold": int(fold_idx),
                            "metric": metric,
                            "arm": arm,
                            "delta": float(val)
                        })

    delta_df = pd.DataFrame(delta_rows)
    delta_df.to_csv(outdir / "paired_deltas_long.csv", index=False)

    for metric in metrics:
        for arm in arms_order:
            if arm == baseline_arm:
                continue
            sub = delta_df[(delta_df["metric"] == metric) & (delta_df["arm"] == arm)]
            if sub.empty:
                continue

            sub.to_csv(outdir / f"pooled_deltas_{metric}_{arm}.csv", index=False)

            plt.figure(figsize=(8, 5))
            if HAS_SEABORN:
                sns.histplot(sub["delta"].values, bins=12, kde=True)
            else:
                plt.hist(sub["delta"].values, bins=12, alpha=0.85)
            plt.axvline(0.0, linestyle="--", linewidth=1.2)
            plt.xlabel(f"delta = {arm} - {baseline_arm} ({metric})")
            plt.ylabel("count")
            plt.title(f"Pooled Delta ({metric}) - {arm}", fontweight="bold")
            plt.tight_layout()

            fig2_path = outdir / f"figure2_pooled_delta_{metric}_{arm}.png"
            plt.savefig(fig2_path, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Wrote Figure 2: {fig2_path}")


def run_comparative_experiment_matrix():
    if not HAS_SCIPY:
        raise RuntimeError("scipy is required for comparative experiment mode (paired tests).")

    if not os.path.exists("gene_family_matrix.csv"):
        raise RuntimeError("gene_family_matrix.csv not found. Run Step 3 first.")
    if not os.path.exists("complete_genomes_with_proteins.csv"):
        raise RuntimeError("complete_genomes_with_proteins.csv not found. Run Step 3 first.")

    # Reproducibility
    random.seed(config.cv_seed)
    np.random.seed(config.cv_seed)

    X_full = pd.read_csv("gene_family_matrix.csv", index_col=0)
    meta = pd.read_csv("complete_genomes_with_proteins.csv")

    if meta.duplicated("assembly_accession").any():
        meta = meta.drop_duplicates("assembly_accession", keep="first")

    common = sorted(set(meta["assembly_accession"]) & set(X_full.index))
    meta = (meta[meta["assembly_accession"].isin(common)]
            .copy()
            .set_index("assembly_accession")
            .loc[common])
    X_full = X_full.loc[common]

    if "genus" not in meta.columns:
        raise RuntimeError("metadata missing genus column")
    meta["genus"] = meta["genus"].fillna("Unknown")

    y_full = meta["is_diazotroph"].astype(int)

    if config.experiment_outdir is None:
        config.experiment_outdir = str(Path(config.output_root) / "comparative_experiment_matrix")
    outdir = Path(config.experiment_outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if config.taxonomy_cache_path is None:
        config.taxonomy_cache_path = str(Path(config.output_root) / "taxonomy_cache.json")

    # Arms requested by user
    arms = [
        {"arm_id": "baseline", "min_genera": 1, "min_families": 1, "min_orders": 1},
        {"arm_id": "arm1_genus4", "min_genera": 4, "min_families": 1, "min_orders": 1},
        {"arm_id": "arm2_genus4_family3", "min_genera": 4, "min_families": 3, "min_orders": 1},
        {"arm_id": "arm3_genus4_family3_order2", "min_genera": 4, "min_families": 3, "min_orders": 2},
    ]
    arms_order = [a["arm_id"] for a in arms]

    print("\n" + "=" * 80)
    print("COMPARATIVE EXPERIMENT MATRIX (MULTI-ARM, MATCHED FOLDS)")
    print("=" * 80)
    print(f"Model: {config.experiment_model}")
    print(f"Train-fold class purity threshold: {config.purity_threshold}")
    print("Arms:")
    for a in arms:
        print(f"  {a['arm_id']}: min_genera={a['min_genera']} min_families={a['min_families']} min_orders={a['min_orders']}")
    print(f"Output dir: {outdir.resolve()}")
    print("=" * 80)

    df_order = _run_multiarm_experiment(
        X_full=X_full,
        y_full=y_full,
        meta=meta,
        group_col="order",
        experiment_name="Order CV",
        arms=arms,
        outdir=outdir / "A_order"
    )

    df_family = _run_multiarm_experiment(
        X_full=X_full,
        y_full=y_full,
        meta=meta,
        group_col="family",
        experiment_name="Family CV",
        arms=arms,
        outdir=outdir / "B_family"
    )

    df_genus = _run_multiarm_experiment(
        X_full=X_full,
        y_full=y_full,
        meta=meta,
        group_col="genus",
        experiment_name="Genus CV",
        arms=arms,
        outdir=outdir / "C_genus"
    )

    all_df = pd.concat([df_order, df_family, df_genus], ignore_index=True)
    all_df.to_csv(outdir / "paired_fold_metrics_all.csv", index=False)

    # Summary table (means and stds per experiment and arm)
    metrics = ["accuracy", "roc_auc", "f1", "precision", "recall", "n_features"]
    summary = (all_df.groupby(["experiment", "arm"])[metrics]
               .agg(["mean", "std", "count"])
               .reset_index())
    summary.to_csv(outdir / "summary_by_experiment_arm.csv", index=False)

    # Figures and pooled deltas
    plot_comparative_inference_multiarm(
        paired_df=all_df,
        metrics=config.experiment_metrics,
        outdir=outdir,
        paired_test=config.paired_test,
        arms_order=arms_order,
        baseline_arm="baseline"
    )

    with open(outdir / "experiment_manifest.json", "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "model": config.experiment_model,
            "paired_test": config.paired_test,
            "metrics": config.experiment_metrics,
            "class_purity_threshold_train_only": config.purity_threshold,
            "cv_seed": config.cv_seed,
            "arms": arms,
            "taxonomy_cache_path": config.taxonomy_cache_path
        }, f, indent=2)

    print("Comparative experiment matrix complete.")
    print(f"Results: {outdir.resolve()}")


# ============================================================================
# PIPELINE STEP 5 and OPTIONAL WRAPPERS
# (kept minimal here; experiment mode is the main focus)
# ============================================================================

def step5_feature_directionality():
    print("Step 5 directionality is not used in experiment mode in this script build.")
    print("If you need it, run your pipeline mode Step 5 from prior pipeline builds.")
    return True


def step6_merge_annotations():
    if not os.path.exists("gf_reps_vs_sprot.tsv"):
        return True
    if os.path.exists("06_merge_annotations_selected.py"):
        _ = subprocess.run(["python3", "06_merge_annotations_selected.py"])
    return True


def step7_expand_gene_families():
    if os.path.exists("07_expand_gene_families.py"):
        _ = subprocess.run(["python3", "07_expand_gene_families.py"])
    return True


def step8_build_narrative_table():
    if os.path.exists("08_build_narrative_table.py"):
        env = os.environ.copy()
        env["TOP_N"] = str(config.top_n_features)
        _ = subprocess.run(["python3", "08_build_narrative_table.py"], env=env)
    return True


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Consolidated Pangenome Pipeline + Multi-Arm Experiment Matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    def _parse_mode(value: str) -> str:
        mode = (value or "").strip().lower()
        if mode == "full":  # backwards-compatible alias from older wrappers
            return "pipeline"
        if mode in {"pipeline", "comparative_experiment_matrix"}:
            return mode
        raise argparse.ArgumentTypeError(
            "Mode must be 'pipeline' (alias 'full') or 'comparative_experiment_matrix'"
        )

    parser.add_argument("--mode", type=_parse_mode, default="pipeline",
                        help="Run standard pipeline steps (pipeline, alias full) or the comparative experiment matrix")

    parser.add_argument("--input", type=str, default="nif_hdk_hits_enriched_with_quality_checkm.csv",
                        help="Input CSV file")
    parser.add_argument("--threads", type=int, default=8, help="Threads for MMseqs2")
    parser.add_argument("--min-genomes", type=int, default=40, help="Min genomes per gene family for Step 3")
    parser.add_argument("--mmseqs-identity", type=float, default=config.mmseqs_identity,
                        help="MMseqs2 minimum sequence identity")
    parser.add_argument("--mmseqs-coverage", type=float, default=config.mmseqs_coverage,
                        help="MMseqs2 minimum coverage fraction")
    parser.add_argument("--purity-threshold", type=float, default=config.purity_threshold,
                        help="Class purity threshold for blocked CV")
    parser.add_argument("--cv-folds", type=int, default=config.cv_folds, help="Number of CV folds")
    parser.add_argument("--cv-seed", type=int, default=config.cv_seed, help="Random seed for CV split")
    parser.add_argument("--cv-group-col", type=str, default=config.cv_group_col,
                        choices=["genus", "family", "order"],
                        help="Blocked CV grouping column")
    parser.add_argument("--min-genera", type=int, default=config.min_genera_per_family,
                        help="Minimum genera per family (breadth filter)")
    parser.add_argument("--min-families", type=int, default=config.min_families_per_family,
                        help="Minimum families per gene family (breadth filter)")
    parser.add_argument("--min-orders", type=int, default=config.min_orders_per_family,
                        help="Minimum orders per gene family (breadth filter)")
    parser.add_argument("--max-genus-purity", type=float, default=config.max_genus_purity,
                        help="Maximum genus purity allowed in genus variant")
    parser.add_argument("--top-n-features", type=int, default=config.top_n_features,
                        help="Top-N features exported in Step 8")
    parser.add_argument("--start-step", type=int, default=1, help="Start step (1-8)")
    parser.add_argument("--end-step", type=int, default=8, help="End step (1-8)")
    parser.add_argument("--skip-download", action="store_true", help="Skip Step 2 download")
    parser.add_argument("--non-interactive", action="store_true", help="Non-interactive mode")

    parser.add_argument("--output-root", type=str, default="results", help="Output root directory")

    parser.add_argument("--download-jobs", type=int, default=config.download_jobs,
                        help="Parallel downloads for Step 2")

    # Taxonomy
    parser.add_argument("--taxonomy-cache", type=str, default=None, help="Path to taxonomy cache JSON")
    parser.add_argument("--taxonomy-sleep", type=float, default=0.34, help="Sleep between Entrez calls")

    # Experiment mode args
    parser.add_argument("--experiment-model", type=str, default="xgboost",
                        choices=["xgboost", "random_forest", "gradient_boosting", "logistic_regression"],
                        help="Model used for experiment mode")
    parser.add_argument("--experiment-outdir", type=str, default=None,
                        help="Output dir for experiment mode (default: <output-root>/comparative_experiment_matrix)")
    parser.add_argument("--paired-test", type=str, default="wilcoxon",
                        choices=["wilcoxon", "ttest"],
                        help="Paired test for arm-vs-baseline p-values")
    parser.add_argument("--experiment-metrics", type=str, default="accuracy,roc_auc,f1,precision,recall",
                        help="Comma-separated metrics to plot")

    args = parser.parse_args()

    config.input_csv = args.input
    config.threads = args.threads
    config.min_genomes_per_family = args.min_genomes
    config.mmseqs_identity = args.mmseqs_identity
    config.mmseqs_coverage = args.mmseqs_coverage
    config.purity_threshold = args.purity_threshold
    config.cv_folds = args.cv_folds
    config.cv_seed = args.cv_seed
    config.cv_group_col = args.cv_group_col
    config.min_genera_per_family = args.min_genera
    config.min_families_per_family = args.min_families
    config.min_orders_per_family = args.min_orders
    config.max_genus_purity = args.max_genus_purity
    config.top_n_features = args.top_n_features
    config.download_jobs = args.download_jobs
    config.non_interactive = args.non_interactive

    config.output_root = args.output_root
    Path(config.output_root).mkdir(parents=True, exist_ok=True)

    config.taxonomy_sleep = args.taxonomy_sleep
    if args.taxonomy_cache is not None:
        config.taxonomy_cache_path = args.taxonomy_cache
    else:
        config.taxonomy_cache_path = str(Path(config.output_root) / "taxonomy_cache.json")

    config.experiment_model = args.experiment_model
    config.experiment_outdir = args.experiment_outdir
    config.paired_test = args.paired_test
    config.experiment_metrics = [m.strip() for m in args.experiment_metrics.split(",") if m.strip()]

    if args.mode == "comparative_experiment_matrix":
        run_comparative_experiment_matrix()
        return

    # Standard pipeline execution (steps 1-3 minimal here)
    steps = [
        (1, "Prepare Data", step1_prepare_data),
        (2, "Download Proteins", step2_download_proteins),
        (3, "Create Gene Families", step3_create_gene_families),
        (5, "Feature Directionality", step5_feature_directionality),
        (6, "Merge Annotations", step6_merge_annotations),
        (7, "Expand Gene Families", step7_expand_gene_families),
        (8, "Build Narrative Table", step8_build_narrative_table),
    ]

    for step_num, step_name, step_func in steps:
        if step_num < args.start_step or step_num > args.end_step:
            continue
        if args.skip_download and step_num == 2:
            print(f"Skipping Step {step_num}: {step_name}")
            continue
        print(f"\n--- Running Step {step_num}: {step_name} ---")
        ok = step_func()
        if not ok:
            print(f"ERROR: Step {step_num} failed.")
            sys.exit(1)


if __name__ == "__main__":
    main()
