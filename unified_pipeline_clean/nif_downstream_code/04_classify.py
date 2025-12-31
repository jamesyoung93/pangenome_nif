#!/usr/bin/env python3
"""
Step 4: Classification with genus-level cross-validation
- Remove highly pure gene families (>90% or <10% diazotroph)
- Run multiple classification models
- Use genus-level cross-validation
- Report performance metrics
"""

import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')  # For HPC without display
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve
)
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

# Optional XGBoost
try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("Warning: XGBoost not available, will skip XGBoost models")


# ------------------------- Utilities -------------------------

def _safe_select(df_or_ser, idx):
    """
    Select rows by either integer positions or index labels.
    If idx items are ints -> use iloc, else assume labels and use loc.
    """
    if len(idx) == 0:
        return df_or_ser.iloc[[]]
    first = idx[0]
    if isinstance(first, (int, np.integer)):
        return df_or_ser.iloc[idx]
    else:
        return df_or_ser.loc[idx]


def filter_pure_gene_families(X, y, threshold=0.9):
    """
    Remove gene families that are highly pure (>threshold diazotroph or <1-threshold diazotroph)
    Returns filtered X and a stats dataframe that is also written to CSV.
    """
    print(f"\nFiltering gene families (removing >={threshold*100:.0f}% or <={100-threshold*100:.0f}% diazotroph)...")

    family_stats = []
    for col in X.columns:
        present = X[col] == 1
        total_count = int(present.sum())
        if total_count == 0:
            continue
        diaz_count = int(y[present].sum())
        diaz_pct = diaz_count / total_count
        family_stats.append({
            'gene_family': col,
            'total_genomes': total_count,
            'diazotroph_count': diaz_count,
            'diazotroph_pct': diaz_pct,
            'is_pure': (diaz_pct >= threshold) or (diaz_pct <= (1 - threshold))
        })

    stats_df = pd.DataFrame(family_stats).sort_values('gene_family')
    pure_families = stats_df.loc[stats_df['is_pure'], 'gene_family'].tolist()
    informative_families = stats_df.loc[~stats_df['is_pure'], 'gene_family'].tolist()

    print(f"  Total gene families: {len(stats_df)}")
    print(f"  Pure families (removed): {len(pure_families)}")
    print(f"  Informative families (kept): {len(informative_families)}")

    # Save purity stats
    stats_df.to_csv('gene_family_purity_stats.csv', index=False)
    print("  Saved purity statistics to: gene_family_purity_stats.csv")

    # Filter X
    X_filtered = X[informative_families]

    return X_filtered, stats_df


def genus_cross_validation_split(metadata_df, n_folds=5, seed=42):
    """
    Create genus-level cross-validation folds.
    Ensures all genomes from the same genus are in the same fold.
    Returns a list of lists containing *index labels* of metadata_df (assembly_accession).
    """
    rng = np.random.default_rng(seed)

    meta = metadata_df.copy()
    # unique genera
    unique_genera = meta['genus'].dropna().unique().tolist()
    rng.shuffle(unique_genera)

    # round-robin to balance
    folds_genera = [[] for _ in range(n_folds)]
    for i, g in enumerate(unique_genera):
        folds_genera[i % n_folds].append(g)

    folds = []
    for fold_genus in folds_genera:
        # use index labels; ensure uniqueness
        idx_labels = meta.index[meta['genus'].isin(fold_genus)].unique().tolist()
        folds.append(idx_labels)

    return folds


def evaluate_model(y_true, y_pred, y_prob):
    """Calculate evaluation metrics."""
    metrics = {
        'accuracy': accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall': recall_score(y_true, y_pred, zero_division=0),
        'f1': f1_score(y_true, y_pred, zero_division=0),
        'roc_auc': roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.0
    }
    return metrics


def run_genus_cv(X, y, metadata_df, model, model_name, n_folds=5, seed=42):
    """
    Run genus-level cross-validation for a given model.
    Uses label-based selection when folds contain accession IDs.
    """
    print(f"\n{'='*60}")
    print(f"Running {model_name} with Genus-level CV")
    print(f"{'='*60}")

    # Build folds (each is a list of accession labels)
    folds = genus_cross_validation_split(metadata_df, n_folds=n_folds, seed=seed)

    fold_metrics = []
    all_y_true = []
    all_y_prob = []
    feature_importances = []

    # Use the accession labels present in X
    all_labels = X.index.tolist()

    for fold_idx, test_labels in enumerate(folds, start=1):
        # keep only labels that actually exist in X
        test_labels = [lab for lab in test_labels if lab in X.index]
        test_set = set(test_labels)
        train_labels = [lab for lab in all_labels if lab not in test_set]

        # Slice safely (labels) and align y on X indices
        X_train = _safe_select(X, train_labels)
        X_test  = _safe_select(X, test_labels)
        y_train = y.loc[X_train.index]
        y_test  = y.loc[X_test.index]

        print(f"\nFold {fold_idx}/{n_folds}")
        print(f"  Train: {len(X_train)} genomes, Test: {len(X_test)} genomes")
        print(f"  Train diazotroph: {int(y_train.sum())}/{len(y_train)} ({y_train.mean()*100:.1f}%)")
        print(f"  Test  diazotroph: {int(y_test.sum())}/{len(y_test)} ({y_test.mean()*100:.1f}%)")

        # Scale only for LR/SVM-like models
        if model_name in ['Logistic Regression', 'SVM']:
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
        else:
            X_train_scaled = X_train
            X_test_scaled = X_test

        # Train & predict
        model.fit(X_train_scaled, y_train)
        y_pred = model.predict(X_test_scaled)

        if hasattr(model, 'predict_proba'):
            y_prob = model.predict_proba(X_test_scaled)[:, 1]
        else:
            raw = model.decision_function(X_test_scaled)
            y_prob = (raw - raw.min()) / (raw.max() - raw.min() + 1e-12)

        # Metrics
        metrics = evaluate_model(y_test, y_pred, y_prob)
        fold_metrics.append(metrics)

        print(f"  Accuracy: {metrics['accuracy']:.4f}")
        print(f"  Precision: {metrics['precision']:.4f}")
        print(f"  Recall: {metrics['recall']:.4f}")
        print(f"  F1: {metrics['f1']:.4f}")
        print(f"  ROC-AUC: {metrics['roc_auc']:.4f}")

        # Aggregate for global ROC
        all_y_true.extend(y_test.tolist())
        all_y_prob.extend(y_prob.tolist())

        # Per-model feature importance
        if hasattr(model, 'feature_importances_'):
            feature_importances.append(model.feature_importances_)
        elif hasattr(model, 'coef_'):
            feature_importances.append(np.abs(model.coef_[0]))

    # Aggregate fold metrics
    metrics_df = pd.DataFrame(fold_metrics)
    mean_metrics = metrics_df.mean()
    std_metrics = metrics_df.std()

    print(f"\n{'='*60}")
    print(f"{model_name} - AGGREGATE RESULTS")
    print(f"{'='*60}")
    for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']:
        print(f"{metric.upper():12s}: {mean_metrics[metric]:.4f} +/- {std_metrics[metric]:.4f}")

    # Aggregate feature importance across folds if present
    if feature_importances:
        avg_importance = np.mean(feature_importances, axis=0)
        feature_importance_df = pd.DataFrame({
            'gene_family': X.columns,
            'importance': avg_importance
        }).sort_values('importance', ascending=False)
    else:
        feature_importance_df = None

    return {
        'fold_metrics': metrics_df,
        'mean_metrics': mean_metrics,
        'std_metrics': std_metrics,
        'all_y_true': np.array(all_y_true),
        'all_y_prob': np.array(all_y_prob),
        'feature_importance': feature_importance_df
    }


def plot_roc_curves(results_dict, output_file='roc_curves.png'):
    """Plot ROC curves for all models (using pooled predictions)."""
    plt.figure(figsize=(10, 8))

    for model_name, results in results_dict.items():
        fpr, tpr, _ = roc_curve(results['all_y_true'], results['all_y_prob'])
        auc = results['mean_metrics']['roc_auc']
        auc_std = results['std_metrics']['roc_auc']
        plt.plot(fpr, tpr, label=f"{model_name} (AUC = {auc:.3f} +/- {auc_std:.3f})", linewidth=2)

    plt.plot([0, 1], [0, 1], 'k--', label='Random Classifier', linewidth=1)
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curves - Genus-Level Cross-Validation', fontsize=14, fontweight='bold')
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved ROC curves to: {output_file}")


def analyze_feature_importance(results_dict, top_n=20):
    """Analyze and save feature importance across models."""
    all_importances = []

    for model_name, results in results_dict.items():
        fi = results['feature_importance']
        if fi is not None:
            tmp = fi.copy()
            tmp['model'] = model_name
            all_importances.append(tmp)

    if not all_importances:
        print("No feature importances available")
        return

    combined = pd.concat(all_importances, ignore_index=True)

    # Average across models
    avg_importance = (combined.groupby('gene_family')['importance']
                      .mean().sort_values(ascending=False))
    top_features = avg_importance.head(top_n).index

    # Pivot per model for the top features
    pivot = (combined[combined['gene_family'].isin(top_features)]
             .pivot(index='gene_family', columns='model', values='importance')
             .fillna(0.0)
             .loc[top_features])

    pivot.to_csv('feature_importance_top.csv')
    print(f"\nSaved top {top_n} feature importances to: feature_importance_top.csv")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))
    pivot.plot(kind='barh', ax=ax)
    ax.set_xlabel('Feature Importance', fontsize=12)
    ax.set_ylabel('Gene Family', fontsize=12)
    ax.set_title(f'Top {top_n} Most Important Gene Families', fontsize=14, fontweight='bold')
    ax.legend(title='Model', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('feature_importance_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved feature importance plot to: feature_importance_plot.png")


# ------------------------- Main -------------------------

def main():
    print("="*80)
    print("STEP 4: CLASSIFICATION WITH GENUS-LEVEL CROSS-VALIDATION")
    print("="*80)

    # Load presence/absence
    print("\nLoading data...")
    gene_family_matrix = pd.read_csv('gene_family_matrix.csv', index_col=0)

    # Metadata: prefer filtered set with proteins
    if os.path.exists('complete_genomes_with_proteins.csv'):
        metadata = pd.read_csv('complete_genomes_with_proteins.csv')
        print("  Using genomes with downloaded proteins")
    else:
        metadata = pd.read_csv('complete_genomes_labeled.csv')
        print("  Using all labeled genomes (assuming all have proteins)")

    # --- DEDUPLICATE by assembly_accession BEFORE intersecting ---
    if metadata.duplicated('assembly_accession').any():
        dup_count = int(metadata.duplicated('assembly_accession').sum())
        print(f"  Note: de-duplicating metadata on assembly_accession (dropping {dup_count} duplicates)")
        metadata = metadata.drop_duplicates('assembly_accession', keep='first')

    # Determine shared accessions and apply a deterministic order
    common = sorted(set(metadata['assembly_accession']) & set(gene_family_matrix.index))
    print(f"  Common genomes in both datasets: {len(common)}")

    if len(common) < 50:
        print("\nWARNING: Very few genomes available for classification!")
        print("Consider downloading more proteins or checking data quality.")

    # Align both X and metadata to the same order & index
    metadata = (metadata[metadata['assembly_accession'].isin(common)]
                .copy()
                .set_index('assembly_accession')
                .loc[common])
    # enforce unique index
    metadata = metadata[~metadata.index.duplicated(keep='first')]
    X_full = gene_family_matrix.loc[common]
    X_full = X_full[~X_full.index.duplicated(keep='first')]

    # If any index mismatch remains, realign to intersection again
    final_common = X_full.index.intersection(metadata.index)
    if len(final_common) != len(common):
        print(f"  Adjusted for uniqueness: {len(final_common)} genomes")
    X_full = X_full.loc[final_common]
    metadata = metadata.loc[final_common]

    print(f"  Gene families: {X_full.shape[1]}")
    print(f"  Genomes: {X_full.shape[0]}")
    print(f"  Diazotrophs: {int(metadata['is_diazotroph'].sum())} ({metadata['is_diazotroph'].mean()*100:.1f}%)")
    print(f"  Unique genera: {metadata['genus'].nunique()}")

    if metadata['genus'].nunique() < 5:
        print("\nWARNING: Very few genera available!")
        print("Genus-level cross-validation may not work well. Consider fewer folds.")

    # Remove highly pure families
    X, purity_stats = filter_pure_gene_families(
        X_full,
        metadata['is_diazotroph'],
        threshold=0.9
    )

    y = metadata['is_diazotroph'].astype(int)

    def filter_min_genera(X, meta, k=3):
        keep = []
        for gf in X.columns:
            g = meta.loc[X[gf] == 1, "genus"].dropna().nunique()
            if g >= k:
                keep.append(gf)
        return X[keep]

    X = filter_min_genera(X, metadata, k=3)

    print(f"\nFinal feature matrix: {X.shape}")

    if X.shape[1] < 10:
        print("\nWARNING: Very few informative gene families!")
        print("Consider lowering purity threshold or checking data quality.")

    # Choose CV folds based on number of genera
    n_genera = metadata['genus'].nunique()
    n_folds = min(5, max(3, n_genera // 3))
    print(f"\nUsing {n_folds}-fold cross-validation ({n_genera} genera)")

    # Define models
    models = {
        'Random Forest': RandomForestClassifier(
            n_estimators=100, max_depth=10, min_samples_leaf=5,
            random_state=42, n_jobs=-1
        ),
        'Gradient Boosting': GradientBoostingClassifier(
            n_estimators=100, max_depth=5, learning_rate=0.1,
            random_state=42
        ),
        'Logistic Regression': LogisticRegression(
            max_iter=1000, random_state=42
        ),
    }
    if HAS_XGB:
        models['XGBoost'] = xgb.XGBClassifier(
            n_estimators=200, max_depth=5, learning_rate=0.1,
            subsample=0.8, colsample_bytree=0.8,
            random_state=42, n_jobs=-1, use_label_encoder=False,
            eval_metric='logloss'
        )

    # Run CV
    results_dict = {}
    for model_name, model in models.items():
        results = run_genus_cv(X, y, metadata, model, model_name, n_folds=n_folds, seed=42)
        results_dict[model_name] = results

    # Summary table
    print("\n" + "="*80)
    print("SUMMARY OF ALL MODELS")
    print("="*80)

    summary_rows = []
    for model_name, results in results_dict.items():
        row = {'Model': model_name}
        for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']:
            mean = results['mean_metrics'][metric]
            std = results['std_metrics'][metric]
            row[metric] = f"{mean:.4f} +/- {std:.4f}"
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))
    summary_df.to_csv('classification_summary.csv', index=False)
    print("\nSaved summary to: classification_summary.csv")

    # Plot ROC
    plot_roc_curves(results_dict, output_file='roc_curves.png')

    # Feature importance overview
    analyze_feature_importance(results_dict, top_n=30)

    # Save per-model fold metrics and importances
    for model_name, results in results_dict.items():
        fold_df = results['fold_metrics'].copy()
        fold_df['model'] = model_name
        fold_df.to_csv(f"fold_metrics_{model_name.replace(' ', '_')}.csv", index=False)

        if results['feature_importance'] is not None:
            results['feature_importance'].to_csv(
                f"feature_importance_{model_name.replace(' ', '_')}.csv", index=False
            )

    print("\n" + "="*80)
    print("STEP 4 COMPLETE")
    print("="*80)
    print("\nAll results saved in current directory.")


if __name__ == '__main__':
    main()
