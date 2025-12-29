#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_merge_annotations_selected.py
- Robustly merge Swiss-Prot hits to UniProt idmapping_selected.tab(.gz)
- Supports: headerless "selected" table, WIDE headered, and LONG 3-col formats
- Backfills GeneName/ProteinName from DIAMOND deflines and entry names
"""

import gzip, re, csv, sys, os
import pandas as pd
import numpy as np

try:
    csv.field_size_limit(sys.maxsize)
except Exception:
    pass

HITS      = "gf_reps_vs_sprot.tsv"          # query, target, evalue, pident, alnlen
INFO_CSV  = "gene_family_info.csv"          # gene_family, representative
DIRCSV    = "feature_directionality_full.csv"
MATRIX    = "gene_family_matrix.csv"
MAPGZ     = "idmapping_selected.tab.gz"
OUT       = "llm_context_features.csv"

ACC_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{4}[0-9]{3})$")
GO_RE  = re.compile(r"\bGO:\d{7}\b", re.I)
EC_RE  = re.compile(r"\bEC[:=]?\s*\d+\.\d+\.\d+\.\d+\b")

def normalize_isoform(tok: str) -> str:
    return re.sub(r"-\d+$", "", (tok or "").strip())

def extract_acc_from_target(raw: str) -> str:
    s = (raw or "").strip()
    if "|" in s:
        parts = [p.strip() for p in s.split("|") if p.strip()]
        for p in parts:
            tok = normalize_isoform(p.split()[0])
            if ACC_RE.match(tok):
                return tok
        if len(parts) > 1:
            return normalize_isoform(parts[1].split()[0])
    tok = normalize_isoform(s.split()[0])
    return tok

def parse_defline_name_and_gene(defline: str):
    s = defline or ""
    entry = ""
    if "|" in s:
        parts = [p for p in s.split("|") if p]
        if len(parts) >= 3:
            entry = parts[2].split()[0]  # e.g., NIFH_AZOVI
    # GN=
    m = re.search(r"\bGN=([A-Za-z0-9_.-]+)", s)
    gene = m.group(1) if m else ""
    if not gene and entry and "_" in entry:
        gene = entry.split("_", 1)[0]
    # recommended protein-like description after entry, up to OS=
    protein = ""
    if "|" in s and len(s.split("|")) >= 3:
        after = s.split("|", 2)[2]
        if " " in after:
            tail = after.split(" ", 1)[1]
            protein = tail.split(" OS=", 1)[0].strip()
    return protein, gene

def load_hits(path: str) -> pd.DataFrame:
    cols = ["query","target","evalue","pident","alnlen"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    df["acc"] = df["target"].map(extract_acc_from_target)
    # derive fallback names from deflines
    protn, genen = [], []
    for s in df["target"].astype(str):
        p, g = parse_defline_name_and_gene(s)
        protn.append(p); genen.append(g)
    df["ProteinName_from_defline"] = protn
    df["GeneName_from_defline"]    = genen
    # best per query
    df = (df.sort_values(["query","evalue","pident"], ascending=[True, True, False])
            .groupby("query", as_index=False).first())
    return df

def rep_to_gf(info_csv: str) -> pd.DataFrame:
    info = pd.read_csv(info_csv)
    need = {"gene_family","representative"}
    if not need.issubset(info.columns):
        raise RuntimeError(f"{info_csv} missing {need}")
    return info.rename(columns={"representative":"rep"})[["rep","gene_family"]]

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore", newline="")
    return open(path, "rt", encoding="utf-8", errors="ignore", newline="")

def stream_uniprot_selected_headerless(path: str, wanted_accs: set) -> pd.DataFrame:
    out = {}
    def ensure(acc):
        return out.setdefault(acc, {"acc": acc, "GO": set(), "EC": set(), "GeneName": set(), "ProteinName": set()})
    with open_maybe_gzip(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for parts in reader:
            if not parts:
                continue
            acc = normalize_isoform((parts[0] if len(parts) > 0 else "").strip())
            if not acc or acc not in wanted_accs:
                continue
            d = ensure(acc)
            entry = (parts[1] if len(parts) > 1 else "").strip()
            if entry:
                if "_" in entry:
                    d["GeneName"].add(entry.split("_", 1)[0].lower())
                d["ProteinName"].add(entry)
            for v in parts[1:]:
                if not v:
                    continue
                for go in GO_RE.findall(v):
                    d["GO"].add(go)
                for ec in EC_RE.findall(v):
                    d["EC"].add(ec.split()[-1])
    if not out:
        return pd.DataFrame(columns=["acc","GO","EC","GeneName","ProteinName"])
    rows = []
    for acc, d in out.items():
        rows.append({
            "acc": acc,
            "GO": ";".join(sorted(d["GO"])) if d["GO"] else "",
            "EC": ";".join(sorted(d["EC"])) if d["EC"] else "",
            "GeneName": ";".join(sorted(d["GeneName"])) if d["GeneName"] else "",
            "ProteinName": ";".join(sorted(d["ProteinName"])) if d["ProteinName"] else "",
        })
    return pd.DataFrame(rows, columns=["acc","GO","EC","GeneName","ProteinName"])

def stream_idmapping_selected(path: str, wanted_accs: set) -> pd.DataFrame:
    with open_maybe_gzip(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        first = next(reader, None)
    if not first:
        return pd.DataFrame(columns=["acc","GO","EC","GeneName","ProteinName"])
    looks_like_ac = bool(first and len(first) and ACC_RE.match(first[0].strip()))
    looks_like_header = any(k in (first or [])
                            for k in ("From","UniProtKB-AC","Entry","Gene Names","Protein names","GO","EC","Type","Database"))
    if looks_like_ac and not looks_like_header:
        return stream_uniprot_selected_headerless(path, wanted_accs)

    # headered WIDE/LONG
    out = {}
    def ensure(acc):
        return out.setdefault(acc, {"acc": acc, "GO": set(), "EC": set(), "GeneName": set(), "ProteinName": set()})
    def add(acc, field, val):
        if not val:
            return
        d = ensure(acc)
        if field in ("GO","EC"):
            for tok in re.split(r"[;,\s]+", val):
                if tok:
                    d[field].add(tok)
        else:
            d[field].add(val)
    def kind_to_field(kind: str):
        k = (kind or "").strip().lower()
        if k == "go" or "gene ontology" in k: return "GO"
        if k.startswith("ec") or "ec number" in k: return "EC"
        if "protein name" in k: return "ProteinName"
        if "gene name" in k or "gene names" in k: return "GeneName"
        if k in ("gene_name","protein_name"): return "GeneName" if k=="gene_name" else "ProteinName"
        return None

    with open_maybe_gzip(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        name = {c.strip(): i for i, c in enumerate(header or [])}
        idx_from = name.get("From") or name.get("Entry") or name.get("UniProtKB-AC")
        idx_to   = name.get("To")
        idx_type = name.get("Type") or name.get("Database")
        wide_markers = (
            "Gene Names","Gene names","Gene names (primary)","Gene_Name","GeneName",
            "Protein names","Protein Names","Protein names (recommended)","Protein name","ProteinName",
            "GO","Gene Ontology (GO)","EC","EC number"
        )
        has_wide = any(k in name for k in wide_markers)

        if idx_from is not None and idx_to is not None and idx_type is not None and not has_wide:
            # LONG 3-col
            for parts in reader:
                if not parts: continue
                if max(idx_from, idx_to, idx_type) >= len(parts): continue
                acc = normalize_isoform(parts[idx_from].strip())
                if acc not in wanted_accs: continue
                fld = kind_to_field(parts[idx_type])
                if not fld: continue
                add(acc, fld, parts[idx_to].strip())
        else:
            # WIDE with headers
            idx_acc  = idx_from if idx_from is not None else (name.get("UniProtKB-AC") or name.get("Entry") or 0)
            idx_go   = name.get("GO") or name.get("Gene Ontology (GO)")
            idx_ec   = name.get("EC") or name.get("EC number")
            idx_gene = next((name[k] for k in ("Gene Names","Gene names","Gene names (primary)","Gene_Name","GeneName") if k in name), None)
            idx_prot = next((name[k] for k in ("Protein names","Protein Names","Protein names (recommended)","Protein name","ProteinName") if k in name), None)
            for parts in reader:
                if not parts: continue
                if idx_acc is None or idx_acc >= len(parts): continue
                acc = normalize_isoform(parts[idx_acc].strip())
                if acc not in wanted_accs: continue
                if idx_go   is not None and idx_go   < len(parts): add(acc,"GO",parts[idx_go].strip())
                if idx_ec   is not None and idx_ec   < len(parts): add(acc,"EC",parts[idx_ec].strip())
                if idx_gene is not None and idx_gene < len(parts): add(acc,"GeneName",parts[idx_gene].strip())
                if idx_prot is not None and idx_prot < len(parts): add(acc,"ProteinName",parts[idx_prot].strip())

    if not out:
        return pd.DataFrame(columns=["acc","GO","EC","GeneName","ProteinName"])
    rows = []
    for acc, d in out.items():
        rows.append({
            "acc": acc,
            "GO": ";".join(sorted(d["GO"])) if d["GO"] else "",
            "EC": ";".join(sorted(d["EC"])) if d["EC"] else "",
            "GeneName": ";".join(sorted(d["GeneName"])) if d["GeneName"] else "",
            "ProteinName": ";".join(sorted(d["ProteinName"])) if d["ProteinName"] else "",
        })
    return pd.DataFrame(rows, columns=["acc","GO","EC","GeneName","ProteinName"])

def main():
    repmap = rep_to_gf(INFO_CSV)
    hits   = load_hits(HITS).merge(repmap, left_on="query", right_on="rep", how="left")

    wanted = set(hits["acc"].dropna().astype(str))
    print(f"[06] Swiss-Prot hits: {len(hits)} ; target accessions extracted: {len(wanted)}")

    ann_small = stream_idmapping_selected(MAPGZ, wanted)
    matched = set(ann_small["acc"]) if len(ann_small) else set()
    print(f"[06] UniProt mapping matches: {len(matched)} of {len(wanted)}")

    ann = hits.merge(ann_small, on="acc", how="left")

    # Backfill from deflines if mapping lacks names
    for col, src in (("ProteinName","ProteinName_from_defline"),
                     ("GeneName","GeneName_from_defline")):
        if col not in ann:
            ann[col] = ""
        need = ann[col].isna() | ann[col].astype(str).str.strip().eq("")
        ann.loc[need, col] = ann.loc[need, src].fillna("")

    # Final fallback: if still empty, use entry name from the defline
    if ann["ProteinName"].astype(str).str.strip().eq("").any():
        entry_from_target = hits["target"].astype(str).map(
            lambda s: s.split("|")[2].split()[0] if ("|" in s and len(s.split("|"))>=3) else ""
        )
        need_pn = ann["ProteinName"].astype(str).str.strip().eq("")
        ann.loc[need_pn, "ProteinName"] = entry_from_target[need_pn].fillna("")
    if ann["GeneName"].astype(str).str.strip().eq("").any():
        entry_from_target = hits["target"].astype(str).map(
            lambda s: (s.split("|")[2].split()[0].split("_",1)[0].lower()
                       if ("|" in s and len(s.split("|"))>=3 and "_" in s.split("|")[2]) else "")
        )
        need_gn = ann["GeneName"].astype(str).str.strip().eq("")
        ann.loc[need_gn, "GeneName"] = entry_from_target[need_gn].fillna("")

    # Directionality + presence rate
    dirfull = pd.read_csv(DIRCSV)
    mat     = pd.read_csv(MATRIX, index_col=0)
    freq    = (mat.sum(axis=0) / mat.shape[0]).rename("presence_rate").reset_index().rename(columns={"index":"gene_family"})

    merged = (dirfull
              .merge(ann[["gene_family","target","acc","pident","evalue","GeneName","ProteinName","GO","EC"]],
                     on="gene_family", how="left")
              .merge(freq, on="gene_family", how="left"))

    if "importance" not in merged.columns:
        merged["importance"] = 0.0
    merged["signed_importance"] = merged["importance"] * np.sign(merged["effect_size"].fillna(0.0))

    cols = [
        "gene_family","signed_importance","importance","direction","effect_size","p_value",
        "presence_rate","n_genomes_with_feature","diaz_rate_present","diaz_rate_absent",
        "target","acc","pident","evalue","GeneName","ProteinName","EC","GO"
    ]
    for c in cols:
        if c not in merged.columns:
            merged[c] = ""
    merged = merged[cols].sort_values(["signed_importance","importance","effect_size"], ascending=[False, False, False])

    merged.to_csv(OUT, index=False)
    print(f"Wrote {OUT} with {len(merged)} rows.")

if __name__ == "__main__":
    main()
