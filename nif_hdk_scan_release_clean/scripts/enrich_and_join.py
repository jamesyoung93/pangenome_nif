#!/usr/bin/env python3
import argparse, csv
ap=argparse.ArgumentParser()
ap.add_argument("--hits", required=True)
ap.add_argument("--missing", default="results/proteomes/missing_proteomes.tsv")
ap.add_argument("--quality", default="")
ap.add_argument("--out", required=True)
a=ap.parse_args()
missing=set()
try:
  with open(a.missing,encoding="utf-8") as f:
    for ln in f:
      ln=ln.strip()
      if ln and not ln.lower().startswith(("assembly","accession")):
        missing.add(ln.split()[0])
except FileNotFoundError:
  pass
rows=[]
with open(a.hits,encoding="utf-8") as f:
  r=csv.DictReader(f)
  for row in r:
    acc=row["assembly_accession"]; dl_ok=acc not in missing
    row["download_ok"]="1" if dl_ok else "0"
    if not dl_ok:
      for k in ("nifH_present","nifD_present","nifK_present"):
        if k in row and row[k] in ("0","0.0"): row[k]=""
    core=(row.get("nifH_present")=="1" and row.get("nifD_present")=="1" and row.get("nifK_present")=="1")
    row["nifHDK_core_call"]="1" if core else "0"
    rows.append(row)
if a.quality:
  qual={}
  with open(a.quality,encoding="utf-8") as f:
    rq=csv.DictReader(f,delimiter="\t")
    # Find an accession-like column name; older files used "accession" while
    # newer NCBI outputs sometimes use "assembly_accession".
    acc_col=None
    for cand in ("accession","assembly_accession","assembly", "assembly accession"):
      if cand in (rq.fieldnames or []):
        acc_col=cand; break
    if acc_col is None:
      raise SystemExit(f"No accession column found in {a.quality}; columns={rq.fieldnames}")
    qcols=[c for c in rq.fieldnames if c!=acc_col]
    for q in rq:
      acc_val=q.get(acc_col)
      if not acc_val: continue
      qual[acc_val]=q
  for row in rows:
    q=qual.get(row["assembly_accession"],{})
    for c in qcols: row[c]=q.get(c,"")
cols=list({k for row in rows for k in row.keys()})
with open(a.out,"w",newline="",encoding="utf-8") as o:
  w=csv.DictWriter(o,fieldnames=cols); w.writeheader(); w.writerows(rows)
print("Wrote ->", a.out)
