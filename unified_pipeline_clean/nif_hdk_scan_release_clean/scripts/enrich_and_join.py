#!/usr/bin/env python3
import argparse, csv, re


def detect_delimiter(path):
  """Return a delimiter from {",", "\t"} using a small sample."""
  with open(path, encoding="utf-8") as f:
    first_line = f.readline()
    sample = first_line + f.read(1024)
  if "\t" in first_line:
    return "\t", first_line
  if "," in first_line:
    return ",", first_line
  try:
    sniffed = csv.Sniffer().sniff(sample, delimiters=",\t")
    return sniffed.delimiter, first_line
  except Exception:
    return ",", first_line


def normalize_header(name):
  return re.sub(r"[\s_-]+", "", (name or "").strip().lower())


def find_accession_key(fieldnames):
  wanted = {"accession", "assemblyaccession"}
  for fn in fieldnames or []:
    if normalize_header(fn) in wanted:
      return fn
  raise ValueError(f"Could not find accession column in headers: {fieldnames}")
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
  delimiter, header_line = detect_delimiter(a.quality)
  with open(a.quality,encoding="utf-8") as f:
    rq=csv.DictReader(f,delimiter=delimiter)
    if rq.fieldnames and len(rq.fieldnames)==1:
      raise ValueError(f"Only one column found when parsing {a.quality} (delimiter '{delimiter}'). Header: {rq.fieldnames[0]}")
    acc_key=find_accession_key(rq.fieldnames)
    qcols=[c for c in rq.fieldnames if c!=acc_key]
    for q in rq:
      acc_val=q.get(acc_key,"")
      if not acc_val: continue
      qual[acc_val]=q
  print(f"Detected delimiter '{delimiter}' for {a.quality} (header starts: {header_line.strip()})")
  print(f"Using accession column '{acc_key}' from quality file")
  for row in rows:
    q=qual.get(row["assembly_accession"],{})
    for c in qcols: row[c]=q.get(c,"")
cols=list({k for row in rows for k in row.keys()})
with open(a.out,"w",newline="",encoding="utf-8") as o:
  w=csv.DictWriter(o,fieldnames=cols); w.writeheader(); w.writerows(rows)
print("Wrote ->", a.out)
