#!/usr/bin/env python3
import argparse
import io
import json
import shutil
import sys
import time
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import pandas as pd

API = "https://rest.uniprot.org"

def chunked(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def http_post_form(url, data_dict, timeout=60):
    data = urlencode(data_dict).encode("utf-8")
    req = Request(url, data=data, method="POST")
    with urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))

def http_get_json(url, timeout=60):
    with urlopen(url, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))

def download_to_file(url, out_path, timeout=300):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=timeout) as r, open(out_path, "wb") as f:
        shutil.copyfileobj(r, f)

def submit_and_wait_idmapping(from_db, to_db, ids, poll_seconds=3, timeout_seconds=3600):
    submit = http_post_form(f"{API}/idmapping/run", {"from": from_db, "to": to_db, "ids": ",".join(ids)})
    job_id = submit.get("jobId")
    if not job_id:
        raise RuntimeError(f"Unexpected submit response (no jobId): {submit}")

    t0 = time.time()
    while True:
        st = http_get_json(f"{API}/idmapping/status/{job_id}")
        job_status = st.get("jobStatus") or st.get("status")

        # Common statuses: RUNNING / FINISHED (varies slightly by API version)
        if job_status in (None, "", "FINISHED", "COMPLETED", "SUCCESS"):
            return job_id

        if job_status in ("RUNNING", "NEW"):
            if time.time() - t0 > timeout_seconds:
                raise TimeoutError(f"Timeout waiting for job {job_id} (last status: {job_status})")
            time.sleep(poll_seconds)
            continue

        if job_status in ("ERROR", "FAILED"):
            raise RuntimeError(f"Job {job_id} failed: {st}")

        # Fallback: if we cannot interpret status, try waiting a bit
        if time.time() - t0 > timeout_seconds:
            raise TimeoutError(f"Timeout waiting for job {job_id} (last payload: {st})")
        time.sleep(poll_seconds)

def normalize_uniprot_go_df(df):
    # UniProt TSV headers are human-readable; normalize key columns.
    col_l = {c.lower(): c for c in df.columns}

    acc_col = None
    for cand in ("entry", "accession"):
        if cand in col_l:
            acc_col = col_l[cand]
            break
    if acc_col is None:
        raise RuntimeError(f"Could not find UniProt accession column in headers: {list(df.columns)}")

    def find_col(substr):
        for c in df.columns:
            if substr in c.lower():
                return c
        return None

    out = pd.DataFrame()
    out["uniprot_accession"] = df[acc_col].astype(str)

    go_ids_col = find_col("gene ontology ids")
    if go_ids_col:
        out["go_ids"] = df[go_ids_col].fillna("").astype(str)
    else:
        out["go_ids"] = ""

    # Optional breakdowns if present
    bp_col = find_col("biological process")
    cc_col = find_col("cellular component")
    mf_col = find_col("molecular function")
    go_col = None
    # Sometimes there is a generic "Gene Ontology (GO)" or just "Gene ontology"
    for c in df.columns:
        cl = c.lower()
        if ("gene ontology" in cl) and ("ids" not in cl) and ("biological process" not in cl) and ("cellular component" not in cl) and ("molecular function" not in cl):
            go_col = c
            break

    if bp_col: out["go_bp"] = df[bp_col].fillna("").astype(str)
    if cc_col: out["go_cc"] = df[cc_col].fillna("").astype(str)
    if mf_col: out["go_mf"] = df[mf_col].fillna("").astype(str)
    if go_col: out["go_all"] = df[go_col].fillna("").astype(str)

    reviewed_col = None
    for cand in ("reviewed",):
        if cand in col_l:
            reviewed_col = col_l[cand]
            break
    if reviewed_col:
        out["reviewed"] = df[reviewed_col].fillna("").astype(str)

    return out

def fetch_uniprot_go(uniprot_accessions, fields, chunk_size=200, out_tsv="uniprot_to_go.tsv"):
    rows = []
    for chunk in chunked(uniprot_accessions, chunk_size):
        q = "(" + " OR ".join([f"accession:{a}" for a in chunk]) + ")"
        url = (
            f"{API}/uniprotkb/stream"
            f"?format=tsv"
            f"&fields={quote(fields, safe=',')}"
            f"&query={quote(q)}"
        )
        tmp = io.StringIO()
        with urlopen(url, timeout=300) as r:
            tmp.write(r.read().decode("utf-8"))
        tmp.seek(0)
        df = pd.read_csv(tmp, sep="\t")
        rows.append(normalize_uniprot_go_df(df))

    go_df = pd.concat(rows, ignore_index=True).drop_duplicates(subset=["uniprot_accession"])
    go_df.to_csv(out_tsv, sep="\t", index=False)
    return go_df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", default="protein_family_cds.tsv")
    ap.add_argument("--out", dest="outfile", default="protein_family_cds_uniprot_go.tsv")
    ap.add_argument("--from-db", default="RefSeq_Protein")
    ap.add_argument("--to-db", default="UniProtKB")
    ap.add_argument("--poll-seconds", type=int, default=3)
    ap.add_argument("--map-chunk", type=int, default=50000, help="IDs per idmapping job (<=100000)")
    ap.add_argument("--go-chunk", type=int, default=200, help="UniProt accessions per GO query")
    ap.add_argument("--collapse-by-protein", action="store_true",
                    help="Collapse multiple UniProt mappings per RefSeq protein into one row (semicolon-joined).")
    args = ap.parse_args()

    inp = pd.read_csv(args.infile, sep="\t")

    # Ensure protein_accession exists
    if "protein_accession" not in inp.columns:
        if "protein_id" in inp.columns and inp["protein_id"].astype(str).str.contains("|", regex=False).any():
            inp["protein_accession"] = inp["protein_id"].astype(str).str.split("|", n=1, expand=True)[1]
        else:
            raise RuntimeError("Input must contain 'protein_accession' or 'protein_id' with genome|protein format.")

    refseq_ids = sorted(set(inp["protein_accession"].dropna().astype(str).tolist()))
    if not refseq_ids:
        raise RuntimeError("No protein_accession values found.")

    # 1) ID mapping: RefSeq -> UniProt
    map_parts = []
    map_out = Path("refseq_to_uniprot.tsv")
    for chunk in chunked(refseq_ids, args.map_chunk):
        job_id = submit_and_wait_idmapping(args.from_db, args.to_db, chunk,
                                           poll_seconds=args.poll_seconds, timeout_seconds=7200)
        part_path = Path("idmapping_parts") / f"{job_id}.tsv"
        download_to_file(f"{API}/idmapping/stream/{job_id}?format=tsv", part_path)
        df = pd.read_csv(part_path, sep="\t")
        # Expect columns like From / To
        if "From" not in df.columns or "To" not in df.columns:
            # try lowercase fallback
            cols = {c.lower(): c for c in df.columns}
            df = df.rename(columns={
                cols.get("from", "From"): "From",
                cols.get("to", "To"): "To",
            })
        map_parts.append(df[["From", "To"]])

    map_df = pd.concat(map_parts, ignore_index=True).dropna()
    map_df.to_csv(map_out, sep="\t", index=False)

    # 2) GO retrieval for UniProt accessions
    uniprot_ids = sorted(set(map_df["To"].astype(str).tolist()))
    if not uniprot_ids:
        raise RuntimeError("No UniProt mappings returned. Check --from-db token and input IDs.")

    # These field names are known to work with uniprotkb/stream TSV output.
    # You can trim if you want fewer columns.
    fields = "accession,reviewed,id,protein_name,gene_names,organism_name,go_id,go_p,go_c,go_f,go"
    go_df = fetch_uniprot_go(uniprot_ids, fields, chunk_size=args.go_chunk, out_tsv="uniprot_to_go.tsv")

    # 3) Join back: RefSeq protein -> UniProt -> GO
    map_go = map_df.merge(go_df, left_on="To", right_on="uniprot_accession", how="left")

    if args.collapse_by_protein:
        # collapse multiple UniProt hits per RefSeq protein
        def uniq_join(s):
            vals = [v for v in s.dropna().astype(str).tolist() if v and v != "nan"]
            return ";".join(sorted(set(vals)))

        map_go = (map_go.groupby("From", as_index=False)
                        .agg({
                            "To": uniq_join,
                            "go_ids": uniq_join,
                            "go_bp": uniq_join if "go_bp" in map_go.columns else "first",
                            "go_cc": uniq_join if "go_cc" in map_go.columns else "first",
                            "go_mf": uniq_join if "go_mf" in map_go.columns else "first",
                            "go_all": uniq_join if "go_all" in map_go.columns else "first",
                        }))
        map_go = map_go.rename(columns={"To": "uniprot_accessions"})
    else:
        map_go = map_go.rename(columns={"To": "uniprot_accession"})

    out = inp.merge(map_go, left_on="protein_accession", right_on="From", how="left")
    out = out.drop(columns=["From"], errors="ignore")

    out.to_csv(args.outfile, sep="\t", index=False)
    print(f"Wrote: {args.outfile}")
    print(f"Wrote: {map_out} (RefSeq -> UniProt)")
    print("Wrote: uniprot_to_go.tsv (UniProt -> GO)")

if __name__ == "__main__":
    try:
        main()
    except HTTPError as e:
        print("HTTPError:", e, file=sys.stderr)
        try:
            print(e.read().decode("utf-8", errors="replace"), file=sys.stderr)
        except Exception:
            pass
        sys.exit(2)
    except URLError as e:
        print("URLError:", e, file=sys.stderr)
        sys.exit(2)
