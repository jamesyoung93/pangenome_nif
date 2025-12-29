#!/usr/bin/env python3
import argparse, csv, io, gzip, re, sys
from pathlib import Path

def open_auto(p: Path):
    return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="ignore") if str(p).endswith(".gz") \
           else open(p, "r", encoding="utf-8", errors="ignore")

def guess_subunit(path: Path):
    n = path.name.lower(); d = path.parent.name.lower()
    for sub in ("nifh","nifd","nifk"):
        if sub in n or sub in d: return sub
    return None

def token_variants(t: str):
    out, seen = [], set()
    def emit(x): 
        if x and x not in seen: seen.add(x); out.append(x)
    emit(t)
    if "|" in t:
        parts = t.split("|"); emit(parts[-1]); emit(parts[0])
    if "." in t:
        base, _, tail = t.rpartition(".")
        if base and tail.isdigit(): emit(base)
    if "/" in t:
        emit(t.split("/")[-1])
    if t and t[0] in "'\"":
        emit(t.strip('"\'')) 
    return out

def parse_domtblouts(dom_root: Path):
    best = {}  # (token, 'nifh'|'nifd'|'nifk') -> (iE, hmmname)
    dom_files = list(Path(dom_root).rglob("*.domtblout"))
    if not dom_files:
        print(f"[warn] no *.domtblout under {dom_root}", file=sys.stderr)
        return best
    rows = 0
    for p in dom_files:
        sub = guess_subunit(p) or "n/a"
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line or line[0] == '#': continue
                parts = line.rstrip("\n").split()
                if len(parts) < 13: continue
                token   = parts[0]
                hmmname = parts[3]            # query name
                try:
                    iE = float(parts[12])     # i-Evalue
                except Exception:
                    try: iE = float(parts[12].replace('e','E'))
                    except Exception: continue
                k = (token, sub)
                if k not in best or iE < best[k][0]:
                    best[k] = (iE, hmmname)
                rows += 1
    print(f"[diag] domtblout rows: {rows}, unique tokens: {len({t for t,_ in best})}", file=sys.stderr)
    return best

def build_token_to_assemblies(downloads_root: Path, tokens):
    asm_pat = re.compile(r'(GC[AF]_\d+\.\d+)')
    token2asms = {}
    targets = set(tokens)
    for faa in Path(downloads_root).rglob("*_protein.faa*"):
        m = asm_pat.search(faa.name) or asm_pat.search(faa.parent.name)
        if not m: continue
        acc = m.group(1)
        with open_auto(faa) as f:
            for ln in f:
                if not ln or ln[0] != '>': continue
                tok = ln[1:].strip().split()[0]
                if tok in targets:
                    token2asms.setdefault(tok, set()).add(acc)
    print(f"[diag] tokens mapped to >=1 assemblies: {sum(1 for t in targets if t in token2asms)} / {len(targets)}", file=sys.stderr)
    return token2asms

def load_assemblies_order(assemblies_tsv: Path):
    order=[]
    with open(assemblies_tsv, "r", encoding="utf-8") as f:
        header=f.readline().rstrip("\n").split("\t")
        idx=header.index("assembly_accession") if "assembly_accession" in header else 0
        for ln in f:
            parts=ln.rstrip("\n").split("\t")
            if len(parts)>idx: order.append(parts[idx])
    seen=set(); out=[]
    for a in order:
        if a not in seen: seen.add(a); out.append(a)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domdir", required=True)
    ap.add_argument("--downloads", required=True)
    ap.add_argument("--assemblies", required=True)
    ap.add_argument("--combined", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--thr", type=float, default=1e-5)
    a=ap.parse_args()

    best = parse_domtblouts(Path(a.domdir))
    tokens = {t for (t,_) in best.keys()}
    token2asms = build_token_to_assemblies(Path(a.downloads), tokens)
    assemblies = load_assemblies_order(Path(a.assemblies))

    subs = ("nifH","nifD","nifK")
    present_tot = {s:0 for s in subs}
    with open(a.out, "w", newline="", encoding="utf-8") as o:
        w=csv.writer(o)
        w.writerow(["assembly_accession","nifH_present","nifH_best_evalue","nifH_best_ref",
                    "nifD_present","nifD_best_evalue","nifD_best_ref",
                    "nifK_present","nifK_best_evalue","nifK_best_ref","best_hit_ids"])
        for acc in assemblies:
            best_e = {s:(None,None,None) for s in subs}  # (iE, ref, token)
            # collect tokens tied to this accession
            tokset=set(t for t in tokens if t in token2asms and acc in token2asms[t])
            for tok in tokset:
                for sub in subs:
                    k=(tok, sub.lower())
                    if k in best:
                        iE, ref = best[k]
                        cur = best_e[sub][0]
                        if cur is None or iE < cur:
                            best_e[sub] = (iE, ref, tok)
            row=[acc]; ids=[]
            for sub in subs:
                iE, ref, tok = best_e[sub]
                present = 1 if (iE is not None and iE <= a.thr) else 0
                if present: present_tot[sub]+=1
                row += [present, f"{iE:.3e}" if iE is not None else "", ref or ""]
                ids.append(f"{sub[3].upper()}={tok or ''}")
            row.append(";".join(ids))
            w.writerow(row)
    print(f"Wrote -> {a.out}  |  H={present_tot['nifH']} D={present_tot['nifD']} K={present_tot['nifK']}", file=sys.stderr)

if __name__ == "__main__":
    main()
