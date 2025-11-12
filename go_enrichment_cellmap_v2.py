#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# v2: wraps labels, nudges to reduce overlap, splits positive/negative panels,
# uses go-basic.obo names if present, keeps matplotlib-only and a single plot per fig.

import argparse, os, re, math, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_go_obo(path="go-basic.obo"):
    if not os.path.exists(path): return {}, {}
    name, ns, cur = {}, {}, None
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("[Term]"): cur = {"id":None, "name":None, "ns":None}; continue
            if line.startswith("id: GO:"): cur["id"]=line.split("id: ")[1].strip()
            elif line.startswith("name: "): cur["name"]=line.split("name: ",1)[1].strip()
            elif line.startswith("namespace: "): cur["ns"]=line.split("namespace: ",1)[1].strip()
            elif not line.strip() and cur and cur["id"]:
                name[cur["id"]]=cur["name"] or cur["id"]; ns[cur["id"]]=cur["ns"] or ""
                cur=None
    if cur and cur.get("id"): name[cur["id"]]=cur.get("name",""); ns[cur["id"]]=cur.get("ns","")
    return name, ns

def parse_go(s):
    if pd.isna(s) or not str(s).strip(): return set()
    return {t for t in re.split(r"[;, \t]+", str(s)) if t.startswith("GO:")}

def bh_fdr(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=np.empty(n,int); r[o]=np.arange(1,n+1)
    q=p*n/np.maximum(r,1); qs=np.minimum.accumulate(q[o][::-1])[::-1]; out=np.empty(n,float); out[o]=np.clip(qs,0,1)
    return out

def fisher_right_tail(a,b,c,d):
    from math import comb
    N=a+b+c+d; K=a+c; n=a+b
    def pmf(x): return comb(K,x)*comb(N-K,n-x)/comb(N,n)
    return min(1.0,sum(pmf(x) for x in range(a, min(K,n)+1)))

def load_table(path):
    df=pd.read_csv(path)
    if "direction" not in df.columns and "effect_size" in df.columns:
        df["direction"]=np.where(df["effect_size"]>0,"Positive",np.where(df["effect_size"]<0,"Negative","Neutral"))
    return df

def gf_go(df):
    gg=df[["gene_family","GO"]].copy()
    gg["GOset"]=gg["GO"].map(parse_go)
    return gg.groupby("gene_family",as_index=False).agg(GOset=("GOset",lambda s:set().union(*s)))

def enrich(df, gg, label="Positive", min_count=2, q_cut=0.1):
    P=set(df.loc[df["direction"]==label,"gene_family"])
    U=gg[gg["GOset"].map(len)>0]; allgf=set(U["gene_family"]); P=P & allgf
    go2gf={}
    for r in U.itertuples():
        for go in r.GOset: go2gf.setdefault(go,set()).add(r.gene_family)
    rows=[]; m=len(P); M=len(allgf)
    for go, members in go2gf.items():
        a=len(members & P); 
        if a<min_count: continue
        b=m-a; c=len(members)-a; d=(M-m)-c
        p=fisher_right_tail(a,b,c,d); rows.append((go,a,len(members),p))
    if not rows: return pd.DataFrame(columns=["GO","pos_hits","term_size","p","q","direction"])
    out=pd.DataFrame(rows,columns=["GO","pos_hits","term_size","p"])
    out["q"]=bh_fdr(out["p"]); out["direction"]=label
    return out[(out["q"]<=q_cut)| (out["pos_hits"]>=3)].sort_values(["q","p","pos_hits"]).reset_index(drop=True)

COMPARTMENT_RULES=[
    ("membrane",  ["membrane","transmembrane","secretion","transport","porin","flagell"]),
    ("dna",       ["dna","replication","recombination","repair","transcription"]),
    ("ribosome",  ["ribosom","translation","rrna","trna","elongation"]),
    ("redox",     ["oxidoreduct","electron","respirat","dehydrogenase","cytochrome","ferredoxin","flavodoxin","quinone"]),
    ("metabolism",["metabolic","biosynth","catabolic","cofactor","sulfur","iron","molybden","heme"]),
    ("stress",    ["stress","oxygen","reactive","peroxid","superoxide","catalase"]),
    ("other",     [""])
]

def compartment(name):
    s=(name or "").lower()
    for comp,keys in COMPARTMENT_RULES:
        if any(k in s for k in keys): return comp
    return "other"

def wrap(txt, width=28):
    w=[]; line=""
    for t in txt.split():
        if len(line)+1+len(t) > width:
            w.append(line); line=t
        else:
            line=(t if not line else line+" "+t)
    if line: w.append(line)
    return "\n".join(w)

def draw_cell(ax):
    ax.set_aspect('equal'); ax.axis('off')
    outer=plt.Circle((0,0), 4.5, fill=False, linewidth=2); ax.add_artist(outer)
    inner=plt.Circle((0,0), 3.8, fill=False, linewidth=1, linestyle="--"); ax.add_artist(inner)
    nuc=plt.Circle((0,0), 2.0, fill=False, linewidth=1); ax.add_artist(nuc)

def place_terms(ax, df, title):
    rng=np.random.default_rng(42)
    draw_cell(ax)
    boxes={
        "membrane":  (-4.5, -4.5, 9.0, 1.4),
        "redox":     (-3.5,  1.2,  7.0, 1.4),
        "metabolism":(-3.2, -1.0,  6.4, 1.6),
        "dna":       (-2.0, -0.6,  4.0, 1.2),
        "ribosome":  (-3.4, -2.8,  6.8, 1.4),
        "stress":    (-3.8,  2.8,  7.6, 1.2),
        "other":     (-1.5, -4.0,  3.0, 1.0),
    }
    coords=[]
    for comp, box in boxes.items():
        sub=df[df["compartment"]==comp]
        if not len(sub): continue
        x0,y0,w,h=box
        xs=rng.uniform(x0, x0+w, size=len(sub))
        ys=rng.uniform(y0, y0+h, size=len(sub))
        coords.append((sub, xs, ys))
    # simple de-overlap nudging
    def nudge(xy):
        X=xy.copy()
        for _ in range(200):
            moved=False
            for i in range(len(X)):
                for j in range(i+1,len(X)):
                    dx=X[i][0]-X[j][0]; dy=X[i][1]-X[j][1]
                    if dx*dx+dy*dy < 0.05:  # close
                        X[i][0]+=0.04; X[i][1]+=0.02; X[j][0]-=0.04; X[j][1]-=0.02; moved=True
            if not moved: break
        return X
    for sub, xs, ys in coords:
        pts=nudge(np.vstack([xs,ys]).T)
        sizes = (-np.log10(np.clip(sub["q"].values, 1e-300, 1.0)) * 150).clip(30, 800)
        ax.scatter(pts[:,0], pts[:,1], s=sizes, marker=("o" if sub["direction"].iloc[0]=="Positive" else "s"), alpha=0.7)
        for (x,y), name in zip(pts, sub["label"].values):
            ax.text(x, y, wrap(name, 24), fontsize=8, ha="center", va="center")
    ax.set_title(title, fontsize=11, fontweight="bold")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--use", default="llm_context_features.csv")
    ap.add_argument("--q", type=float, default=0.10)
    ap.add_argument("--min-count", type=int, default=2)
    args=ap.parse_args()

    df=load_table(args.use)
    gg=gf_go(df)
    pos=enrich(df, gg, "Positive", args.min_count, args.q)
    neg=enrich(df, gg, "Negative", args.min_count, args.q)

    go_name, go_ns = read_go_obo("go-basic.obo")
    for d in (pos,neg):
        d["name"]=[go_name.get(g,g) for g in d["GO"]]
        d["namespace"]=[go_ns.get(g,"") for g in d["GO"]]
        d["compartment"]=[compartment(n) for n in d["name"]]
        d["label"]=d["name"]  # already wrapped downstream

    pos.to_csv("enriched_go_positive.csv", index=False)
    neg.to_csv("enriched_go_negative.csv", index=False)

    # separate panels
    if len(pos):
        fig,ax=plt.subplots(figsize=(10,7)); place_terms(ax, pos, "Positive-enriched terms")
        plt.tight_layout(); fig.savefig("go_cell_map_positive.svg", bbox_inches="tight"); plt.close(fig)
    if len(neg):
        fig,ax=plt.subplots(figsize=(10,7)); place_terms(ax, neg, "Negative-enriched terms")
        plt.tight_layout(); fig.savefig("go_cell_map_negative.svg", bbox_inches="tight"); plt.close(fig)

    # overview (side-by-side if both exist)
    if len(pos) and len(neg):
        fig,axs=plt.subplots(1,2,figsize=(14,6))
        place_terms(axs[0], pos, "Positive"); place_terms(axs[1], neg, "Negative")
        plt.tight_layout(); fig.savefig("go_cell_map_overview.svg", bbox_inches="tight"); plt.close(fig)

if __name__=="__main__":
    main()
