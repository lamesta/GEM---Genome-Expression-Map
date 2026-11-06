"""Regression: genomemap core (Node) vs data embedded in the original R10.79 viewer."""
import json, sys, numpy as np, pandas as pd
orig_html = sys.argv[1]; variant = sys.argv[2] if len(sys.argv) > 2 else "zero_tolerant"
s = open(orig_html, encoding="utf-8").read(); i0 = s.index("const DATA = ") + 13; i1 = s.index(";\n", i0)
D = json.loads(s[i0:i1].replace("<\\/", "</"))
O = json.load(open(f"tests/out/out_r1079_{variant}.json")); ds = O["ds"]; B = O["bins"]["reps"][0]
rows = []
def add(item, value, tol, ok, note=""): rows.append(dict(check=item, value=value, tolerance=tol, pass_=bool(ok), note=note))
f = ds["features"]; g = D["genes"]
add("feature ids & order identical", f["id"] == g["id"], "exact", f["id"] == g["id"])
for k, kk in [("s", "s"), ("e", "e"), ("st", "st"), ("type", "ft"), ("gene", "g"), ("cog", "cog"), ("cogcat", "cc"), ("group", "sg")]:
    same = sum(a == b for a, b in zip(f[k], g[kk])); add(f"feature field {k}", f"{same}/{len(g[kk])} identical", "exact", same == len(g[kk]))
from urllib.parse import unquote
prod_same = sum(a == b for a, b in zip(f["product"], g["p"])); prod_dec = sum(a == unquote(b) for a, b in zip(f["product"], g["p"]))
add("feature product (after URL-decoding original)", f"{prod_dec}/{len(g['p'])} identical; {len(g['p'])-prod_same} originals showed undecoded %2C", "exact", prod_dec == len(g["p"]), "original viewer displayed GFF escapes such as %2C; new parser decodes them")
add("expression filter", f"{sum(a == b for a, b in zip(ds['expressed'], g['x']))}/{len(g['x'])} identical", "exact", ds["expressed"] == g["x"])
libs = [l["id"] for l in ds["libraries"]]; order = [libs.index(k) for k in D["libkeys"]]
N = np.array([ds["normalized"][j] for j in order], float); G0 = np.array(g["n"], float); dN = np.abs(N - G0); okN = bool((dN <= 0.051 + 2e-6 * np.abs(G0)).all())
add("normalised counts vs original (rounded to 0.1)", f"max |Δ| {np.nanmax(dN):.4f}", "0.05 rounding + 0.001 + 2e-6·x", okN, "residual: implied size factors differ from the original by ≤1.6e-6 (relative)")
cids = [c["id"] for c in ds["contrasts"]]
for c, ck in enumerate(D["ck"]):
    a = np.array(ds["log2fc"][cids.index(ck)], float); b = np.array(g["fc"][c], float); d = np.nanmax(np.abs(a - b))
    add(f"log2FC {ck} max |Δ| (orig rounded to 1e-4)", float(d), 5e-5 + 1e-9, d <= 5e-5 + 1e-9)
DB = D["bins"]; add("bin count", B["n"], "exact", B["n"] == DB["n"])
chg = set(np.where((np.array(B["ng"]) != np.array(DB["ng"])) | (np.array(B["nx"]) != np.array(DB["nx"])))[0].tolist())
add("bins with different membership", sorted(chg), "≤2 boundary bins", len(chg) <= 2, "new convention: position p belongs to bin floor((p-1)/bw), matching how bins are drawn; original used floor(mid/bw)")
for c, ck in enumerate(D["ck"]):
    ci = cids.index(ck)
    for key in ["med", "fu", "fd"]:
        a = np.array([np.nan if v is None else v for v in B[key][ci]]); b = np.array([np.nan if v is None else v for v in DB[key][c]])
        keep = np.array([b_ not in chg for b_ in range(len(a))])
        same_null = bool((np.isnan(a[keep]) == np.isnan(b[keep])).all()); d = float(np.nanmax(np.abs(a[keep] - b[keep])))
        add(f"bins {key} {ck} (excluding boundary bins)", f"max |Δ| {d:.2e}; null pattern identical {same_null}", 1e-4, same_null and d <= 1e-4)
lm_ok = [(L["s"], L["e"]) for L in ds["landmarks"]] == [(L["s"], L["e"]) for L in D["lm"]]
add("landmarks coordinates", lm_ok, "exact", lm_ok)
df = pd.DataFrame(rows); df.to_csv(f"tests/regression_R10.79_{variant}.tsv", sep="\t", index=False)
print(df.to_string()); print("ALL PASS" if df.pass_.all() else "FAILURES: %d" % (~df.pass_).sum())
