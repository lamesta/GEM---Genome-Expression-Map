"""Cross-check: 'standard' size factors from genomemap core vs pydeseq2 (DESeq2 port) on the R10.79 CDS counts."""
import json, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet
import pydeseq2
fc = pd.read_csv("example_R10.79/featureCounts_counts_all.txt", sep="\t", comment="#")
fc.columns = [c.split("/")[-1].replace(".bam", "") for c in fc.columns]
O = json.load(open("tests/out/out_r1079_standard.json")); core = O["ds"]["meta"]["size_factors"]
types = dict(zip(O["ds"]["features"]["id"], O["ds"]["features"]["type"]))
# CDS rows per annotation (core uses annotation type CDS for size factors)
import re
cds_ids = set()
for l in open("example_R10.79/R10.79_bakta.gff3"):
    if l.startswith("##FASTA"): break
    f = l.split("\t")
    if len(f) > 8 and f[2] == "CDS": cds_ids.add(re.search(r"locus_tag=([^;]+)", f[8]).group(1))
libs = list(core.keys())
M = fc[fc.Geneid.isin(cds_ids)].set_index("Geneid")[libs].T
meta = pd.DataFrame({"condition": ["A"] * len(libs)}, index=libs)
dds = DeseqDataSet(counts=M, metadata=meta, design="~1", quiet=True)
dds.fit_size_factors()
py = dict(zip(libs, dds.obs["size_factors"].values if "size_factors" in dds.obs else dds.obsm["size_factors"]))
rows = [dict(library=l, genomemap_standard=core[l], pydeseq2=float(py[l]), rel_diff=abs(core[l] / float(py[l]) - 1)) for l in libs]
df = pd.DataFrame(rows); df.to_csv("tests/crosscheck_size_factors_pydeseq2.tsv", sep="\t", index=False)
print("pydeseq2", pydeseq2.__version__); print(df.to_string()); print("MAX REL DIFF", df.rel_diff.max())
