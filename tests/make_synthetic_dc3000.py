"""SYNTHETIC test data (not biological): negative-binomial counts on the DC3000 NCBI annotation (chromosome + pDC3000A + pDC3000B),
3 MC + 3 AER replicates. Planted effects: pDC3000A genes +2 log2 in AER; chromosome 1.00–1.10 Mb −1.5 log2. Used only to test layout and Mode B."""
import re, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261001)
rows = []
for l in open("tests/fixtures/DC3000_NCBI.gff"):
    if l.startswith("##FASTA"): break
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS": continue
    m = re.search(r"Parent=gene-([^;]+)", f[8]); lt = re.search(r"locus_tag=([^;]+)", f[8])
    rows.append((f[0], int(f[3]), int(f[4]), (lt.group(1) if lt else m.group(1))))
g = pd.DataFrame(rows, columns=["rep", "s", "e", "id"]).drop_duplicates("id")
lt_gene = {}
for l in open("tests/fixtures/DC3000_NCBI.gff"):
    f = l.split("\t")
    if len(f) > 8 and f[2] == "gene":
        a = dict(x.split("=", 1) for x in f[8].strip().split(";") if "=" in x)
        if "locus_tag" in a: lt_gene["gene-" + a["locus_tag"]] = a["locus_tag"]
mu = np.exp(rng.normal(5.5, 1.5, len(g)))
lfc = np.zeros(len(g)); lfc[g.rep.eq("NC_004633.1").values] = 2.0; lfc[(g.rep.eq("NC_004578.1") & g.s.between(1_000_000, 1_100_000)).values] = -1.5
sf = {"MC_1": 0.9, "MC_2": 1.1, "MC_3": 1.0, "AER_1": 0.8, "AER_2": 1.2, "AER_3": 1.05}
disp = 0.05
C = {}
for s, k in sf.items():
    m = mu * k * (2 ** lfc if s.startswith("AER") else 1)
    p = 1 / (1 + m * disp); C[s] = rng.negative_binomial(1 / disp, p)
M = pd.DataFrame(C, index=g.id); M.index.name = "gene_id"
M.to_csv("tests/synthetic_DC3000/SYNTHETIC_counts_matrix.tsv", sep="\t")
pd.DataFrame([dict(sample_id=s, label=s.replace("_", " rep "), condition=s.split("_")[0], replicate=s.split("_")[1], batch="SYN") for s in sf]).to_csv("tests/synthetic_DC3000/samplesheet.tsv", sep="\t", index=False)
pd.DataFrame([dict(contrast_id="AER_vs_MC", label="AER vs MC (synthetic, 3 vs 3)", short_label="AER vs MC", treatment="condition=AER", control="condition=MC", role="test", de_file="SYNTHETIC_DESeq2_AER_vs_MC.csv"),
              dict(contrast_id="MC1_vs_MC2", label="MC rep 1 vs MC rep 2 (noise)", short_label="noise MC1/MC2", treatment="MC_1", control="MC_2", role="reference", de_file="")]).to_csv("tests/synthetic_DC3000/contrasts.tsv", sep="\t", index=False)
pd.DataFrame([dict(name="planted −1.5 region", replicon="NC_004578.1", start=1000000, end=1100000, type="synthetic effect", color="", label="planted −1.5"),
              dict(name="pDC3000A (planted +2)", replicon="NC_004633.1", start=1, end=73661, type="synthetic effect", color="", label="planted +2")]).to_csv("tests/synthetic_DC3000/landmarks.tsv", sep="\t", index=False)
# Mode B: pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
inf = DefaultInference(n_cpus=1)  # sandbox: no multiprocessing
meta = pd.DataFrame({"condition": [s.split("_")[0] for s in sf]}, index=list(sf))
dds = DeseqDataSet(counts=M.T, metadata=meta, design="~condition", quiet=True, inference=inf); dds.deseq2()
res = DeseqStats(dds, contrast=["condition", "AER", "MC"], quiet=True, inference=inf); res.summary()
R = res.results_df; R.to_csv("tests/synthetic_DC3000/SYNTHETIC_DESeq2_AER_vs_MC.csv")
pd.DataFrame([dict(contrast_id="AER_vs_MC", label="AER vs MC (synthetic, DESeq2 Wald)", short_label="AER vs MC", treatment="", control="", role="test", de_file="SYNTHETIC_DESeq2_AER_vs_MC.csv")]).to_csv("tests/synthetic_DC3000/contrasts_modeB.tsv", sep="\t", index=False)
pl = g.rep.eq("NC_004633.1").values; reg = (g.rep.eq("NC_004578.1") & g.s.between(1_000_000, 1_100_000)).values
print("genes", len(g), "plasmid A genes", pl.sum(), "region genes", reg.sum())
print("DESeq2 mean LFC plasmid A", round(R.log2FoldChange.values[pl].mean(), 3), "region", round(R.log2FoldChange.values[reg].mean(), 3), "padj<0.05:", int((R.padj < 0.05).sum()))
