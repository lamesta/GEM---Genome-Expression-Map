# genomemap — data model (dataset format v1)

## 1. Audit of the current R10.79 viewer (artifact `R10.79_genome_map_interactive (1).html`)

**Structure.** One self-contained HTML file:
- embedded fonts (254 KB);
- CSS (3.7 KB);
- body markup (4.3 KB);
- D3 v7 (inline);
- `const DATA = {…}` (about 1 MB);
- the app (37 KB, one IIFE).

**What `DATA` contains.**

| key | content | generic? |
|---|---|---|
| `G` | genome length (single replicon) | no: one replicon only |
| `libs`, `libkeys` | 7 library labels / ids | yes |
| `ck`, `cl` | 5 contrast ids / labels; **the last one is hard-wired as the noise contrast** | partly |
| `genes.{id,ft,s,e,st,g,p,cog,cc,sg,x}` | locus tag, feature type, start, end, strand, gene symbol, product, COG id, COG letters, COG super-group, expression-filter flag | yes |
| `genes.n[lib][gene]` | DESeq2-normalised counts | yes |
| `genes.fc[contrast][gene]` | log2((A+1)/(B+1)) | yes |
| `bins.{bw,n,ng,nx,med,fu,fd}` | 10-kb bins: n genes, n passing filter, median log2FC, fraction >+1 / <−1 | derivable in the viewer |
| `lm` | 7 landmarks with types `ina`/`t3ss`/`rrna` (colours hard-coded) | no |
| `meta.{norm,filter,operons}` | method text | yes |

**R10.79-specific items hard-coded in the app.**
- Title, circle centre text ("R10.79", "5.87 Mb") and the notes paragraph.
- `SHORT` contrast labels.
- `OVR = {MGGFHP_01089: "ina"}`.
- Landmark types and colours.
- The initial window centred on the INA gene, and the self-test gene list.

Everything else (circle + linear rendering, zoom/pan, tooltips, search, CSV export, gene-level / bin-median / bin-fraction modes, expression-filter display) is generic and is kept.

## 2. Concepts

- **Replicon:** a chromosome, plasmid or contig, with a name, a length and a topology (circular or linear).
- **Feature:** a gene-level row (CDS, tRNA, rRNA, ncRNA, tmRNA, …) on one replicon.
- **Library:** one sequenced sample, described by sample_id, condition, replicate and batch.
- **Contrast:** treatment group vs control group (e.g. AER vs MC), each group being one or more libraries.
  - `kind = descriptive`: values computed in the viewer from counts (Mode A); no p-values.
  - `kind = statistical`: imported from DESeq2/edgeR (Mode B); carries padj.
  - `role = test` is a normal contrast. `role = reference` is a noise contrast, e.g. control vs control from another experiment; it is shown on the inner ring and as an optional linear track.
- **Landmark:** a named region (gene cluster, island, gene of interest) drawn outside the ruler.
- **Label override:** a display label for a locus (e.g. MGGFHP_01089 → *ina*), with an optional note.

## 3. Dataset JSON (what the viewer renders)

Arrays are column-oriented (one array per attribute, indexed by feature) to keep files small.

```json
{
  "format": "genomemap-dataset", "version": 1,
  "title": "Pseudomonas syringae R10.79 — transcriptional responses",
  "organism": "Pseudomonas syringae", "strain": "R10.79",
  "replicons": [{"id": "contig_1", "name": "chromosome", "length": 5870538, "topology": "circular"}],
  "features": {
    "id": ["MGGFHP_00001", "..."], "rep": [0, "..."], "s": [2, "..."], "e": [10078, "..."], "st": [-1, "..."],
    "type": ["CDS", "..."], "gene": ["solG", "..."], "product": ["...", "..."],
    "cog": ["", "..."], "cogcat": ["", "..."], "group": ["none", "..."]
  },
  "libraries": [{"id": "EXP3_MC", "label": "EXP3 MC", "condition": "MC", "replicate": 1, "batch": "EXP3"}],
  "normalized": [[2264.1, "..."]],
  "expressed": [true, "..."],
  "contrasts": [{"id": "E3_AER_vs_MC", "label": "AER vs MC · EXP3", "short": "AER EXP3",
                 "treatment": ["EXP3_AER"], "control": ["EXP3_MC"],
                 "kind": "descriptive", "role": "test"}],
  "log2fc": [[0.1045, "..."]],
  "padj": [null],
  "basemean": [null],
  "landmarks": [{"name": "ina (MGGFHP_01089)", "rep": 0, "s": 1271643, "e": 1275671, "type": "gene", "color": "#9C4B22", "label": "ina"}],
  "labels": {"MGGFHP_01089": {"label": "ina", "note": "identity from project domain analysis"}},
  "meta": {"mode": "A", "normalization": "...", "filter": "...", "pseudocount": 1,
           "created": "ISO-8601", "software": "genomemap 1.0",
           "inputs": [{"file": "counts.txt", "sha256": "..."}], "warnings": []}
}
```

**Rules.**
- Features are sorted by (replicon order, start). Coordinates are 1-based and inclusive.
- `normalized[l][i]` is present only in Mode A; it may be `null` in Mode B.
- `log2fc[c][i]` is a number, or `null` if the gene is absent from that contrast's table.
- `padj[c]` is either `null` (descriptive contrast) or an array.
- `expressed[i]` is the expression-filter result (Mode A). In Mode B it is true when baseMean > 0 and padj is not NA in at least one contrast.
- Bin statistics (10-kb median, fraction |log2FC| > 1) are computed by the viewer from `log2fc` and `expressed`, per replicon. Genes are assigned to a bin by their midpoint.

The formal schema is in `genomemap_dataset_v1.schema.json`.

## 4. Computation in Mode A (counts → dataset)

The rules below were checked against the current R10.79 viewer; full regression is in step 6.

1. **Size factors (median-of-ratios; Anders & Huber 2010).** Computed on rows of type CDS only (configurable). Two variants:
   - **`standard`** (default; identical to DESeq2 `estimateSizeFactors`): use only genes with a count > 0 in every library, $s_j = \mathrm{median}_i\, \big(k_{ij} / (\prod_j k_{ij})^{1/m}\big)$.
   - **`zero_tolerant`** (the variant actually used in the current R10.79 viewer, audited 1 Oct 2026): zero counts are treated as missing. The geometric mean is taken over the non-zero libraries of each gene, and the median over the non-zero ratios of each library. On the R10.79 counts it reproduces the embedded normalised counts to rounding (max |Δlog2FC| = 5e-5).
   - The two variants differ by ≤ 0.09% in size factors on the R10.79 data. The text "DESeq2 median-of-ratios" in the current viewer is therefore approximately, but not exactly, correct; the new viewer states which variant was used.
2. **Normalised counts:** $k_{ij} / s_j$ for all feature types.
3. **Expression filter:** CPM > 1 in at least *n* libraries (default *n* = 2). CPM is computed against the CDS library totals.
4. **Contrast value:** $\log_2\big((\bar A + 1)/(\bar B + 1)\big)$, where $\bar A$ and $\bar B$ are the means of normalised counts over the treatment and control libraries. With one library per group this is the single-library ratio used now.
5. **No test statistic is computed.** With replicates, users who need inference should run DESeq2/edgeR and load the result tables (Mode B). The viewer then shows padj.
