# genomemap — input files

The viewer (`genome_map_viewer.html`) and the command-line builder (`genomemap_build.py`) accept the same files. The format is detected automatically from the file content, not from its name.

| # | File | Required | Purpose |
|---|---|---|---|
| 1 | Genome annotation: GFF3 or GenBank | **yes** | replicons, gene coordinates, names, products, COG |
| 2a | Mode A: count table (featureCounts, generic matrix, or HTSeq-count files) | one of 2a/2b | raw read counts per gene and sample |
| 2b | Mode B: DE result tables (DESeq2 / edgeR / limma), one per contrast | one of 2a/2b | log2FC and adjusted p-values computed elsewhere |
| 3 | `samplesheet.tsv` | Mode A: yes | describes each library |
| 4 | `contrasts.tsv` | **yes** | which pairs or groups to compare (several allowed) |
| 5 | `landmarks.tsv` | no | named regions drawn around the circle |
| 6 | `labels.tsv` | no | display names for selected loci (e.g. *ina*) |

All tables are tab- or comma-separated text, UTF-8, with a header row. Lines starting with `#` are ignored.

---

## 1. Genome annotation

**Accepted inputs.**
- **GFF3** from Bakta, Prokka, PGAP/NCBI (`*_genomic.gff`) or RAST. A `##FASTA` section, if present, is skipped.
- **GenBank** flat file (`.gbk`, `.gbff`, `.gb`) with one `LOCUS` entry per replicon.

**What is read.**

| Element | GFF3 | GenBank |
|---|---|---|
| Replicon id and length | `##sequence-region` lines, or `region` features | `LOCUS` line |
| Topology | `Is_circular=true` on `region` | `circular` / `linear` on `LOCUS` |
| Gene-level features: CDS, tRNA, rRNA, ncRNA, tmRNA (optionally others) | yes | yes |
| Feature id: `locus_tag`, else `ID` (Prokka/Bakta), else `Name`; for PGAP the `locus_tag` of the parent `gene` | yes | `/locus_tag` |
| `gene`, `product` | yes | yes |
| COG id / COG category (`Dbxref=COG:COG2197,COG:TK` in Bakta; `/note`/`/db_xref` COG where present) | yes | yes |

**Rule:** the feature ids must be the same ids used in the count table. If you counted with featureCounts on a GTF/GFF made from the same annotation, this is automatic. The viewer reports the match rate and stops if it is below 90%.

## 2a. Mode A — counts (computed in the viewer)

1. **featureCounts output** (recommended). Use the file as written by featureCounts, including the `# Program:` comment line. Columns are `Geneid Chr Start End Strand Length <sample1> <sample2> …`. Sample columns that are BAM paths are shortened to the file name without `.bam`.
2. **Generic matrix**: first column `gene_id` (or `Geneid`, `locus_tag`), then one column per sample. Only integer counts are accepted; TPM or FPKM are rejected.
3. **HTSeq-count**: one two-column file per sample (`gene_id<TAB>count`). Summary rows starting with `__` are removed. The file name without its extensions (`.gz`, `.txt`, `.tsv`, `.csv`, `.count(s)`, `.htseq`; e.g. `MC_1.htseq.txt` → `MC_1`) is the sample id, unless `count_column` in the sample sheet says otherwise. Load all per-sample files together.

Producing counts, e.g. featureCounts for paired-end, reverse-stranded libraries (adapt `-s` to your kit):
```
featureCounts -a annotation.gtf -F GTF -t exon -g gene_id -p --countReadPairs -s 2 -T 8 -o counts.txt *.bam
```
or directly on the Bakta GFF3:
```
featureCounts -a annotation.gff3 -F GFF -t CDS,tRNA,rRNA,ncRNA,tmRNA -g locus_tag -p --countReadPairs -s 2 -o counts.txt *.bam
```

**Computed in the viewer** (documented in `data_model.md` §4):
- median-of-ratios size factors on CDS rows, `standard` (DESeq2) or `zero_tolerant` variant;
- normalised counts;
- expression filter: CPM > 1 in ≥ *n* libraries, default *n* = 2;
- log2((mean treatment + 1) / (mean control + 1)).

**No p-values are computed in Mode A.**

## 2b. Mode B — DE result tables (imported)

One table per contrast, named in the `de_file` column of `contrasts.tsv`. Columns are recognised by name:

| Quantity | DESeq2 `results()` | edgeR `topTags()` | limma `topTable()` |
|---|---|---|---|
| gene id | first column / row names, or `gene_id`, `Geneid`, `locus_tag` | same | same |
| log2 fold change | `log2FoldChange` | `logFC` | `logFC` |
| adjusted p | `padj` | `FDR` | `adj.P.Val` |
| abundance (optional; shown in the tooltip and CSV export under its own column name) | `baseMean` | `logCPM` | `AveExpr` |

Write the tables with e.g. `write.csv(as.data.frame(res), "AER_vs_MC.csv")` in R. If you used `lfcShrink`, the shrunken values are shown; say so in the contrast label. Genes missing from a table are shown as "not tested".

## 3. `samplesheet.tsv` (Mode A)

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | must match a count column (after BAM-path shortening) or `count_column` |
| `label` | no | short display name (default: `sample_id`) |
| `condition` | yes | e.g. `MC`, `AER` |
| `replicate` | no | 1, 2, 3 … |
| `batch` | no | experiment / run, e.g. `EXP3` |
| `count_column` | no | count-table column name, if it differs from `sample_id` |

## 4. `contrasts.tsv`

One row per comparison. You can list as many as you like; the viewer switches between them, and up to 6 can be shown as rings at the same time.

| column | required | meaning |
|---|---|---|
| `contrast_id` | yes | unique id, no spaces |
| `label` | yes | full display label, e.g. `AER vs MC · EXP3` |
| `short_label` | no | ≤ 16 characters, used for linear-view track names |
| `treatment` | Mode A | sample ids separated by `;`, **or** `condition=AER` (optionally `condition=AER&batch=EXP3`) |
| `control` | Mode A | as above, e.g. `EXP3_MC_19626` or `condition=MC&batch=EXP3` |
| `role` | no | `test` (default) or `reference`. A reference contrast is a noise baseline (e.g. MC vs MC from two experiments), drawn on the inner ring and as an optional linear track. |
| `de_file` | Mode B | file name of the DE table for this contrast |

**Pairs.** For "control + sample" pairs, write one row per pair, e.g. `treatment = EXP3_AER_19626`, `control = EXP3_MC_19626`. With one library per side the value is a descriptive ratio, and the viewer shows a "no replication" banner.

## 5. `landmarks.tsv` (optional)

| column | required | meaning |
|---|---|---|
| `name` | yes | tooltip / list name |
| `replicon` | yes | replicon id as in the annotation |
| `start`, `end` | yes | 1-based, inclusive |
| `type` | no | free text (gene, cluster, island…) |
| `color` | no | hex colour, default by type |
| `label` | no | short text printed on the map (default: `name`) |

A landmark may also be given as a single `feature_id` instead of `replicon` / `start` / `end`; it then takes that feature's coordinates. Keep the `name` column in either case: the viewer recognises a landmark table by `name` plus `start`/`end` or `feature_id`.

## 6. `labels.tsv` (optional)

`feature_id`, `label`, `note`. Use this to show a curated name for a locus whose automatic annotation is uninformative (e.g. MGGFHP_01089 → *ina*).

---

### Worked example
`example_R10.79/` contains the full R10.79 input set: Bakta GFF3, featureCounts output, `samplesheet.tsv`, `contrasts.tsv` (four test pairs vs their own MC plus the MC EXP4 / MC EXP3 noise reference), `landmarks.tsv` and `labels.tsv`. Empty templates are in `templates/`.
