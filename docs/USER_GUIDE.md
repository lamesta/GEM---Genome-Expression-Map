# genomemap user guide

genomemap draws bacterial RNA-seq results on a circular and a linear genome map. You give it a genome annotation and either raw counts or finished differential-expression tables. You choose which treatment-vs-control pairs to compare, switch between them, and save the result as a single HTML report that anyone can open.

It is one HTML file: `dist/genome_map_viewer.html`. There is nothing to install. It runs in a current Chrome, Edge, Firefox or Safari, and works offline. Your data are read inside the browser tab and are not uploaded anywhere.

Contents:

1. Quick start
2. Which files do I need?
3. Mode A and Mode B: what the numbers mean
4. Preparing the input files
5. Describing samples and pairs
6. Reading the map
7. Saving and sharing
8. Command-line builder (optional)
9. Limitations and interpretation
10. Troubleshooting
11. References

## 1. Quick start

1. Open `dist/genome_map_viewer.html` in the browser.
2. Drag all input files onto the drop zone at once (or use *choose files*). The table shows what each file was detected as. If a detection is wrong, change it in the *detected as* column.
3. If you did not load a `contrasts.tsv`, define the pairs in *Samples and contrasts* (treatment vs control, *test* or *reference*, label). Then use *add pair* for each one.
4. Check the *Options* (defaults are explained in §3), then click **Build map**.
5. Tick the contrasts you want as rings, explore, and use **Save report (.html)** to keep or share the result.

To try it, load the six files in `example_R10.79/` (annotation, featureCounts output, sample sheet, contrasts, landmarks, labels).

## 2. Which files do I need?

| What you have | Mode | Files to load |
|---|---|---|
| BAM files | not yet | first count reads with featureCounts or HTSeq-count (§4.2), then Mode A |
| Raw read counts per gene and sample | **A** | annotation + count table + `samplesheet.tsv` + `contrasts.tsv` (or define pairs in the page) |
| DESeq2 / edgeR / limma result tables, one per comparison | **B** | annotation + DE tables + `contrasts.tsv` (with a `de_file` column) |
| A previously saved dataset (`*_genomemap_dataset.json`) | as saved | just that file |
| A saved report (`*_genome_map_report.html`) | as saved | open it directly; nothing to load |

In every case you may also add:

- `landmarks.tsv`: named regions drawn around the circle, e.g. *ina*, the T3SS cluster or rRNA operons;
- `labels.tsv`: curated names for loci whose automatic annotation is uninformative.

The exact column specifications are in `docs/INPUT_FORMATS.md`, and empty templates are in `templates/`.

**The one rule that matters most:** the gene ids in the count or DE table must be the ids of the annotation you load (normally the `locus_tag`). The viewer reports the match rate and refuses to build below 90%.

## 3. Mode A and Mode B: what the numbers mean

### Mode A: counts, computed in the viewer
For each contrast the viewer computes:

1. **Size factors** by median-of-ratios (Anders & Huber 2010) on CDS rows.
    - *DESeq2 standard* (default): only genes with a non-zero count in every library are used. These factors are identical to DESeq2's `estimateSizeFactors` as implemented in pydeseq2 (tested, see `TEST_REPORT.md`).
    - *zero-tolerant*: zeros are treated as missing. This is the variant that the 2026 R10.79 viewer used; keep it only if you need to reproduce that viewer's numbers exactly. On R10.79 the two variants differ by ≤ 0.09%.
2. **Normalised counts** = raw counts / size factor.
3. **Expression filter**: a gene is "expressed" if CPM > 1 in at least 2 libraries (both adjustable). Genes failing the filter are drawn in grey or hidden.
4. **log2 fold change** = log2((mean normalised count in treatment + 1) / (mean normalised count in control + 1)).

**Mode A computes no p-values.** The values are descriptive.

- With one library per side (a "pair"), the log2FC is a single-library ratio, and the viewer shows a *no replication* banner.
- With replicates, Mode A still only compares means. For statistics, run DESeq2 or edgeR and use Mode B.

The pseudocount of 1 shrinks fold changes of lowly expressed genes towards 0. This is intentional: it keeps the map from being dominated by noisy low-count genes.

**Reference contrasts.** Give a contrast `role = reference` when it compares two libraries that should be biologically equivalent. An example is the control from experiment 4 vs the control from experiment 3. The viewer draws the first ticked reference as the inner ring *N*, and optionally as a linear track. It shows how large fold changes get from batch and library noise alone. **Read a test contrast against it**: a gene at log2FC +1.2 in a test pair means little if the noise ring routinely reaches ±1.2 in the same region.

### Mode B: imported DE tables
log2FC and adjusted p-values are taken as they are from your DESeq2 / edgeR / limma tables. Bars with padj ≥ the threshold (default 0.05) are faded, and the switch for this is in the controls. The abundance column (`baseMean`, `logCPM` or `AveExpr`) is shown in the tooltip and in the CSV export. If you used `lfcShrink`, the shrunken values are what you see; say so in the contrast label.

## 4. Preparing the input files

### 4.1 Annotation
Use the annotation the reads were counted against, as GFF3 (Bakta, Prokka, NCBI/PGAP) or GenBank. Multi-contig assemblies and plasmids are supported. Each replicon becomes an arc on the circle, and *Replicon* lets you show one at a time.

### 4.2 Counting reads
**featureCounts** (Liao et al. 2014) is recommended. Use the output file unchanged.

Count directly on a Bakta GFF3, using locus tags as ids:
```
featureCounts -a genome.gff3 -F GFF -t CDS,tRNA,rRNA,ncRNA,tmRNA -g locus_tag \
              -p --countReadPairs -s 2 -T 8 -o counts.txt *.bam
```
Or count on a GTF derived from the same annotation:
```
featureCounts -a genome.gtf -F GTF -t exon -g gene_id -p --countReadPairs -s 2 -T 8 -o counts.txt *.bam
```
This second form is how the R10.79 example counts were produced (featureCounts v2.1.1; see the first line of `example_R10.79/featureCounts_counts_all.txt`).

Options:

- `-p --countReadPairs` is for paired-end reads; omit both for single-end reads.
- `-s` sets strandedness: `0` unstranded, `1` stranded, `2` reverse-stranded. It depends on the library kit, so check it rather than copying it. One practical check is to run featureCounts with `-s 1` and with `-s 2` on one BAM and keep the setting that assigns far more reads. RSeQC `infer_experiment.py` does the same more formally.

**HTSeq-count** (Anders et al. 2015) output is accepted as one file per sample; load all of them together. A generic integer matrix (gene × sample) is also accepted. TPM or FPKM tables are rejected, because size factors must be computed on raw counts.

**tRNA and rRNA** are off by default in *Features shown*. After rRNA depletion their counts mostly reflect depletion efficiency, not expression.

### 4.3 DE tables (Mode B)
In R, for example:
```r
res <- results(dds, contrast = c("condition", "AER", "MC"))
write.csv(as.data.frame(res), "AER_vs_MC.csv")
```
Put the file name in the `de_file` column of `contrasts.tsv`.

## 5. Describing samples and pairs

`samplesheet.tsv` has one row per library. `sample_id` must equal the count-column name; featureCounts BAM paths are shortened automatically to the file name without `.bam`.

`contrasts.tsv` has one row per comparison. `treatment` and `control` take either:

- explicit sample ids separated by `;`, e.g. `EXP3_AER_19626`; or
- a selection rule, e.g. `condition=AER&batch=EXP3`.

Example: the R10.79 set (`example_R10.79/contrasts.tsv`). Each treated library is compared with the control (MC) of its own experiment, plus one noise reference:

| contrast_id | treatment | control | role |
|---|---|---|---|
| E3_AER_vs_MC | EXP3_AER_19626 | EXP3_MC_19626 | test |
| E4_AER_vs_MC | EXP4_AER_23626 | EXP4_MC_23626 | test |
| E4_1h_vs_MC | EXP4_NaCl1hr_23626 | EXP4_MC_23626 | test |
| C_bubbled5h_vs_MC_Tsh | Troubleshoot_BBL_NORM_09626 | BBL_Troubleshoot_MC_NORM_09626 | test |
| noise_MC4_vs_MC3 | EXP4_MC_23626 | EXP3_MC_19626 | reference |

With replicates, the same file becomes, for example, `treatment = condition=AER&batch=EXP5` and `control = condition=MC&batch=EXP5`.

If you define pairs in the page instead, *download samplesheet.tsv + contrasts.tsv* saves them for next time.

## 6. Reading the map

- **Circle, outside → in**:
    - ruler with coordinates per replicon;
    - landmarks;
    - ring *S*: genes by strand (+ strand outer half, dark grey; − strand inner half, light grey);
    - ring *C*: COG functional category of each gene, drawn only when the annotation carries COG assignments (e.g. Bakta);
    - one ring per ticked test contrast (up to 6);
    - the reference ring *N*.
- **Bars** point outward for up-regulation and inward for down-regulation. Colour saturates at the chosen range (±2, ±3 or ±5 log2).
- **Circle bars** can show individual genes (*gene level*) or 10-kb bins: the bin median log2FC, or the fraction of genes with |log2FC| > 1.
- **Navigation**:
    - wheel/pinch zooms the circle and drag pans;
    - shift + drag selects an arc, which opens it in the linear view;
    - a click centres a 30-kb window;
    - *Region* jumps to coordinates;
    - *Search* finds locus tags, gene names, labels and products.
- **Tooltips** show the locus tag, gene, product, COG, every contrast's log2FC, padj (Mode B) and normalised counts per library (Mode A). Duplicated gene symbols (e.g. one name assigned to several loci) are flagged: rely on the locus tag.
- **Exports**:
    - *Export genes in window (CSV)* writes all values for the genes in the linear window;
    - *save view: PNG / SVG / PDF* (above the circle and above the linear view) saves what you are looking at: the circle at its current zoom, or the linear window with its tracks. With *title + legend* ticked, a block underneath gives the title, the mode, the region shown, the ring key (circle) and the colour scale;
    - PNG is 4× the screen resolution. SVG is vector graphics, editable in Inkscape or Illustrator. PDF is one page at the on-screen size holding the same high-resolution image (about 380 dpi); it is not vector, so for a vector PDF open the SVG in Inkscape or Illustrator and save it as PDF.

## 7. Saving and sharing

- **Save report (.html)** writes one self-contained HTML file containing the viewer and your data. Colleagues open it directly in a browser; nothing else is needed. This is the recommended way to share results.
- **Save dataset (.json)** writes only the computed data. Load it into `genome_map_viewer.html` later; it is smaller and can be archived next to the analysis.
- **Load other data** returns to the loader.

Each dataset records how it was made: mode, size-factor variant and values, filter, input file names with SHA-256 checksums, software version, date and warnings. This record appears in the *Data and methods* panel; the full checksums are in the saved dataset JSON (`meta.inputs`) and on hover.

## 8. Command-line builder (optional)

For repeated or batch builds, `cli/genomemap_build.py` produces the same report without a browser. `dist/genomemap.py` is the same tool as one self-contained file (all viewer source is stored inside it), convenient for sending by e-mail: use `python genomemap.py` with the same commands below, plus `python genomemap.py extract DIR` to unpack the source. It needs Python ≥ 3.9 and Node.js ≥ 18; no Python packages are required for building.
```
python cli/genomemap_build.py viewer -o genome_map_viewer.html          # empty viewer
python cli/genomemap_build.py build --annotation genome.gff3 --counts counts.txt \
       --samplesheet samplesheet.tsv --contrasts contrasts.tsv \
       --landmarks landmarks.tsv --labels labels.tsv \
       --title "My strain, RNA-seq" -o my_report.html                  # Mode A
python cli/genomemap_build.py build --annotation genome.gbff --contrasts contrasts.tsv \
       --de-dir de_tables/ -o my_report_modeB.html                       # Mode B
```
`python cli/genomemap_build.py build -h` lists all options (`--sf`, `--cpm`, `--min-libs`, `--types`, `--alpha`, `--bin-width`, `--dataset-json`). The browser and the CLI share the same computation code (`src/genomemap-core.js`) and give identical results (tested).

## 9. Limitations and interpretation

- **Single-library pairs cannot separate treatment from noise.** Without replicates, judge a change only against the reference (noise) ring, and confirm key genes independently, e.g. by RT-qPCR.
- **Median-of-ratios assumes most genes do not change.** If a treatment changes transcription globally, every fold change is relative to the median gene, not absolute. Examples are a general shutdown after aerosolisation, or a large loss of mRNA before fixation. Spike-ins would be needed to measure absolute changes.
- **Sample handling is part of the treatment.** If conditions differ in time-to-fixation or collection medium, those differences are in the fold changes too. Compare each condition with the control processed in the same way and in the same experiment where possible. In the R10.79 design, this is why each treated library is paired with its own experiment's MC.
- **Mode A fold changes are not shrunken or tested.** Use Mode B for gene lists you intend to report as differentially expressed.
- Real-browser rendering was tested only indirectly (see `TEST_REPORT.md` §5).

## 10. Troubleshooting

| Message | Cause | Fix |
|---|---|---|
| *only N% of annotation features have a row in the count table* | ids differ, e.g. counts by `gene_id` = `gene-XYZ` but annotation `locus_tag` = `XYZ`, or the counts belong to another annotation | count with `-g locus_tag`, or load the annotation that was used for counting |
| *sample 'X' not found among count columns* | `sample_id` ≠ count column name | fix `sample_id`, or add `count_column` to the sample sheet |
| *contains non-integer values* | TPM/FPKM or normalised values | provide raw counts |
| *contrast 'X': sample 'Y' is not in the samplesheet* | typo in `contrasts.tsv` | check ids; selection rules use sample-sheet columns |
| *N annotation features have no count row; shown as missing* | a few features (often ncRNA/tmRNA) were not counted | harmless if N is small |
| a `.gz` file is not read | very old browser (needs `DecompressionStream`) | update the browser or decompress first |
| a file is detected as the wrong type | ambiguous table layout | change the type in the *detected as* column |

## 11. References

- Anders S, Huber W (2010). Differential expression analysis for sequence count data. *Genome Biology* 11:R106. (median-of-ratios size factors)
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39:btad547.
- Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30:923-930.
- Anders S, Pyl PT, Huber W (2015). HTSeq: a Python framework to work with high-throughput sequencing data. *Bioinformatics* 31:166-169.
- Schwengers O et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics* 7:000685.
- Bostock M, Ogievetsky V, Heer J (2011). D³ data-driven documents. *IEEE Trans. Vis. Comput. Graph.* 17:2301-2309.
