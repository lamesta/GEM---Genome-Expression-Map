# genomemap 1.0.0 — test report

Date: 1 October 2026. Reproduce with `bash tests/run_all_tests.sh [original_R10.79_viewer.html]`; the run log is written to `tests/out/test_log.txt`.

Software used for testing: Node v26.10.0, Python 3.12.14, numpy 2.5.3, pandas 3.0.6, pytest 9.1.1, pydeseq2 0.5.4, jsdom 24.1.3, cairosvg 2.9.1. The viewer embeds D3 v7.9.0.

**Result of the final run: 17 test groups passed, 0 failed, 0 skipped.** The pytest suite inside that run has 16 tests.

## 1. What was tested and against what

| # | Test | Reference it is compared to | Result |
|---|---|---|---|
| 1 | R10.79 rebuilt from the raw inputs (Bakta GFF3 + featureCounts) with `--sf zero_tolerant` | data embedded in the original R10.79 viewer (artifact `b0b7214b…`) | ALL PASS (35 checks, `tests/regression_R10.79_zero_tolerant.tsv`) |
| 2 | `standard` size factors on the R10.79 CDS counts | pydeseq2 0.5.4 (Python port of DESeq2) | identical (relative difference 0.0 for all 7 libraries; `tests/crosscheck_size_factors_pydeseq2.tsv`) |
| 3 | Mode A arithmetic on a 6-library synthetic data set, both size-factor variants | an independent numpy re-implementation in `tests/test_cli.py` | size factors: max relative difference ≤ 8e-16; log2FC: max abs. difference 5.0e-5 (= rounding to 4 decimals); expression filter identical |
| 4 | Recovery of planted effects (synthetic DC3000 layout) | the values planted by `tests/make_synthetic_dc3000.py` | see §3 |
| 5 | Mode B import (DESeq2 csv; hand-written edgeR and limma tables) | the tables themselves | log2FC identical within 4-decimal rounding; padj within 3-significant-digit rounding; NA padj stay NA; edgeR `FDR`/`logCPM` and limma `adj.P.Val`/`AveExpr` columns recognised |
| 6 | Annotation parsers | the same genome in another format | DC3000 NCBI GFF3 vs GBFF: 5933/5933 features with identical coordinates. R10.79 Bakta GFF3 vs GBK: 5216/5216 identical coordinates, 5215/5216 identical products, 5216/5216 identical COG categories. Prokka-style GFF (2 contigs) parsed |
| 7 | Compressed input | the same files uncompressed | `.gz` annotation + counts give an identical dataset |
| 8 | HTSeq-count input (6 per-sample files, with `__` summary rows) | the same counts as one matrix | identical size factors and log2FC (CLI and browser) |
| 9 | Browser loader (drag-and-drop path) | CLI build of the same files | R10.79: identical size factors; HTSeq and Mode B loads render without errors |
| 10 | Viewer self-test (jsdom) on R10.79, synthetic Mode A (3 replicons), synthetic Mode B | built-in `#selftest` assertions (rings, linear tracks, replicon switching, search, zoom labels, SVG export, CSV export) | pass, no script errors |
| 11 | In-browser "Save report" round trip | — | the saved HTML reopens directly into the map and passes the self-test |
| 12 | Error handling | expected error messages | ids that do not match the annotation, non-integer (TPM-like) counts, and a contrast that names an unknown sample all stop the build with a specific message |
| 13 | File-type auto-detection | known types of 11 example files | all correct (including a landmark table given by `feature_id`) |
| 14 | Figure export (PNG / SVG / PDF, both views) | exported SVG rendered independently (cairosvg); PDF writer output rendered by pdfium | legend, ring key and colour scale present; PDF page size equals on-screen size (1 px = 0.75 pt); rendered PDF image matches the source pixels (mean abs. difference < 3 of 255 for Flate, < 8 for JPEG). The browser's canvas rasterisation step itself was not exercised (jsdom) |

## 2. Expected, documented differences from the original R10.79 viewer

These are deliberate changes, not regressions.

1. **Size factors.** The original viewer used a *zero-tolerant* median-of-ratios (zeros treated as missing). The new default is the standard DESeq2 estimator (genes with a zero in any library are excluded; Anders & Huber 2010). On R10.79 the two variants differ by ≤ 0.09%. The original variant remains available as `zero_tolerant`, and reproduces the original log2FC to the 4th decimal.
2. **Product names are URL-decoded.** 405 products in the original viewer showed a literal `%2C` instead of a comma.
3. **Bin boundaries.** Genes are now assigned to 10-kb bins by start coordinate, `floor((start − 1) / width)`; the original used the gene midpoint. Two bins (230 and 231) change membership; all other 586 bins are identical within 1e-4.
4. **Normalised counts** agree with the original within its 0.1 rounding.

## 3. Synthetic multi-replicon test (not biological data)

`tests/synthetic_DC3000/` holds negative-binomial counts (dispersion 0.05) generated on the real NCBI DC3000 annotation (GCF_000007805.1: chromosome NC_004578.1, plasmids pDC3000A NC_004633.1 and pDC3000B NC_004632.1). There are 5850 CDS rows and 3 MC + 3 AER libraries. The true library size factors were 0.9, 1.1, 1.0, 0.8, 1.2 and 1.05. Planted effects in AER: all 77 pDC3000A genes +2 log2; the 87 chromosome genes between 1.00 and 1.10 Mb −1.5 log2. Nothing else was planted.

| Quantity | Planted | Mode A (viewer, standard SF) | Mode B (pydeseq2 Wald) |
|---|---|---|---|
| size factors MC_1…AER_3 | 0.9, 1.1, 1.0, 0.8, 1.2, 1.05 | 0.910, 1.112, 1.010, 0.807, 1.203, 1.059 | — |
| pDC3000A median log2FC | +2.0 | +2.03 | mean +2.04 |
| chromosome 1.00–1.10 Mb median log2FC | −1.5 | −1.49 | mean −1.49 |
| rest of chromosome, median | 0 | 0.00 | — |
| pDC3000B (no effect planted) | 0 | +0.10 (median) | mean +0.08 (s.e. 0.03), 0 genes with padj < 0.05 |
| genes with padj < 0.05 | 164 planted | — (Mode A has no test) | 167 |

The small positive shift on pDC3000B (78 genes) is not planted. No gene there reaches padj < 0.05. We read it as sampling noise in a small replicon, and it is a reminder that a replicon-wide median of a few dozen genes is not by itself evidence of an effect.

The circle and linear renders (`docs/img/synth_A_circle.png`, `synth_B_circle.png`, `synth_A_linear.png`; R10.79: `docs/img/r1079_report_circle.png`, `r1079_report_linear.png`) show three things:
- the three replicons, separated by gaps, with per-replicon coordinates;
- the +2 block on pDC3000A and the −1.5 block at 1.0–1.1 Mb;
- in Mode B, faded bars for padj ≥ 0.05.

Replicons too short for a printed name are identified by hovering over their ruler arc; clicking the arc shows that replicon alone.

## 4. Bugs found and fixed by these tests

| Bug | Found by | Fix |
|---|---|---|
| HTSeq files named `sample.htseq.txt` were given the sample id `sample.htseq` | pytest HTSeq test | all known extensions are stripped |
| A large single HTSeq file could be classified as a count matrix, because type detection looked at a 20-kb slice whose last line was cut | browser-loader HTSeq test | detection ignores the truncated last line |
| A landmark table given only by `feature_id` was classified as a label table | review while writing this report | detection now keys on the `name` column |
| Labels of small replicons and their "0" ticks overlapped | visual check of the synthetic render | names are printed only where the arc is long enough; otherwise shown on hover |
| `pydeseq2` failed in the sandbox (multiprocessing blocked) | synthetic data generation | `DefaultInference(n_cpus=1)`; irrelevant for normal use |

## 5. Not tested / limitations

- **Real browsers.** Headless Chrome could not run in the build environment, so the UI tests use jsdom with a stub canvas. Canvas drawing itself was not exercised; the SVG export path, which draws through the same painter interface, was. Opening `dist/genome_map_viewer.html` once in Chrome, Firefox and Safari before distribution is still recommended.
- **edgeR and limma tables** were tested only on small hand-written examples in the layout written by `write.csv(topTags(...))` and `write.csv(topTable(...))` (row names in the first, unnamed column). They were not tested on outputs from real analyses.
- **RAST and PGAP GFF3 variants** beyond the NCBI RefSeq DC3000 file were not tested.
- **Very large genomes or many libraries** (> 10 Mb, > 50 libraries) were not benchmarked; the largest data sets tested are R10.79 (5.9 Mb, 7 libraries) and the synthetic DC3000 set (6.5 Mb, 3 replicons, 6 libraries).
