#!/usr/bin/env bash
# Runs the complete genome-map test battery and writes tests/out/test_log.txt.
# Usage (from the package root):  bash tests/run_all_tests.sh [ORIGINAL_R10.79_VIEWER.html]
# Required: node >= 18, python >= 3.9 with numpy, pandas, pytest.
# Optional: pydeseq2 (size-factor cross-check), jsdom in tests/jsdom_env (UI tests:
#           cd tests/jsdom_env && npm install jsdom@24), cairosvg (PNG renders of the SVG exports),
#           the original R10.79 viewer HTML as $1 (regression against the previous viewer).
set -u
cd "$(dirname "$0")/.."
mkdir -p tests/out
LOG=tests/out/test_log.txt; : > "$LOG"
PASS=0; FAIL=0; SKIP=0
run() { local name="$1"; shift; echo "=== $name" | tee -a "$LOG"
  if "$@" >> "$LOG" 2>&1; then echo "    PASS" | tee -a "$LOG"; PASS=$((PASS+1)); else echo "    FAIL" | tee -a "$LOG"; FAIL=$((FAIL+1)); fi; }
skip() { echo "=== $1"$'\n'"    SKIP ($2)" | tee -a "$LOG"; SKIP=$((SKIP+1)); }
S=tests/synthetic_DC3000

run "build empty viewer (dist/genome_map_viewer.html)" python cli/genomemap_build.py viewer -o dist/genome_map_viewer.html
run "build R10.79 example report" python cli/genomemap_build.py build --annotation example_R10.79/R10.79_bakta.gff3 \
  --counts example_R10.79/featureCounts_counts_all.txt --samplesheet example_R10.79/samplesheet.tsv --contrasts example_R10.79/contrasts.tsv \
  --landmarks example_R10.79/landmarks.tsv --labels example_R10.79/labels.tsv --title "P. syringae R10.79 RNA-seq genome map" --strain R10.79 \
  --dataset-json example_R10.79/R10.79_dataset.json -o example_R10.79/R10.79_genome_map_report.html
run "build synthetic DC3000 Mode A report" python cli/genomemap_build.py build --annotation tests/fixtures/DC3000_NCBI.gff \
  --counts $S/SYNTHETIC_counts_matrix.tsv --samplesheet $S/samplesheet.tsv --contrasts $S/contrasts.tsv --landmarks $S/landmarks.tsv \
  --title "SYNTHETIC test: P. syringae pv. tomato DC3000 layout" --strain "DC3000 (synthetic)" -o tests/out/SYNTH_DC3000_modeA.html
run "build synthetic DC3000 Mode B report (GenBank + DESeq2 table)" python cli/genomemap_build.py build --annotation tests/fixtures/DC3000_NCBI.gbff \
  --contrasts $S/contrasts_modeB.tsv --de-dir $S --landmarks $S/landmarks.tsv \
  --title "SYNTHETIC test: DC3000, Mode B (pydeseq2)" --strain "DC3000 (synthetic)" -o tests/out/SYNTH_DC3000_modeB.html
run "parser tests (GFF3 Bakta/NCBI/Prokka, GenBank, DE tables)" bash -c "node tests/test_parsers.js > tests/out/parsers.json && tail -1 tests/out/parsers.json | grep -q 'PARSER TESTS PASSED'"
run "pytest CLI suite (independent numpy Mode A, gz, HTSeq, Mode B, errors)" python -m pytest -q -p no:cacheprovider tests/test_cli.py
run "core build of R10.79 (standard + zero_tolerant) for regression" bash -c "node tests/build_r1079_dataset.js standard >/dev/null && node tests/build_r1079_dataset.js zero_tolerant >/dev/null"
if python -c "import pydeseq2" 2>/dev/null; then run "size factors vs pydeseq2" python tests/crosscheck_pydeseq2.py; else skip "size factors vs pydeseq2" "pydeseq2 not installed"; fi
if [ $# -ge 1 ] && [ -f "$1" ]; then run "regression vs original R10.79 viewer" bash -c "python tests/regression_r1079.py '$1' zero_tolerant | tail -1 | grep -q 'ALL PASS'"
else skip "regression vs original R10.79 viewer" "original viewer HTML not given"; fi
if [ -d tests/jsdom_env/node_modules/jsdom ]; then
  chk() { python -c "import json,sys; a=json.load(open(sys.argv[1])); st=a.get('selftest') or {}; sys.exit(0 if (not a['errors'] and st.get('ok')) else 1)" "$1"; }
  run "UI self-test: R10.79 report" bash -c "node tests/ui_test.js example_R10.79/R10.79_genome_map_report.html tests/out/r1079_report >/dev/null && $(declare -f chk); chk tests/out/r1079_report_result.json"
  run "UI self-test: synthetic Mode A (3 replicons)" bash -c "node tests/ui_test.js tests/out/SYNTH_DC3000_modeA.html tests/out/synth_A >/dev/null && $(declare -f chk); chk tests/out/synth_A_result.json"
  run "UI self-test: synthetic Mode B (padj fade)" bash -c "node tests/ui_test.js tests/out/SYNTH_DC3000_modeB.html tests/out/synth_B >/dev/null && $(declare -f chk); chk tests/out/synth_B_result.json"
  run "browser loader: featureCounts R10.79 files" bash -c "node tests/ui_test.js dist/genome_map_viewer.html tests/out/loader_r1079 example_R10.79/R10.79_bakta.gff3 example_R10.79/featureCounts_counts_all.txt example_R10.79/samplesheet.tsv example_R10.79/contrasts.tsv example_R10.79/landmarks.tsv example_R10.79/labels.tsv >/dev/null && $(declare -f chk); chk tests/out/loader_r1079_result.json"
  run "browser loader: 6 HTSeq files + samplesheet" bash -c "node tests/ui_test.js dist/genome_map_viewer.html tests/out/loader_htseq tests/fixtures/DC3000_NCBI.gff $S/htseq/*.txt $S/samplesheet.tsv $S/contrasts.tsv $S/landmarks.tsv >/dev/null && $(declare -f chk); chk tests/out/loader_htseq_result.json"
  run "browser loader: Mode B (GenBank + DESeq2 csv)" bash -c "node tests/ui_test.js dist/genome_map_viewer.html tests/out/loader_modeB tests/fixtures/DC3000_NCBI.gbff $S/SYNTHETIC_DESeq2_AER_vs_MC.csv $S/contrasts_modeB.tsv $S/landmarks.tsv >/dev/null && $(declare -f chk); chk tests/out/loader_modeB_result.json"
  run "browser 'Save report' round trip" node tests/save_report_test.js
  if python -c "import cairosvg" 2>/dev/null; then
    run "rasterise SVG exports to PNG" python -c "
import cairosvg
for p in ['r1079_report_circle','r1079_report_linear','synth_A_circle','synth_B_circle','synth_A_linear']:
    cairosvg.svg2png(url=f'tests/out/{p}.svg', write_to=f'tests/out/{p}.png', output_width=1400)"
  fi
else skip "UI / browser-loader tests" "jsdom not installed in tests/jsdom_env"; fi
echo "SUMMARY: $PASS passed, $FAIL failed, $SKIP skipped" | tee -a "$LOG"
[ $FAIL -eq 0 ]
