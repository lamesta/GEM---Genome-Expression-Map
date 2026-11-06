#!/usr/bin/env node
/* Node entry used by genomemap_build.py: reads a JSON config (argv[2]), builds the dataset with genomemap-core.js, writes JSON (argv[3]). */
const fs = require("fs"); const path = require("path"); const zlib = require("zlib"); const crypto = require("crypto");
const C = require(path.join(__dirname, "..", "src", "genomemap-core.js"));
const cfg = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const read = f => { const b = fs.readFileSync(f); return (/\.gz$/i.test(f) ? zlib.gunzipSync(b) : b).toString("utf8"); };
const inputs = [];
const R = f => { const t = read(f); inputs.push({ file: path.basename(f), bytes: t.length, sha256: crypto.createHash("sha256").update(t).digest("hex") }); return t; };
try {
  const annText = R(cfg.annotation); const annType = C.detectType(cfg.annotation, annText);
  const ann = annType === "genbank" ? C.parseGenBank(annText) : C.parseGFF3(annText);
  const showTypes = cfg.types || ["CDS", "ncRNA", "tmRNA"];
  let V;
  if (cfg.counts) {
    let cnt;
    if (cfg.counts.length > 1) cnt = C.parseHTSeqFiles(cfg.counts.map(f => ({ name: path.basename(f), text: R(f) })));
    else { const t = R(cfg.counts[0]); const ty = C.detectType(cfg.counts[0], t); cnt = ty === "featurecounts" ? C.parseFeatureCounts(t) : (ty === "htseq" ? C.parseHTSeqFiles([{ name: path.basename(cfg.counts[0]), text: t }]) : C.parseMatrix(t)); }
    const samples = cfg.samplesheet ? C.parseSampleSheet(R(cfg.samplesheet)) : cnt.samples.map(s => ({ id: s, label: s, condition: s, replicate: "", batch: "", count_column: "" }));
    V = C.computeModeA(ann, cnt, samples, C.parseContrasts(R(cfg.contrasts)), { sfVariant: cfg.sf || "standard", cpm: cfg.cpm ?? 1, minLibs: cfg.min_libs ?? 2, showTypes, sfTypes: cfg.sf_types || ["CDS"] });
  } else {
    const contrasts = C.parseContrasts(R(cfg.contrasts)); const tables = {};
    for (const c of contrasts) tables[c.id] = C.parseDETable(R(path.join(cfg.de_dir || path.dirname(cfg.contrasts), c.de_file)));
    V = C.computeModeB(ann, tables, contrasts, { showTypes, alpha: cfg.alpha ?? 0.05 });
  }
  const ds = C.buildDataset({ annotation: ann, values: V, landmarks: cfg.landmarks ? C.parseLandmarks(R(cfg.landmarks)) : [], labels: cfg.labels ? C.parseLabels(R(cfg.labels)) : {},
    title: cfg.title || ann.organism || "", organism: ann.organism || "", strain: cfg.strain || "", inputs });
  ds.meta.bin_width = cfg.bin_width || 10000; if (V.meta.mode === "B") ds.meta.alpha = cfg.alpha ?? 0.05;
  ds.meta.build = { tool: "genomemap_build.py", config: Object.fromEntries(Object.entries(cfg).map(([k, v]) => [k, Array.isArray(v) ? v.map(x => path.basename(String(x))) : (typeof v === "string" && fs.existsSync(v) ? path.basename(v) : v)])) };
  fs.writeFileSync(process.argv[3], JSON.stringify(C.compactDataset(ds)));
  process.stdout.write(JSON.stringify({ ok: true, features: ds.features.id.length, contrasts: ds.contrasts.length, warnings: ds.meta.warnings, report: V.report, size_factors: V.sizeFactors || null }) + "\n");
} catch (e) { process.stdout.write(JSON.stringify({ ok: false, error: e.message }) + "\n"); process.exit(2); }
