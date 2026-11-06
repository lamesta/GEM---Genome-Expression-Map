// Build the R10.79 dataset with the core (legacy zero-tolerant variant) and write JSON for regression.
const fs = require("fs"); const C = require("../src/genomemap-core.js");
const ex = __dirname + "/../example_R10.79/";
const rd = f => fs.readFileSync(ex + f, "utf8");
const variant = process.argv[2] || "zero_tolerant";
const t0 = Date.now();
const types = ["R10.79_bakta.gff3", "featureCounts_counts_all.txt", "samplesheet.tsv", "contrasts.tsv", "landmarks.tsv", "labels.tsv"].map(f => [f, C.detectType(f, rd(f))]);
const ann = C.parseGFF3(rd("R10.79_bakta.gff3"));
const cnt = C.parseFeatureCounts(rd("featureCounts_counts_all.txt"));
const V = C.computeModeA(ann, cnt, C.parseSampleSheet(rd("samplesheet.tsv")), C.parseContrasts(rd("contrasts.tsv")), { sfVariant: variant });
const ds = C.buildDataset({ annotation: ann, values: V, landmarks: C.parseLandmarks(rd("landmarks.tsv")), labels: C.parseLabels(rd("labels.tsv")), title: "Pseudomonas syringae R10.79", organism: "Pseudomonas syringae", strain: "R10.79" });
const bins = C.computeBins(ds, 10000);
fs.mkdirSync(__dirname + "/out", { recursive: true }); fs.writeFileSync(__dirname + `/out/out_r1079_${variant}.json`, JSON.stringify({ types, ds, bins, ms: Date.now() - t0 }));
console.log(JSON.stringify({ types, nRep: ann.replicons.length, rep: ann.replicons, nAnn: ann.features.length, nShown: ds.features.id.length, sf: V.sizeFactors, report: V.report, warnings: ds.meta.warnings, ms: Date.now() - t0 }, null, 1));
