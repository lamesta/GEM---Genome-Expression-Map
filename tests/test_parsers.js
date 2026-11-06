const fs = require("fs"); const C = require("../src/genomemap-core.js"); const assert = require("assert");
const fx = f => fs.readFileSync(__dirname + "/fixtures/" + f, "utf8");
const res = {};
// NCBI PGAP GFF (DC3000: chromosome + 2 plasmids)
let t = fx("DC3000_NCBI.gff"); res.dc_gff_type = C.detectType("x", t);
let a = C.parseGFF3(t); res.dc_gff = { reps: a.replicons, n: a.features.length, types: count(a.features.map(f => f.type)), ex: a.features.slice(0, 2), withProduct: a.features.filter(f => f.product).length, warnings: a.warnings };
// NCBI GenBank (same genome)
t = fx("DC3000_NCBI.gbff"); res.dc_gbk_type = C.detectType("x", t);
let b = C.parseGenBank(t); res.dc_gbk = { reps: b.replicons, n: b.features.length, types: count(b.features.map(f => f.type)), ex: b.features.slice(0, 2) };
const ga = new Map(a.features.map(f => [f.id, f])); let same = 0, coord = 0; for (const f of b.features) { const g = ga.get(f.id); if (g) { same++; if (g.s === f.s && g.e === f.e && g.st === f.st) coord++; } }
res.dc_gff_vs_gbk = { gbkFeatures: b.features.length, idsInGff: same, identicalCoords: coord };
// Bakta GenBank vs Bakta GFF (R10.79)
const rg = C.parseGFF3(fs.readFileSync(__dirname + "/../example_R10.79/R10.79_bakta.gff3", "utf8")); const rb = C.parseGenBank(fx("R10.79_bakta.gbff"));
const m = new Map(rg.features.map(f => [f.id, f])); let s2 = 0, c2 = 0, p2 = 0, cog2 = 0; for (const f of rb.features) { const g = m.get(f.id); if (g) { s2++; if (g.s === f.s && g.e === f.e && g.st === f.st) c2++; if (g.product === f.product) p2++; if (g.cogcat === f.cogcat) cog2++; } }
res.r1079_gff_vs_gbk = { gff: rg.features.length, gbk: rb.features.length, shared: s2, identicalCoords: c2, identicalProduct: p2, identicalCOGcat: cog2, gbkReps: rb.replicons };
// Prokka-style synthetic
const p = C.parseGFF3(fx("prokka_style_synthetic.gff")); res.prokka = { reps: p.replicons, feats: p.features.map(f => [f.id, f.rep, f.type, f.gene, f.product, f.cog, f.group]) };
assert.strictEqual(p.features.length, 5); assert.strictEqual(p.replicons.length, 2); assert.strictEqual(p.features[4].product, "Replication protein, RepA family");
// group resolution
const ss = [{ id: "a", condition: "MC", batch: "E3" }, { id: "b", condition: "AER", batch: "E3" }, { id: "c", condition: "MC", batch: "E4" }];
assert.deepStrictEqual(C.resolveGroup("condition=MC", ss), ["a", "c"]); assert.deepStrictEqual(C.resolveGroup("condition=MC&batch=E4", ss), ["c"]); assert.deepStrictEqual(C.resolveGroup("a;b", ss), ["a", "b"]);
// size factors: standard variant on a toy matrix (hand-computable)
const toy = [[10, 20, 0, 5], [20, 40, 3, 10], [40, 80, 6, 20]];
const sfS = C.sizeFactors(toy, null, "standard"); res.toy_sf_standard = sfS;
// rows 0,1,3 used; geo means; ratios lib0: 10/20,20/40,5/10 = 0.5 → exactly 0.5,1,2
assert(Math.abs(sfS[0] - 0.5) < 1e-12 && Math.abs(sfS[1] - 1) < 1e-12 && Math.abs(sfS[2] - 2) < 1e-12);
// DE table parsing (DESeq2 write.csv with row names; edgeR)
const de1 = C.parseDETable('"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"\n"G1",10.5,1.2,0.3,4,0.0001,0.001\n"G2",3,NA,NA,NA,NA,NA\n');
const de2 = C.parseDETable("gene_id\tlogFC\tlogCPM\tPValue\tFDR\nG1\t-0.5\t4.2\t0.2\t0.4\n");
res.de = { deseq2: de1, edger: de2 };
assert.deepStrictEqual(de1.genes, ["G1", "G2"]); assert.strictEqual(de1.log2fc[1], null); assert.strictEqual(de2.padj[0], 0.4);
function count(a) { const o = {}; for (const x of a) o[x] = (o[x] || 0) + 1; return o; }
fs.mkdirSync(__dirname + "/out", { recursive: true }); fs.writeFileSync(__dirname + "/out/out_parsers.json", JSON.stringify(res, null, 1)); console.log(JSON.stringify(res, (k, v) => (k === "ex" ? undefined : v), 1).slice(0, 4000));
console.log("PARSER TESTS PASSED");
