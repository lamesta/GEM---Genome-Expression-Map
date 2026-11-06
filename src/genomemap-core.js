/* genomemap-core.js — parsers and computations shared by the HTML viewer and the CLI builder.
 * Pure JavaScript, no dependencies. Works in the browser (window.GenomeMapCore) and in Node (module.exports).
 * Coordinates: 1-based, inclusive. All tables: UTF-8 text with a header row.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.GenomeMapCore = factory();
}(typeof self !== "undefined" ? self : this, function () {
  "use strict";
  const VERSION = "1.0.0";

  // ------------------------------------------------------------------ text utilities
  function lines(text) { return text.replace(/\r\n?/g, "\n").split("\n"); }
  function stripBOM(t) { return t.charCodeAt(0) === 0xFEFF ? t.slice(1) : t; }
  function splitCSVLine(line, delim) {
    if (delim === "\t" || line.indexOf('"') < 0) return line.split(delim);
    const out = []; let cur = "", q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
      else if (c === '"') q = true; else if (c === delim) { out.push(cur); cur = ""; } else cur += c;
    }
    out.push(cur); return out;
  }
  /** Parse a delimited table. Lines starting with '#' are skipped. Returns {header, rows, delim}. */
  function parseTable(text) {
    const L = lines(stripBOM(text)).filter(l => l.length && l[0] !== "#");
    if (!L.length) return { header: [], rows: [], delim: "\t" };
    const first = L[0];
    const delim = (first.split("\t").length >= first.split(",").length) ? "\t" : ",";
    const header = splitCSVLine(first, delim).map(h => h.trim());
    const rows = [];
    for (let i = 1; i < L.length; i++) { const r = splitCSVLine(L[i], delim); if (r.length === 1 && r[0].trim() === "") continue; rows.push(r.map(x => x.trim())); }
    return { header, rows, delim };
  }
  const lc = s => String(s || "").toLowerCase();
  function findCol(header, names) {
    const H = header.map(lc);
    for (const n of names) { const i = H.indexOf(lc(n)); if (i >= 0) return i; }
    return -1;
  }
  function basenameSample(s) { return String(s).split(/[\\/]/).pop().replace(/\.(bam|sam|cram)$/i, "").replace(/(\.sorted|\.dedup)$/i, ""); }

  // ------------------------------------------------------------------ file type detection
  /** Returns one of: gff3, genbank, featurecounts, htseq, matrix, de_table, samplesheet, contrasts, landmarks, labels, dataset_json, unknown */
  function detectType(name, text) {
    const full = stripBOM(text); let t = full.slice(0, 20000);
    if (full.length > 20000 && t.lastIndexOf("\n") > 0) t = t.slice(0, t.lastIndexOf("\n")); // never classify on a truncated last line
    const head = t.trimStart();
    if (head.startsWith("{")) { try { const j = JSON.parse(text); if (j && j.format === "genomemap-dataset") return "dataset_json"; } catch (e) { } return "unknown"; }
    if (/^##gff-version\s+3/m.test(t.slice(0, 500))) return "gff3";
    if (/^LOCUS\s+\S+/.test(head)) return "genbank";
    if (/^#\s*Program:featureCounts/m.test(t) || /^Geneid\tChr\tStart\tEnd\tStrand\tLength/m.test(t)) return "featurecounts";
    const tab = parseTable(t.split("\n").slice(0, 200).join("\n"));
    const H = tab.header.map(lc);
    if (H.length === 0) return "unknown";
    const has = n => H.indexOf(n) >= 0;
    if (has("contrast_id")) return "contrasts";
    if (has("sample_id") && (has("condition") || has("label"))) return "samplesheet";
    if (has("name") && ((has("start") && has("end")) || has("feature_id"))) return "landmarks";
    if (has("feature_id") && has("label") && !has("start")) return "labels";
    if (has("log2foldchange") || has("logfc") || has("log2fc")) return "de_table";
    // HTSeq: two columns, no header-like first line, typically ends with __ rows
    const raw = lines(t).filter(l => l.length);
    if (raw.length && raw.every(l => l.split("\t").length === 2) && /\t\d+$/.test(raw[0]) ) return "htseq";
    if (tab.header.length >= 2 && tab.rows.length && tab.rows.slice(0, 50).every(r => r.slice(1).every(v => v === "" || /^-?\d+(\.0+)?$/.test(v)))) return "matrix";
    return "unknown";
  }

  // ------------------------------------------------------------------ annotation parsers
  const DEFAULT_TYPES = ["CDS", "ncRNA", "tmRNA", "tRNA", "rRNA"];
  // Sequence Ontology / INSDC RNA subtypes reported as ncRNA (NCBI GFF uses SO terms such as RNase_P_RNA, SRP_RNA, antisense_RNA)
  const NCRNA_SUBTYPES = new Set(["ncRNA", "RNase_P_RNA", "SRP_RNA", "antisense_RNA", "sRNA", "misc_RNA", "RNase_MRP_RNA", "riboswitch_RNA"]);
  function normType(t) { return NCRNA_SUBTYPES.has(t) ? "ncRNA" : t; }
  function decodeAttr(v) { try { return decodeURIComponent(v); } catch (e) { return v; } }
  function parseAttrs(s) {
    const o = {};
    for (const kv of s.split(";")) { if (!kv) continue; const i = kv.indexOf("="); if (i < 0) continue; o[kv.slice(0, i).trim()] = decodeAttr(kv.slice(i + 1).trim()); }
    return o;
  }
  function cogFromXrefs(x) {
    let cog = "", cat = "";
    if (!x) return { cog, cat };
    for (const d of x.split(",")) {
      const m = d.trim().match(/^COG:(.+)$/i); if (!m) continue;
      const v = m[1].trim();
      if (/^COG\d+$/i.test(v)) cog = cog ? cog + "," + v.toUpperCase() : v.toUpperCase();
      else if (/^[A-Z]{1,6}$/.test(v)) cat = v;
    }
    return { cog, cat };
  }
  /** GFF3 → {replicons:[{id,name,length,topology}], features:[{id,rep,s,e,st,type,gene,product,cog,cogcat}]} */
  function parseGFF3(text, opts) {
    opts = opts || {}; const types = new Set(opts.types || DEFAULT_TYPES);
    const reps = new Map(); const repOrder = [];
    const addRep = (id, len, top) => { if (!reps.has(id)) { reps.set(id, { id, name: id, length: len || 0, topology: top || "linear" }); repOrder.push(id); } const r = reps.get(id); if (len && len > r.length) r.length = len; if (top) r.topology = top; };
    const genes = new Map(); const kids = []; let organism = "";
    for (const l of lines(stripBOM(text))) {
      if (l.startsWith("##FASTA")) break;
      if (l.startsWith("##sequence-region")) { const p = l.split(/\s+/); addRep(p[1], +p[3]); continue; }
      if (l.startsWith("# organism ")) { organism = l.slice(11).trim(); continue; }
      if (!l || l[0] === "#") continue;
      const f = l.split("\t"); if (f.length < 9) continue;
      const [seq, , type, s, e, , strand, , attr] = f; const a = parseAttrs(attr);
      if (type === "region" || type === "chromosome" || type === "plasmid" || type === "contig" || type === "source") {
        addRep(seq, +e, (a.Is_circular === "true" || a.is_circular === "true") ? "circular" : (reps.has(seq) ? null : "linear"));
        if (a.Name && a.Name !== seq && type !== "region") reps.get(seq).name = a.Name;
        if (type === "region" && a.genome === "chromosome" && reps.get(seq).name === seq) reps.get(seq).name = "chromosome";
        if (type === "region" && a.genome === "plasmid" && a.plasmid_name) reps.get(seq).name = a["plasmid-name"] || a.plasmid_name;
        if (type === "plasmid" || (a.plasmid_name || a["plasmid-name"])) reps.get(seq).name = a.plasmid_name || a["plasmid-name"] || reps.get(seq).name;
        continue;
      }
      if (!reps.has(seq)) addRep(seq, 0, null);
      const rec = { seq, type, s: +s, e: +e, st: strand === "+" ? 1 : (strand === "-" ? -1 : 0), a };
      if (type === "gene" || type === "pseudogene") { if (a.ID) genes.set(a.ID, rec); }
      else if (types.has(normType(type, a))) { rec.type = normType(type, a); kids.push(rec); }
    }
    const out = []; const seen = new Map(); const usedGenes = new Set();
    for (const k of kids) {
      const par = k.a.Parent ? genes.get(k.a.Parent.split(",")[0]) : null; if (par) usedGenes.add(k.a.Parent.split(",")[0]);
      const id = k.a.locus_tag || (par && par.a.locus_tag) || k.a.old_locus_tag || k.a.ID || k.a.Name;
      if (!id) continue;
      const x = cogFromXrefs(k.a.Dbxref || k.a.db_xref || (par && (par.a.Dbxref || par.a.db_xref)));
      const feat = { id, seq: k.seq, s: Math.min(k.s, k.e), e: Math.max(k.s, k.e), st: k.st, type: k.type,
        gene: k.a.gene || (par && (par.a.gene || par.a.Name)) || "", product: k.a.product || (k.a.Name && k.a.Name !== id ? k.a.Name : "") || "",
        cog: x.cog, cogcat: x.cat, pseudo: (k.a.pseudo === "True" || k.a.pseudo === "true" || !!(par && par.type === "pseudogene")) };
      if (seen.has(id)) { const o = seen.get(id); o.s = Math.min(o.s, feat.s); o.e = Math.max(o.e, feat.e); continue; } // multi-segment CDS
      seen.set(id, feat); out.push(feat);
    }
    if (opts.includeGeneOnly) for (const [gid, g] of genes) { if (usedGenes.has(gid)) continue; const id = g.a.locus_tag || gid; if (seen.has(id)) continue; const f = { id, seq: g.seq, s: g.s, e: g.e, st: g.st, type: g.type, gene: g.a.gene || g.a.Name || "", product: g.a.product || g.a.description || "", cog: "", cogcat: "", pseudo: g.type === "pseudogene" }; seen.set(id, f); out.push(f); }
    const res = finishAnnotation(repOrder.map(id => reps.get(id)), out, "GFF3"); res.organism = organism; return res;
  }
  function parseGBLocation(loc) {
    const st = /complement\(/.test(loc) ? -1 : 1;
    const nums = (loc.match(/\d+/g) || []).map(Number);
    if (!nums.length) return null;
    return { s: Math.min(...nums), e: Math.max(...nums), st, joinWraps: /join\(/.test(loc) && nums.length >= 4 && nums[0] > nums[nums.length - 1] };
  }
  /** GenBank flat file (one or more LOCUS entries) → same structure as parseGFF3 */
  function parseGenBank(text, opts) {
    opts = opts || {}; const types = new Set(opts.types || DEFAULT_TYPES);
    const reps = []; const feats = []; const geneQ = new Map();
    const L = lines(stripBOM(text)); let i = 0; let rep = null; let inFeat = false; let cur = null; let organism = "";
    const flush = () => { if (!cur) return; const loc = parseGBLocation(cur.loc); if (loc) { cur.s = loc.s; cur.e = loc.e; cur.st = loc.st; cur.wraps = loc.joinWraps; cur.seq = rep.id;
          if (cur.key === "source") { if (cur.q.plasmid) rep.name = cur.q.plasmid; else if (cur.q.chromosome) rep.name = "chromosome " + cur.q.chromosome; }
          else if (cur.key === "gene") { const lt = cur.q.locus_tag; if (lt) geneQ.set(rep.id + "|" + lt, cur); } else if (types.has(normType(cur.key, cur.q))) { cur.key = normType(cur.key, cur.q); feats.push(cur); } } cur = null; };
    for (; i < L.length; i++) {
      const l = L[i];
      if (l.startsWith("LOCUS")) { const p = l.split(/\s+/); const len = +(l.match(/(\d+)\s+bp/) || [0, 0])[1]; rep = { id: p[1], name: p[1], length: len, topology: /\bcircular\b/i.test(l) ? "circular" : "linear" }; reps.push(rep); inFeat = false; continue; }
      if (!rep) continue;
      if (/^ {2}ORGANISM/.test(l) && !organism) { organism = l.slice(12).trim(); continue; }
      if (l.startsWith("DEFINITION")) { const d = l.slice(12); const pm = d.match(/plasmid\s+([^\s,]+)/i); if (pm) rep.name = pm[1]; else if (/chromosome/i.test(d)) rep.name = "chromosome"; continue; }
      if (l.startsWith("VERSION")) { const v = l.split(/\s+/)[1]; if (v) { rep.accession = v; rep.locusName = rep.id; if (rep.name === rep.id) rep.name = v; rep.id = v; } continue; }
      if (l.startsWith("FEATURES")) { inFeat = true; continue; }
      if (l.startsWith("ORIGIN") || l.startsWith("//") || l.startsWith("CONTIG")) { flush(); inFeat = false; continue; }
      if (!inFeat) continue;
      const fm = l.match(/^ {5}(\S+)\s+(\S.*)$/);
      if (fm) { flush(); cur = { key: fm[1], loc: fm[2].trim(), q: {}, last: null }; continue; }
      if (!cur) continue;
      const t = l.trim();
      if (t.startsWith("/")) { const m = t.match(/^\/([^=]+)(=(.*))?$/); if (m) { const k = m[1]; let v = m[3] === undefined ? "true" : m[3]; v = v.replace(/^"/, ""); cur.q[k] = cur.q[k] && k === "db_xref" ? cur.q[k] + "," + v.replace(/"$/, "") : v.replace(/"$/, ""); cur.last = k; cur.open = /^"/.test(m[3] || "") && !/"$/.test(m[3] || "\""); } }
      else if (cur.last && cur.open) { cur.q[cur.last] += (cur.last === "translation" ? "" : " ") + t.trim().replace(/"$/, ""); if (/"$/.test(t.trim())) cur.open = false; }
      else if (!cur.last) cur.loc += t.trim();
    }
    flush();
    for (const r of reps) { /* source /plasmid qualifier is not tracked here; replicon name = LOCUS name */ }
    const out = []; const seen = new Map();
    for (const f of feats) {
      const g = f.q.locus_tag ? geneQ.get(f.seq + "|" + f.q.locus_tag) : null;
      const id = f.q.locus_tag || f.q.old_locus_tag || (g && g.q.locus_tag) || f.q.protein_id; if (!id) continue;
      const x = cogFromXrefs((f.q.db_xref || "").split(",").map(s => s.trim()).join(","));
      const feat = { id, seq: f.seq, s: f.s, e: f.e, st: f.st, type: f.key, gene: f.q.gene || (g && g.q.gene) || "", product: f.q.product || "", cog: x.cog, cogcat: x.cat, pseudo: !!f.q.pseudo };
      if (seen.has(id)) continue; seen.set(id, feat); out.push(feat);
    }
    const res = finishAnnotation(reps, out, "GenBank"); res.organism = organism; return res;
  }
  const COG_GROUP = {};
  "JAKLB".split("").forEach(c => COG_GROUP[c] = "info");
  "DYVTMNZWUOX".split("").forEach(c => COG_GROUP[c] = "cell");
  "CGEFHIPQ".split("").forEach(c => COG_GROUP[c] = "metab");
  "RS".split("").forEach(c => COG_GROUP[c] = "poor");
  function cogGroup(cat) { return cat ? (COG_GROUP[cat[0]] || "none") : "none"; }
  function finishAnnotation(reps, feats, source) {
    const ridx = new Map(reps.map((r, i) => [r.id, i]));
    const warnings = [];
    for (const r of reps) if (!r.length) { const mx = Math.max(0, ...feats.filter(f => f.seq === r.id).map(f => f.e)); r.length = mx; warnings.push(`replicon ${r.id}: length not given in file; using last feature end (${mx} bp)`); }
    feats.forEach(f => { f.rep = ridx.get(f.seq); f.group = cogGroup(f.cogcat); });
    feats.sort((a, b) => a.rep - b.rep || a.s - b.s || a.e - b.e);
    const usedReps = new Set(feats.map(f => f.rep));
    return { source, replicons: reps, features: feats, warnings, nRepliconsWithFeatures: usedReps.size };
  }

  // ------------------------------------------------------------------ count parsers
  /** → {genes:[id], samples:[name], counts:[[sampleCol] of Float64Array], meta:{...}} */
  function parseFeatureCounts(text) {
    const tab = parseTable(text);
    const H = tab.header; if (H[0] !== "Geneid") throw new Error("featureCounts: first column must be 'Geneid'");
    const sIdx = []; for (let j = 6; j < H.length; j++) sIdx.push(j);
    const samples = sIdx.map(j => basenameSample(H[j]));
    const genes = tab.rows.map(r => r[0]);
    const counts = sIdx.map(j => Float64Array.from(tab.rows, r => +r[j]));
    const cmd = (lines(text).find(l => /^#\s*Program:featureCounts/.test(l)) || "").slice(0, 2000);
    return { genes, samples, counts, meta: { format: "featureCounts", command: cmd } };
  }
  function parseMatrix(text) {
    const tab = parseTable(text);
    const samples = tab.header.slice(1).map(basenameSample);
    const genes = tab.rows.map(r => r[0]);
    const counts = samples.map((_, j) => Float64Array.from(tab.rows, r => r[j + 1] === "" ? NaN : +r[j + 1]));
    for (let j = 0; j < counts.length; j++) for (let i = 0; i < Math.min(counts[j].length, 5000); i++) {
      const v = counts[j][i]; if (!Number.isNaN(v) && (v < 0 || Math.abs(v - Math.round(v)) > 1e-6)) throw new Error(`count matrix column '${samples[j]}' contains non-integer values (e.g. ${v}); provide raw read counts, not TPM/FPKM/normalised values`);
    }
    return { genes, samples, counts, meta: { format: "matrix" } };
  }
  function parseHTSeqFiles(files) { // files: [{name, text}]
    const all = new Map(); const samples = [];
    files.forEach((f, j) => { let nm = f.name; while (/\.(gz|txt|tsv|csv|counts?|htseq|htseq-count)$/i.test(nm)) nm = nm.replace(/\.(gz|txt|tsv|csv|counts?|htseq|htseq-count)$/i, ""); samples.push(basenameSample(nm)); for (const l of lines(f.text)) { if (!l || l.startsWith("__")) continue; const [g, v] = l.split("\t"); if (!all.has(g)) all.set(g, new Array(files.length).fill(0)); all.get(g)[j] = +v; } });
    const genes = [...all.keys()]; const counts = samples.map((_, j) => Float64Array.from(genes, g => all.get(g)[j]));
    return { genes, samples, counts, meta: { format: "HTSeq-count" } };
  }
  /** DE table → {genes, log2fc, padj, basemean, columns} */
  function parseDETable(text) {
    const tab = parseTable(text); let H = tab.header; let rows = tab.rows;
    let idc = findCol(H, ["gene_id", "geneid", "locus_tag", "gene", "id", "feature_id"]);
    if (idc < 0) { idc = 0; if (H[0] === "" || rows.length && rows[0].length === H.length + 1) { /* R row names */ if (rows.length && rows[0].length === H.length + 1) H = [""].concat(H); } }
    const fcc = findCol(H, ["log2FoldChange", "logFC", "log2FC", "log2_fold_change"]);
    const pc = findCol(H, ["padj", "FDR", "adj.P.Val", "qvalue", "q_value", "p_adj"]);
    const bc = findCol(H, ["baseMean", "logCPM", "AveExpr"]);
    if (fcc < 0) throw new Error("DE table: no log2 fold-change column (log2FoldChange / logFC) found");
    const num = v => (v === "" || v === "NA" || v === "NaN" || v === "nan") ? null : +v;
    return { genes: rows.map(r => r[idc]), log2fc: rows.map(r => num(r[fcc])), padj: pc >= 0 ? rows.map(r => num(r[pc])) : null, basemean: bc >= 0 ? rows.map(r => num(r[bc])) : null,
      columns: { id: H[idc], log2fc: H[fcc], padj: pc >= 0 ? H[pc] : null, abundance: bc >= 0 ? H[bc] : null } };
  }

  // ------------------------------------------------------------------ design tables
  function rowsToObjects(tab) { return tab.rows.map(r => Object.fromEntries(tab.header.map((h, i) => [h.trim(), (r[i] === undefined ? "" : r[i])]))); }
  function parseSampleSheet(text) {
    const o = rowsToObjects(parseTable(text));
    if (!o.length || !("sample_id" in o[0])) throw new Error("samplesheet: column 'sample_id' is required");
    return o.map(r => ({ id: r.sample_id, label: r.label || r.sample_id, condition: r.condition || "", replicate: r.replicate || "", batch: r.batch || "", count_column: r.count_column || "" }));
  }
  function parseContrasts(text) {
    const o = rowsToObjects(parseTable(text));
    if (!o.length || !("contrast_id" in o[0])) throw new Error("contrasts: column 'contrast_id' is required");
    return o.map(r => ({ id: r.contrast_id, label: r.label || r.contrast_id, short: r.short_label || r.label || r.contrast_id, treatment: r.treatment || "", control: r.control || "", role: (lc(r.role) === "reference" ? "reference" : "test"), de_file: r.de_file || "" }));
  }
  function parseLandmarks(text) {
    return rowsToObjects(parseTable(text)).map(r => ({ name: r.name, replicon: r.replicon || "", s: r.start ? +r.start : null, e: r.end ? +r.end : null, feature_id: r.feature_id || "", type: r.type || "", color: r.color || "", label: r.label || "" }));
  }
  function parseLabels(text) { const out = {}; for (const r of rowsToObjects(parseTable(text))) if (r.feature_id) out[r.feature_id] = { label: r.label || "", note: r.note || "" }; return out; }
  /** group spec: "id1;id2" or "condition=AER&batch=EXP3" → array of sample ids */
  function resolveGroup(spec, samples) {
    spec = String(spec || "").trim(); if (!spec) return [];
    if (/=/.test(spec)) {
      const crit = spec.split("&").map(x => x.split("=").map(y => y.trim()));
      return samples.filter(s => crit.every(([k, v]) => String(s[k] === undefined ? "" : s[k]) === v)).map(s => s.id);
    }
    return spec.split(/[;,]/).map(x => x.trim()).filter(Boolean);
  }

  // ------------------------------------------------------------------ Mode A computation
  function median(arr) { const a = arr.filter(v => v !== null && !Number.isNaN(v)).sort((x, y) => x - y); const n = a.length; if (!n) return null; return n % 2 ? a[(n - 1) >> 1] : (a[n / 2 - 1] + a[n / 2]) / 2; }
  /** counts: array (libraries) of arrays (rows). rowMask: rows used for size factors. */
  function sizeFactors(counts, rowMask, variant) {
    const m = counts.length, n = counts[0].length; const ratios = counts.map(() => []);
    for (let i = 0; i < n; i++) {
      if (rowMask && !rowMask[i]) continue;
      let sum = 0, k = 0, allPos = true;
      for (let j = 0; j < m; j++) { const v = counts[j][i]; if (v > 0) { sum += Math.log(v); k++; } else allPos = false; }
      if (variant === "standard") { if (!allPos) continue; }
      else if (k === 0) continue;
      const lg = sum / k;
      for (let j = 0; j < m; j++) { const v = counts[j][i]; if (v > 0) ratios[j].push(Math.log(v) - lg); }
    }
    return ratios.map(r => Math.exp(median(r)));
  }
  /** Mode A: annotation + counts + samplesheet + contrasts → partial dataset (values only). */
  function computeModeA(ann, cnt, samples, contrasts, opts) {
    opts = Object.assign({ sfVariant: "standard", sfTypes: ["CDS"], showTypes: ["CDS", "ncRNA", "tmRNA"], cpm: 1, minLibs: 2, pseudo: 1 }, opts || {});
    const warnings = [];
    // map sample sheet to count columns
    const colIdx = new Map(cnt.samples.map((s, j) => [s, j]));
    const libs = samples.map(s => { const key = s.count_column || s.id; if (!colIdx.has(key)) throw new Error(`sample '${s.id}' not found among count columns (${cnt.samples.join(", ")})`); return Object.assign({}, s, { col: colIdx.get(key) }); });
    const unused = cnt.samples.filter(s => !libs.some(l => cnt.samples[l.col] === s)); if (unused.length) warnings.push(`count columns not in samplesheet (ignored): ${unused.join(", ")}`);
    // feature ↔ count rows
    const rowOf = new Map(cnt.genes.map((g, i) => [g, i]));
    const showT = new Set(opts.showTypes), sfT = new Set(opts.sfTypes);
    const feats = ann.features.filter(f => showT.has(f.type));
    const annIds = new Set(ann.features.map(f => f.id));
    const matched = feats.filter(f => rowOf.has(f.id)).length;
    const report = { annotationFeatures: ann.features.length, shownFeatures: feats.length, matchedShown: matched, countRows: cnt.genes.length, countRowsInAnnotation: cnt.genes.filter(g => annIds.has(g)).length };
    report.matchRate = feats.length ? matched / feats.length : 0;
    if (report.matchRate < 0.9) throw new Error(`only ${(100 * report.matchRate).toFixed(1)}% of annotation features (${matched}/${feats.length}) have a row in the count table — check that the count table uses the same ids (locus_tag) as the annotation`);
    if (matched < feats.length) warnings.push(`${feats.length - matched} annotation features have no count row; shown as missing`);
    // size-factor rows: count rows whose annotation type is in sfTypes
    const typeOf = new Map(ann.features.map(f => [f.id, f.type]));
    const sfMask = cnt.genes.map(g => sfT.has(typeOf.get(g)));
    const C = libs.map(l => cnt.counts[l.col]);
    const sf = sizeFactors(C, sfMask, opts.sfVariant);
    const libTot = C.map(col => { let t = 0; for (let i = 0; i < col.length; i++) if (sfMask[i]) t += col[i]; return t; }); // CPM vs CDS totals
    const nF = feats.length;
    const normalized = libs.map(() => new Array(nF).fill(null));
    const expressed = new Array(nF).fill(false);
    for (let i = 0; i < nF; i++) {
      const r = rowOf.get(feats[i].id); if (r === undefined) continue;
      let k = 0;
      for (let j = 0; j < libs.length; j++) { const v = C[j][r]; normalized[j][i] = v / sf[j]; if (v / libTot[j] * 1e6 > opts.cpm) k++; }
      expressed[i] = k >= opts.minLibs;
    }
    const libIndex = new Map(libs.map((l, j) => [l.id, j]));
    const outC = [], log2fc = [];
    for (const c of contrasts) {
      const T = resolveGroup(c.treatment, libs), K = resolveGroup(c.control, libs);
      for (const x of T.concat(K)) if (!libIndex.has(x)) throw new Error(`contrast '${c.id}': sample '${x}' is not in the samplesheet`);
      if (!T.length || !K.length) throw new Error(`contrast '${c.id}': treatment or control group is empty`);
      const ti = T.map(x => libIndex.get(x)), ki = K.map(x => libIndex.get(x));
      const v = new Array(nF);
      for (let i = 0; i < nF; i++) {
        if (normalized[0][i] === null) { v[i] = null; continue; }
        let a = 0, b = 0; for (const j of ti) a += normalized[j][i]; for (const j of ki) b += normalized[j][i];
        v[i] = Math.log2((a / ti.length + opts.pseudo) / (b / ki.length + opts.pseudo));
      }
      log2fc.push(v);
      outC.push({ id: c.id, label: c.label, short: c.short, treatment: T, control: K, kind: "descriptive", role: c.role, nTreatment: T.length, nControl: K.length });
      if (T.length < 2 || K.length < 2) c._noRep = true;
    }
    const noRep = outC.filter(c => c.nTreatment < 2 || c.nControl < 2).map(c => c.id);
    if (noRep.length) warnings.push(`no replication in ${noRep.length} contrast(s) (${noRep.join(", ")}): values are single-library descriptive ratios, not tested differences`);
    return { features: feats, libraries: libs.map(l => ({ id: l.id, label: l.label, condition: l.condition, replicate: l.replicate, batch: l.batch })),
      normalized, expressed, contrasts: outC, log2fc, padj: outC.map(() => null), basemean: outC.map(() => null),
      sizeFactors: Object.fromEntries(libs.map((l, j) => [l.id, sf[j]])), libTotalsSF: Object.fromEntries(libs.map((l, j) => [l.id, libTot[j]])),
      report, warnings,
      meta: { mode: "A", size_factor_variant: opts.sfVariant,
        normalization: `median-of-ratios size factors (${opts.sfVariant === "standard" ? "DESeq2 standard: genes with non-zero counts in all libraries" : "zero-tolerant: zeros treated as missing"}) on ${[...sfT].join("/")} rows; log2FC = log2((mean normalised treatment + ${opts.pseudo}) / (mean normalised control + ${opts.pseudo}))`,
        filter: `CPM > ${opts.cpm} in ≥ ${opts.minLibs} of ${libs.length} libraries (CPM relative to ${[...sfT].join("/")} library totals)`, pseudocount: opts.pseudo } };
  }
  // ------------------------------------------------------------------ Mode B import
  function computeModeB(ann, deTables, contrasts, opts) {
    opts = Object.assign({ showTypes: ["CDS", "ncRNA", "tmRNA"], alpha: 0.05 }, opts || {});
    const showT = new Set(opts.showTypes); const feats = ann.features.filter(f => showT.has(f.type));
    const idx = new Map(feats.map((f, i) => [f.id, i])); const warnings = [];
    const outC = [], log2fc = [], padj = [], basemean = []; const tested = new Array(feats.length).fill(false);
    const report = { shownFeatures: feats.length, perContrast: {} };
    for (const c of contrasts) {
      const t = deTables[c.id]; if (!t) throw new Error(`contrast '${c.id}': DE table '${c.de_file}' not loaded`);
      const fc = new Array(feats.length).fill(null), pa = t.padj ? new Array(feats.length).fill(null) : null, bm = t.basemean ? new Array(feats.length).fill(null) : null;
      let hit = 0;
      t.genes.forEach((g, r) => { const i = idx.get(g); if (i === undefined) return; hit++; fc[i] = t.log2fc[r]; if (pa) pa[i] = t.padj[r]; if (bm) bm[i] = t.basemean[r]; if (t.log2fc[r] !== null) tested[i] = true; });
      const rate = feats.length ? hit / feats.length : 0; report.perContrast[c.id] = { rows: t.genes.length, matched: hit, matchRate: rate, columns: t.columns };
      if (hit / Math.max(1, t.genes.length) < 0.5) throw new Error(`contrast '${c.id}': only ${hit} of ${t.genes.length} DE-table ids match annotation ids`);
      if (!t.padj) warnings.push(`contrast '${c.id}': no adjusted p-value column; significance overlay disabled`);
      log2fc.push(fc); padj.push(pa); basemean.push(bm);
      outC.push({ id: c.id, label: c.label, short: c.short, treatment: [], control: [], kind: "statistical", role: c.role, abundance: t.columns.abundance || null, source: `${c.de_file} (${t.columns.log2fc}${t.columns.padj ? ", " + t.columns.padj : ""})` });
    }
    return { features: feats, libraries: [], normalized: null, expressed: tested, contrasts: outC, log2fc, padj, basemean, report, warnings,
      meta: { mode: "B", size_factor_variant: "none", normalization: "imported from DE result tables (computed outside the viewer)", filter: "genes present with a non-missing log2FC in ≥1 table", pseudocount: null, alpha: opts.alpha } };
  }

  // ------------------------------------------------------------------ assemble dataset
  function resolveLandmarks(lm, ann, feats) {
    const rid = new Map(ann.replicons.map((r, i) => [r.id, i])); const fid = new Map(feats.map(f => [f.id, f])); const out = []; const warn = [];
    for (const L of lm || []) {
      if (L.feature_id && fid.has(L.feature_id)) { const f = fid.get(L.feature_id); out.push({ name: L.name || L.feature_id, rep: f.rep, s: f.s, e: f.e, type: L.type || "gene", color: L.color || "", label: L.label || L.name || L.feature_id }); continue; }
      const r = rid.has(L.replicon) ? rid.get(L.replicon) : (ann.replicons.length === 1 ? 0 : -1);
      if (r < 0 || !(L.s > 0) || !(L.e > 0)) { warn.push(`landmark '${L.name}' skipped (unknown replicon or coordinates)`); continue; }
      out.push({ name: L.name, rep: r, s: Math.min(L.s, L.e), e: Math.max(L.s, L.e), type: L.type || "", color: L.color || "", label: L.label || L.name });
    }
    return { landmarks: out, warnings: warn };
  }
  function buildDataset(parts) {
    const { annotation: ann, values: V, landmarks, labels, title, organism, strain, inputs } = parts;
    const F = V.features; const lmR = resolveLandmarks(landmarks, ann, F);
    const ds = {
      format: "genomemap-dataset", version: 1, title: title || "", organism: organism || "", strain: strain || "",
      replicons: ann.replicons.map(r => ({ id: r.id, name: r.name || r.id, length: r.length, topology: r.topology || "linear" })),
      features: { id: F.map(f => f.id), rep: F.map(f => f.rep), s: F.map(f => f.s), e: F.map(f => f.e), st: F.map(f => f.st), type: F.map(f => f.type),
        gene: F.map(f => f.gene || ""), product: F.map(f => f.product || ""), cog: F.map(f => f.cog || ""), cogcat: F.map(f => f.cogcat || ""), group: F.map(f => f.group || "none") },
      libraries: V.libraries, normalized: V.normalized, expressed: V.expressed, contrasts: V.contrasts,
      log2fc: V.log2fc, padj: V.padj, basemean: V.basemean, landmarks: lmR.landmarks, labels: labels || {},
      meta: Object.assign({}, V.meta, { created: new Date().toISOString(), software: "genomemap " + VERSION, annotation_source: ann.source,
        size_factors: V.sizeFactors || null, inputs: inputs || [], report: V.report, warnings: (ann.warnings || []).concat(V.warnings || [], lmR.warnings) })
    };
    return ds;
  }
  /** Round numeric arrays for compact embedding. */
  function compactDataset(ds, digits) {
    digits = digits || { fc: 4, norm: 1, p: 6 };
    const r = (v, d) => (v === null || v === undefined || Number.isNaN(v)) ? null : +v.toFixed(d);
    const rp = v => (v === null || v === undefined || Number.isNaN(v)) ? null : (v === 0 ? 0 : +v.toPrecision(3));
    const o = Object.assign({}, ds);
    o.log2fc = ds.log2fc.map(a => a.map(v => r(v, digits.fc)));
    if (ds.normalized) o.normalized = ds.normalized.map(a => a.map(v => r(v, digits.norm)));
    o.padj = ds.padj.map(a => a ? a.map(rp) : null);
    o.basemean = ds.basemean.map(a => a ? a.map(v => r(v, 1)) : null);
    return o;
  }

  // ------------------------------------------------------------------ bins (viewer-side)
  /** Per replicon bins of width bw. Returns {bw, reps:[{rep, n, start:[...], ng:[], nx:[], med:[[c][b]], fu:[[c][b]], fd:[[c][b]]}]} */
  function computeBins(ds, bw) {
    bw = bw || 10000; const F = ds.features, NC = ds.contrasts.length; const out = { bw, reps: [] };
    ds.replicons.forEach((R, ri) => {
      const nb = Math.max(1, Math.ceil(R.length / bw)); const members = Array.from({ length: nb }, () => []); const ng = new Array(nb).fill(0), nx = new Array(nb).fill(0);
      for (let i = 0; i < F.id.length; i++) { if (F.rep[i] !== ri) continue; const mid = (F.s[i] + F.e[i]) / 2; const b = Math.min(Math.floor((mid - 1) / bw), nb - 1); ng[b]++; if (ds.expressed[i]) { nx[b]++; members[b].push(i); } }
      const med = [], fu = [], fd = [];
      for (let c = 0; c < NC; c++) {
        const m = new Array(nb), u = new Array(nb), d = new Array(nb); const v = ds.log2fc[c];
        for (let b = 0; b < nb; b++) { const vals = members[b].map(i => v[i]).filter(x => x !== null); if (!vals.length) { m[b] = null; u[b] = null; d[b] = null; continue; }
          m[b] = median(vals); u[b] = vals.filter(x => x > 1).length / vals.length; d[b] = vals.filter(x => x < -1).length / vals.length; }
        med.push(m); fu.push(u); fd.push(d);
      }
      out.reps.push({ rep: ri, n: nb, ng, nx, med, fu, fd });
    });
    return out;
  }

  /** One-page PDF holding one image. o = {data: Uint8Array (zlib-compressed RGB rows for FlateDecode, or JPEG bytes for DCTDecode),
   *  filter: "FlateDecode" | "DCTDecode", imgW, imgH (pixels), pageW, pageH (points), title}. Returns Uint8Array. */
  function buildPdf(o) {
    const enc = s => { const b = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) b[i] = s.charCodeAt(i) & 255; return b; };
    const hex16 = t => "<FEFF" + Array.from(String(t || "")).map(ch => { const c = ch.codePointAt(0); const u = c > 0xFFFF ? [0xD800 + ((c - 0x10000) >> 10), 0xDC00 + ((c - 0x10000) & 1023)] : [c]; return u.map(x => x.toString(16).padStart(4, "0")).join(""); }).join("").toUpperCase() + ">";
    const f = v => (+v).toFixed(2).replace(/\.?0+$/, "");
    const d = new Date(), pad = n => String(n).padStart(2, "0");
    const date = `D:${d.getUTCFullYear()}${pad(d.getUTCMonth() + 1)}${pad(d.getUTCDate())}${pad(d.getUTCHours())}${pad(d.getUTCMinutes())}${pad(d.getUTCSeconds())}Z`;
    const content = `q ${f(o.pageW)} 0 0 ${f(o.pageH)} 0 0 cm /Im0 Do Q\n`;
    const parts = [], offs = []; let pos = 0;
    const push = b => { parts.push(b); pos += b.length; };
    const obj = (n, body) => { offs[n] = pos; push(enc(`${n} 0 obj\n${body}\nendobj\n`)); };
    push(new Uint8Array([0x25, 0x50, 0x44, 0x46, 0x2D, 0x31, 0x2E, 0x34, 0x0A, 0x25, 0xE2, 0xE3, 0xCF, 0xD3, 0x0A])); // %PDF-1.4 + binary marker
    obj(1, "<< /Type /Catalog /Pages 2 0 R >>");
    obj(2, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    obj(3, `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${f(o.pageW)} ${f(o.pageH)}] /Resources << /XObject << /Im0 5 0 R >> >> /Contents 4 0 R >>`);
    obj(4, `<< /Length ${content.length} >>\nstream\n${content}endstream`);
    offs[5] = pos; push(enc(`5 0 obj\n<< /Type /XObject /Subtype /Image /Width ${o.imgW} /Height ${o.imgH} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /${o.filter} /Length ${o.data.length} >>\nstream\n`));
    push(o.data); push(enc("\nendstream\nendobj\n"));
    obj(6, `<< /Title ${hex16(o.title)} /Producer (genomemap ${VERSION}) /CreationDate (${date}) >>`);
    const xref = pos; let x = "xref\n0 7\n0000000000 65535 f \n"; for (let n = 1; n <= 6; n++) x += String(offs[n]).padStart(10, "0") + " 00000 n \n";
    push(enc(x + `trailer\n<< /Size 7 /Root 1 0 R /Info 6 0 R >>\nstartxref\n${xref}\n%%EOF\n`));
    const out = new Uint8Array(pos); let k = 0; for (const b of parts) { out.set(b, k); k += b.length; } return out;
  }

  return { VERSION, parseTable, detectType, parseGFF3, parseGenBank, parseFeatureCounts, parseMatrix, parseHTSeqFiles, parseDETable,
    parseSampleSheet, parseContrasts, parseLandmarks, parseLabels, resolveGroup, sizeFactors, computeModeA, computeModeB, buildDataset,
    compactDataset, computeBins, cogGroup, median, basenameSample, buildPdf };
}));
