/* genomemap loader — reads user files in the browser, detects their type, builds the dataset with GenomeMapCore, then calls GenomeMapViewer.start(ds). */
(function () {
  "use strict";
  const C = window.GenomeMapCore;
  const $ = id => document.getElementById(id);
  // pristine copy of this page (loader.js is the last script, so the whole document is already parsed); used by "Save report"
  const PRISTINE = "<!DOCTYPE html>\n" + document.documentElement.outerHTML;
  const TYPES = ["gff3", "genbank", "featurecounts", "matrix", "htseq", "de_table", "samplesheet", "contrasts", "landmarks", "labels", "dataset_json", "ignore"];
  const TYPE_LABEL = { gff3: "annotation (GFF3)", genbank: "annotation (GenBank)", featurecounts: "counts (featureCounts)", matrix: "counts (matrix)", htseq: "counts (HTSeq, one sample)", de_table: "DE result table (Mode B)",
    samplesheet: "sample sheet", contrasts: "contrasts", landmarks: "landmarks", labels: "locus labels", dataset_json: "saved genomemap dataset", unknown: "unknown — choose", ignore: "ignore" };
  const files = []; let design = { samples: null, contrasts: null, auto: false };
  function esc(s) { return String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }
  function msg(kind, text) { const d = document.createElement("div"); d.className = "msg " + kind; d.innerHTML = text; $("msgs").appendChild(d); }
  function clearMsgs() { $("msgs").innerHTML = ""; }
  async function readFile(f) {
    if (/\.gz$/i.test(f.name)) {
      if (typeof DecompressionStream === "undefined") throw new Error(`${f.name}: this browser cannot decompress .gz; please decompress the file first`);
      const ds = f.stream().pipeThrough(new DecompressionStream("gzip")); return await new Response(ds).text();
    }
    if (typeof f.text === "function") return await f.text();
    return await new Promise((res, rej) => { const r = new FileReader(); r.onload = () => res(r.result); r.onerror = () => rej(r.error); r.readAsText(f); });
  }
  async function addFiles(list) {
    clearMsgs();
    for (const f of list) {
      try {
        const text = await readFile(f); const name = f.name.replace(/\.gz$/i, "");
        const i = files.findIndex(x => x.name === name); const rec = { name, text, size: f.size, type: C.detectType(name, text) };
        if (i >= 0) files[i] = rec; else files.push(rec);
      } catch (e) { msg("err", esc(e.message)); }
    }
    refresh();
  }
  function byType(t) { return files.filter(f => f.type === t); }
  function refresh() {
    const tb = $("fileTable").querySelector("tbody");
    tb.innerHTML = files.map((f, i) => `<tr><td>${esc(f.name)}</td><td>${(f.size / 1024).toFixed(0)} kB</td><td><select data-i="${i}">${TYPES.concat(f.type === "unknown" ? ["unknown"] : []).map(t => `<option value="${t}" ${t === f.type ? "selected" : ""}>${TYPE_LABEL[t]}</option>`).join("")}</select></td><td><a href="#" data-del="${i}">remove</a></td></tr>`).join("");
    tb.querySelectorAll("select").forEach(s => s.onchange = () => { files[+s.dataset.i].type = s.value; refresh(); });
    tb.querySelectorAll("a[data-del]").forEach(a => a.onclick = ev => { ev.preventDefault(); files.splice(+a.dataset.del, 1); refresh(); });
    const saved = byType("dataset_json");
    const ann = byType("gff3").concat(byType("genbank"));
    const modeA = byType("featurecounts").length + byType("matrix").length + byType("htseq").length > 0;
    const modeB = byType("de_table").length > 0;
    $("optPanel").style.display = ann.length || saved.length ? "block" : "none";
    document.querySelectorAll(".modeA").forEach(e => e.style.display = modeA ? "" : "none");
    document.querySelectorAll(".modeB").forEach(e => e.style.display = modeB ? "" : "none");
    if (ann.length && !$("optTitle").value) { try { const a = parseAnnotation(ann[0]); $("optTitle").value = a.organism || ann[0].name; fillTypes(a); } catch (e) { } }
    else if (ann.length && !$("optTypes").innerHTML) { try { fillTypes(parseAnnotation(ann[0])); } catch (e) { } }
    updateDesign(modeA, modeB);
    let ready = false, why = "";
    if (saved.length) { ready = true; why = "a saved dataset will be opened"; }
    else if (!ann.length) why = "an annotation (GFF3 or GenBank) is required";
    else if (ann.length > 1) why = "load only one annotation file";
    else if (modeA && modeB) why = "load either counts (Mode A) or DE tables (Mode B), not both";
    else if (!modeA && !modeB) why = "add counts (Mode A) or DE result tables (Mode B)";
    else if (!(design.contrasts && design.contrasts.length)) why = "define at least one contrast";
    else { ready = true; why = modeA ? "Mode A: values computed here from counts" : "Mode B: values imported from DE tables"; }
    $("buildBtn").disabled = !ready; $("buildMsg").textContent = why;
  }
  let annCache = null;
  function parseAnnotation(f) { if (annCache && annCache.f === f && annCache.types === typesKey()) return annCache.a; const opt = { types: ["CDS", "ncRNA", "tmRNA", "tRNA", "rRNA"] }; const a = f.type === "gff3" ? C.parseGFF3(f.text, opt) : C.parseGenBank(f.text, opt); annCache = { f, a, types: typesKey() }; return a; }
  function typesKey() { return "all"; }
  function fillTypes(a) {
    const n = {}; a.features.forEach(f => n[f.type] = (n[f.type] || 0) + 1);
    const def = new Set(["CDS", "ncRNA", "tmRNA"]);
    $("optTypes").innerHTML = Object.keys(n).map(t => `<label style="margin-right:8px"><input type="checkbox" class="ftype" value="${t}" ${def.has(t) ? "checked" : ""}> ${t} (${n[t]})</label>`).join("") +
      `<span class="small"> — tRNA/rRNA are off by default: rRNA-depleted libraries make their counts unreliable</span>`;
  }
  function countsObject() {
    if (byType("featurecounts").length) return C.parseFeatureCounts(byType("featurecounts")[0].text);
    if (byType("matrix").length) return C.parseMatrix(byType("matrix")[0].text);
    if (byType("htseq").length) return C.parseHTSeqFiles(byType("htseq"));
    return null;
  }
  function updateDesign(modeA, modeB) {
    const show = modeA || modeB; $("designPanel").style.display = show ? "block" : "none"; if (!show) return;
    const ssF = byType("samplesheet")[0], ctF = byType("contrasts")[0];
    let info = [];
    try {
      if (modeA) {
        const cnt = countsObject();
        if (ssF) { design.samples = C.parseSampleSheet(ssF.text); info.push(`sample sheet: ${design.samples.length} libraries`); }
        else { design.samples = cnt.samples.map(s => ({ id: s, label: s, condition: s, replicate: "", batch: "", count_column: "" })); info.push(`no sample sheet: ${cnt.samples.length} count columns used as libraries (condition = column name)`); }
        const opts = design.samples.map(s => `<option value="${esc(s.id)}">${esc(s.label)}</option>`).join("");
        $("peT").innerHTML = opts; $("peC").innerHTML = opts;
      } else { design.samples = []; }
      if (ctF) { design.contrasts = C.parseContrasts(ctF.text); design.auto = false; info.push(`contrasts file: ${design.contrasts.length} contrasts`); }
      else if (modeB) { design.contrasts = byType("de_table").map(f => { const id = f.name.replace(/\.(csv|tsv|txt)$/i, ""); return { id, label: id, short: id.slice(0, 16), treatment: "", control: "", role: "test", de_file: f.name }; }); design.auto = true; info.push("no contrasts file: one contrast per DE table (label = file name)"); }
      else if (!design.contrasts || !design.auto) { design.contrasts = design.contrasts && design.auto ? design.contrasts : []; design.auto = true; }
    } catch (e) { info.push(`<span style="color:#9C4B22">${esc(e.message)}</span>`); }
    $("pairEditor").style.display = (modeA && !ctF) ? "block" : "none";
    $("designInfo").innerHTML = info.join(" · ");
    const tb = $("contrastTable").querySelector("tbody");
    tb.innerHTML = (design.contrasts || []).map((c, i) => `<tr><td>${esc(c.id)}</td><td>${esc(c.label)}</td><td>${esc(c.treatment || "—")}</td><td>${esc(c.control || "—")}</td><td>${c.role}</td><td>${design.auto && modeA ? `<a href="#" data-cdel="${i}">remove</a>` : (c.de_file ? esc(c.de_file) : "")}</td></tr>`).join("");
    tb.querySelectorAll("a[data-cdel]").forEach(a => a.onclick = ev => { ev.preventDefault(); design.contrasts.splice(+a.dataset.cdel, 1); refresh(); });
  }
  $("peAdd").onclick = () => {
    const t = $("peT").value, c = $("peC").value; if (!t || !c || t === c) { msg("warn", "choose two different libraries"); return; }
    const S = id => design.samples.find(s => s.id === id).label;
    const label = $("peLabel").value.trim() || `${S(t)} vs ${S(c)}`;
    let id = label.replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "").slice(0, 40) || ("c" + (design.contrasts.length + 1));
    while (design.contrasts.some(x => x.id === id)) id += "_2";
    design.contrasts.push({ id, label, short: label.slice(0, 16), treatment: t, control: c, role: $("peRole").value, de_file: "" }); design.auto = true; $("peLabel").value = ""; refresh();
  };
  $("dlSheets").onclick = ev => {
    ev.preventDefault();
    const ss = ["sample_id\tlabel\tcondition\treplicate\tbatch\tcount_column"].concat((design.samples || []).map(s => [s.id, s.label, s.condition, s.replicate, s.batch, s.count_column].join("\t"))).join("\n") + "\n";
    const ct = ["contrast_id\tlabel\tshort_label\ttreatment\tcontrol\trole\tde_file"].concat((design.contrasts || []).map(c => [c.id, c.label, c.short, c.treatment, c.control, c.role, c.de_file].join("\t"))).join("\n") + "\n";
    download("samplesheet.tsv", ss, "text/tab-separated-values"); setTimeout(() => download("contrasts.tsv", ct, "text/tab-separated-values"), 300);
  };
  function download(name, text, type) { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; document.body.appendChild(a); a.click(); a.remove(); }
  window.GenomeMapDownload = download;
  async function sha256(text) { if (!(window.crypto && crypto.subtle)) return ""; const h = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text)); return [...new Uint8Array(h)].map(b => b.toString(16).padStart(2, "0")).join(""); }
  $("buildBtn").onclick = async () => {
    clearMsgs();
    try {
      const saved = byType("dataset_json")[0];
      if (saved) { const ds = JSON.parse(saved.text); msg("ok", `opened saved dataset “${esc(ds.title || saved.name)}” (${ds.features.id.length} features, ${ds.contrasts.length} contrasts)`); return go(ds, +$("optBin").value); }
      const annF = byType("gff3").concat(byType("genbank"))[0]; const ann = parseAnnotation(annF);
      const showTypes = [...document.querySelectorAll(".ftype:checked")].map(x => x.value);
      if (!showTypes.length) throw new Error("tick at least one feature type to show");
      let V;
      if (byType("de_table").length) {
        const tables = {}; for (const c of design.contrasts) { const f = files.find(x => x.name === c.de_file || x.name.replace(/\.(csv|tsv|txt)$/i, "") === c.de_file.replace(/\.(csv|tsv|txt)$/i, "")); if (!f) throw new Error(`contrast '${c.id}': DE file '${c.de_file}' not loaded`); tables[c.id] = C.parseDETable(f.text); }
        V = C.computeModeB(ann, tables, design.contrasts, { showTypes, alpha: +$("optAlpha").value });
      } else {
        V = C.computeModeA(ann, countsObject(), design.samples, design.contrasts, { sfVariant: $("optSF").value, cpm: +$("optCPM").value, minLibs: +$("optMinLibs").value, showTypes });
      }
      const lm = byType("landmarks").length ? C.parseLandmarks(byType("landmarks")[0].text) : [];
      const labels = byType("labels").length ? C.parseLabels(byType("labels")[0].text) : {};
      const inputs = []; for (const f of files) if (f.type !== "ignore") inputs.push({ file: f.name, type: f.type, bytes: f.text.length, sha256: await sha256(f.text) });
      const ds = C.buildDataset({ annotation: ann, values: V, landmarks: lm, labels, title: $("optTitle").value, organism: ann.organism || "", inputs });
      if (ds.meta.mode === "B") ds.meta.alpha = +$("optAlpha").value;
      for (const w of ds.meta.warnings) msg("warn", esc(w));
      const R = V.report; msg("ok", ds.meta.mode === "A" ? `annotation: ${R.annotationFeatures} features, ${R.shownFeatures} shown; ${R.matchedShown} matched to count rows (${(100 * R.matchRate).toFixed(1)}%).` :
        `annotation: ${R.shownFeatures} features shown; ` + Object.entries(R.perContrast).map(([k, v]) => `${esc(k)}: ${v.matched}/${v.rows} table rows matched`).join("; "));
      go(ds, +$("optBin").value);
    } catch (e) { msg("err", esc(e.message)); console.error(e); }
  };
  function go(ds, bw) { $("loader").style.display = "none"; window.GenomeMapViewer.start(ds, { binWidth: bw }); }
  // drag & drop
  const drop = $("drop");
  ["dragenter", "dragover"].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add("over"); }));
  ["dragleave", "drop"].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove("over"); }));
  drop.addEventListener("drop", e => addFiles([...e.dataTransfer.files]));
  $("fileIn").onchange = e => addFiles([...e.target.files]);
  // start: embedded dataset (report built by the CLI) or loader
  window.addEventListener("DOMContentLoaded", () => {
    const emb = document.getElementById("gm-dataset");
    let skip = false; try { skip = sessionStorage.getItem("gm-skip-embedded") === "1"; sessionStorage.removeItem("gm-skip-embedded"); } catch (e) { /* storage unavailable */ }
    if (emb && emb.textContent.trim() && !skip) {
      const ds = JSON.parse(emb.textContent); const bw = ds.meta && ds.meta.bin_width ? ds.meta.bin_width : 10000;
      $("loader").style.display = "none"; window.GenomeMapViewer.start(ds, { binWidth: bw });
    } else $("loader").style.display = "block";
  });
  /** Write a self-contained HTML report: this viewer with the dataset embedded (opens directly in the map). */
  function saveReport(ds, name) {
    const json = JSON.stringify(C.compactDataset(ds)).replace(/<\//g, "<\\/");
    const d = new DOMParser().parseFromString(PRISTINE, "text/html"); // scripts are not executed by DOMParser
    const old = d.getElementById("gm-dataset"); if (old) old.remove();
    const tag = d.createElement("script"); tag.id = "gm-dataset"; tag.type = "application/json"; tag.textContent = json;
    const first = d.body.querySelector("script"); d.body.insertBefore(tag, first);
    download(name, "<!DOCTYPE html>\n" + d.documentElement.outerHTML, "text/html");
  }
  window.GenomeMapSaveReport = saveReport;
  window.GenomeMapLoader = { addFiles, files };
})();
