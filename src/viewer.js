/* genomemap viewer — circular + linear genome map of per-gene log2 fold changes.
 * Generalised from the 2026 R10.79 viewer: any number of replicons, contrasts, libraries; Mode A (descriptive) or Mode B (DE tables).
 * Entry point: GenomeMapViewer.start(dataset, {binWidth}) */
(function () {
"use strict";
function start(ds, opt) {
opt = opt || {};
const C = window.GenomeMapCore;
const $ = id => document.getElementById(id);
$("viewer").style.display = "block";
const F = ds.features, N = F.id.length, REP = ds.replicons, CT = ds.contrasts, NC = CT.length;
const FCV = ds.log2fc, PADJ = ds.padj || CT.map(() => null), BM = ds.basemean || CT.map(() => null), anyBM = BM.some(b => b), NORM = ds.normalized, X = ds.expressed, LIBS = ds.libraries || [];
const LAB = ds.labels || {}; const MODE = (ds.meta && ds.meta.mode) || "A"; const ALPHA = (ds.meta && ds.meta.alpha) || 0.05;
const isRef = c => CT[c].role === "reference";
const TEST = [], REFS = []; for (let c = 0; c < NC; c++) (isRef(c) ? REFS : TEST).push(c);
const hasCOG = F.group && F.group.some(gp => gp && gp !== "none");
const anyPadj = PADJ.some(p => p);
const BW = opt.binWidth || 10000; const BINS = C.computeBins(ds, BW);
// feature index ranges per replicon (features sorted by rep, start)
const repLo = REP.map(() => -1), repHi = REP.map(() => -2);
for (let i = 0; i < N; i++) { const r = F.rep[i]; if (repLo[r] < 0) repLo[r] = i; repHi[r] = i; }
const repsWithFeat = REP.map((_, r) => r).filter(r => repLo[r] >= 0);
// ---------- texts ----------
const esc = s => String(s === undefined || s === null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmt = d3.format(","); const f2 = v => v === null || v === undefined ? "n.a." : (v >= 0 ? "+" : "") + v.toFixed(2);
const fmtP = p => p === null || p === undefined ? "n.a." : (p < 1e-3 ? p.toExponential(1) : p.toFixed(3));
const totLen = REP.reduce((a, r) => a + r.length, 0);
$("vTitle").innerHTML = esc(ds.title || ds.organism || "genome map") + " — genome map of transcriptional responses";
const typesShown = [...new Set(F.type)].join(", ");
$("vSub").innerHTML = `${REP.length} replicon${REP.length > 1 ? "s" : ""} (${REP.map(r => `${esc(r.name || r.id)}${r.name && r.name !== r.id ? " [" + esc(r.id) + "]" : ""}, ${fmt(r.length)} bp${r.topology === "circular" ? ", circular" : ""}`).join("; ")}). ` +
  `Features shown: ${esc(typesShown)} (n = ${fmt(N)}). ${MODE === "A" ? "Values computed from counts (descriptive)." : "Values imported from DE result tables."} Colour range ±<span class="vmTxt">3</span> log2, values beyond are clipped.`;
const noRep = CT.filter(c => c.kind === "descriptive" && (c.treatment.length < 2 || c.control.length < 2));
let bh = "";
if (MODE === "A" && noRep.length) bh += `<div class="note">${noRep.length === CT.length ? "single-library descriptive ratios — no replication in any contrast" : `${noRep.length} of ${CT.length} contrasts are single-library descriptive ratios (no replication)`}</div>`;
if (MODE === "A" && noRep.length < CT.length) bh += `<div class="note" style="font-weight:400">contrasts with replicates: ratio of mean normalised counts — descriptive, not a statistical test</div>`;
if (MODE === "B") bh += `<div class="note" style="font-weight:400">log2FC and padj imported from DE tables; padj threshold ${ALPHA}</div>`;
$("banner").innerHTML = bh;
$("filtWord").textContent = MODE === "A" ? "genes failing expression filter" : "genes not tested";
if (anyPadj) { $("sigCtl").style.display = ""; $("alphaTxt").textContent = ALPHA; }
// ---------- colour ----------
const STOPS = ["#2F4682", "#3D5A9E", "#5B7FC4", "#A9BEE0", "#F2F2F7", "#F7DCC6", "#F0B48A", "#C86A3C", "#9C4B22"];
let VM = 3; let col;
function setVM(v) { VM = v; col = d3.scaleLinear().domain(STOPS.map((_, i) => -VM + 2 * VM * i / 8)).range(STOPS).interpolate(d3.interpolateRgb).clamp(true);
  document.querySelectorAll(".vmTxt").forEach(e => e.textContent = VM);
  const tk = [-VM, -2 * VM / 3, -VM / 3, 0, VM / 3, 2 * VM / 3, VM].map(v => +v.toFixed(2));
  $("lgTicks").innerHTML = tk.map((v, k) => `<span style="left:${(k / 6 * 100).toFixed(3)}%">${v > 0 ? "+" : (v < 0 ? "−" : "")}${Math.abs(v)}${k === 0 || k === 6 ? "*" : ""}</span>`).join(""); }
setVM(3);
const GREY = "#C4C4CC", INK = "#3A3A3A", INK2 = "#6B6B6B", NEUT = "#B8B8C0", GRID = "#DCDCE6";
const SG = { info: ["#3D5A9E", "Information storage & processing (J,A,K,L,B)", "COG: information"], cell: ["#C86A3C", "Cellular processes & signalling (D,Y,V,T,M,N,Z,W,U,O,X)", "COG: cellular"], metab: ["#A9BEE0", "Metabolism (C,G,E,F,H,I,P,Q)", "COG: metabolism"], poor: ["#F0B48A", "Poorly characterised (R,S)", "COG: poorly characterised"], none: ["#E6E6EC", "No COG assignment", "no COG"] };
const LM_PAL = ["#9C4B22", "#2F4682", "#6B6B6B", "#117A65", "#7D3C98", "#B7950B"]; const lmTypeCol = {};
const lmColor = L => L.color || (lmTypeCol[L.type || ""] = lmTypeCol[L.type || ""] || LM_PAL[Object.keys(lmTypeCol).length % LM_PAL.length]);
const st = { mode: "gene", ring: "all", showFilt: true, showNoise: false, sig: true, sel: TEST.slice(0, Math.min(4, TEST.length)), refSel: REFS.slice(0, 1),
  reps: repsWithFeat.slice(), win: { r: repsWithFeat[0], s: 1, e: 30000 }, hi: -1, dmode: "pan" };
if (!st.sel.length && REFS.length) { st.sel = [REFS[0]]; }
const mod = (a, m) => ((a % m) + m) % m;
// ---------- labels ----------
const symLoci = {}; for (let i = 0; i < N; i++) { if (F.gene[i]) (symLoci[F.gene[i]] = symLoci[F.gene[i]] || []).push(F.id[i]); }
const OVR = id => LAB[id] && LAB[id].label;
function baseLabel(i) { return OVR(F.id[i]) || F.gene[i] || F.id[i]; }
function isDup(i) { return F.gene[i] && symLoci[F.gene[i]].length > 1; }
function linLabel(i) { return OVR(F.id[i]) ? OVR(F.id[i]) : (isDup(i) ? `${F.gene[i]} (${F.id[i].slice(-5)})` : baseLabel(i)); }
const repName = r => REP[r].name || REP[r].id;
function sigAlpha(c, i) { const P = PADJ[c]; if (!P || !st.sig) return 1; const p = P[i]; return (p !== null && p < ALPHA) ? 1 : 0.3; }
function tipGene(i) {
  const L = LAB[F.id[i]];
  let h = `<b>${esc(F.id[i])}</b>${F.gene[i] ? ` · <i>${esc(F.gene[i])}</i>` : ""}${L && L.label ? ` · <i>${esc(L.label)}</i>${L.note ? ` (${esc(L.note)})` : ""}` : ""} <span class="small">(${esc(F.type[i])})</span><br>${esc(F.product[i] || "—")}<br>`;
  h += `<span class="small">${esc(repName(F.rep[i]))}: ${fmt(F.s[i])}–${fmt(F.e[i])} (${F.st[i] > 0 ? "+" : "−"}) · ${fmt(F.e[i] - F.s[i] + 1)} bp${hasCOG ? ` · COG ${esc(F.cog[i] || "—")} ${F.cogcat[i] ? "[" + esc(F.cogcat[i]) + "]" : ""}` : ""}</span><br>`;
  if (isDup(i)) { const Lc = symLoci[F.gene[i]]; h += `<span class="dup">symbol “${esc(F.gene[i])}” is assigned to ${Lc.length} loci (${esc(Lc.length > 6 ? Lc.slice(0, 6).join(", ") + ", …" : Lc.join(", "))}); rely on locus tag and product.</span><br>`; }
  if (!X[i]) h += `<span class="flag">${MODE === "A" ? "fails expression filter" : "not tested in any table"}</span><br>`;
  const showP = anyPadj;
  h += `<table><tr><td><b>contrast</b></td><td><b>log2FC</b></td>${showP ? "<td><b>padj</b></td>" : ""}${anyBM ? "<td><b>abundance</b></td>" : ""}${NORM ? '<td></td><td colspan="2"><b>normalised counts</b></td>' : ""}</tr>`;
  const rows = Math.max(NC, NORM ? LIBS.length : 0);
  for (let r = 0; r < rows; r++) {
    h += "<tr>";
    if (r < NC) { const v = FCV[r][i]; const big = v !== null && Math.abs(v) >= 1; h += `<td>${esc(CT[r].label)}${isRef(r) ? " <span class='small'>(ref)</span>" : ""}</td><td style="color:${big ? "#9C4B22" : INK};font-weight:${big ? 700 : 400}">${f2(v)}</td>`; if (showP) h += `<td>${PADJ[r] ? fmtP(PADJ[r][i]) : ""}</td>`; if (anyBM) h += `<td class="small">${BM[r] && BM[r][i] !== null ? esc(CT[r].abundance || "") + " " + BM[r][i].toFixed(1) : ""}</td>`; }
    else h += "<td></td><td></td>" + (showP ? "<td></td>" : "") + (anyBM ? "<td></td>" : "");
    if (NORM) h += r < LIBS.length ? `<td style="width:8px"></td><td>${esc(LIBS[r].label || LIBS[r].id)}</td><td>${NORM[r][i] === null ? "n.a." : NORM[r][i].toFixed(1)}</td>` : "<td></td><td></td><td></td>";
    h += "</tr>";
  }
  return h + "</table>";
}
const tip = $("tip");
function showTip(html, ev) { tip.innerHTML = html; tip.style.display = "block"; const w = tip.offsetWidth, hgt = tip.offsetHeight; let x = ev.clientX + 14, y = ev.clientY + 14;
  if (x + w > window.innerWidth - 6) x = ev.clientX - w - 14; if (y + hgt > window.innerHeight - 6) y = Math.max(4, ev.clientY - hgt - 14); tip.style.left = x + "px"; tip.style.top = y + "px"; }
function hideTip() { tip.style.display = "none"; }
function geneAt(r, pos) { let lo = repLo[r], hi = repHi[r], best = -1; if (lo < 0) return -1;
  while (lo <= hi) { const m = (lo + hi) >> 1; if (F.s[m] <= pos) { best = m; lo = m + 1; } else hi = m - 1; }
  for (let k = best; k >= Math.max(repLo[r], best - 6); k--) { if (k >= 0 && F.s[k] <= pos && F.e[k] >= pos) return k; } return -1; }
function nearestGene(r, pos, tolBp) { const b = geneAt(r, pos); if (b >= 0) return b; let lo = repLo[r], hi = repHi[r]; if (lo < 0) return -1;
  while (lo < hi) { const m = (lo + hi) >> 1; if (F.s[m] < pos) lo = m + 1; else hi = m; }
  let best = -1, bd = Infinity; for (const k of [lo - 1, lo]) { if (k < repLo[r] || k > repHi[r]) continue; const d = Math.min(Math.abs(F.s[k] - pos), Math.abs(F.e[k] - pos)); if (d < bd) { bd = d; best = k; } }
  return bd <= tolBp ? best : -1; }
// ---------- controls: replicons, rings, contrasts ----------
const repSel = $("repSel");
repSel.innerHTML = (repsWithFeat.length > 1 ? `<option value="all">all (${repsWithFeat.length})</option>` : "") + repsWithFeat.map(r => `<option value="${r}">${esc(repName(r))} (${(REP[r].length / 1e6).toFixed(REP[r].length < 1e6 ? 3 : 2)} Mb)</option>`).join("");
const ringSel = $("ringSel");
function fillRingSel() { ringSel.innerHTML = `<option value="all">selected contrasts</option>` + TEST.concat(REFS).map(c => `<option value="${c}">enlarge: ${esc(CT[c].label)}</option>`).join(""); ringSel.value = st.ring; }
fillRingSel();
function fillCtList() {
  $("ctList").innerHTML = TEST.map(c => `<div><label><input type="checkbox" class="ctT" value="${c}" ${st.sel.includes(c) ? "checked" : ""}>${esc(CT[c].label)}</label>${CT[c].kind === "statistical" ? '<span class="badge">padj</span>' : (CT[c].treatment.length > 1 || CT[c].control.length > 1 ? `<span class="badge">${CT[c].treatment.length} vs ${CT[c].control.length}</span>` : '<span class="badge">1 vs 1</span>')}</div>`).join("") +
    REFS.map(c => `<div class="ref"><label><input type="checkbox" class="ctR" value="${c}" ${st.refSel.includes(c) ? "checked" : ""}>${esc(CT[c].label)}</label><span class="badge">reference</span></div>`).join("");
  document.querySelectorAll(".ctT").forEach(el => el.onchange = () => { const v = [...document.querySelectorAll(".ctT:checked")].map(e => +e.value); if (v.length > 6) { el.checked = false; return; } st.sel = v.length ? v : st.sel; if (!v.length) el.checked = true; redrawAll(); });
  document.querySelectorAll(".ctR").forEach(el => el.onchange = () => { st.refSel = [...document.querySelectorAll(".ctR:checked")].map(e => +e.value); redrawAll(); });
  $("refLab").style.display = REFS.length ? "" : "none";
}
fillCtList();
// ---------- circular geometry ----------
const S = 900, cx = S / 2, cy = S / 2;
const GEO = { start: {}, span: {}, gap: 0 };
const A0 = new Float64Array(N), A1 = new Float64Array(N), SM = new Float64Array(N), CM = new Float64Array(N), S0 = new Float64Array(N), C0 = new Float64Array(N), S1 = new Float64Array(N), C1 = new Float64Array(N), VIS = new Uint8Array(N);
function computeGeo() {
  const reps = st.reps; const tot = reps.reduce((a, r) => a + REP[r].length, 0);
  const gap = reps.length === 1 ? 6 * Math.PI / 180 : Math.max(0.6, Math.min(3, 40 / reps.length)) * Math.PI / 180;
  const avail = 2 * Math.PI - gap * reps.length; let a = gap / 2; GEO.start = {}; GEO.span = {}; GEO.gap = gap;
  for (const r of reps) { GEO.start[r] = a; GEO.span[r] = avail * REP[r].length / tot; a += GEO.span[r] + gap; }
  for (let i = 0; i < N; i++) { const r = F.rep[i]; if (!(r in GEO.start)) { VIS[i] = 0; continue; } VIS[i] = 1; A0[i] = ang(r, F.s[i]); A1[i] = ang(r, F.e[i] + 1) ; if (A1[i] < A0[i]) A1[i] = A0[i];
    const m = (A0[i] + A1[i]) / 2; SM[i] = Math.sin(m); CM[i] = Math.cos(m); S0[i] = Math.sin(A0[i]); C0[i] = Math.cos(A0[i]); S1[i] = Math.sin(A1[i]); C1[i] = Math.cos(A1[i]); }
}
function ang(r, p) { const L = REP[r].length; return GEO.start[r] + (Math.min(Math.max(p - 1, 0), L) / L) * GEO.span[r]; }
function posFromAng(a) { for (const r of st.reps) { const s0 = GEO.start[r], sp = GEO.span[r]; if (a >= s0 && a <= s0 + sp) return { r, p: Math.round((a - s0) / sp * REP[r].length) + 1 }; } return null; }
computeGeo();
const R_TICK = 338, R_STR = [318, 332], R_COG = [305, 314];
let Z = d3.zoomIdentity; const CXs = () => Z.x + Z.k * cx, CYs = () => Z.y + Z.k * cy;
function layout() {
  const ref = st.refSel.length ? st.refSel[0] : null; const rings = [];
  if (st.ring === "all") { const rc = st.sel; const n = rc.length; const inner = ref !== null ? 100 : 60; const th = Math.min(48, (300 - inner - 6 * (n - 1)) / Math.max(n, 1)); let top = 300;
    rc.forEach((c, k) => { rings.push({ c, r0: top - th, r1: top, lab: String(k + 1) }); top -= th + 6; });
    if (ref !== null && !rc.includes(ref)) rings.push({ c: ref, r0: 48, r1: 84, lab: "N" }); }
  else { const c = +st.ring; rings.push({ c, r0: 120, r1: 300, lab: isRef(c) ? "N" : String(Math.max(1, st.sel.indexOf(c) + 1)) }); if (ref !== null && ref !== c) rings.push({ c: ref, r0: 60, r1: 110, lab: "N" }); }
  return rings;
}
const MARG = 60;
function inView(sx, sy, m) { m = m || MARG; return sx > -m && sx < S + m && sy > -m && sy < S + m; }
function arcVisible(sa, ca, sb, cb, smid, cmid, r) { const Xc = CXs(), Yc = CYs(), rr = r * Z.k; return inView(Xc + rr * smid, Yc - rr * cmid) || inView(Xc + rr * sa, Yc - rr * ca) || inView(Xc + rr * sb, Yc - rr * cb); }
function geneVisible(i, r) { return VIS[i] && (Z.k < 1.5 || arcVisible(S0[i], C0[i], S1[i], C1[i], SM[i], CM[i], r)); }
const canvas = $("circCanvas"), dpr = window.devicePixelRatio || 1;
canvas.width = S * dpr; canvas.height = S * dpr; canvas.style.width = S + "px"; canvas.style.height = S + "px";
const ctx = canvas.getContext("2d");
// painter abstraction: canvas (screen) or SVG string (export)
let P = null;
const canvasPainter = { begin() { ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, S, S); },
  wedge(a0, a1, R0, R1, fill, alpha) { ctx.globalAlpha = alpha; ctx.beginPath(); ctx.arc(CXs(), CYs(), R1, a0 - Math.PI / 2, a1 - Math.PI / 2, false); ctx.arc(CXs(), CYs(), R0, a1 - Math.PI / 2, a0 - Math.PI / 2, true); ctx.closePath(); ctx.fillStyle = fill; ctx.fill(); ctx.globalAlpha = 1; },
  annulus(R0, R1, fill) { ctx.beginPath(); ctx.arc(CXs(), CYs(), R1, 0, 2 * Math.PI); ctx.arc(CXs(), CYs(), R0, 0, 2 * Math.PI, true); ctx.fillStyle = fill; ctx.fill(); },
  arcLine(a0, a1, R, color, w) { ctx.lineWidth = w; ctx.strokeStyle = color; ctx.beginPath(); ctx.arc(CXs(), CYs(), R, a0 - Math.PI / 2, a1 - Math.PI / 2); ctx.stroke(); } };
function svgPainter() { const parts = []; const Xc = () => CXs(), Yc = () => CYs(); const pt = (a, R) => [(Xc() + R * Math.sin(a)).toFixed(2), (Yc() - R * Math.cos(a)).toFixed(2)];
  return { parts, begin() { }, wedge(a0, a1, R0, R1, fill, alpha) { const lg = (a1 - a0) > Math.PI ? 1 : 0; const [x1, y1] = pt(a0, R1), [x2, y2] = pt(a1, R1), [x3, y3] = pt(a1, R0), [x4, y4] = pt(a0, R0);
      parts.push(`<path d="M${x1},${y1}A${R1.toFixed(2)},${R1.toFixed(2)} 0 ${lg} 1 ${x2},${y2}L${x3},${y3}A${R0.toFixed(2)},${R0.toFixed(2)} 0 ${lg} 0 ${x4},${y4}Z" fill="${fill}"${alpha < 1 ? ` fill-opacity="${alpha}"` : ""}/>`); },
    annulus(R0, R1, fill) { parts.push(`<circle cx="${Xc().toFixed(2)}" cy="${Yc().toFixed(2)}" r="${((R0 + R1) / 2).toFixed(2)}" fill="none" stroke="${fill}" stroke-width="${(R1 - R0).toFixed(2)}"/>`); },
    arcLine(a0, a1, R, color, w) { const [x1, y1] = pt(a0, R), [x2, y2] = pt(a1, R); parts.push(`<path d="M${x1},${y1}A${R.toFixed(2)},${R.toFixed(2)} 0 ${(a1 - a0) > Math.PI ? 1 : 0} 1 ${x2},${y2}" fill="none" stroke="${color}" stroke-width="${w}"/>`); } }; }
function wedge(a0, a1, r0, r1, fill, alpha) { const k = Z.k, R0 = r0 * k, R1 = r1 * k; const minA = 0.7 / Math.max(R1, 1);
  if (a1 - a0 < minA) { const m = (a0 + a1) / 2; a0 = m - minA / 2; a1 = m + minA / 2; } P.wedge(a0, a1, R0, R1, fill, alpha === undefined ? 1 : alpha); }
function drawCircle(painter) {
  P = painter || canvasPainter; P.begin(); const k = Z.k;
  const rm = (R_STR[0] + R_STR[1]) / 2;
  for (let i = 0; i < N; i++) { if (!st.showFilt && !X[i]) continue; if (!geneVisible(i, rm)) continue; const al = X[i] ? 1 : 0.35;
    if (F.st[i] > 0) wedge(A0[i], A1[i], rm + 0.5 / k, R_STR[1], INK2, al); else wedge(A0[i], A1[i], R_STR[0], rm - 0.5 / k, NEUT, al); }
  if (hasCOG) for (let i = 0; i < N; i++) { if (!st.showFilt && !X[i]) continue; if (!geneVisible(i, R_COG[1])) continue; wedge(A0[i], A1[i], R_COG[0], R_COG[1], SG[F.group[i]] ? SG[F.group[i]][0] : SG.none[0], X[i] ? 1 : 0.35); }
  for (const R of layout()) {
    const base = (R.r0 + R.r1) / 2, half = (R.r1 - R.r0) / 2 - 1; P.annulus(R.r0 * k, R.r1 * k, "#FBFBFD"); const V = FCV[R.c];
    if (st.mode === "gene") {
      for (let i = 0; i < N; i++) { const v = V[i]; const pass = X[i] && v !== null; if (!pass && !st.showFilt) continue; if (v === null) continue; if (!geneVisible(i, base)) continue;
        const Lh = Math.min(Math.abs(v), VM) / VM * half; if (Lh * k < 0.3) continue; const al = pass ? sigAlpha(R.c, i) : 0.45;
        if (v >= 0) wedge(A0[i], A1[i], base, base + Lh, pass ? col(v) : GREY, al); else wedge(A0[i], A1[i], base - Lh, base, pass ? col(v) : GREY, al); }
    } else {
      for (const r of st.reps) { const B = BINS.reps[r], bw = BINS.bw, Lr = REP[r].length;
        for (let b = 0; b < B.n; b++) { const a0 = ang(r, b * bw + 1), a1 = ang(r, Math.min((b + 1) * bw, Lr) + 1);
          if (Z.k >= 1.5 && !arcVisible(Math.sin(a0), Math.cos(a0), Math.sin(a1), Math.cos(a1), Math.sin((a0 + a1) / 2), Math.cos((a0 + a1) / 2), base)) continue;
          if (st.mode === "med") { const v = B.med[R.c][b]; if (v === null) continue; const Lh = Math.min(Math.abs(v), VM) / VM * half; if (Lh * k < 0.3) continue;
            if (v >= 0) wedge(a0, a1, base, base + Lh, col(v)); else wedge(a0, a1, base - Lh, base, col(v)); }
          else { const fu = B.fu[R.c][b], fd = B.fd[R.c][b]; if (fu === null) continue; const Lu = Math.min(fu, 0.5) / 0.5 * half, Ld = Math.min(fd, 0.5) / 0.5 * half;
            if (Lu * k > 0.3) wedge(a0, a1, base, base + Lu, "#C86A3C"); if (Ld * k > 0.3) wedge(a0, a1, base - Ld, base, "#3D5A9E"); } } }
    }
    for (const r of st.reps) P.arcLine(GEO.start[r], GEO.start[r] + GEO.span[r], base * k, INK, 0.6);
  }
}
// ---------- SVG overlay ----------
const svg = d3.select("#circSvg").attr("width", S).attr("height", S).style("font-family", "inherit");
const gOver = svg.append("g"), gGL = svg.append("g"), gSel = svg.append("g"), gLm = svg.append("g");
const arcGen = d3.arc(); function arcPathPx(a0, a1, R0, R1) { return arcGen({ innerRadius: R0, outerRadius: R1, startAngle: a0, endAngle: a1 }); }
const TICK_STEPS = [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 2e4, 1e4, 5e3, 2e3, 1e3, 500, 200]; const LABEL_MODE_K = 20;
function tickFmt(p, step) { if (p === 0) return "0"; if (step >= 1e5) return (p / 1e6).toFixed(1) + " Mb"; if (step >= 1e4) return (p / 1e6).toFixed(2) + " Mb"; if (step >= 1e3) return (p / 1e6).toFixed(3) + " Mb"; return fmt(p); }
function tangentRot(a) { let d = a * 180 / Math.PI; if (d > 90 && d < 270) d += 180; return d; }
function drawOverlay() {
  gOver.selectAll("*").remove(); gGL.selectAll("*").remove();
  const k = Z.k, Xc = CXs(), Yc = CYs(), Rt = R_TICK * k, labelMode = k >= LABEL_MODE_K;
  const T = gOver.append("g").attr("transform", `translate(${Xc},${Yc})`);
  const r0 = st.reps[0]; const bpPx = REP[r0].length / (GEO.span[r0] * Rt); let step = TICK_STEPS[0]; for (const s of TICK_STEPS) { if (s / bpPx >= 90) step = s; }
  const minor = step / 5 >= 100 ? step / 5 : step;
  for (const r of st.reps) {
    const Lr = REP[r].length; T.append("path").attr("d", arcPathPx(GEO.start[r], GEO.start[r] + GEO.span[r], Rt, Rt + 0.8)).attr("fill", INK);
    if (st.reps.length > 1) T.append("path").attr("d", arcPathPx(GEO.start[r], GEO.start[r] + GEO.span[r], Rt - 2, Rt + 9)).attr("fill", "transparent").style("cursor", "pointer")
      .on("mousemove", ev => showTip(`<b>${esc(repName(r))}</b> <span class="small">${esc(REP[r].id)}</span><br>${fmt(Lr)} bp · ${REP[r].topology}<br><span class="small">click: show this replicon alone</span>`, ev)).on("mouseleave", hideTip)
      .on("mousedown", ev => ev.stopPropagation()).on("click", ev => { ev.stopPropagation(); hideTip(); setReps(String(r)); });
    const nT = Lr / minor; if (nT > 4000) continue;
    for (let p = 0; p <= Lr; p += minor) { const a = ang(r, p + 1), sa = Math.sin(a), ca = Math.cos(a); if (!inView(Xc + Rt * sa, Yc - Rt * ca, 80)) continue;
      const major = (Math.round(p) % step === 0); const Lt = major ? 7 : 3.5;
      T.append("line").attr("x1", Rt * sa).attr("y1", -Rt * ca).attr("x2", (Rt + Lt) * sa).attr("y2", -(Rt + Lt) * ca).attr("stroke", INK).attr("stroke-width", major ? 1 : 0.6);
      if (major && GEO.span[r] * Rt > 60) { const txt = tickFmt(p, step);
        if (!labelMode) { const tx = (Rt + 11) * sa, ty = -(Rt + 11) * ca; T.append("text").attr("x", tx).attr("y", ty).attr("text-anchor", sa > 0.3 ? "start" : (sa < -0.3 ? "end" : "middle")).attr("dominant-baseline", "middle").attr("font-size", 10.5).attr("fill", INK).text(txt); }
        else { const tx = (Rt + 16) * sa, ty = -(Rt + 16) * ca; T.append("text").attr("transform", `translate(${tx},${ty}) rotate(${tangentRot(a)})`).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-size", 10.5).attr("fill", INK).text(txt); } } }
    if (st.reps.length > 1 && GEO.span[r] * (Rt + 50) > repName(r).length * 7 + 10) { const am = GEO.start[r] + GEO.span[r] / 2, rr = Rt + 50;
      if (inView(Xc + rr * Math.sin(am), Yc - rr * Math.cos(am), 100)) T.append("text").attr("transform", `translate(${rr * Math.sin(am)},${-rr * Math.cos(am)}) rotate(${tangentRot(am)})`).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-size", 11).attr("font-weight", 700).attr("fill", "#3D5A9E").text(repName(r)); }
  }
  const rings = layout(); const lab = [];
  lab.push({ r: (R_STR[0] + R_STR[1]) / 2, t: "S" }); if (hasCOG) lab.push({ r: (R_COG[0] + R_COG[1]) / 2, t: "C" });
  rings.forEach(R => lab.push({ r: (R.r0 + R.r1) / 2, t: R.lab }));
  const la = st.reps.length > 1 ? 0 : 0; // label column at 12 o'clock (inside the start gap)
  for (const L of lab) { if (!inView(Xc, Yc - L.r * k, 10)) continue; T.append("text").attr("x", 0).attr("y", -L.r * k).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-size", 9.5).attr("font-weight", 700).attr("fill", "#3D5A9E").text(L.t); }
  for (const R of rings) for (const r of st.reps) T.append("path").attr("d", arcPathPx(GEO.start[r], GEO.start[r] + GEO.span[r], R.r0 * k, R.r1 * k)).attr("fill", "none").attr("stroke", GRID).attr("stroke-width", 0.6);
  if (inView(Xc, Yc, 40)) { T.append("text").attr("y", -8).attr("text-anchor", "middle").attr("font-size", 13).attr("font-weight", 600).attr("fill", INK).text(ds.strain || (st.reps.length === 1 ? repName(st.reps[0]) : "")).style("font-style", ds.strain ? "normal" : "normal");
    const shownLen = st.reps.reduce((a, r) => a + REP[r].length, 0); T.append("text").attr("y", 9).attr("text-anchor", "middle").attr("font-size", 10.5).attr("fill", INK2).text((shownLen / 1e6).toFixed(shownLen < 1e6 ? 3 : 2) + " Mb" + (st.reps.length > 1 ? ` · ${st.reps.length} replicons` : "")); }
  let nLab = 0;
  if (labelMode) { const TG = gGL.append("g").attr("transform", `translate(${Xc},${Yc})`);
    const rC = Math.hypot(S / 2 - Xc, S / 2 - Yc) / k; const cands = [{ r: R_TICK + 52 / k }].concat(rings.map(R => ({ r: R.r1 + 3 / k })));
    let best = cands[0]; for (const c of cands) { if (Math.abs(c.r - rC) < Math.abs(best.r - rC)) best = c; }
    const rr0 = best.r * k; svg.node().dataset.labelR = best.r.toFixed(1);
    for (let i = 0; i < N; i++) { if (!VIS[i]) continue; if (!st.showFilt && !X[i]) continue; const w = (A1[i] - A0[i]) * rr0; if (w < 10.5) continue; if (!inView(Xc + rr0 * SM[i], Yc - rr0 * CM[i], 20)) continue;
      const am = (A0[i] + A1[i]) / 2, deg = am * 180 / Math.PI, flip = deg > 180;
      TG.append("text").attr("transform", `translate(${rr0 * SM[i]},${-rr0 * CM[i]}) rotate(${flip ? deg + 90 : deg - 90})`).attr("text-anchor", flip ? "end" : "start").attr("dominant-baseline", "middle")
        .attr("font-size", Math.min(10.5, w - 1)).attr("font-style", (F.gene[i] || OVR(F.id[i])) ? "italic" : "normal").attr("fill", X[i] ? INK : INK2).style("cursor", "pointer")
        .text(linLabel(i)).datum(i).on("mousemove", (ev, d) => showTip(tipGene(d), ev)).on("mouseleave", hideTip).on("mousedown", ev => ev.stopPropagation())
        .on("click", (ev, d) => { ev.stopPropagation(); st.hi = d; setWin(F.rep[d], F.s[d] - 10000, F.e[d] + 10000); });
      nLab++; } }
  svg.node().dataset.nlabels = nLab;
  let cpos = ""; if (k > 1.001) { let ac = Math.atan2(S / 2 - Xc, -(S / 2 - Yc)); if (ac < 0) ac += 2 * Math.PI; const pp = posFromAng(ac); if (pp) cpos = ` · view centre ≈ ${st.reps.length > 1 ? repName(pp.r) + " " : ""}${(pp.p / 1e6).toFixed(3)} Mb`; }
  $("zoomInfo").textContent = k > 1.001 ? `zoom ×${k.toFixed(1)}${cpos}${labelMode ? " · gene labels on" : ` · labels from ×${LABEL_MODE_K}`}` : "";
  let h = `<div><b style="color:#3D5A9E">S</b> strand (outer half +, inner half −)</div>`; if (hasCOG) h += `<div><b style="color:#3D5A9E">C</b> COG functional group</div>`;
  for (const R of rings) h += `<div><b style="color:#3D5A9E">${R.lab}</b> ${esc(CT[R.c].label)}${isRef(R.c) ? " (reference / noise baseline)" : ""}</div>`;
  h += `<div class="small">Bars point outward for up, inward for down, from the baseline circle.${anyPadj ? ` Faded bars: padj ≥ ${ALPHA} (when the fade option is on).` : ""}</div>`;
  $("ringKey").innerHTML = h;
  drawLandmarks(); drawSelection();
}
function lmZoom(L) { const pad = Math.max(8000, (L.e - L.s) * 0.3); setWin(L.rep, L.s - pad, L.e + pad); }
const LMS = ds.landmarks || [];
function drawLandmarks() {
  gLm.selectAll("*").remove(); const k = Z.k, Xc = CXs(), Yc = CYs(), Rt = R_TICK * k, labelMode = k >= LABEL_MODE_K;
  const T = gLm.append("g").attr("transform", `translate(${Xc},${Yc})`);
  LMS.forEach(L => { if (!(L.rep in GEO.start)) return;
    let a0 = ang(L.rep, L.s), a1 = ang(L.rep, L.e + 1); const minA = 2 / Rt; if (a1 - a0 < minA) { const m = (a0 + a1) / 2; a0 = m - minA / 2; a1 = m + minA / 2; }
    const am = (a0 + a1) / 2; if (!(inView(Xc + (Rt + 30) * Math.sin(am), Yc - (Rt + 30) * Math.cos(am), 150) || inView(Xc + (Rt + 30) * Math.sin(a0), Yc - (Rt + 30) * Math.cos(a0)) || inView(Xc + (Rt + 30) * Math.sin(a1), Yc - (Rt + 30) * Math.cos(a1)))) return;
    const cl = lmColor(L);
    const gg = T.append("g").style("cursor", "pointer").on("mousedown", ev => ev.stopPropagation()).on("click", ev => { ev.stopPropagation(); lmZoom(L); })
      .on("mousemove", ev => showTip(`<b>${esc(L.name)}</b>${L.type ? ` <span class="small">(${esc(L.type)})</span>` : ""}<br>${esc(repName(L.rep))}: ${fmt(L.s)}–${fmt(L.e)}<br><span class="small">click to zoom</span>`, ev)).on("mouseleave", hideTip);
    gg.append("path").attr("d", arcPathPx(a0, a1, Rt + 26, Rt + 33)).attr("fill", cl);
    const label = L.label || L.name; const deg = am * 180 / Math.PI, flip = deg > 180; const ital = /^[a-z]{3}[A-Z]?$/.test(label);
    if (!labelMode) { const r = Rt + 37; gg.append("text").attr("transform", `translate(${r * Math.sin(am)},${-r * Math.cos(am)}) rotate(${flip ? deg + 90 : deg - 90})`).attr("text-anchor", flip ? "end" : "start").attr("dominant-baseline", "middle").attr("font-size", 10.5).attr("font-weight", 600).attr("font-style", ital ? "italic" : "normal").attr("fill", cl).text(label); }
    else { const r = Rt + 42; gg.append("text").attr("transform", `translate(${r * Math.sin(am)},${-r * Math.cos(am)}) rotate(${tangentRot(am)})`).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-size", 10.5).attr("font-weight", 600).attr("font-style", ital ? "italic" : "normal").attr("fill", cl).text(label); }
  });
}
function fillLmList() { if (!LMS.length) { $("lmPanel").style.display = "none"; return; }
  $("lmList").innerHTML = LMS.map((L, k) => `<div><span class="lm" data-k="${k}">${esc(L.name)}</span> <span class="small">${REP.length > 1 ? esc(repName(L.rep)) + ": " : ""}${fmt(L.s)}–${fmt(L.e)}</span></div>`).join("");
  $("lmList").querySelectorAll(".lm").forEach(el => el.addEventListener("click", () => { const L = LMS[+el.dataset.k]; if (!(L.rep in GEO.start)) { setReps(repsWithFeat.length > 1 ? "all" : String(L.rep)); } lmZoom(L); })); }
fillLmList();
function drawSelection() {
  gSel.selectAll("*").remove(); const { r, s, e } = st.win; if (!(r in GEO.start)) return; const G = REP[r].length; const len = e - s; if (len >= G) return;
  const s1 = mod(s - 1, G) + 1; const a0 = ang(r, s1); let a1 = a0 + (len / G) * GEO.span[r];
  gSel.append("path").attr("transform", `translate(${CXs()},${CYs()})`).attr("d", arcPathPx(a0, Math.max(a1, a0 + 1.5 / (R_TICK * Z.k)), 44 * Z.k, (R_TICK + 22) * Z.k))
    .attr("fill", "#3D5A9E").attr("fill-opacity", 0.10).attr("stroke", "#3D5A9E").attr("stroke-width", 1).style("pointer-events", "none");
}
// ---------- zoom ----------
let raf = 0; function scheduleDraw() { if (raf) return; raf = requestAnimationFrame(() => { raf = 0; drawCircle(); drawOverlay(); }); }
const zoom = d3.zoom().extent([[0, 0], [S, S]]).scaleExtent([1, 300]).translateExtent([[0, 0], [S, S]]).clickDistance(4)
  .filter(ev => (ev.type === "wheel" || ev.type === "dblclick" || (!ev.shiftKey && st.dmode === "pan")) && !ev.button && (!ev.ctrlKey || ev.type === "wheel"))
  .on("zoom", ev => { Z = ev.transform; hideTip(); scheduleDraw(); });
svg.call(zoom);
function zoomBy(f) { svg.transition().duration(250).call(zoom.scaleBy, f, [S / 2, S / 2]); }
$("czin").onclick = () => zoomBy(2); $("czout").onclick = () => zoomBy(0.5);
$("creset").onclick = () => svg.transition().duration(300).call(zoom.transform, d3.zoomIdentity);
function zoomToWindow(animate) { const { r, s, e } = st.win; if (!(r in GEO.start)) return; const G = REP[r].length, len = Math.max(e - s, 500); const am = ang(r, mod(s - 1, G) + 1) + (len / G) * GEO.span[r] / 2;
  const arc = (len / G) * GEO.span[r] * R_TICK; const k = Math.max(1, Math.min(300, 520 / arc)); const rr = 292, px = cx + rr * Math.sin(am), py = cy - rr * Math.cos(am);
  const t = d3.zoomIdentity.translate(S / 2 - k * px, S / 2 - k * py).scale(k); if (animate) svg.transition().duration(400).call(zoom.transform, t); else svg.call(zoom.transform, t); }
$("czwin").onclick = () => zoomToWindow(true);
document.querySelectorAll("input[name=dmode]").forEach(el => el.addEventListener("change", () => { st.dmode = el.value; svg.style("cursor", st.dmode === "select" ? "crosshair" : "grab"); }));
svg.style("cursor", "grab");
function polar(ev) { const rect = svg.node().getBoundingClientRect(); const x = ev.clientX - rect.left - CXs(), y = ev.clientY - rect.top - CYs(); let a = Math.atan2(x, -y); if (a < 0) a += 2 * Math.PI; return { a, r: Math.hypot(x, y) / Z.k }; }
let drag = null;
svg.on("mousedown.sel", ev => { if (!(ev.shiftKey || st.dmode === "select") || ev.button) return; const Pp = polar(ev); if (Pp.r < 44 || Pp.r > R_TICK + 22) return; drag = { a0: Pp.a, last: Pp.a, cum: 0, moved: false }; ev.preventDefault(); });
window.addEventListener("mousemove", ev => { if (!drag) return; const Pp = polar(ev); let d = Pp.a - drag.last; if (d > Math.PI) d -= 2 * Math.PI; if (d < -Math.PI) d += 2 * Math.PI;
  drag.cum += d; drag.last = Pp.a; if (Math.abs(drag.cum) * R_TICK * Z.k > 4) drag.moved = true;
  if (drag.moved) { let a0 = drag.a0, a1 = drag.a0 + drag.cum; if (a1 < a0) { const t = a0; a0 = a1; a1 = t; }
    gSel.selectAll("*").remove(); gSel.append("path").attr("transform", `translate(${CXs()},${CYs()})`).attr("d", arcPathPx(a0, a1, 44 * Z.k, (R_TICK + 22) * Z.k)).attr("fill", "#3D5A9E").attr("fill-opacity", 0.12).attr("stroke", "#3D5A9E").style("pointer-events", "none"); } });
window.addEventListener("mouseup", () => { if (!drag) return; const D = drag; drag = null;
  if (!D.moved) { const pp = posFromAng(D.a0); if (pp) setWin(pp.r, pp.p - 15000, pp.p + 15000); return; }
  let a0 = D.a0, a1 = D.a0 + D.cum; if (a1 < a0) { const t = a0; a0 = a1; a1 = t; }
  const p0 = posFromAng(mod(a0, 2 * Math.PI)), p1 = posFromAng(mod(a1, 2 * Math.PI)); if (!p0) { drawSelection(); return; }
  const r = p0.r; const span = (a1 - a0) / GEO.span[r] * REP[r].length; setWin(r, p0.p, p0.p + Math.round(span)); });
svg.on("click.centre", ev => { if (ev.defaultPrevented || ev.shiftKey || st.dmode === "select") return; const Pp = polar(ev); if (Pp.r < 44 || Pp.r > R_TICK + 22) return; const pp = posFromAng(Pp.a); if (!pp) return; setWin(pp.r, pp.p - 15000, pp.p + 15000); });
svg.on("mousemove.tip", ev => { if (drag || ev.buttons) return; const Pp = polar(ev); if (Pp.r > R_TICK + 22 || Pp.r < 44) { hideTip(); return; }
  const pp = posFromAng(Pp.a); if (!pp) { hideTip(); return; }
  const rings = layout(); let ring = null; for (const R of rings) if (Pp.r >= R.r0 && Pp.r <= R.r1) ring = R;
  const inStrand = Pp.r >= R_STR[0] && Pp.r <= R_STR[1], inCog = hasCOG && Pp.r >= R_COG[0] && Pp.r <= R_COG[1];
  if (!ring && !inStrand && !inCog) { hideTip(); return; }
  if (ring && st.mode !== "gene") { const B = BINS.reps[pp.r], b = Math.min(Math.floor((pp.p - 1) / BINS.bw), B.n - 1); const m = B.med[ring.c][b];
    showTip(`<b>${esc(CT[ring.c].label)}</b><br>${esc(repName(pp.r))} bin ${fmt(b * BINS.bw + 1)}–${fmt(Math.min((b + 1) * BINS.bw, REP[pp.r].length))}<br>${B.ng[b]} features (${B.nx[b]} ${MODE === "A" ? "pass filter" : "tested"})<br>` +
      (m === null ? "no expressed genes" : `median log2FC ${f2(m)}<br>fraction &gt; +1: ${(B.fu[ring.c][b] * 100).toFixed(0)}% · &lt; −1: ${(B.fd[ring.c][b] * 100).toFixed(0)}%`) + `<br><span class="small">bins use ${MODE === "A" ? "expression-filtered" : "tested"} genes, assigned by midpoint</span>`, ev); return; }
  const i = nearestGene(pp.r, pp.p, Math.max(50, 800 / Z.k)); if (i < 0) { hideTip(); return; } if (!st.showFilt && !X[i]) { hideTip(); return; }
  showTip((ring ? `<span class="small">ring: ${esc(CT[ring.c].label)}</span><br>` : "") + tipGene(i), ev); });
svg.on("mouseleave.tip", hideTip);
// ---------- legend + notes ----------
$("lgCog").innerHTML = hasCOG ? Object.keys(SG).map(k => `<span title="${SG[k][1]}"><span class="sw" style="background:${SG[k][0]};margin-top:0"></span>${SG[k][2]}</span>`).join("") : "";
$("lgMode").innerHTML = MODE === "A" ? "log2FC (descriptive ratio of normalised counts); * clipped" : "log2FC from DE tables; * clipped";
$("lgFilt").textContent = MODE === "A" ? `fails expression filter (${(ds.meta && ds.meta.filter) || ""})` : "not tested / missing in DE table";
(function notes() { const M = ds.meta || {}; let h = "";
  h += `<b>Mode ${MODE}</b>. ${esc(M.normalization || "")}.<br>`;
  if (M.filter) h += `Expression filter: ${esc(M.filter)}.<br>`;
  if (M.size_factors) h += `Size factors: ${Object.entries(M.size_factors).map(([k, v]) => `${esc(k)} ${(+v).toFixed(3)}`).join(", ")}.<br>`;
  h += `Contrasts: ${CT.map(c => `<i>${esc(c.label)}</i>${c.kind === "descriptive" ? ` = [${esc(c.treatment.join(", "))}] vs [${esc(c.control.join(", "))}]` : (c.source ? ` (${esc(c.source)})` : "")}${isRef(CT.indexOf(c)) ? " — reference" : ""}`).join("; ")}.<br>`;
  if (Object.keys(LAB).length) h += `Curated labels: ${Object.entries(LAB).map(([k, v]) => `${esc(k)} → <i>${esc(v.label)}</i>${v.note ? " (" + esc(v.note) + ")" : ""}`).join("; ")}.<br>`;
  h += `Gene symbols come from the annotation; where one symbol is assigned to several loci the linear-view label carries the locus suffix.<br>`;
  if ((M.warnings || []).length) h += `<b>Warnings:</b> ${M.warnings.map(esc).join(" · ")}<br>`;
  h += `<span class="small">${esc(M.software || "")} · built ${esc((M.created || "").slice(0, 19).replace("T", " "))} · annotation: ${esc(M.annotation_source || "")}${(M.inputs || []).length ? " · inputs: " + M.inputs.map(x => esc(x.file) + (x.sha256 ? ` <span title="SHA-256 ${esc(x.sha256)}">(sha256 ${esc(x.sha256.slice(0, 12))}…)</span>` : "")).join(", ") : ""}</span>`;
  $("notes").innerHTML = h; })();
// ---------- linear view ----------
const LW = Math.max(1100, $("linWrap").clientWidth - 20), ML = 170, MR = 36;
const GENE_Y0 = 36, GENE_H = 15, LANE_GAP = 26, TR_Y0 = GENE_Y0 + 2 * LANE_GAP + 16, TRK_H = 54, TRK_GAP = 8;
function linTracks() { const t = [];
  if (st.ring === "all") { for (const c of st.sel) t.push({ c, h: TRK_H, tall: false }); }
  else t.push({ c: +st.ring, h: Math.max(4, st.sel.length) * TRK_H + 3 * TRK_GAP, tall: true });
  if (st.showNoise) for (const c of st.refSel) if (!t.some(x => x.c === c)) t.push({ c, h: TRK_H, tall: false });
  let y = TR_Y0; for (const Tk of t) { Tk.y0 = y; y += Tk.h + TRK_GAP; } return t; }
const lsvg = d3.select("#linSvg").attr("width", LW).style("font-family", "inherit");
const lclip = lsvg.append("defs").append("clipPath").attr("id", "lclip").append("rect").attr("x", ML).attr("y", 0).attr("width", LW - ML - MR);
const gAxis = lsvg.append("g"), gBrush = lsvg.append("g"), gGenes = lsvg.append("g").attr("clip-path", "url(#lclip)"), gTr = lsvg.append("g");
function winGenes() { const { r, s, e } = st.win, G = REP[r].length, L = e - s, out = []; const circ = REP[r].topology === "circular";
  for (let i = repLo[r]; i <= repHi[r]; i++) { if (circ) { let o0 = mod(F.s[i] - s, G), o1 = mod(F.e[i] - s, G); if (o0 <= L || o1 <= L) { if (o0 > o1) o0 -= G; out.push([i, o0, o1]); } }
    else { if (F.e[i] >= s && F.s[i] <= e) out.push([i, F.s[i] - s, F.e[i] - s]); } }
  return out; }
function shortLab(c) { return CT[c].short || CT[c].label; }
function drawLinear() {
  const { r, s, e } = st.win, G = REP[r].length, L = e - s; const tracks = linTracks();
  const LH = (tracks.length ? tracks[tracks.length - 1].y0 + tracks[tracks.length - 1].h : TR_Y0) + 16; lsvg.attr("height", LH); lclip.attr("height", LH);
  $("rs").value = mod(s - 1, G) + 1; $("re").value = mod(e - 1, G) + 1;
  $("linTitle").textContent = `${repName(r)} : ${fmt(mod(s - 1, G) + 1)} – ${fmt(mod(e - 1, G) + 1)} (${fmt(L)} bp)` + (s < 1 || e > G ? " · wraps the origin" : "");
  const x = d3.scaleLinear().domain([0, L]).range([ML, LW - MR]);
  gAxis.selectAll("*").remove(); gGenes.selectAll("*").remove(); gTr.selectAll("*").remove();
  const ax = d3.axisTop(x).ticks(10).tickFormat(o => { const p = mod(s + o - 1, G) + 1; return L > 200000 ? (p / 1e6).toFixed(3) + " Mb" : fmt(Math.round(p)); });
  gAxis.append("g").attr("transform", "translate(0,22)").call(ax).call(a => a.selectAll("text").attr("font-size", 10).attr("fill", INK)).call(a => a.selectAll("line,path").attr("stroke", INK));
  gAxis.append("text").attr("x", ML - 8).attr("y", GENE_Y0 + GENE_H / 2 + 2).attr("text-anchor", "end").attr("font-size", 11).attr("font-weight", 600).attr("fill", INK).text("genes (+ strand)");
  gAxis.append("text").attr("x", ML - 8).attr("y", GENE_Y0 + LANE_GAP + GENE_H / 2 + 2).attr("text-anchor", "end").attr("font-size", 11).attr("font-weight", 600).attr("fill", INK).text("genes (− strand)");
  const W = winGenes(); const pxbp = (LW - ML - MR) / L;
  const gg = gGenes.selectAll("g.gene").data(W).join("g").attr("class", "gene").style("cursor", "pointer");
  gg.each(function (d) { const [i, o0, o1] = d; const x0 = x(o0), x1 = x(o1 + 1), w = Math.max(x1 - x0, 0.8), y = GENE_Y0 + (F.st[i] > 0 ? 0 : LANE_GAP);
    const head = Math.min(7, w * 0.45); let pts;
    if (F.st[i] > 0) pts = [[x0, y], [x0 + w - head, y], [x0 + w, y + GENE_H / 2], [x0 + w - head, y + GENE_H], [x0, y + GENE_H]];
    else pts = [[x0 + w, y], [x0 + head, y], [x0, y + GENE_H / 2], [x0 + head, y + GENE_H], [x0 + w, y + GENE_H]];
    const sel = d3.select(this);
    sel.append("polygon").attr("points", pts.map(p => p.join(",")).join(" ")).attr("fill", i === st.hi ? "#F0B48A" : (F.type[i] === "CDS" ? "#E6E6EC" : "#B8B8C0")).attr("fill-opacity", X[i] ? 1 : 0.45)
      .attr("stroke", i === st.hi ? "#9C4B22" : INK2).attr("stroke-width", i === st.hi ? 1.6 : 0.7).attr("stroke-dasharray", X[i] ? null : "2,1.5");
    const lab = linLabel(i); const need = lab.length * 5.6 + 6;
    if (w > need) sel.append("text").attr("x", x0 + w / 2).attr("y", y + GENE_H / 2 + 0.5).attr("text-anchor", "middle").attr("dominant-baseline", "middle").attr("font-size", 9.5)
      .attr("font-style", (F.gene[i] || OVR(F.id[i])) ? "italic" : "normal").attr("fill", INK).style("pointer-events", "none").text(lab);
  }).on("mousemove", (ev, d) => showTip(tipGene(d[0]), ev)).on("mouseleave", hideTip).on("click", (ev, d) => { st.hi = d[0]; drawLinear(); });
  for (const TK of tracks) { const c = TK.c, y0 = TK.y0, H = TK.h, base = y0 + H / 2, half = H / 2 - 2; const V = FCV[c];
    const yv = v => base - Math.max(-VM, Math.min(VM, v)) / VM * half;
    const T = gTr.append("g").attr("class", "track").attr("data-c", c);
    T.append("rect").attr("x", ML).attr("y", y0).attr("width", LW - ML - MR).attr("height", H).attr("fill", "#FBFBFD").attr("stroke", GRID);
    const gs = TK.tall ? [-2 * VM / 3, -VM / 3, VM / 3, 2 * VM / 3] : [-VM / 3, VM / 3];
    gs.forEach(v => T.append("line").attr("x1", ML).attr("x2", LW - MR).attr("y1", yv(v)).attr("y2", yv(v)).attr("stroke", GRID).attr("stroke-dasharray", "3,3"));
    if (TK.tall) { const lines = CT[c].label.split(" · ");
      T.append("text").attr("x", ML - 8).attr("y", base - 14).attr("text-anchor", "end").attr("font-size", 12).attr("font-weight", 600).attr("fill", INK).text(lines[0]);
      if (lines[1]) T.append("text").attr("x", ML - 8).attr("y", base + 1).attr("text-anchor", "end").attr("font-size", 12).attr("font-weight", 600).attr("fill", INK).text(lines[1]);
      T.append("text").attr("x", ML - 8).attr("y", base + 17).attr("text-anchor", "end").attr("font-size", 9.5).attr("fill", INK2).text(`log2FC, ±${VM} range`);
      T.append("text").attr("x", ML - 8).attr("y", y0 + 12).attr("text-anchor", "end").attr("font-size", 9.5).attr("fill", "#9C4B22").text("▲ up");
      T.append("text").attr("x", ML - 8).attr("y", y0 + H - 5).attr("text-anchor", "end").attr("font-size", 9.5).attr("fill", "#3D5A9E").text("▼ down"); }
    else { T.append("text").attr("x", ML - 8).attr("y", base - 4).attr("text-anchor", "end").attr("font-size", 11).attr("font-weight", 600).attr("fill", isRef(c) ? INK2 : INK).text(shortLab(c));
      T.append("text").attr("x", ML - 8).attr("y", base + 10).attr("text-anchor", "end").attr("font-size", 9).attr("fill", INK2).text(isRef(c) ? "reference" : `log2FC, ±${VM} range`); }
    const ticks = TK.tall ? [-VM, -2 * VM / 3, -VM / 3, 0, VM / 3, 2 * VM / 3, VM] : [-VM, -VM / 3, VM / 3, VM];
    ticks.forEach(v => { const yy = Math.min(Math.max(yv(v), y0 + 5), y0 + H - 5); const vv = +v.toFixed(2);
      T.append("line").attr("x1", LW - MR).attr("x2", LW - MR + 3).attr("y1", yv(v)).attr("y2", yv(v)).attr("stroke", INK2).attr("stroke-width", 0.6);
      T.append("text").attr("class", "ylab").attr("x", LW - MR + 5).attr("y", yy).attr("text-anchor", "start").attr("dominant-baseline", "middle").attr("font-size", 8.5).attr("fill", INK2).text(vv === 0 ? "0" : (vv > 0 ? "+" : "−") + Math.abs(vv)); });
    const Bc = T.append("g").attr("clip-path", "url(#lclip)");
    Bc.selectAll("rect.b").data(W.filter(d => V[d[0]] !== null && (st.showFilt || X[d[0]]))).join("rect").attr("class", "b")
      .attr("x", d => x(d[1])).attr("width", d => Math.max(x(d[2] + 1) - x(d[1]) - (pxbp > 0.02 ? 0.5 : 0), 0.8))
      .attr("y", d => Math.min(yv(V[d[0]]), base)).attr("height", d => Math.max(Math.abs(yv(V[d[0]]) - base), 0.6))
      .attr("fill", d => X[d[0]] ? col(V[d[0]]) : GREY).attr("fill-opacity", d => X[d[0]] ? sigAlpha(c, d[0]) : 0.5)
      .attr("stroke", d => d[0] === st.hi ? "#9C4B22" : null).attr("stroke-width", 1.2)
      .on("mousemove", (ev, d) => showTip(`<span class="small">track: ${esc(CT[c].label)}</span><br>` + tipGene(d[0]), ev)).on("mouseleave", hideTip);
    T.append("line").attr("x1", ML).attr("x2", LW - MR).attr("y1", base).attr("y2", base).attr("stroke", INK).attr("stroke-width", 0.7);
  }
  gBrush.selectAll("*").remove();
  const br = d3.brushX().extent([[ML, GENE_Y0 - 6], [LW - MR, GENE_Y0 + LANE_GAP + GENE_H + 4]]).on("end", ev => { if (!ev.selection) return; const [a, b] = ev.selection.map(x.invert); if (b - a < 200) return; setWin(r, Math.round(s + a), Math.round(s + b)); });
  gBrush.call(br); gBrush.select(".overlay").style("cursor", "crosshair");
  const n = $("linSvg"); n.dataset.ngenes = W.length; n.dataset.ntracks = tracks.length; n.dataset.tracks = tracks.map(t => t.c + ":" + t.h).join(",");
}
function setWin(r, s, e) { const G = REP[r].length; const circ = REP[r].topology === "circular"; s = Math.round(s); e = Math.round(e);
  if (e - s < 500) { const m = (s + e) / 2; s = Math.round(m - 250); e = Math.round(m + 250); }
  if (e - s > G) { s = 1; e = G; }
  if (circ) { if (s < 1 - G) { s += G; e += G; } if (s > G) { s -= G; e -= G; } }
  else { if (s < 1) { e += 1 - s; s = 1; } if (e > G) { s -= e - G; e = G; } s = Math.max(1, s); }
  st.win = { r, s, e }; drawSelection(); drawLinear(); }
function setReps(v) { st.reps = v === "all" ? repsWithFeat.slice() : [+v]; repSel.value = v; computeGeo(); svg.call(zoom.transform, d3.zoomIdentity); drawCircle(); drawOverlay(); if (!(st.win.r in GEO.start)) setWin(st.reps[0], 1, 30000); }
$("zin").onclick = () => { const { r, s, e } = st.win, m = (s + e) / 2, h = (e - s) / 4; setWin(r, m - h, m + h); };
$("zout").onclick = () => { const { r, s, e } = st.win, m = (s + e) / 2, h = (e - s); setWin(r, m - h, m + h); };
$("pl").onclick = () => { const { r, s, e } = st.win, d = (e - s) / 2; setWin(r, s - d, e - d); };
$("pr").onclick = () => { const { r, s, e } = st.win, d = (e - s) / 2; setWin(r, s + d, e + d); };
$("goBtn").onclick = () => { let s = +$("rs").value, e = +$("re").value; if (!(s > 0) || !(e > 0)) return; const G = REP[st.win.r].length; if (e < s) e += G; setWin(st.win.r, s, e); };
$("modeSel").onchange = ev => { st.mode = ev.target.value; drawCircle(); };
$("vmSel").onchange = ev => { setVM(+ev.target.value); drawCircle(); drawLinear(); };
repSel.onchange = ev => setReps(ev.target.value);
function redrawAll() { drawCircle(); drawOverlay(); drawLinear(); }
ringSel.onchange = ev => { st.ring = ev.target.value; redrawAll(); };
$("showFilt").onchange = ev => { st.showFilt = ev.target.checked; redrawAll(); };
$("showNoise").onchange = ev => { st.showNoise = ev.target.checked; drawLinear(); };
$("sigOnly").onchange = ev => { st.sig = ev.target.checked; drawCircle(); drawLinear(); };
const sIn = $("search"), sRes = $("searchRes");
function doSearch(q) { q = q.trim().toLowerCase(); if (q.length < 2) return [];
  const hits = []; for (let i = 0; i < N; i++) { const id = F.id[i].toLowerCase(), sy = (F.gene[i] || "").toLowerCase(), ov = (OVR(F.id[i]) || "").toLowerCase(), p = (F.product[i] || "").toLowerCase();
    let rank = -1; if (id === q || sy === q || ov === q) rank = 0; else if (sy.startsWith(q) || ov.startsWith(q) || id.includes(q)) rank = 1; else if (sy.includes(q) || p.includes(q)) rank = 2; if (rank >= 0) hits.push([rank, i]); }
  hits.sort((a, b) => a[0] - b[0] || a[1] - b[1]); return hits.slice(0, 40).map(h => h[1]); }
function jump(i) { st.hi = i; if (!(F.rep[i] in GEO.start)) setReps(repsWithFeat.length > 1 ? "all" : String(F.rep[i])); setWin(F.rep[i], F.s[i] - 10000, F.e[i] + 10000); sRes.style.display = "none"; }
sIn.addEventListener("input", () => { const r = doSearch(sIn.value); if (!r.length) { sRes.style.display = "none"; return; }
  sRes.innerHTML = r.map(i => `<div data-i="${i}"><b>${esc(F.id[i])}</b> ${F.gene[i] ? "<i>" + esc(F.gene[i]) + "</i> " : ""}${OVR(F.id[i]) ? "<i>" + esc(OVR(F.id[i])) + "</i> " : ""}<span class="small">${esc(F.product[i])}</span></div>`).join("") + (r.length >= 40 ? '<div class="small">… first 40 matches</div>' : "");
  sRes.style.display = "block"; sRes.querySelectorAll("div[data-i]").forEach(el => el.onclick = () => jump(+el.dataset.i)); });
sIn.addEventListener("keydown", ev => { if (ev.key === "Enter") { const r = doSearch(sIn.value); if (r.length) jump(r[0]); } });
document.addEventListener("click", ev => { if (!sRes.contains(ev.target) && ev.target !== sIn) sRes.style.display = "none"; });
const fileStem = (ds.strain || ds.title || "genomemap").replace(/[^A-Za-z0-9._-]+/g, "_").slice(0, 40);
function csvRows() { const W = winGenes(); const escv = v => { v = v === null || v === undefined ? "" : String(v); return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; };
  const head = ["feature_id", "replicon", "feature_type", "start", "end", "strand", "gene", "label", "product", "cog", "cog_category", MODE === "A" ? "expr_filter_pass" : "tested"]
    .concat(NORM ? LIBS.map(l => "norm_" + l.id) : []).concat(CT.map(c => "log2fc_" + c.id)).concat(CT.filter((c, k) => PADJ[k]).map(c => "padj_" + c.id)).concat(CT.filter((c, k) => BM[k]).map(c => (c.abundance || "abundance") + "_" + c.id));
  const lines = [head.join(",")];
  for (const [i] of W) { const row = [F.id[i], REP[F.rep[i]].id, F.type[i], F.s[i], F.e[i], F.st[i] > 0 ? "+" : "-", F.gene[i], OVR(F.id[i]) || "", F.product[i], F.cog[i], F.cogcat[i], X[i] ? "TRUE" : "FALSE"]
      .concat(NORM ? NORM.map(a => a[i]) : []).concat(FCV.map(a => a[i])).concat(PADJ.filter(p => p).map(p => p[i])).concat(BM.filter(b => b).map(b => b[i])); lines.push(row.map(escv).join(",")); }
  return lines; }
$("csvBtn").onclick = () => { const { r, s, e } = st.win, G = REP[r].length; window.GenomeMapDownload(`${fileStem}_${REP[r].id}_${mod(s - 1, G) + 1}-${mod(e - 1, G) + 1}.csv`, csvRows().join("\n") + "\n", "text/csv"); };
$("jsonBtn").onclick = () => { const o = Object.assign({}, ds); o.meta = Object.assign({}, ds.meta, { bin_width: BW }); window.GenomeMapDownload(`${fileStem}_genomemap_dataset.json`, JSON.stringify(C.compactDataset(o)), "application/json"); };
$("htmlBtn").onclick = () => { const o = Object.assign({}, ds); o.meta = Object.assign({}, ds.meta, { bin_width: BW }); window.GenomeMapSaveReport(o, `${fileStem}_genome_map_report.html`); };
$("newBtn").onclick = () => { try { sessionStorage.setItem("gm-skip-embedded", "1"); } catch (e) { } location.reload(); };
function circleSVG() { const sp = svgPainter(); drawCircle(sp); drawCircle(); const over = $("circSvg").innerHTML;
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${S}" height="${S}" viewBox="0 0 ${S} ${S}" font-family="Open Sauce Sans, Helvetica Neue, Arial, sans-serif"><rect width="${S}" height="${S}" fill="#fff"/>${sp.parts.join("")}${over}</svg>`; }
// ---------- figure export: PNG / SVG / PDF of the circle or the linear view ----------
const FIG_FONT = "Open Sauce Sans, Helvetica Neue, Arial, sans-serif";
let _ff = null; // @font-face rules of the page, embedded so that exported figures keep the viewer font
function fontFaces() { if (_ff === null) { const t = [...document.querySelectorAll("style")].map(x => x.textContent).join("\n"); _ff = (t.match(/@font-face\s*\{[^}]*\}/g) || []).join("\n"); } return _ff; }
function legendSVG(W, kind) { // title, settings and colour scale drawn under the figure
  const P = 20, L = []; let y = 26;
  const txt = (x, yy, t, o = {}) => `<text x="${x}" y="${yy}" font-size="${o.fs || 11}" font-weight="${o.fw || 400}" fill="${o.c || INK}"${o.anchor ? ` text-anchor="${o.anchor}"` : ""}>${esc(t)}</text>`;
  L.push(txt(P, y, ds.title || ds.strain || "genome map", { fs: 15, fw: 700 })); y += 18;
  const sig = anyPadj && st.sig ? ` · faded bars: padj ≥ ${ALPHA}` : "";
  const mode = MODE === "A" ? "Mode A: descriptive log2 fold change of normalised counts (no p-values)" : "Mode B: log2FC and padj from DE tables" + sig;
  const where = kind === "linear" ? $("linTitle").textContent : (st.reps.length === REP.length ? (REP.length > 1 ? `all ${REP.length} replicons` : repName(st.reps[0])) : st.reps.map(repName).join(", ")) + ($("zoomInfo").textContent ? " · " + $("zoomInfo").textContent : "");
  L.push(txt(P, y, `${mode} · ${where}`, { fs: 10.5, c: INK2 })); y += 20;
  if (kind === "circle") {
    const rows = [["S", "genes by strand (outer half +, inner half −)"]]; if (hasCOG) rows.push(["C", "COG functional group"]);
    for (const R of layout()) rows.push([R.lab, CT[R.c].label + (isRef(R.c) ? " (reference / noise baseline)" : "")]);
    const cols = W > 820 && rows.length > 3 ? 2 : 1, per = Math.ceil(rows.length / cols), cw = (W - 2 * P) / cols;
    rows.forEach((r, i) => { const cx = P + Math.floor(i / per) * cw, cy = y + (i % per) * 15; L.push(txt(cx, cy, r[0], { fw: 700, c: "#3D5A9E" }) + txt(cx + 18, cy, r[1])); });
    y += per * 15 + 6; L.push(txt(P, y, "Bars point outward for higher, inward for lower values in the treatment.", { fs: 10, c: INK2 })); y += 16;
  }
  const bw = Math.min(300, W - 2 * P - 220), bx = P, by = y; const gid = "gmLg" + kind;
  L.push(`<defs><linearGradient id="${gid}" x1="0" x2="1" y1="0" y2="0">${STOPS.map((c, i) => `<stop offset="${(i / 8).toFixed(4)}" stop-color="${c}"/>`).join("")}</linearGradient></defs>`);
  L.push(`<rect x="${bx}" y="${by}" width="${bw}" height="10" fill="url(#${gid})" stroke="${INK}" stroke-width="0.5"/>`);
  [-VM, -VM / 2, 0, VM / 2, VM].forEach((v, i) => { const xx = bx + bw * i / 4; L.push(`<line x1="${xx}" x2="${xx}" y1="${by + 10}" y2="${by + 13}" stroke="${INK}" stroke-width="0.6"/>` + txt(xx, by + 24, (v > 0 ? "+" : v < 0 ? "−" : "") + +Math.abs(v).toFixed(2), { fs: 9.5, anchor: "middle" })); });
  L.push(txt(bx + bw + 12, by + 9, `log2 fold change (treatment / control); colour saturates at ±${VM}`, { fs: 10, c: INK2 }));
  if (st.showFilt) L.push(`<rect x="${bx + bw + 12}" y="${by + 15}" width="10" height="8" fill="#C4C4CC" opacity="0.6"/>` + txt(bx + bw + 26, by + 22, MODE === "A" ? "genes failing the expression filter" : "genes not tested", { fs: 10, c: INK2 }));
  y = by + 34; L.push(txt(W - P, y, `genomemap ${C.VERSION} · ${new Date().toISOString().slice(0, 10)}`, { fs: 9, c: INK2, anchor: "end" }));
  return { g: L.join(""), h: y + 10 };
}
function figureSVG(kind, withLegend) {
  let W, H, inner;
  if (kind === "circle") { const sp = svgPainter(); drawCircle(sp); drawCircle(); W = S; H = S; inner = sp.parts.join("") + $("circSvg").innerHTML; }
  else { const n = $("linSvg"); W = +n.getAttribute("width"); H = +n.getAttribute("height"); inner = n.innerHTML; }
  const lg = withLegend ? legendSVG(W, kind) : { g: "", h: 0 };
  const svgText = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H + lg.h}" viewBox="0 0 ${W} ${H + lg.h}" font-family="${FIG_FONT}"><style>${fontFaces()}</style><rect width="${W}" height="${H + lg.h}" fill="#fff"/><g>${inner}</g>${lg.h ? `<g transform="translate(0,${H})">${lg.g}</g>` : ""}</svg>`;
  return { svg: svgText, w: W, h: H + lg.h };
}
function figureName(kind, ext) { if (kind === "circle") return `${fileStem}_circle.${ext}`; const { r, s, e } = st.win, G = REP[r].length; return `${fileStem}_linear_${REP[r].id}_${mod(s - 1, G) + 1}-${mod(e - 1, G) + 1}.${ext}`; }
function rasterize(fig, scale) { return new Promise((res, rej) => { const img = new Image(); const url = URL.createObjectURL(new Blob([fig.svg], { type: "image/svg+xml" }));
  img.onload = () => { const cv = document.createElement("canvas"); cv.width = Math.round(fig.w * scale); cv.height = Math.round(fig.h * scale); const c2 = cv.getContext("2d");
    c2.fillStyle = "#fff"; c2.fillRect(0, 0, cv.width, cv.height); c2.scale(scale, scale); c2.drawImage(img, 0, 0, fig.w, fig.h); URL.revokeObjectURL(url); res(cv); };
  img.onerror = () => { URL.revokeObjectURL(url); rej(new Error("could not render the figure")); }; img.src = url; }); }
async function deflate(bytes) { const cs = new CompressionStream("deflate"); const w = cs.writable.getWriter(); w.write(bytes); w.close(); return new Uint8Array(await new Response(cs.readable).arrayBuffer()); }
async function pdfFromCanvas(cv, fig, title) { // one page, physical size = on-screen size (1 px = 0.75 pt), image at the canvas resolution
  let data, filter;
  if (typeof CompressionStream !== "undefined") { const px = cv.getContext("2d").getImageData(0, 0, cv.width, cv.height).data, rgb = new Uint8Array(cv.width * cv.height * 3);
    for (let i = 0, j = 0; i < px.length; i += 4) { rgb[j++] = px[i]; rgb[j++] = px[i + 1]; rgb[j++] = px[i + 2]; } data = await deflate(rgb); filter = "FlateDecode"; }
  else { const b = await new Promise(r => cv.toBlob(r, "image/jpeg", 0.95)); data = new Uint8Array(await b.arrayBuffer()); filter = "DCTDecode"; }
  return C.buildPdf({ data, filter, imgW: cv.width, imgH: cv.height, pageW: fig.w * 0.75, pageH: fig.h * 0.75, title });
}
const EXPORT_SCALE = 4; // 4 px per screen px = 384 dpi at the on-screen size
async function saveFigure(kind, fmt) {
  const withLegend = $(kind === "circle" ? "cLeg" : "lLeg").checked, fig = figureSVG(kind, withLegend), name = figureName(kind, fmt);
  const btn = $(kind === "circle" ? "cExpMsg" : "lExpMsg"); btn.textContent = fmt === "svg" ? "" : "rendering…";
  try {
    if (fmt === "svg") { window.GenomeMapDownload(name, fig.svg, "image/svg+xml"); return; }
    const cv = await rasterize(fig, EXPORT_SCALE);
    if (fmt === "png") { const b = await new Promise(r => cv.toBlob(r, "image/png")); window.GenomeMapDownload(name, b, "image/png"); }
    else { const t = `${ds.title || ds.strain || "genome map"} · ${kind === "circle" ? "circular map" : "linear view " + $("linTitle").textContent}`; window.GenomeMapDownload(name, await pdfFromCanvas(cv, fig, t), "application/pdf"); }
    btn.textContent = "";
  } catch (e) { btn.textContent = "export failed: " + e.message; }
}
for (const [kind, pre] of [["circle", "c"], ["linear", "l"]]) for (const fmt of ["png", "svg", "pdf"]) $(pre + "Exp_" + fmt).onclick = () => saveFigure(kind, fmt);
// ---------- init ----------
drawCircle(); drawOverlay();
(function initWin() { const L0 = LMS.find(L => L.rep in GEO.start); if (L0) { const gi = F.id.findIndex((id, i) => F.rep[i] === L0.rep && F.s[i] === L0.s && F.e[i] === L0.e); st.hi = gi; setWin(L0.rep, L0.s - 15000, L0.e + 15000); } else setWin(st.reps[0], 1, 30000); })();
// ---------- self-test (#selftest) ----------
if (location.hash === "#selftest") {
  const out = { fonts: document.fonts ? document.fonts.check('12px "Open Sauce Sans"') : null, d3: d3.version, nGenes: N, nReps: REP.length, nContrasts: NC, mode: MODE };
  const lin = $("linSvg"); out.linGenes = +lin.dataset.ngenes; out.tracks_default = lin.dataset.tracks; out.csvLines = csvRows().length; out.csvHead = csvRows()[0];
  st.showNoise = true; drawLinear(); out.tracks_ref = lin.dataset.tracks; st.showNoise = false;
  if (TEST.length) { st.ring = String(TEST[0]); redrawAll(); out.tracks_enlarged = lin.dataset.tracks; st.ring = "all"; redrawAll(); }
  st.mode = "med"; drawCircle(); st.mode = "frac"; drawCircle(); st.mode = "gene"; drawCircle();
  if (repsWithFeat.length > 1) { setReps(String(repsWithFeat[1])); out.rep2_linTitle = $("linTitle").textContent; out.rep2_genes = +lin.dataset.ngenes; setReps("all"); out.repsAll = st.reps.length; }
  const q = (LMS[0] && LMS[0].label) || F.gene.find(x => x) || F.id[0]; out.search = [q, doSearch(q).slice(0, 3).map(i => F.id[i])];
  const svgTxt = circleSVG(); out.svgBytes = svgTxt.length; out.svgPaths = (svgTxt.match(/<path/g) || []).length;
  const fc = figureSVG("circle", true), fl = figureSVG("linear", true); out.figCircle = [fc.w, fc.h, fc.svg.includes("colour saturates")]; out.figLinear = [fl.w, fl.h, fl.svg.includes(ds.title ? esc(ds.title) : "genome map")];
  out.figNames = [figureName("circle", "pdf"), figureName("linear", "png")];
  const pdf = C.buildPdf({ data: new Uint8Array([1, 2, 3]), filter: "FlateDecode", imgW: 1, imgH: 1, pageW: 10, pageH: 10, title: "t" }); out.pdfHead = String.fromCharCode(...pdf.slice(0, 8));
  svg.call(zoom.transform, d3.zoomIdentity); drawCircle(); drawOverlay(); out.labels_k1 = +svg.node().dataset.nlabels;
  zoomToWindow(false); drawCircle(); drawOverlay(); out.zoomK = +Z.k.toFixed(2); out.labels_zoomed = +svg.node().dataset.nlabels; svg.call(zoom.transform, d3.zoomIdentity);
  out.ok = true; $("selftest").textContent = JSON.stringify(out);
}
window.GenomeMapViewer.current = { ds, st, setWin, setReps, doSearch, csvRows, circleSVG, figureSVG, saveFigure };
}
window.GenomeMapViewer = { start };
})();
