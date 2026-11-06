// jsdom: open a report, call the figure exporter for both views (with legend), write the SVGs.
const fs = require("fs"), path = require("path");
const { JSDOM, VirtualConsole } = require(path.join(__dirname, "jsdom_env/node_modules/jsdom"));
const [, , htmlPath, outPrefix] = process.argv; const errors = []; const vc = new VirtualConsole(); vc.on("jsdomError", e => errors.push(String(e.message || e)));
const dom = new JSDOM(fs.readFileSync(htmlPath, "utf8"), { runScripts: "dangerously", pretendToBeVisual: true, url: "http://localhost/r.html#selftest", virtualConsole: vc,
  beforeParse(w) { const noop = () => {}; const ctx = new Proxy({}, { get: (t, k) => (k in t ? t[k] : noop), set: (t, k, v) => { t[k] = v; return true; } });
    w.HTMLCanvasElement.prototype.getContext = () => ctx; w.HTMLCanvasElement.prototype.toBlob = noop; w.URL.createObjectURL = () => "blob:x"; w.SVGElement.prototype.getBBox = () => ({ x: 0, y: 0, width: 10, height: 10 }); } });
dom.window.addEventListener("load", () => setTimeout(() => { const V = dom.window.GenomeMapViewer.current;
  for (const k of ["circle", "linear"]) { const f = V.figureSVG(k, true); fs.writeFileSync(`${outPrefix}_${k}_legend.svg`, f.svg); }
  const st = JSON.parse(dom.window.document.getElementById("selftest").textContent || "null");
  console.log(JSON.stringify({ errors, ok: st && st.ok, figCircle: st && st.figCircle, figLinear: st && st.figLinear, figNames: st && st.figNames, pdfHead: st && st.pdfHead })); process.exit(0); }, 300));
