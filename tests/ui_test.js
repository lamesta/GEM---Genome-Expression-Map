// Headless UI test with jsdom: runs the viewer self-test, exports circle SVG + linear SVG, optionally drives the loader with files.
const fs = require("fs"); const path = require("path");
const { JSDOM, VirtualConsole } = require(path.join(__dirname, "jsdom_env/node_modules/jsdom"));
const [, , htmlPath, outPrefix, ...loadFiles] = process.argv;
const errors = []; const vc = new VirtualConsole(); vc.on("jsdomError", e => errors.push(String(e.message || e))); vc.on("error", e => errors.push(String(e)));
const html = fs.readFileSync(htmlPath, "utf8");
const dom = new JSDOM(html, { runScripts: "dangerously", pretendToBeVisual: true, url: "http://localhost/viewer.html#selftest", virtualConsole: vc,
  beforeParse(w) {
    const noop = () => {}; const ctx = new Proxy({}, { get: (t, k) => (k in t ? t[k] : noop), set: (t, k, v) => { t[k] = v; return true; } });
    w.HTMLCanvasElement.prototype.getContext = () => ctx; w.HTMLCanvasElement.prototype.toBlob = noop;
    w.URL.createObjectURL = () => "blob:x"; w.SVGElement.prototype.getBBox = () => ({ x: 0, y: 0, width: 10, height: 10 });
  } });
const w = dom.window;
function finish() {
  const out = { errors };
  const stEl = w.document.getElementById("selftest"); out.selftest = stEl && stEl.textContent ? JSON.parse(stEl.textContent) : null;
  const V = w.GenomeMapViewer && w.GenomeMapViewer.current;
  if (V) { fs.writeFileSync(outPrefix + "_circle.svg", V.circleSVG());
    const lin = w.document.getElementById("linSvg"); fs.writeFileSync(outPrefix + "_linear.svg", `<svg xmlns="http://www.w3.org/2000/svg" width="${lin.getAttribute("width")}" height="${lin.getAttribute("height")}" font-family="Helvetica, Arial, sans-serif"><rect width="100%" height="100%" fill="#fff"/>${lin.innerHTML}</svg>`);
    out.title = w.document.getElementById("vTitle").textContent; out.banner = w.document.getElementById("banner").textContent; out.sub = w.document.getElementById("vSub").textContent;
    out.ringKey = w.document.getElementById("ringKey").textContent; out.notes = w.document.getElementById("notes").textContent.slice(0, 1500);
    out.datasetMeta = { mode: V.ds.meta.mode, sf: V.ds.meta.size_factors, warnings: V.ds.meta.warnings }; out.fc_sample = V.ds.log2fc.map(a => a.slice(0, 3)); }
  out.loaderMsgs = [...w.document.querySelectorAll("#msgs .msg")].map(e => e.className + ": " + e.textContent);
  fs.writeFileSync(outPrefix + "_result.json", JSON.stringify(out, null, 1)); console.log(JSON.stringify(out, null, 1).slice(0, 6000)); process.exit(0);
}
w.addEventListener("load", async () => {
  if (loadFiles.length) { // drive the loader
    const files = loadFiles.map(f => new w.File([fs.readFileSync(f, "utf8")], path.basename(f)));
    await w.GenomeMapLoader.addFiles(files);
    const opt = loadFiles.find(f => f.startsWith("--")); 
    if (process.env.GM_SF) w.document.getElementById("optSF").value = process.env.GM_SF;
    w.document.getElementById("buildBtn").click();
    setTimeout(finish, 4000);
  } else setTimeout(finish, 3000);
});
setTimeout(() => { errors.push("timeout"); finish(); }, 60000);
