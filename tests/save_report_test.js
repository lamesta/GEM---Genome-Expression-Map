// Drive loader -> viewer -> "Save report", then reopen the saved report and run its self-test.
const fs = require("fs"); const path = require("path");
const { JSDOM, VirtualConsole } = require(path.join(__dirname, "jsdom_env/node_modules/jsdom"));
function open(html, hash, files, cb) {
  const errors = []; const vc = new VirtualConsole(); vc.on("jsdomError", e => errors.push(String(e.message || e)));
  let saved = null;
  const dom = new JSDOM(html, { runScripts: "dangerously", pretendToBeVisual: true, url: "http://localhost/v.html" + hash, virtualConsole: vc,
    beforeParse(w) { const noop = () => {}; const ctx = new Proxy({}, { get: (t, k) => (k in t ? t[k] : noop), set: (t, k, v) => { t[k] = v; return true; } });
      w.HTMLCanvasElement.prototype.getContext = () => ctx; w.URL.createObjectURL = b => { w.__lastBlob = b; return "blob:x"; };
      w.HTMLAnchorElement.prototype.click = function () { saved = { name: this.download, blob: w.__lastBlob }; }; } });
  const w = dom.window;
  w.addEventListener("load", async () => {
    if (files) { await w.GenomeMapLoader.addFiles(files.map(f => new w.File([fs.readFileSync(f, "utf8")], path.basename(f)))); w.document.getElementById("optSF").value = "zero_tolerant"; w.document.getElementById("buildBtn").click(); }
    setTimeout(async () => { cb({ w, errors, getSaved: () => saved }); }, 3500);
  });
}
const E = path.join(__dirname, "../example_R10.79/");
open(fs.readFileSync(path.join(__dirname, "../dist/genome_map_viewer.html"), "utf8"), "", ["R10.79_bakta.gff3", "featureCounts_counts_all.txt", "samplesheet.tsv", "contrasts.tsv", "landmarks.tsv", "labels.tsv"].map(f => E + f), async ({ w, errors, getSaved }) => {
  w.document.getElementById("htmlBtn").click();
  const s = getSaved(); const text = await new Promise(res => { const r = new w.FileReader(); r.onload = () => res(r.result); r.readAsText(s.blob); });
  fs.writeFileSync(path.join(__dirname, "out", s.name), text);
  console.log("saved", s.name, (text.length / 1e6).toFixed(2), "MB; errors in session 1:", errors);
  open(text, "#selftest", null, ({ w: w2, errors: e2 }) => { const st = w2.document.getElementById("selftest").textContent;
    console.log("reopened report: errors", e2, "loader visible:", w2.document.getElementById("loader").style.display, "selftest:", st.slice(0, 300)); process.exit(0); });
});
setTimeout(() => { console.log("timeout"); process.exit(1); }, 60000);
