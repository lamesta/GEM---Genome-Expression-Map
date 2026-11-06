#!/usr/bin/env python3
"""genomemap command-line builder.

  python genomemap_build.py viewer -o genome_map_viewer.html
      writes the empty, self-contained viewer (drag & drop loader).

  python genomemap_build.py build --annotation ann.gff3 --counts counts.txt --samplesheet samplesheet.tsv \\
      --contrasts contrasts.tsv [--landmarks landmarks.tsv] [--labels labels.tsv] -o report.html
      Mode A: computes the dataset with the same JavaScript core as the viewer (needs Node.js >= 18)
      and writes a shareable HTML report with the dataset embedded.

  python genomemap_build.py build --annotation ann.gff3 --contrasts contrasts.tsv --de-dir de_tables/ -o report.html
      Mode B: DE result tables named in the de_file column of contrasts.tsv.

Without Node.js, open the viewer, load the files, and use "Save report (.html)" in the browser:
the output is identical.
"""
import argparse, json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, "..", "src")
VERSION = "1.0.0"


def read(p): return open(os.path.join(SRC, p), encoding="utf-8").read()


def assemble(dataset_json=None, title="genomemap: genome map of transcriptional responses"):
    """Assemble the self-contained HTML. Order matters: loader.js must be the LAST script (it captures the page source)."""
    emb = ""
    if dataset_json is not None:
        emb = '<script id="gm-dataset" type="application/json">' + dataset_json.replace("</", "<\\/") + "</script>\n"
    html = ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"generator\" content=\"genomemap " + VERSION + "\">\n"
            "<title>" + title + "</title>\n<style>\n" + read("fonts.css") + "\n" + read("style.css") + "\n</style>\n</head>\n"
            + read("body.html") + "\n" + emb
            + "<script>\n" + read("d3.min.js") + "\n</script>\n"
            + "<script>\n" + read("genomemap-core.js") + "\n</script>\n"
            + "<script>\n" + read("viewer.js") + "\n</script>\n"
            + "<script>\n" + read("loader.js") + "\n</script>\n</body>\n</html>\n")
    return html


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("viewer", help="write the empty viewer"); v.add_argument("-o", "--out", default="genome_map_viewer.html")
    b = sub.add_parser("build", help="build a report with an embedded dataset")
    b.add_argument("--annotation", required=True, help="GFF3 or GenBank (.gz accepted)")
    b.add_argument("--counts", nargs="+", help="featureCounts output, count matrix, or several HTSeq-count files (Mode A)")
    b.add_argument("--samplesheet"); b.add_argument("--contrasts", required=True)
    b.add_argument("--de-dir", help="folder with DE tables named in contrasts.tsv de_file (Mode B)")
    b.add_argument("--landmarks"); b.add_argument("--labels")
    b.add_argument("--types", default="CDS,ncRNA,tmRNA", help="feature types shown (comma-separated)")
    b.add_argument("--sf", choices=["standard", "zero_tolerant"], default="standard", help="size-factor variant (Mode A)")
    b.add_argument("--sf-types", default="CDS", help="feature types used for size factors and CPM totals")
    b.add_argument("--cpm", type=float, default=1.0); b.add_argument("--min-libs", type=int, default=2)
    b.add_argument("--alpha", type=float, default=0.05); b.add_argument("--bin-width", type=int, default=10000)
    b.add_argument("--title", default=""); b.add_argument("--strain", default="")
    b.add_argument("--dataset-json", help="also write the dataset JSON here")
    b.add_argument("--node", default=shutil.which("node") or "node"); b.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    if a.cmd == "viewer":
        open(a.out, "w", encoding="utf-8").write(assemble()); print(f"wrote {a.out} ({os.path.getsize(a.out)/1e6:.2f} MB)"); return
    if not (a.counts or a.de_dir): ap.error("give --counts (Mode A) or --de-dir (Mode B)")
    if a.counts and not a.samplesheet: print("note: no --samplesheet; each count column becomes its own library", file=sys.stderr)
    cfg = dict(annotation=os.path.abspath(a.annotation), counts=[os.path.abspath(x) for x in a.counts] if a.counts else None,
               samplesheet=os.path.abspath(a.samplesheet) if a.samplesheet else None, contrasts=os.path.abspath(a.contrasts),
               de_dir=os.path.abspath(a.de_dir) if a.de_dir else None, landmarks=os.path.abspath(a.landmarks) if a.landmarks else None,
               labels=os.path.abspath(a.labels) if a.labels else None, types=a.types.split(","), sf=a.sf, sf_types=a.sf_types.split(","),
               cpm=a.cpm, min_libs=a.min_libs, alpha=a.alpha, bin_width=a.bin_width, title=a.title, strain=a.strain)
    cfg = {k: v for k, v in cfg.items() if v is not None}
    with tempfile.TemporaryDirectory() as td:
        cp, dp = os.path.join(td, "cfg.json"), os.path.join(td, "ds.json"); json.dump(cfg, open(cp, "w"))
        try: r = subprocess.run([a.node, os.path.join(HERE, "build_dataset.js"), cp, dp], capture_output=True, text=True)
        except FileNotFoundError: sys.exit("Node.js not found (install it, e.g. `conda install -c conda-forge nodejs`), or use the browser viewer and 'Save report (.html)'.")
        res = json.loads(r.stdout.strip().splitlines()[-1]) if r.stdout.strip() else {"ok": False, "error": r.stderr}
        if not res.get("ok"): sys.exit("error: " + res.get("error", "unknown"))
        dsj = open(dp, encoding="utf-8").read()
    if a.dataset_json: open(a.dataset_json, "w", encoding="utf-8").write(dsj)
    title = (a.title or json.loads(dsj).get("title") or "genomemap") + ": genome map"
    open(a.out, "w", encoding="utf-8").write(assemble(dsj, title))
    for w in res.get("warnings", []): print("warning:", w, file=sys.stderr)
    print(f"wrote {a.out} ({os.path.getsize(a.out)/1e6:.2f} MB): {res['features']} features, {res['contrasts']} contrasts")


if __name__ == "__main__":
    main()
