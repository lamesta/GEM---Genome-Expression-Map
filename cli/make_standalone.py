"""Write dist/genomemap.py: ONE self-contained Python file that carries the whole viewer source
(src/*, cli/build_dataset.js) compressed inside it. Same commands as cli/genomemap_build.py, plus `extract`.

Usage (from the package root):  python cli/make_standalone.py
"""
import base64, re, zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = {f"src/{n}": ROOT / "src" / n for n in ["fonts.css", "style.css", "body.html", "d3.min.js", "genomemap-core.js", "viewer.js", "loader.js"]}
FILES["cli/build_dataset.js"] = ROOT / "cli" / "build_dataset.js"

src = (ROOT / "cli" / "genomemap_build.py").read_text(encoding="utf-8")
blobs = "\n".join(f'    "{k}": """\n{base64.encodebytes(zlib.compress(p.read_bytes(), 9)).decode()}""",' for k, p in FILES.items())

head = '''#!/usr/bin/env python3
"""genomemap 1.0.0, single-file edition (viewer + command-line builder in one file).

  python genomemap.py viewer -o genome_map_viewer.html
      writes the empty, self-contained viewer (open it in a browser and drag in your files).

  python genomemap.py build --annotation ann.gff3 --counts counts.txt --samplesheet samplesheet.tsv \\\\
      --contrasts contrasts.tsv [--landmarks landmarks.tsv] [--labels labels.tsv] -o report.html
      Mode A. Needs Node.js >= 18 (the computation runs in the same JavaScript core as the viewer).

  python genomemap.py build --annotation ann.gff3 --contrasts contrasts.tsv --de-dir de_tables/ -o report.html
      Mode B: DE result tables named in the de_file column of contrasts.tsv.

  python genomemap.py extract DIR
      writes the embedded source files (src/*.js, *.css, body.html, cli/build_dataset.js) to DIR.

Only the Python standard library is needed. Without Node.js, use `viewer`, open the HTML,
load the files and press "Save report (.html)": the result is identical.
The source is stored below as zlib-compressed base64 text.
"""
import argparse, base64, json, os, shutil, subprocess, sys, tempfile, zlib

_EMBED = {
''' + blobs + '''
}


def _get(name): return zlib.decompress(base64.b64decode(_EMBED[name])).decode("utf-8")
'''

body = src.split("import argparse, json, os, shutil, subprocess, sys, tempfile\n", 1)[1]
body = re.sub(r'HERE = .*?\n', "", body, count=1)
body = body.replace('def read(p): return open(os.path.join(SRC, p), encoding="utf-8").read()', 'def read(p): return _get("src/" + p)')
# node entry: write core + build_dataset.js into the temp dir and run it there
old_run = 'try: r = subprocess.run([a.node, os.path.join(HERE, "build_dataset.js"), cp, dp], capture_output=True, text=True)'
assert old_run in body
new_run = ('open(os.path.join(td, "package.json"), "w").write(\'{"type": "commonjs"}\')  # stops Node searching parent folders\n'
           '        open(os.path.join(td, "genomemap-core.js"), "w", encoding="utf-8").write(_get("src/genomemap-core.js"))\n'
           '        open(os.path.join(td, "build_dataset.js"), "w", encoding="utf-8").write(_get("cli/build_dataset.js").replace(\'path.join(__dirname, "..", "src", "genomemap-core.js")\', \'path.join(__dirname, "genomemap-core.js")\'))\n'
           '        try: r = subprocess.run([a.node, os.path.join(td, "build_dataset.js"), cp, dp], capture_output=True, text=True)')
body = body.replace(old_run, new_run)
# extract command
old_sub = 'b = sub.add_parser("build", help="build a report with an embedded dataset")'
assert old_sub in body
body = body.replace(old_sub, 'x = sub.add_parser("extract", help="write the embedded source files to a folder"); x.add_argument("dir")\n    ' + old_sub)
old_v = '    if a.cmd == "viewer":'
body = body.replace(old_v, '''    if a.cmd == "extract":
        for k in _EMBED:
            p = os.path.join(a.dir, *k.split("/")); os.makedirs(os.path.dirname(p), exist_ok=True)
            open(p, "w", encoding="utf-8").write(_get(k)); print("wrote", p)
        return
''' + old_v, 1)
assert 'path.join(__dirname, "..", "src", "genomemap-core.js")' in FILES["cli/build_dataset.js"].read_text()
out = ROOT / "dist" / "genomemap.py"
out.write_text(head + body, encoding="utf-8")
print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
