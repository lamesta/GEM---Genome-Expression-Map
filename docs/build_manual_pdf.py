"""Build docs/USER_GUIDE.pdf from docs/USER_GUIDE.md.

Usage (from the package root):  python docs/build_manual_pdf.py
Requires: markdown, weasyprint (conda-forge). Uses the viewer's own font (src/fonts.css) and Menlo/DejaVu Sans Mono for code.
The Markdown file stays the single source of text; this script adds a title page, a table of contents with
page numbers, and two figures (docs/img/r1079_report_circle.png, docs/img/r1079_report_linear.png) in section 6.
"""
import re
from pathlib import Path

import markdown
from weasyprint import HTML

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
VERSION = "1.0.0"
DATE = "October 2026"

md_text = (DOCS / "USER_GUIDE.md").read_text(encoding="utf-8")
md = markdown.Markdown(extensions=["tables", "fenced_code", "toc", "sane_lists"], extension_configs={"toc": {"toc_depth": "2"}})
body = md.convert(md_text)

# drop the H1 (it becomes the title page) and the hand-written contents list (replaced by a paged TOC)
body = re.sub(r"<h1[^>]*>.*?</h1>\s*", "", body, count=1, flags=re.S)
body = re.sub(r"<p>Contents:</p>\s*<ol>.*?</ol>\s*", "<!--TOC-->", body, count=1, flags=re.S)

# paged table of contents from the H2 headings
h2 = re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body)
toc = '<nav class="toc"><h2 class="toc-title">Contents</h2><ul>' + "".join(
    f'<li><a href="#{i}">{t}</a></li>' for i, t in h2) + "</ul></nav>"
body = body.replace("<!--TOC-->", toc)

# figures for section 6, inserted before section 7
fig = f"""
<figure class="fig circle">
  <img src="{(DOCS / 'img' / 'r1079_report_circle.png').as_uri()}" alt="Circular map of P. syringae R10.79">
  <figcaption><b>Figure 1.</b> Circular map of the example data set (<i>P. syringae</i> R10.79, <code>example_R10.79/</code>, Mode A).
  From outside in: ruler; landmarks (<i>ina</i>, T3SS <i>hrp/hrc</i> cluster, rRNA operons 1 to 5); ring S, genes by strand;
  ring C, COG category; rings 1 to 4, the test contrasts (1 AER vs MC EXP3, 2 AER vs MC EXP4, 3 AER + 1 h NaCl vs MC EXP4,
  4 bubbled 5 h vs MC Tsh); ring N, the noise reference (MC EXP4 vs MC EXP3). Orange bars point outward (higher in the
  treatment), blue bars inward (lower). The blue wedge marks the window open in the linear view (Figure 2). Every contrast
  here compares one library with one library, so the values are descriptive ratios without p-values.</figcaption>
</figure>
<figure class="fig linear">
  <img src="{(DOCS / 'img' / 'r1079_report_linear.png').as_uri()}" alt="Linear view around ina">
  <figcaption><b>Figure 2.</b> Linear view of the 34-kb window around <i>ina</i> (MGGFHP_01089, outlined), with one track per
  selected test contrast on a ±3 log2 scale. Colour intensity follows the size of the log2 fold change. Hovering a gene
  shows its locus tag, product, every contrast's value and the normalised counts per library.</figcaption>
</figure>
"""
body = re.sub(r'(<h2 id="7-saving-and-sharing">)', fig + r"\1", body, count=1)
assert "Figure 1." in body and "Figure 2." in body and 'class="toc"' in body

fonts_css = (ROOT / "src" / "fonts.css").read_text(encoding="utf-8")
css = fonts_css + """
@page { size: A4; margin: 22mm 20mm 20mm 20mm;
  @top-right { content: "genomemap """ + VERSION + """ · User guide"; font: 8pt 'Open Sauce Sans'; color: #7a7f8a; }
  @bottom-center { content: "page " counter(page) " of " counter(pages); font: 8pt 'Open Sauce Sans'; color: #7a7f8a; } }
@page :first { @top-right { content: none; } @bottom-center { content: none; } }
html { font-family: 'Open Sauce Sans', 'Helvetica Neue', Arial, sans-serif; font-size: 9.8pt; line-height: 1.45; color: #1f2329; }
.title-page { page: first; height: 230mm; display: flex; flex-direction: column; justify-content: center; break-after: page; }
.title-page .kicker { color: #9C4B22; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase; font-size: 10pt; }
.title-page h1 { font-size: 34pt; font-weight: 700; color: #2F4682; margin: 4mm 0 2mm 0; line-height: 1.1; }
.title-page .sub { font-size: 14pt; color: #3a3f48; margin: 0 0 14mm 0; }
.title-page .meta { font-size: 10pt; color: #5b616c; border-top: 1.2pt solid #2F4682; padding-top: 4mm; width: 120mm; }
.title-page img { width: 92mm; margin: 0 0 10mm 0; }
h2 { font-size: 15pt; font-weight: 700; color: #2F4682; margin: 9mm 0 3mm 0; break-after: avoid; border-bottom: 0.6pt solid #c9cfdb; padding-bottom: 1.5mm; }
h3 { font-size: 11.5pt; font-weight: 600; color: #2F4682; margin: 6mm 0 2mm 0; break-after: avoid; }
p { margin: 0 0 2.6mm 0; } ul, ol { margin: 0 0 2.6mm 0; padding-left: 6mm; } li { margin: 0.6mm 0; }
a { color: #2F4682; text-decoration: none; }
code { font-family: Menlo, 'DejaVu Sans Mono', monospace; font-size: 8.4pt; background: #f1f3f6; padding: 0 1.2pt; border-radius: 1.5pt; }
pre { font-family: Menlo, 'DejaVu Sans Mono', monospace; font-size: 7.9pt; line-height: 1.4; background: #f5f6f8; border-left: 2pt solid #2F4682;
  padding: 2.5mm 3mm; white-space: pre-wrap; break-inside: avoid; margin: 0 0 3mm 0; }
pre code { background: none; padding: 0; font-size: inherit; }
table { border-collapse: collapse; width: 100%; margin: 1mm 0 4mm 0; font-size: 8.5pt; line-height: 1.35; }
th { background: #e9edf5; color: #1f2329; font-weight: 600; text-align: left; }
th, td { border: 0.5pt solid #c9cfdb; padding: 1.4mm 2mm; vertical-align: top; }
tr { break-inside: avoid; }
strong, b { font-weight: 700; }
.toc { break-after: page; } .toc-title { border: none; }
.toc ul { list-style: none; padding: 0; } .toc li { margin: 1.8mm 0; font-size: 10.5pt; }
.toc a { color: #1f2329; } .toc a::after { content: leader('.') target-counter(attr(href), page); color: #5b616c; }
figure.fig { margin: 5mm 0; break-inside: avoid; text-align: center; }
figure.circle img { width: 135mm; } figure.linear img { width: 170mm; }
figcaption { text-align: left; font-size: 8.4pt; line-height: 1.4; color: #3a3f48; margin-top: 2mm; }
"""

title = f"""
<section class="title-page">
  <div class="kicker">Bacterial RNA-seq · genome maps</div>
  <h1>genomemap<br>User guide</h1>
  <div class="sub">Circular and linear maps of transcriptional responses for any annotated bacterial genome</div>
  <img src="{(DOCS / 'img' / 'r1079_report_circle.png').as_uri()}" alt="">
  <div class="meta">Version {VERSION} · {DATE}<br>Viewer: <code>dist/genome_map_viewer.html</code> · Text source: <code>docs/USER_GUIDE.md</code></div>
</section>
"""
html = f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><title>genomemap {VERSION} user guide</title><style>{css}</style></head><body>{title}{body}</body></html>"
out = DOCS / "USER_GUIDE.pdf"
HTML(string=html, base_url=str(DOCS)).write_pdf(out)
print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
