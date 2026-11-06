"""pytest suite for the genome-map CLI builder (cli/genomemap_build.py).

Run from the package root:  python -m pytest -q tests/test_cli.py
Requires: python >= 3.9, numpy, pandas, pytest, and node >= 18 on PATH.

The Mode-A arithmetic is re-implemented here independently in numpy, so these
tests check the JavaScript core against a second implementation, not against itself.
"""
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "cli" / "genomemap_build.py"
SYN = ROOT / "tests" / "synthetic_DC3000"
FIX = ROOT / "tests" / "fixtures"
EX = ROOT / "example_R10.79"


def run_build(tmp_path, *args, expect_ok=True, name="out"):
    out = tmp_path / f"{name}.html"
    js = tmp_path / f"{name}.json"
    cmd = [sys.executable, str(CLI), "build", *map(str, args), "--dataset-json", str(js), "-o", str(out)]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if expect_ok:
        assert p.returncode == 0, p.stdout + p.stderr
        return json.loads(js.read_text()), out.read_text(encoding="utf8"), p
    assert p.returncode != 0, "build was expected to fail"
    return None, None, p


def col(ds, cid, key="log2fc"):
    i = [c["id"] for c in ds["contrasts"]].index(cid)
    return pd.Series(ds[key][i], index=ds["features"]["id"], dtype="float64")


# ---------- independent numpy implementation of Mode A ----------
def mode_a_reference(counts, ann_ids, ann_types, treat, ctrl, sf_variant="standard", cpm=1.0, min_libs=2):
    show = [i for i, t in zip(ann_ids, ann_types) if t in ("CDS", "ncRNA", "tmRNA")]
    cds = [i for i, t in zip(ann_ids, ann_types) if t == "CDS"]
    M = counts.reindex(show)
    K = counts.reindex(cds).dropna()
    logs = np.log(K.where(K > 0))
    if sf_variant == "standard":
        keep = (K > 0).all(axis=1)
        lg = logs[keep]
        sf = np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))
    else:  # zero-tolerant: zeros treated as missing
        gm = logs.mean(axis=1, skipna=True)
        sf = np.exp(logs.sub(gm, axis=0).median(axis=0, skipna=True))
    norm = M / sf
    cpm_m = M / K.sum(axis=0) * 1e6
    expressed = (cpm_m > cpm).sum(axis=1) >= min_libs
    l2 = np.log2((norm[treat].mean(axis=1) + 1) / (norm[ctrl].mean(axis=1) + 1))
    return sf, l2, expressed


@pytest.fixture(scope="module")
def dc_ann():
    rows = []
    for l in (FIX / "DC3000_NCBI.gff").read_text().splitlines():
        if l.startswith("#") or not l.strip():
            continue
        f = l.split("\t")
        if len(f) < 9:
            continue
        rows.append(f)
    return rows


# ---------- tests ----------
def test_empty_viewer(tmp_path):
    out = tmp_path / "v.html"
    p = subprocess.run([sys.executable, str(CLI), "viewer", "-o", str(out)], cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    html = out.read_text(encoding="utf8")
    assert 'id="gm-dataset"' not in html
    assert "GenomeMapViewer" in html and "d3" in html
    assert out.stat().st_size < 2_000_000


def test_r1079_example_matches_committed_dataset(tmp_path):
    ds, html, _ = run_build(tmp_path, "--annotation", EX / "R10.79_bakta.gff3", "--counts", EX / "featureCounts_counts_all.txt",
                            "--samplesheet", EX / "samplesheet.tsv", "--contrasts", EX / "contrasts.tsv",
                            "--landmarks", EX / "landmarks.tsv", "--labels", EX / "labels.tsv")
    ref = json.loads((EX / "R10.79_dataset.json").read_text())
    assert 'id="gm-dataset"' in html
    assert ds["features"]["id"] == ref["features"]["id"]
    assert len(ds["contrasts"]) == 5
    for a, b in zip(ds["log2fc"], ref["log2fc"]):
        assert a == b
    assert ds["meta"]["size_factor_variant"] == "standard"
    assert any("no replication" in w for w in ds["meta"]["warnings"])


@pytest.mark.parametrize("variant", ["standard", "zero_tolerant"])
def test_mode_a_against_numpy(tmp_path, dc_ann, variant):
    ds, _, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", SYN / "SYNTHETIC_counts_matrix.tsv",
                         "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", "--sf", variant)
    counts = pd.read_csv(SYN / "SYNTHETIC_counts_matrix.tsv", sep="\t", index_col=0)
    ids, types = ds["features"]["id"], ds["features"]["type"]
    sf, l2, expressed = mode_a_reference(counts, ids, types, ["AER_1", "AER_2", "AER_3"], ["MC_1", "MC_2", "MC_3"], variant)
    sf_js = pd.Series(ds["meta"]["size_factors"])
    assert np.allclose(sf_js[sf.index].values, sf.values, rtol=1e-9)
    js = col(ds, "AER_vs_MC")
    m = js.notna() & l2.reindex(js.index).notna()
    assert m.sum() > 5800
    assert np.max(np.abs(js[m] - l2.reindex(js.index)[m])) <= 5e-5  # JS rounds to 4 decimals
    exp_js = pd.Series(ds["expressed"], index=ids)
    assert (exp_js == expressed.reindex(ids).fillna(False)).all()


def test_planted_signal_recovered(tmp_path):
    ds, _, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", SYN / "SYNTHETIC_counts_matrix.tsv",
                         "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv")
    js = col(ds, "AER_vs_MC")
    reps = [r["id"] for r in ds["replicons"]]
    rep = pd.Series([reps[i] for i in ds["features"]["rep"]], index=js.index)
    s = pd.Series(ds["features"]["s"], index=js.index)
    assert abs(js[rep == "NC_004633.1"].median() - 2.0) < 0.3
    region = (rep == "NC_004578.1") & (s >= 1_000_000) & (s < 1_100_000)
    assert abs(js[region].median() + 1.5) < 0.3
    assert abs(js[(rep == "NC_004578.1") & ~region].median()) < 0.1


def test_gzip_inputs_identical(tmp_path):
    for src in [FIX / "DC3000_NCBI.gff", SYN / "SYNTHETIC_counts_matrix.tsv"]:
        with open(src, "rb") as fi, gzip.open(tmp_path / (src.name + ".gz"), "wb") as fo:
            shutil.copyfileobj(fi, fo)
    a, _, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", SYN / "SYNTHETIC_counts_matrix.tsv",
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", name="plain")
    b, _, _ = run_build(tmp_path, "--annotation", tmp_path / "DC3000_NCBI.gff.gz", "--counts", tmp_path / "SYNTHETIC_counts_matrix.tsv.gz",
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", name="gz")
    assert a["log2fc"] == b["log2fc"] and a["features"] == b["features"]


def test_htseq_files_identical_to_matrix(tmp_path):
    counts = pd.read_csv(SYN / "SYNTHETIC_counts_matrix.tsv", sep="\t", index_col=0)
    files = []
    for s in counts.columns:
        f = tmp_path / f"{s}.htseq.txt"
        body = "\n".join(f"{g}\t{int(v)}" for g, v in counts[s].items())
        f.write_text(body + "\n__no_feature\t1234\n__ambiguous\t56\n__too_low_aQual\t0\n__not_aligned\t0\n__alignment_not_unique\t789\n")
        files.append(f)
    a, _, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", SYN / "SYNTHETIC_counts_matrix.tsv",
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", name="matrix")
    b, _, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", *files,
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", name="htseq")
    assert a["log2fc"] == b["log2fc"]
    assert a["meta"]["size_factors"] == b["meta"]["size_factors"]


def test_mode_b_passthrough_and_genbank(tmp_path):
    ds, html, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gbff", "--contrasts", SYN / "contrasts_modeB.tsv",
                            "--de-dir", SYN, "--alpha", "0.05")
    de = pd.read_csv(SYN / "SYNTHETIC_DESeq2_AER_vs_MC.csv", index_col=0)
    js, pj = col(ds, "AER_vs_MC"), col(ds, "AER_vs_MC", "padj")
    common = js.index.intersection(de.index)
    assert len(common) == len(de)
    assert np.nanmax(np.abs(js[common] - de.loc[common, "log2FoldChange"])) <= 5e-5
    pj = pj.reindex(de.index)
    m = de["padj"].notna()
    assert pj[~m].isna().all()  # DESeq2 NA padj (independent filtering) stays missing
    assert np.allclose(pj[m], de.loc[m, "padj"], rtol=5.01e-3, atol=1e-300)  # stored with 3 significant digits (max rel. rounding error 0.5%)
    assert ds["meta"]["mode"] == "B" and ds["meta"]["alpha"] == 0.05
    assert [r["name"] for r in ds["replicons"]][1:] == ["pDC3000A", "pDC3000B"]


def test_unmatched_ids_fail_clearly(tmp_path):
    counts = pd.read_csv(SYN / "SYNTHETIC_counts_matrix.tsv", sep="\t", index_col=0)
    counts.index = ["X_" + i for i in counts.index]
    bad = tmp_path / "bad_ids.tsv"
    counts.to_csv(bad, sep="\t")
    _, _, p = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", bad,
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", expect_ok=False)
    assert "have a row in the count table" in (p.stdout + p.stderr)


def test_non_integer_counts_rejected(tmp_path):
    counts = pd.read_csv(SYN / "SYNTHETIC_counts_matrix.tsv", sep="\t", index_col=0).astype(float)
    counts.iloc[0, 0] += 0.37
    bad = tmp_path / "tpm_like.tsv"
    counts.to_csv(bad, sep="\t")
    _, _, p = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", bad,
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", expect_ok=False)
    assert "non-integer" in (p.stdout + p.stderr)


def test_unknown_sample_in_contrast_rejected(tmp_path):
    c = tmp_path / "contrasts_bad.tsv"
    c.write_text("contrast_id\ttreatment\tcontrol\trole\nX\tAER_9\tMC_1\ttest\n")
    _, _, p = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", SYN / "SYNTHETIC_counts_matrix.tsv",
                        "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", c, expect_ok=False)
    out = p.stdout + p.stderr
    assert "AER_9" in out


def test_prokka_style_gff(tmp_path):
    p = subprocess.run(["node", "-e", f"""
const C=require({json.dumps(str(ROOT / 'src' / 'genomemap-core.js'))});
const fs=require('fs'); const a=C.parseGFF3(fs.readFileSync({json.dumps(str(FIX / 'prokka_style_synthetic.gff'))},'utf8'));
console.log(JSON.stringify({{n:a.features.length, reps:a.replicons.length, ids:a.features.slice(0,3).map(f=>f.id)}}));"""],
                       cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    r = json.loads(p.stdout)
    assert r["n"] > 0 and r["reps"] >= 1


def test_detect_types():
    files = {SYN / "htseq" / "MC_1.htseq.txt": "htseq", SYN / "SYNTHETIC_counts_matrix.tsv": "matrix",
             SYN / "samplesheet.tsv": "samplesheet", SYN / "contrasts.tsv": "contrasts", SYN / "landmarks.tsv": "landmarks",
             SYN / "SYNTHETIC_DESeq2_AER_vs_MC.csv": "de_table", EX / "featureCounts_counts_all.txt": "featurecounts",
             EX / "R10.79_bakta.gff3": "gff3", FIX / "DC3000_NCBI.gbff": "genbank", EX / "labels.tsv": "labels"}
    js = "const C=require(%s),fs=require('fs');const F=%s;console.log(JSON.stringify(F.map(f=>C.detectType(f,fs.readFileSync(f,'utf8')))));" % (
        json.dumps(str(ROOT / "src" / "genomemap-core.js")), json.dumps([str(f) for f in files]))
    p = subprocess.run(["node", "-e", js], cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert json.loads(p.stdout) == list(files.values())


def test_landmark_by_feature_id_detected_and_resolved(tmp_path):
    lm = tmp_path / "landmarks_fid.tsv"
    lm.write_text("name\tfeature_id\ttype\tlabel\nfirst CDS\tPSPTO_RS00005\tgene\tRS00005\n")
    js = "const C=require(%s),fs=require('fs');console.log(C.detectType('x',fs.readFileSync(%s,'utf8')));" % (
        json.dumps(str(ROOT / "src" / "genomemap-core.js")), json.dumps(str(lm)))
    p = subprocess.run(["node", "-e", js], cwd=ROOT, capture_output=True, text=True)
    assert p.stdout.strip() == "landmarks"
    ds, _, _ = run_build(tmp_path, "--annotation", FIX / "DC3000_NCBI.gff", "--counts", SYN / "SYNTHETIC_counts_matrix.tsv",
                         "--samplesheet", SYN / "samplesheet.tsv", "--contrasts", SYN / "contrasts.tsv", "--landmarks", lm)
    i = ds["features"]["id"].index("PSPTO_RS00005")
    L = ds["landmarks"][0]
    assert (L["s"], L["e"], L["rep"]) == (ds["features"]["s"][i], ds["features"]["e"][i], ds["features"]["rep"][i])


def test_limma_and_edger_rowname_tables():
    limma = '"","logFC","AveExpr","t","P.Value","adj.P.Val","B"\n"G1",1.5,7.2,5.1,1e-5,0.002,3.1\n"G2",-0.3,4.4,-0.8,0.4,0.9,-5\n'
    edger = '"","logFC","logCPM","F","PValue","FDR"\n"G1",-2.1,5.5,30,1e-6,1e-4\n'
    js = "const C=require(%s);console.log(JSON.stringify([C.parseDETable(%s),C.parseDETable(%s)]));" % (
        json.dumps(str(ROOT / "src" / "genomemap-core.js")), json.dumps(limma), json.dumps(edger))
    p = subprocess.run(["node", "-e", js], cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    L, E = json.loads(p.stdout)
    assert L["genes"] == ["G1", "G2"] and L["log2fc"] == [1.5, -0.3] and L["padj"] == [0.002, 0.9] and L["basemean"] == [7.2, 4.4]
    assert L["columns"]["padj"] == "adj.P.Val" and L["columns"]["abundance"] == "AveExpr"
    assert E["genes"] == ["G1"] and E["log2fc"] == [-2.1] and E["padj"] == [1e-4] and E["columns"]["abundance"] == "logCPM"


def test_pdf_writer(tmp_path):
    pdfium = pytest.importorskip("pypdfium2")
    Image = pytest.importorskip("PIL.Image")
    pre = tmp_path / "p"
    p = subprocess.run(["node", str(ROOT / "tests" / "test_pdf.js"), str(pre)], cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    rgb = np.frombuffer((tmp_path / "p_rgb.raw").read_bytes(), dtype=np.uint8).reshape(200, 300, 3)
    Image.fromarray(rgb).save(tmp_path / "p.jpg", quality=95)
    p = subprocess.run(["node", str(ROOT / "tests" / "test_pdf.js"), str(pre), str(tmp_path / "p.jpg")], cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    for name, tol in [("p_flate.pdf", 3.0), ("p_dct.pdf", 8.0)]:
        doc = pdfium.PdfDocument(str(tmp_path / name))
        assert len(doc) == 1 and tuple(round(v) for v in doc[0].get_size()) == (225, 150)
        img = np.asarray(doc[0].render(scale=300 / 225).to_pil().convert("RGB").resize((300, 200)), dtype=float)
        assert np.abs(img - rgb).mean() < tol, name
    assert "test figure" in pdfium.PdfDocument(str(tmp_path / "p_flate.pdf")).get_metadata_dict().get("Title", "")
