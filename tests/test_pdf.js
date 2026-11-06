// Writes two one-page PDFs with the core's buildPdf: FlateDecode (RGB test pattern) and DCTDecode (JPEG given as argv[3]).
const fs = require("fs"), zlib = require("zlib"), C = require("../src/genomemap-core.js");
const W = 300, H = 200, rgb = new Uint8Array(W * H * 3);
for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) { const i = (y * W + x) * 3; rgb[i] = Math.round(255 * x / (W - 1)); rgb[i + 1] = Math.round(255 * y / (H - 1)); rgb[i + 2] = (x > 100 && x < 200 && y > 50 && y < 150) ? 0 : 200; }
fs.writeFileSync(process.argv[2] + "_rgb.raw", rgb);
fs.writeFileSync(process.argv[2] + "_flate.pdf", C.buildPdf({ data: zlib.deflateSync(rgb), filter: "FlateDecode", imgW: W, imgH: H, pageW: 225, pageH: 150, title: "test figure · log2FC ≥ 1 ü" }));
if (process.argv[3]) fs.writeFileSync(process.argv[2] + "_dct.pdf", C.buildPdf({ data: fs.readFileSync(process.argv[3]), filter: "DCTDecode", imgW: W, imgH: H, pageW: 225, pageH: 150, title: "jpeg test" }));
console.log("ok");
