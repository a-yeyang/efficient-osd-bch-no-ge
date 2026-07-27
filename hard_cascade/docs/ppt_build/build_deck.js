// Build the SIMPLIFIED 5-page Chinese RS+BCH cascade deck with pptxgenjs.
//
// Design rules (per advisor + user feedback):
//   - Every schematic / flow / principle diagram is built from NATIVE PPT
//     shapes + text boxes so it stays fully editable in PowerPoint.
//   - ONLY the simulation-result figures (BER curve, latency-KPI bars) are
//     kept as images — those are real experiment output, cannot be redrawn
//     with shapes, and the advisor explicitly said "just paste the AI report's
//     result figures in".
//   - 5 pages total: 封面/背景 → 原理与流程 → 创新点 → 结果与分析 → 总结.
//
// Palette: white-dominant, light-green bordered — Beijing Institute of
// Technology official emblem colors (VI A1-08: PANTONE 347C / 349C / 1535C).

const pptxgen = require("pptxgenjs");
const path = require("path");

const ROOT = "/Users/chenshiyang.10/workspace/llosd_reproduction";
const RSFIG = path.join(ROOT, "rs_bch_cascade/figures");
const OUT = path.join(ROOT, "hard_cascade/docs/rs_bch_cascade_deck.pptx");

// ---- palette (BIT official VI colors) ----
const GREEN      = "009A44"; // PANTONE 347C — primary school green
const GREEN_75   = "14AE68"; // 75% tint — softer light-green border
const GREEN_50   = "89C997"; // 50% tint
const GREEN_DARK = "015C31"; // PANTONE 349C — dark green (titles, strong text)
const BROWN      = "A23E0A"; // PANTONE 1535C — brown (caution / error accent)
const INK  = "1A1A1A";
const GRAY = "5C6670";
const GRAYL = "8A929B";
const WHITE = "FFFFFF";
const GREEN_WASH = "EEF8F0";  // very light green card background
const GREEN_TINT = "DCF0E1";  // stronger green card background
const BROWN_TINT = "FBEEE6";  // light brown card background (caveat)
const GRAY_TINT = "F2F3F4";   // neutral (de-emphasized / traditional)

const FONT = "Microsoft YaHei";

const pres = new pptxgen();
pres.defineLayout({ name: "W", width: 13.333, height: 7.5 });
pres.layout = "W";
pres.author = "陈诗阳";
pres.title = "RS+BCH 级联码低时延译码方案";

const W = 13.333, H = 7.5, M = 0.6;
const AR_BER = 1120 / 700;   // n255_scheme_a_ber
const AR_KPI = 1680 / 630;   // n255_kpi_bars

// ---------------------------------------------------------------- helpers
function title(s, txt, opts) {
  opts = opts || {};
  s.addText(txt, {
    x: M, y: 0.34, w: W - 2 * M, h: 0.9, fontFace: FONT, fontSize: opts.size || 30,
    bold: true, color: opts.color || GREEN_DARK, align: "left", valign: "middle", margin: 0,
  });
}
function kicker(s, txt, color) {
  s.addText(txt, {
    x: M, y: 0.06, w: W - 2 * M, h: 0.32, fontFace: FONT, fontSize: 12.5, bold: true,
    color: color || GREEN, align: "left", valign: "middle", margin: 0, charSpacing: 2,
  });
}
function topBand(s, color) {
  s.addShape(pres.ShapeType.rect, {
    x: 0, y: 0, w: W, h: 0.14, fill: { color: color || GREEN }, line: { type: "none" },
  });
}
function card(s, x, y, w, h, fill, line) {
  s.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.09, fill: { color: fill },
    line: line === undefined ? { color: GREEN_75, width: 1.25 }
      : (line ? { color: line, width: 1.25 } : { type: "none" }),
  });
}
function badge(s, x, y, d, r, fill) {
  r = r || 0.32;
  s.addShape(pres.ShapeType.ellipse, {
    x, y, w: r * 2, h: r * 2, fill: { color: fill || GREEN }, line: { type: "none" },
  });
  s.addText(String(d), {
    x, y, w: r * 2, h: r * 2, fontFace: FONT, fontSize: 15, bold: true,
    color: WHITE, align: "center", valign: "middle", margin: 0,
  });
}
// a labeled process box (white fill, colored border + colored bold title text)
function pbox(s, x, y, w, h, txt, color, opts) {
  opts = opts || {};
  s.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.05, fill: { color: opts.fill || WHITE },
    line: { color: color, width: opts.lw || 1.75 },
  });
  s.addText(txt, {
    x: x + 0.03, y, w: w - 0.06, h, fontFace: FONT, fontSize: opts.fs || 12.5,
    bold: opts.bold !== false, color: opts.tc || color, align: "center",
    valign: "middle", margin: 0, lineSpacingMultiple: 0.98,
  });
}
// straight arrow between two points
function arrow(s, x1, y1, x2, y2, color, w) {
  s.addShape(pres.ShapeType.line, {
    x: Math.min(x1, x2), y: Math.min(y1, y2),
    w: Math.abs(x2 - x1), h: Math.abs(y2 - y1),
    line: { color: color || INK, width: w || 1.75, endArrowType: "triangle" },
    flipH: x2 < x1, flipV: y2 < y1,
  });
}

// ============================================================ SLIDE 1 — 封面 / 背景
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  topBand(s, GREEN);

  s.addText("华为横向 · FEC 前向纠错", {
    x: M, y: 0.95, w: 8, h: 0.4, fontFace: FONT, fontSize: 13.5, bold: true,
    color: GREEN, align: "left", margin: 0, charSpacing: 2,
  });
  s.addText("RS + BCH 级联码\n低时延译码方案", {
    x: M, y: 1.45, w: 8.5, h: 1.9, fontFace: FONT, fontSize: 42, bold: true,
    color: GREEN_DARK, align: "left", margin: 0, lineSpacingMultiple: 1.02,
  });
  s.addText([
    { text: "背景：", options: { bold: true, color: GREEN } },
    { text: "PAM4 高速链路误码需要强纠错，但传统 OSD 类软判译码 ", options: {} },
    { text: "延迟高、难并行", options: { bold: true, color: BROWN } },
    { text: "。本方案用 Lagrange 插值取代高斯消元把延迟降下来。", options: {} },
  ], {
    x: M, y: 3.55, w: 11.8, h: 0.7, fontFace: FONT, fontSize: 15.5, color: INK,
    align: "left", margin: 0, lineSpacingMultiple: 1.15, valign: "top",
  });

  // three benefit chips (native)
  const chips = [
    ["纠错更强", "外码 RS 兜底内码残余错误，级联叠加纠错力"],
    ["延迟更低", "Lagrange 插值取代高斯消元，关键路径大幅缩短"],
    ["抗突发", "矩阵交织把突发错误打散到多个 BCH 码字"],
  ];
  const cy = 4.55, cw = (W - 2 * M - 2 * 0.35) / 3;
  chips.forEach((c, i) => {
    const x = M + i * (cw + 0.35);
    card(s, x, cy, cw, 1.55, GREEN_WASH, GREEN_75);
    s.addText(c[0], { x: x + 0.22, y: cy + 0.18, w: cw - 0.44, h: 0.45,
      fontFace: FONT, fontSize: 17, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
    s.addText(c[1], { x: x + 0.22, y: cy + 0.68, w: cw - 0.44, h: 0.78,
      fontFace: FONT, fontSize: 12.5, color: INK, align: "left", margin: 0,
      valign: "top", lineSpacingMultiple: 1.1 });
  });

  s.addText("陈诗阳 · 2026-07-27", {
    x: M, y: 6.75, w: 6, h: 0.35, fontFace: FONT, fontSize: 11.5, color: GRAYL,
    align: "left", margin: 0,
  });
  s.addNotes("开场：为什么做 RS+BCH 级联码。高速 PAM4 链路误码需要强纠错，传统软判译码延迟高、难并行。我们的方案在保证纠错增益的同时，用 Lagrange 插值把译码延迟降下来，并用交织抗突发。三点价值：纠错更强、延迟更低、抗突发。");
})();

// ============================================================ SLIDE 2 — 原理与流程图（全原生）
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "原理 · 完整链路", GREEN);
  title(s, "级联码构造与收发链路");

  // ---- Row A: cascade encode chain ----
  s.addText("① 级联码构造：外码 RS（符号级）⊕ 内码 BCH（比特级）", {
    x: M, y: 1.3, w: W - 2 * M, h: 0.35, fontFace: FONT, fontSize: 14.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });
  const ay = 1.82, ah = 0.95;
  const encChain = [
    ["信息符号\nk 个", GRAY],
    ["外码 RS 编码\n(符号级校验)", GREEN_DARK],
    ["符号 → 比特\n串行化", GRAY],
    ["分块 → 内码\nBCH 编码", GREEN],
    ["编码比特\n→ 信道", GRAY],
  ];
  const n = encChain.length, bw = 2.15, gap = (W - 2 * M - n * bw) / (n - 1);
  encChain.forEach((b, i) => {
    const x = M + i * (bw + gap);
    pbox(s, x, ay, bw, ah, b[0], b[1], { fs: 12 });
    if (i) arrow(s, x - gap + 0.02, ay + ah / 2, x - 0.02, ay + ah / 2, GRAY, 1.75);
  });
  // annotation bands under RS / BCH
  s.addText("抗突发 · 符号级", { x: M + (bw + gap), y: ay + ah + 0.02, w: bw, h: 0.28,
    fontFace: FONT, fontSize: 10.5, color: GREEN_DARK, align: "center", margin: 0 });
  s.addText("纠随机错 · 比特级", { x: M + 3 * (bw + gap), y: ay + ah + 0.02, w: bw, h: 0.28,
    fontFace: FONT, fontSize: 10.5, color: GREEN, align: "center", margin: 0 });

  // ---- Row B: full TX→channel→RX link ----
  s.addText("② 完整链路：编码 → PAM4 调制 → AWGN 信道 → 两级译码", {
    x: M, y: 3.35, w: W - 2 * M, h: 0.35, fontFace: FONT, fontSize: 14.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });
  const ty = 3.9, th = 0.95;
  const tx = [
    ["信息比特", GRAY], ["RS 外码编码", GREEN_DARK], ["BCH 内码编码", GREEN],
    ["PAM4 调制", GREEN_50], ["AWGN 信道", BROWN],
  ];
  const tn = tx.length, tbw = 2.15, tgap = (W - 2 * M - tn * tbw) / (tn - 1);
  tx.forEach((b, i) => {
    const x = M + i * (tbw + tgap);
    pbox(s, x, ty, tbw, th, b[0], b[1], { fs: 12.5, fill: i === 4 ? BROWN : WHITE, tc: i === 4 ? WHITE : b[1] });
    if (i) arrow(s, x - tgap + 0.02, ty + th / 2, x - 0.02, ty + th / 2, GRAY, 1.75);
  });
  // down arrow from AWGN to RX row
  const chX = M + 4 * (tbw + tgap) + tbw / 2;
  arrow(s, chX, ty + th + 0.02, chX, 5.15, BROWN, 2.2);

  const by = 5.15, bh2 = 0.95;
  const rx = [
    ["PAM4 解调\n(软信息 LLR)", GREEN_50], ["BCH 内码译码\nLLOSD", GREEN],
    ["RS 外码译码\nLCC-BR / BM", GREEN_DARK], ["信息比特", GRAY],
  ];
  // right-to-left layout: box0 sits under the channel column
  const rn = rx.length, rbw = 2.15, rgap = (W - 2 * M - rn * rbw) / (rn - 1);
  // position so the first RX box aligns under channel; simplest: full-width even spread, reversed arrows
  rx.forEach((b, i) => {
    const x = M + (rn - 1 - i) * (rbw + rgap);  // rightmost = index 0
    pbox(s, x, by, rbw, bh2, b[0], b[1], { fs: 12 });
    if (i) {
      const xPrev = M + (rn - i) * (rbw + rgap);
      arrow(s, xPrev - 0.02, by + bh2 / 2, x + rbw + 0.02, by + bh2 / 2, GRAY, 1.75);
    }
  });

  // bottom note
  card(s, M, 6.45, W - 2 * M, 0.72, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "两级译码是延迟与开销的核心：", options: { bold: true, color: GREEN_DARK } },
    { text: "内码用 LLOSD、外码用 LCC-BR，二者同在一个 GF(2ᵐ) 上可共享 Lagrange 代数结构。", options: { color: INK } },
  ], { x: M + 0.25, y: 6.5, w: W - 2 * M - 0.5, h: 0.62, fontFace: FONT, fontSize: 12.5,
       align: "left", margin: 0, valign: "middle", lineSpacingMultiple: 1.05 });

  s.addNotes("这页把原理和链路一次讲清。上排：信息符号先过外码 RS 编码(符号级、抗突发)，串行成比特后分块做内码 BCH 编码(比特级、纠随机错)，再进信道。下排完整链路：编码→PAM4调制→AWGN信道，再反向 PAM4解调(算LLR)→BCH内码译码(LLOSD)→RS外码译码(LCC-BR/BM)→信息比特。两级译码是延迟和开销的核心，内外码同域可共享代数结构。");
})();

// ============================================================ SLIDE 3 — 创新点（全原生）
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "核心创新", GREEN);
  title(s, "创新点：Lagrange 插值取代高斯消元");

  // ---- LEFT: traditional OSD serial chain (de-emphasized gray) ----
  const lx = M, lw = 3.7;
  s.addText("传统 OSD：高斯消元", { x: lx, y: 1.4, w: lw, h: 0.4, fontFace: FONT,
    fontSize: 15, bold: true, color: GRAY, align: "center", margin: 0 });
  s.addText("对 k×k 矩阵逐行消元", { x: lx, y: 1.8, w: lw, h: 0.32, fontFace: FONT,
    fontSize: 11.5, color: GRAY, align: "center", margin: 0 });
  const steps = ["消第 1 行", "消第 2 行 (依赖上一行)", "消第 3 行 (依赖上一行)", "消第 4 行 (依赖上一行)"];
  const sy = 2.25, sh = 0.62, sgap = 0.28;
  steps.forEach((t, i) => {
    const y = sy + i * (sh + sgap);
    pbox(s, lx + 0.35, y, lw - 0.7, sh, t, GRAY, { fs: 11, bold: false, tc: INK, fill: GRAY_TINT });
    if (i) arrow(s, lx + lw / 2, y - sgap + 0.02, lx + lw / 2, y - 0.02, GRAYL, 1.75);
  });
  s.addText("关键路径 ∝ k（串行、需选主元）→ 延迟高", {
    x: lx, y: sy + 4 * (sh + sgap) + 0.05, w: lw, h: 0.6, fontFace: FONT, fontSize: 11.5,
    bold: true, color: BROWN, align: "center", margin: 0, valign: "top", lineSpacingMultiple: 1.1 });

  // divider arrow "→"
  s.addText("→", { x: lx + lw + 0.02, y: 3.4, w: 0.7, h: 0.8, fontFace: FONT,
    fontSize: 34, bold: true, color: GREEN, align: "center", valign: "middle", margin: 0 });

  // ---- RIGHT: Lagrange parallel structure (green) ----
  const rx = lx + lw + 0.75, rw = W - rx - M;
  s.addText("本方法 LLOSD：Lagrange 插值", { x: rx, y: 1.4, w: rw, h: 0.4, fontFace: FONT,
    fontSize: 15, bold: true, color: GREEN_DARK, align: "center", margin: 0 });
  s.addText("把重编码看成多项式插值，免高斯消元", { x: rx, y: 1.8, w: rw, h: 0.32,
    fontFace: FONT, fontSize: 11.5, color: GRAY, align: "center", margin: 0 });

  // node box (left) -> 3 parallel basis funcs -> eval box (right)
  const nodeX = rx + 0.1, nodeW = 2.0, nodeY = 2.55, nodeH = 1.9;
  pbox(s, nodeX, nodeY, nodeW, nodeH, "k 个最可靠\n位置的值\n(插值节点)", GREEN_DARK, { fs: 12 });
  const bx = rx + rw / 2 - 1.1, bw2 = 2.2, bh2 = 0.52;
  const basis = ["基函数 L₀(x)", "基函数 L₁(x)", "基函数 L₂(x)"];
  const byTop = 2.5, bgap = 0.28;
  basis.forEach((t, i) => {
    const y = byTop + i * (bh2 + bgap);
    pbox(s, bx, y, bw2, bh2, t, GREEN, { fs: 11.5 });
    arrow(s, nodeX + nodeW + 0.02, nodeY + nodeH / 2, bx - 0.02, y + bh2 / 2, GREEN_50, 1.4);
  });
  const evX = rx + rw - 2.05, evW = 2.0, evY = 2.55, evH = 1.9;
  pbox(s, evX, evY, evW, evH, "并行求值\n得所有位置\n重编码符号", GREEN_50, { fs: 12, tc: GREEN_DARK });
  basis.forEach((t, i) => {
    const y = byTop + i * (bh2 + bgap);
    arrow(s, bx + bw2 + 0.02, y + bh2 / 2, evX - 0.02, evY + evH / 2, GREEN_50, 1.4);
  });
  s.addText("基函数与求值都可并行、无主元依赖 → 关键路径大幅缩短、延迟低", {
    x: rx, y: 4.65, w: rw, h: 0.5, fontFace: FONT, fontSize: 12, bold: true,
    color: GREEN_DARK, align: "center", margin: 0, valign: "top", lineSpacingMultiple: 1.1 });

  // ---- bottom: two insight cards ----
  const yb = 5.5, cwB = (W - 2 * M - 0.4) / 2;
  card(s, M, yb, cwB, 1.55, BROWN_TINT, BROWN);
  s.addText("关键洞察（IEEE TIT 2025）", { x: M + 0.22, y: yb + 0.14, w: cwB - 0.44, h: 0.38,
    fontFace: FONT, fontSize: 13.5, bold: true, color: BROWN, align: "left", margin: 0 });
  s.addText([
    { text: "BCH ⊂ RS", options: { bold: true, color: BROWN } },
    { text: "，故 BCH 可借 RS 的 Lagrange 插值来译码，把 OSD 里最贵的高斯消元换成并行插值求值。", options: { color: INK } },
  ], { x: M + 0.22, y: yb + 0.54, w: cwB - 0.44, h: 0.95, fontFace: FONT, fontSize: 12.5,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.15 });

  card(s, M + cwB + 0.4, yb, cwB, 1.55, GREEN_WASH, GREEN_75);
  const ox = M + cwB + 0.4;
  s.addText("我们的扩展：内外码共享 Lagrange", { x: ox + 0.22, y: yb + 0.14, w: cwB - 0.44, h: 0.38,
    fontFace: FONT, fontSize: 13.5, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  s.addText([
    { text: "内外码同在一个 GF(2ᵐ) 上，可共享 α 幂表、两两差值表、分母积、基函数结构。", options: { color: INK } },
    { text: "免主元、易并行、关键路径短 → 延迟低。", options: { bold: true, color: GREEN_DARK } },
  ], { x: ox + 0.22, y: yb + 0.54, w: cwB - 0.44, h: 0.95, fontFace: FONT, fontSize: 12.5,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.15 });

  s.addNotes("创新点核心。传统 OSD 要对 k×k 矩阵做高斯消元：逐行、依赖上一行、需选主元，是串行关键路径，延迟高难并行。IEEE TIT 2025 的洞察：BCH 是 RS 的子码，可用 RS 的 Lagrange 插值来译码——把重编码看成多项式插值，基函数和求值都能并行、无主元依赖，关键路径大幅缩短。我们的扩展：内码 LLOSD 与外码 LCC-BR 同域，可共享幂表/差值表/分母积/基函数，进一步省算力。");
})();

// ============================================================ SLIDE 4 — 结果图与分析（保留仿真图片）
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "仿真结果 · PAM4-AWGN", GREEN);
  title(s, "结果与分析：纠错增益 vs 延迟代价");

  // LEFT: BER curve figure (kept as image — real simulation output)
  s.addText("① 纠错性能：级联 vs 纯 RS（BER-SNR，n=255）", {
    x: M, y: 1.28, w: 6.2, h: 0.35, fontFace: FONT, fontSize: 13.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0 });
  // BER 1120x700 → fit into 6.1 wide box
  const berW = 6.1, berH = berW / AR_BER;
  s.addImage({ path: path.join(RSFIG, "n255_scheme_a_ber.png"), x: M, y: 1.7, w: berW, h: berH });

  // RIGHT: KPI latency bars (kept as image)
  const rx = 7.0;
  s.addText("② 延迟 KPI：对软/硬两种基线的延迟增幅", {
    x: rx, y: 1.28, w: W - rx - M, h: 0.35, fontFace: FONT, fontSize: 13.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0 });
  const kpiW = W - rx - M, kpiH = kpiW / AR_KPI;
  s.addImage({ path: path.join(RSFIG, "n255_kpi_bars.png"), x: rx, y: 1.9, w: kpiW, h: kpiH });

  // analysis cards under the right figure
  const ay = 1.9 + kpiH + 0.25;
  card(s, rx, ay, W - rx - M, 1.5, GREEN_WASH, GREEN_75);
  s.addText("延迟结论", { x: rx + 0.2, y: ay + 0.12, w: W - rx - M - 0.4, h: 0.34,
    fontFace: FONT, fontSize: 13, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  s.addText([
    { text: "对软判基线 LCC-BR：低 SNR 延迟 −20%~−27%，KPI（≤+10%）达标；", options: { color: INK, breakLine: true } },
    { text: "对硬判基线 BM：软判法必然更贵（+1156%~+1285%），需 Direct 硬判路线达标。", options: { color: BROWN } },
  ], { x: rx + 0.2, y: ay + 0.5, w: W - rx - M - 0.4, h: 0.95, fontFace: FONT, fontSize: 12,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.15 });

  // BER analysis strip (bottom-left, under BER figure)
  const by = 1.7 + berH + 0.15;
  card(s, M, by, 6.1, 1.5, BROWN_TINT, BROWN);
  s.addText("纠错结论（诚实）", { x: M + 0.2, y: by + 0.12, w: 5.7, h: 0.34,
    fontFace: FONT, fontSize: 13, bold: true, color: BROWN, align: "left", margin: 0 });
  s.addText([
    { text: "门限效应：", options: { bold: true, color: BROWN } },
    { text: "低 SNR（6–7.5 dB）级联劣于纯 RS；高 SNR（≥8 dB）反超，9.5 dB 级联 BER=0、纯 RS 仍 2.2×10⁻⁴。BCH t=1 几乎无增益 → 建议 t≥2。", options: { color: INK } },
  ], { x: M + 0.2, y: by + 0.5, w: 5.7, h: 0.95, fontFace: FONT, fontSize: 12,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.15 });

  s.addNotes("结果与分析，两张真实仿真图。左：BER-SNR，级联 vs 纯RS。要诚实讲门限效应——低SNR段级联反而差(内码BCH误纠引额外错)，8dB以上反超，9.5dB级联BER到0而纯RS还有2e-4；BCH t=1几乎无增益，建议t≥2。右：延迟KPI柱状图。跟软判基线LCC-BR比延迟降20-27%达标；跟硬判基线BM比因为软判对硬判本就贵会+1156%~+1285%，这不是bug，要用Direct硬判路线满足硬判KPI。核心是：拿纠错增益换延迟，讲清跟谁比。");
})();

// ============================================================ SLIDE 5 — 总结 & 下一步
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  topBand(s, GREEN);
  s.addText("总结 & 下一步", {
    x: M, y: 0.55, w: W - 2 * M, h: 0.9, fontFace: FONT, fontSize: 32, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });

  const items = [
    ["做了什么", GREEN,
      "在 PAM4-AWGN 全链路上实现 RS+BCH 级联码：外码 RS 抗突发、内码 BCH 纠随机错，编码→调制→信道→两级译码闭环可跑。"],
    ["核心创新", GREEN,
      "用 Lagrange 插值替代 OSD 的高斯消元（免主元、可并行、关键路径短），内外码共享有限域代数结构进一步省算力。"],
    ["抗突发交织", GREEN,
      "复现 APCC-2022 三步矩阵交织，w=1bit 把突发打散到多码字、每字≤1 错，零码率代价的抗突发增益。"],
    ["诚实结论", BROWN,
      "开销约 12%（外码 RS × 内码 BCH 码率）；延迟对软判基线达标、对硬判基线偏贵；BER 低 SNR 有门限、高 SNR 反超；内码建议 t≥2。"],
  ];
  const y0 = 1.65, rh = 1.2;
  items.forEach((it, i) => {
    const y = y0 + i * rh;
    badge(s, M, y + 0.08, i + 1, 0.32, it[1]);
    s.addText(it[0], { x: M + 0.85, y: y, w: 2.5, h: 0.5, fontFace: FONT,
      fontSize: 16, bold: true, color: GREEN_DARK, align: "left", margin: 0, valign: "middle" });
    s.addText(it[2], {
      x: M + 3.45, y: y - 0.04, w: W - M - 3.45 - M, h: rh - 0.12,
      fontFace: FONT, fontSize: 13, color: INK, align: "left", margin: 0,
      valign: "middle", lineSpacingMultiple: 1.12 });
  });

  card(s, M, 6.45, W - 2 * M, 0.72, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "下一步：", options: { bold: true, color: GREEN_DARK } },
    { text: "① 内码升到 BCH t≥2 提升级联增益；② P=64 并行 + Lagrange 共享压延迟；③ 补软判 Chase 交织增益。", options: { color: INK } },
  ], { x: M + 0.25, y: 6.5, w: W - 2 * M - 0.5, h: 0.62, fontFace: FONT,
       fontSize: 12, align: "left", margin: 0, valign: "middle", lineSpacingMultiple: 1.05 });

  s.addNotes("总结四点：做了什么(全链路级联码闭环)、核心创新(Lagrange替代高斯消元+共享)、抗突发交织(APCC2022三步交织)、诚实结论(开销约12%、延迟看基线、BER有门限、建议t≥2)。下一步：升内码t、并行压延迟、补软判交织。结束。");
})();

pres.writeFile({ fileName: OUT }).then((f) => {
  console.log("WROTE", f);
}).catch((e) => { console.error("ERR", e); process.exit(1); });
