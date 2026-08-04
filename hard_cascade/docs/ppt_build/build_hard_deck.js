// Build the 4-page HARD-DECISION-ONLY RS+BCH cascade deck with pptxgenjs.
//
// Scope (per advisor): hard-decision RS+BCH cascade ONLY — no LLOSD/OSD soft-decision
// content here (that lives in rs_bch_cascade_deck.pptx). 4 pages:
//   1. 创新点 (innovation points)
//   2. 级联码构造原理 + 译码原理（构造矩阵、交织、查表）
//   3. 折线图/柱状图仿真结果
//   4. 时延分析、复杂度分析、乘法器数量、LUT 数量
//
// Reuses the 3 existing native drawio diagrams under docs/diagrams/ (as images —
// they are already fully-formed schematics) plus the 5 real result figures under
// figures/. Writes to a NEW file — does NOT touch rs_bch_级联码-Hard.pptx (open in
// PowerPoint) and does NOT touch rs_bch_cascade_deck.pptx (soft-decision deck).
//
// Palette / helpers follow build_deck.js conventions (BIT official VI colors).

const pptxgen = require("pptxgenjs");
const path = require("path");

const ROOT = "/Users/chenshiyang.10/workspace/llosd_reproduction";
const DOCS = path.join(ROOT, "hard_cascade/docs");
const DIAG = path.join(DOCS, "diagrams");
const FIG = path.join(ROOT, "hard_cascade/figures");
const CROPS = path.join(DOCS, "ppt_build/crops");
const OUT = path.join(DOCS, "rs_bch_cascade_hard_deck.pptx");

// ---- palette (BIT official VI colors) ----
const GREEN      = "009A44";
const GREEN_75   = "14AE68";
const GREEN_50   = "89C997";
const GREEN_DARK = "015C31";
const BROWN      = "A23E0A";
const INK  = "1A1A1A";
const GRAY = "5C6670";
const GRAYL = "8A929B";
const WHITE = "FFFFFF";
const GREEN_WASH = "EEF8F0";
const GREEN_TINT = "DCF0E1";
const BROWN_TINT = "FBEEE6";
const GRAY_TINT = "F2F3F4";

const FONT = "Microsoft YaHei";

const pres = new pptxgen();
pres.defineLayout({ name: "W", width: 13.333, height: 7.5 });
pres.layout = "W";
pres.author = "陈诗阳";
pres.title = "RS+BCH 级联码硬判决方案";

const W = 13.333, H = 7.5, M = 0.6;

// ---------------------------------------------------------------- helpers
function title(s, txt, opts) {
  opts = opts || {};
  s.addText(txt, {
    x: M, y: 0.34, w: W - 2 * M, h: 0.9, fontFace: FONT, fontSize: opts.size || 28,
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
function arrow(s, x1, y1, x2, y2, color, w) {
  s.addShape(pres.ShapeType.line, {
    x: Math.min(x1, x2), y: Math.min(y1, y2),
    w: Math.abs(x2 - x1), h: Math.abs(y2 - y1),
    line: { color: color || INK, width: w || 1.75, endArrowType: "triangle" },
    flipH: x2 < x1, flipV: y2 < y1,
  });
}

// ============================================================ SLIDE 1 — 创新点
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  topBand(s, GREEN);
  kicker(s, "硬判决 RS+BCH 级联码 · 创新点", GREEN);
  title(s, "创新点：三级硬件复用把级联开销压到 KPI 之内", { size: 26 });

  s.addText([
    { text: "背景：", options: { bold: true, color: GREEN } },
    { text: "外码 RS 抗突发、内码 BCH 纠随机错的级联码能带来编码增益，但 内码译码是级联链路上 ", options: {} },
    { text: "多出来的时延", options: { bold: true, color: BROWN } },
    { text: "——若不优化，会超过导师给定的 10% 时延 KPI 红线。", options: {} },
  ], {
    x: M, y: 1.18, w: W - 2 * M, h: 0.5, fontFace: FONT, fontSize: 13.5, color: INK,
    align: "left", margin: 0, lineSpacingMultiple: 1.15, valign: "top",
  });

  // three-step innovation chain (native)
  const steps = [
    ["Cascade", GRAY, "引入内码 BCH 兜底信道随机错\n代价：级联本身多出一段内码译码时延\n(Config1 +4cyc / Config2 +8cyc)"],
    ["Direct", GREEN, "内码 BCH 译码从 Conv(BM+Chien)\n换成 Direct 查表求根\n(Config1 省2cyc / Config2 省5cyc，\n最大单项节省)"],
    ["Lagrange 共享", GREEN_DARK, "内外码共享 syndrome 计算单元(v1，省1cyc)\n再共享 GF乘法器阵列(v2，再省1cyc)\n把 BCH↔RS 的握手时延也省掉"],
  ];
  const cy = 1.85, cw = (W - 2 * M - 2 * 0.4) / 3, ch = 2.15;
  steps.forEach((c, i) => {
    const x = M + i * (cw + 0.4);
    card(s, x, cy, cw, ch, i === 2 ? GREEN_WASH : WHITE, c[1]);
    badge(s, x + 0.2, cy + 0.18, i + 1, 0.26, c[1]);
    s.addText(c[0], { x: x + 0.82, y: cy + 0.14, w: cw - 1.02, h: 0.42,
      fontFace: FONT, fontSize: 16.5, bold: true, color: c[1], align: "left", valign: "middle", margin: 0 });
    s.addText(c[2], { x: x + 0.22, y: cy + 0.68, w: cw - 0.44, h: ch - 0.85,
      fontFace: FONT, fontSize: 12, color: INK, align: "left", margin: 0,
      valign: "top", lineSpacingMultiple: 1.2 });
    if (i < 2) arrow(s, x + cw + 0.04, cy + ch / 2, x + cw + 0.36, cy + ch / 2, GRAYL, 2);
  });

  // KPI outcome table (native)
  const ty = 4.28;
  s.addText("三个创新点逐项叠加后的 KPI 结果", { x: M, y: ty, w: 6, h: 0.32,
    fontFace: FONT, fontSize: 13.5, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  s.addTable([
    [
      { text: "Config", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "Pure RS-BM", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "+Cascade(Conv)", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "+Direct", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "+Lagrange v1", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "+Lagrange v2", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
    ],
    [
      { text: "Config 1 (KP4)", options: { bold: true } },
      "33 cyc（基线）",
      { text: "37 cyc  +12.1%  ✗", options: { color: BROWN } },
      { text: "35 cyc  +6.1%  ✓", options: { color: GREEN_DARK, bold: true } },
      { text: "34 cyc  +3.0%  ✓", options: { color: GREEN_DARK } },
      { text: "33 cyc  +0.0%  ✓", options: { color: GREEN_DARK } },
    ],
    [
      { text: "Config 2 (n=255)", options: { bold: true } },
      "19 cyc（基线）",
      { text: "27 cyc  +42.1%  ✗", options: { color: BROWN } },
      { text: "22 cyc  +15.8%  ✗", options: { color: BROWN } },
      { text: "21 cyc  +10.5%  ✗", options: { color: BROWN } },
      { text: "20 cyc  +5.3%  ✓", options: { color: GREEN_DARK, bold: true } },
    ],
  ], {
    x: M, y: ty + 0.36, w: W - 2 * M, h: 1.15, fontFace: FONT, fontSize: 11.5,
    border: { type: "solid", color: GREEN_75, pt: 0.75 }, align: "center", valign: "middle",
    autoPage: false,
  });

  card(s, M, 5.95, W - 2 * M, 1.1, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "两个 config 的达标路径不同：", options: { bold: true, color: GREEN_DARK } },
    { text: "Config 1 外码大(t=15)、内码极小(t=1)，仅 ", options: { color: INK } },
    { text: "Direct 一项", options: { bold: true, color: GREEN_DARK } },
    { text: " 即可达标；Config 2 时延基数小(19cyc)、内码开销占比更重，必须", options: { color: INK } },
    { text: " 三项创新齐全", options: { bold: true, color: GREEN_DARK } },
    { text: " 才达标（KPI ≤10%）。外码越大、内码越小，级联时延开销越可忽略。", options: { color: INK } },
  ], { x: M + 0.25, y: 6.02, w: W - 2 * M - 0.5, h: 0.95, fontFace: FONT, fontSize: 12.5,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.2 });

  s.addNotes("创新点：三级硬件复用。①Cascade引入内码BCH是时延净开销来源；②Direct把BCH译码从Conv(BM+Chien)换成查表求根，是最大单项节省；③Lagrange硬件共享，内外码共享syndrome单元(v1)、再共享GF乘法器阵列(v2)。两个config达标路径不同：Config1外码大内码小，仅Direct即达标；Config2内码开销占比重，必须三项齐全才达标。");
})();

// ============================================================ SLIDE 2 — 构造原理 + 译码原理
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "原理 · 构造 + 译码 + 交织 + 查表", GREEN);
  title(s, "级联码构造原理与译码原理", { size: 26 });

  // ---- Row A: construction chain + config params (reuse fig1 diagram, cropped) ----
  s.addText("① 码构造：外码 RS（符号级）⊕ 内码 BCH（比特级），两组实用码型", {
    x: M, y: 1.1, w: W - 2 * M, h: 0.3, fontFace: FONT, fontSize: 13, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });
  const fig1AR = 1602 / 594, fig1H = 2.55, fig1W = fig1H * fig1AR;
  s.addImage({ path: path.join(CROPS, "fig1_top.png"), x: (W - fig1W) / 2, y: 1.42, w: fig1W, h: fig1H });

  // ---- Row B: BCH decode two paths (reuse fig3 diagram, cropped) ----
  s.addText("② 内码 BCH 译码：Conventional（传统）vs Direct（查表求根，本方案）", {
    x: M, y: 4.08, w: W - 2 * M, h: 0.3, fontFace: FONT, fontSize: 13, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });
  const fig3AR = 1802 / 310, fig3H = 1.15, fig3W = fig3H * fig3AR;
  s.addImage({ path: path.join(CROPS, "fig3_bch_tight.png"), x: (W - fig3W) / 2, y: 4.4, w: fig3W, h: fig3H });

  // ---- Row C: interleaving 3-step (native, compact — not covered by existing diagrams) ----
  s.addText("③ 抗突发：三步矩阵交织（复现 APCC 2022，接入 Config 1）", {
    x: M, y: 5.95, w: W - 2 * M, h: 0.28, fontFace: FONT, fontSize: 13, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });
  const iy = 6.26, ih = 0.58;
  const ilv = [
    ["Step1 预处理：RS码字符号轮转交错", GREEN_DARK],
    ["Step2 矩阵交织：按列写入·按行做BCH编码", GREEN],
    ["Step3 信道传输：矩阵内按行发比特", GREEN_50],
  ];
  const iw = (W - 2 * M - 2 * 0.35) / 3;
  ilv.forEach((b, i) => {
    const x = M + i * (iw + 0.35);
    pbox(s, x, iy, iw, ih, b[0], b[1], { fs: 10.5 });
    if (i) arrow(s, x - 0.35 + 0.02, iy + ih / 2, x - 0.02, iy + ih / 2, GRAY, 1.5);
  });

  s.addNotes("原理页三部分：①码构造总链路+两组实用码型参数(Config1是KP4以太网FEC缩短码，Config2是均衡型本原码，padding不计入开销)，直接复用图1构造流程图；②内码BCH两条译码路径，Conventional传统BM+Chien 4/8周期，Direct本方案查表求根2/3周期，复用图3译码流程图；③三步矩阵交织抗突发，复现APCC2022论文——预处理RS码字轮转交错、按列写入按行BCH编码、信道按行发送，列宽w有1bit/2bits/1symbol三档权衡。");
})();

// ============================================================ SLIDE 3 — Direct 译码原理
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  topBand(s, GREEN);
  kicker(s, "原理 · BCH Direct 译码 · 查表替代迭代求根", GREEN);
  title(s, "Direct 译码原理：闭式代数变换 + 查表求根", { size: 26 });

  s.addText([
    { text: "思路：", options: { bold: true, color: GREEN } },
    { text: "传统 BM 译码需迭代 2t 次求错误定位多项式 Λ(x)，再做 Chien 搜索遍历所有位置。Lagendijk 论文指出：t≤4 时，Λ(x) 的根可以用", options: {} },
    { text: " 闭式代数变换 + LUT 查表", options: { bold: true, color: GREEN_DARK } },
    { text: " 直接求出，省掉迭代和遍历。", options: {} },
  ], {
    x: M, y: 1.14, w: W - 2 * M, h: 0.42, fontFace: FONT, fontSize: 13, color: INK,
    align: "left", margin: 0, lineSpacingMultiple: 1.15, valign: "top",
  });

  // t=1 degenerate case strip
  card(s, M, 1.62, W - 2 * M, 0.6, WHITE, GRAY);
  s.addText([
    { text: "t=1（退化闭式）：  ", options: { bold: true, color: GRAY } },
    { text: "S₁=0 → 无错；否则错误位置 p = log_α(S₁)，直接查 α 幂表后修正 1 个比特。", options: { color: INK } },
    { text: "   共 2 cycles：syndrome + (log 查表 + 修正)。", options: { color: GRAY } },
  ], { x: M + 0.22, y: 1.62, w: W - 2 * M - 0.4, h: 0.6, fontFace: FONT, fontSize: 12.5,
       align: "left", valign: "middle", margin: 0, lineSpacingMultiple: 1.1 });

  // t=2 three-cycle chain (native, mirrors slide1 chain style)
  s.addText("t=2（通用情况）：代数变换分 3 个 cycle 完成求根", {
    x: M, y: 2.42, w: W - 2 * M, h: 0.3, fontFace: FONT, fontSize: 13.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0,
  });
  const steps = [
    ["Cycle 1 · Syndrome", GRAY, "计算 {S₁, S₃}\n用 n 个 GF 乘法器并行求值\n（与外码 RS 的 syndrome\n可共享同一硬件，见下页）"],
    ["Cycle 2 · 变量代换", GREEN, "Λ_monic(x) = x²+S₁x+(S₁³+S₃)/S₁\n令 x = S₁·Y\nA(Y) = Y²+Y+k,  k=(S₁³+S₃)/S₁³\n（3 次 GF 乘法 + 1 次除法）"],
    ["Cycle 3 · LUT 求根+修正", GREEN_DARK, "查 2ᵏ 大小 LUT[k] → (Y₁,Y₂)\nXᵢ = S₁·Yᵢ,  pᵢ = log_α(Xᵢ)\n直接翻转错误位（无需 Chien\n遍历全部 n 个位置）"],
  ];
  const cy = 2.76, cw = (W - 2 * M - 2 * 0.4) / 3, ch = 1.95;
  steps.forEach((c, i) => {
    const x = M + i * (cw + 0.4);
    card(s, x, cy, cw, ch, i === 2 ? GREEN_WASH : WHITE, c[1]);
    badge(s, x + 0.2, cy + 0.18, i + 1, 0.26, c[1]);
    s.addText(c[0], { x: x + 0.82, y: cy + 0.12, w: cw - 1.02, h: 0.46,
      fontFace: FONT, fontSize: 13, bold: true, color: c[1], align: "left", valign: "middle", margin: 0 });
    s.addText(c[2], { x: x + 0.22, y: cy + 0.62, w: cw - 0.44, h: ch - 0.78,
      fontFace: FONT, fontSize: 11.5, color: INK, align: "left", margin: 0,
      valign: "top", lineSpacingMultiple: 1.18 });
    if (i < 2) arrow(s, x + cw + 0.04, cy + ch / 2, x + cw + 0.36, cy + ch / 2, GRAYL, 2);
  });

  // key insight highlight box
  card(s, M, 4.94, W - 2 * M, 1.15, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "关键性质（char-2 域）：", options: { bold: true, color: GREEN_DARK } },
    { text: "在 GF(2ᵏ) 上，若 Y₁ 是 A(Y) 的一个根，则 Y₂ = Y₁+1 也一定是根——这是二次多项式在特征 2 域上的固有性质。因此", options: { color: INK } },
    { text: " 只需一个 2ᵏ 大小的 LUT", options: { bold: true, color: GREEN_DARK } },
    { text: "（而非两个独立 LUT）即可同时得到两个根，LUT 可完全展开为组合逻辑，1 个 clock cycle 内完成。", options: { color: INK } },
  ], { x: M + 0.25, y: 5.02, w: W - 2 * M - 0.5, h: 1.0, fontFace: FONT, fontSize: 12.5,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.2 });

  // bottom comparison note
  s.addText([
    { text: "对比：", options: { bold: true, color: BROWN } },
    { text: " Conventional 需 BM 迭代 2t 次求 Λ(x) + Chien 搜索遍历 n 个位置 → t=2 需 8 cycle；", options: { color: INK } },
    { text: "Direct 用代数变换替代迭代、单次 LUT 替代遍历 → t=2 仅 3 cycle，t=1 仅 2 cycle。", options: { color: INK } },
  ], { x: M, y: 6.24, w: W - 2 * M, h: 0.5, fontFace: FONT, fontSize: 12, align: "left",
       margin: 0, valign: "top", lineSpacingMultiple: 1.2 });

  s.addNotes("Direct译码原理页。t=1退化情况：S1=0无错，否则p=log_α(S1)直接查表，2cycle。t=2通用情况分3个cycle：Cycle1算syndrome{S1,S3}；Cycle2变量代换x=S1·Y把Λ_monic(x)变成标准形A(Y)=Y²+Y+k；Cycle3查一个2^m大小的LUT直接得到两个根——关键是char-2域上二次多项式的性质：若Y1是根则Y1+1也是根，所以只需1个LUT不需要2个，LUT可展开为纯组合逻辑1个cycle内完成，不需要Chien那种遍历全部n个位置的搜索。对比Conventional的BM迭代(2t次)+Chien遍历，Direct的关键优势是用代数变换和查表替代了迭代与遍历这两个时延大头。");
})();

// ============================================================ SLIDE 4 — 拉格朗日硬件共享原理
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "原理 · 硬件级联复用 · Lagrange 共享 v1/v2", GREEN);
  title(s, "拉格朗日共享原理：硬件模块复用逐级压缩级联时延", { size: 25 });

  s.addText([
    { text: "说明：", options: { bold: true, color: GREEN } },
    { text: "这里的“拉格朗日共享”指", options: {} },
    { text: " 硬件模块复用", options: { bold: true, color: GREEN_DARK } },
    { text: "（内外码共享同一套 syndrome 计算单元 / GF 乘法器阵列），不是拉格朗日插值运算本身——BCH-Direct 用查表求根，RS 仍用 BM 迭代，共享的是电路，不是算法。", options: {} },
  ], {
    x: M, y: 1.12, w: W - 2 * M, h: 0.46, fontFace: FONT, fontSize: 12.5, color: INK,
    align: "left", margin: 0, lineSpacingMultiple: 1.15, valign: "top",
  });

  // v0/v1/v2 three-step chain
  const steps2 = [
    ["v0 · 无共享", GRAY, "T = T_BCH + T_RS-BM\n内外码各自独立硬件模块\n（各自的 α 幂表、syndrome\n计算单元、乘法器阵列）"],
    ["v1 · 共享 syndrome 单元", GREEN, "T = T_BCH + T_RS-BM − 1\n内外码共享 α 幂表 +\nsyndrome 计算单元，同一个\npolynomial evaluator 并行\n算出 BCH 和 RS 的 syndrome"],
    ["v2 · 再共享 GF 乘法器阵列", GREEN_DARK, "T = T_BCH + T_RS-BM − 2\nBCH-Direct 的 LUT 求根\n cycle 与 RS-BM 的 Chien\n搜索分时复用同一套物理\n乘法器阵列，握手 cycle 也省"],
  ];
  const cy2 = 1.68, cw2 = (W - 2 * M - 2 * 0.4) / 3, ch2 = 2.05;
  steps2.forEach((c, i) => {
    const x = M + i * (cw2 + 0.4);
    card(s, x, cy2, cw2, ch2, i === 2 ? GREEN_WASH : WHITE, c[1]);
    badge(s, x + 0.2, cy2 + 0.18, i, 0.26, c[1]);
    s.addText(c[0], { x: x + 0.82, y: cy2 + 0.12, w: cw2 - 1.02, h: 0.46,
      fontFace: FONT, fontSize: 12.5, bold: true, color: c[1], align: "left", valign: "middle", margin: 0 });
    s.addText(c[2], { x: x + 0.22, y: cy2 + 0.62, w: cw2 - 0.44, h: ch2 - 0.78,
      fontFace: FONT, fontSize: 11, color: INK, align: "left", margin: 0,
      valign: "top", lineSpacingMultiple: 1.16 });
    if (i < 2) arrow(s, x + cw2 + 0.04, cy2 + ch2 / 2, x + cw2 + 0.36, cy2 + ch2 / 2, GRAYL, 2);
  });

  // v2 core-observation keybox (quoting report_v2 §2 boxed insight)
  card(s, M, 3.98, W - 2 * M, 1.25, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "v2 共享的核心观察：", options: { bold: true, color: GREEN_DARK } },
    { text: "BCH-Direct 的 Cycle 3（LUT 查表+错误修正）与 RS-BM 的 Chien 搜索都需要 n 路并行 GF 乘法器阵列，但", options: { color: INK } },
    { text: " 二者时间上不冲突", options: { bold: true, color: GREEN_DARK } },
    { text: "——BCH Cycle 3 完成后，RS-BM 才开始 BM 迭代，Chien 阶段又在 BM 之后。因此可以让两个阶段分时复用同一套物理乘法器阵列，无需两份硬件；阵列共享后，原本 BCH 输出→RS 输入之间的数据握手 cycle 也可以省掉，再降 1 cycle。", options: { color: INK } },
  ], { x: M + 0.25, y: 4.06, w: W - 2 * M - 0.5, h: 1.1, fontFace: FONT, fontSize: 11.8,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.18 });

  // LUT-size tradeoff note (report_v2 §3)
  s.addText("延伸：LUT 规模也影响可压缩空间", { x: M, y: 5.42, w: W - 2 * M, h: 0.28,
    fontFace: FONT, fontSize: 12.5, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  s.addTable([
    [
      { text: "GF 域", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "LUT 大小", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "查表 MUX 深度", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "Direct 总周期", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
    ],
    ["n=255（GF(2⁸)）", "256 项", "8-input", "3 cycle"],
    [{ text: "n=127（GF(2⁷)）", options: { bold: true, color: GREEN_DARK } }, "128 项（减半）",
      "7-input（更浅）", { text: "2 cycle（precompute+LUT+修正合并1cycle）", options: { bold: true, color: GREEN_DARK } }],
  ], {
    x: M, y: 5.72, w: W - 2 * M, h: 0.75, fontFace: FONT, fontSize: 10.5,
    border: { type: "solid", color: GREEN_75, pt: 0.75 }, align: "center", valign: "middle",
    autoPage: false,
  });

  s.addNotes("拉格朗日硬件共享原理页。先澄清：这里共享的是硬件模块(syndrome计算单元/GF乘法器阵列)，不是拉格朗日插值运算——BCH-Direct用查表，RS仍用BM迭代。v0无共享T=T_BCH+T_RS-BM；v1共享syndrome计算单元(内外码共享α幂表+同一个polynomial evaluator并行算两个syndrome)省1cycle；v2再共享GF乘法器阵列，核心观察是BCH-Direct的Cycle3和RS-BM的Chien搜索时间上不冲突(BCH先完成，RS的BM+Chien在其后)，所以可以分时复用同一套物理乘法器阵列，且阵列共享后握手cycle也能省，再降1cycle。延伸表格：LUT规模和GF域大小的关系，n=127比n=255的LUT小一半、MUX更浅，能把precompute+LUT+修正合并到1cycle，Direct从3cycle降到2cycle。");
})();

// ============================================================ SLIDE 5 — 时延计算依据（表）
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  topBand(s, GREEN);
  kicker(s, "时延 · 计算依据 · 公式逐项代入", GREEN);
  title(s, "时延计算依据：公式代入两组 config 得到逐项 cycle 数", { size: 25 });

  s.addText([
    { text: "依据：", options: { bold: true, color: GREEN } },
    { text: "下表把", options: {} },
    { text: " 前两页原理公式", options: { bold: true, color: GREEN_DARK } },
    { text: "（RS-BM(t)=2t+3、BCH-Direct(t)、v0/v1/v2 共享）逐项代入 Config1(KP4) 与 Config2(n=255) 的实际参数，得到与第 1 页 KPI 表、第 6 页柱状图完全对应的 cycle 数——不是仿真测出来的，是公式算出来的。", options: {} },
  ], {
    x: M, y: 1.14, w: W - 2 * M, h: 0.5, fontFace: FONT, fontSize: 12.5, color: INK,
    align: "left", margin: 0, lineSpacingMultiple: 1.18, valign: "top",
  });

  s.addText("Config 1（KP4）：RS(544,514,t=15) + BCH(144,136,t=1)　　|　　Config 2（n=255）：RS(255,239,t=8) + BCH(255,239,t=2)", {
    x: M, y: 1.72, w: W - 2 * M, h: 0.3, fontFace: FONT, fontSize: 11, italic: true,
    color: GRAY, align: "left", margin: 0,
  });

  s.addTable([
    [
      { text: "方案", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "公式", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "Config1 代入 (t_rs=15, t_bch=1)", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "cyc", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "Config2 代入 (t_rs=8, t_bch=2)", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
      { text: "cyc", options: { bold: true, fill: { color: GREEN_DARK }, color: WHITE } },
    ],
    ["Pure RS-BM（基线）", "T = 2t_rs+3",
      "2×15+3", { text: "33", options: { bold: true } },
      "2×8+3", { text: "19", options: { bold: true } }],
    ["+Cascade（BCH-Conv, v0）", "T = T_BCH-Conv + T_RS-BM",
      "4 + 33", { text: "37", options: { color: BROWN } },
      "8 + 19", { text: "27", options: { color: BROWN } }],
    ["+Direct（BCH-Direct, v0）", "T = T_BCH-Direct + T_RS-BM",
      "2 + 33", { text: "35", options: { color: BROWN } },
      "3 + 19", { text: "22", options: { color: BROWN } }],
    ["+Lagrange v1（共享 syndrome）", "T = T_BCH-Direct + T_RS-BM − 1",
      "35 − 1", { text: "34", options: { color: BROWN } },
      "22 − 1", { text: "21", options: { color: BROWN } }],
    ["+Lagrange v2（再共享乘法器阵列）", "T = T_BCH-Direct + T_RS-BM − 2",
      "35 − 2", { text: "33", options: { bold: true, color: GREEN_DARK } },
      "22 − 2", { text: "20", options: { bold: true, color: GREEN_DARK } }],
  ], {
    x: M, y: 2.1, w: W - 2 * M, h: 2.85, fontFace: FONT, fontSize: 11.5, valign: "middle",
    border: { type: "solid", color: GREEN_75, pt: 0.75 }, align: "center", autoPage: false,
    colW: [2.35, 2.85, 2.0, 0.65, 2.0, 0.65],
  });

  card(s, M, 5.15, W - 2 * M, 1.55, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "两组 config 读法不同：", options: { bold: true, color: GREEN_DARK } },
    { text: " Config1 基线 33cyc 大、内码 t_bch=1 开销小，Direct 一项代入后 35→34→33，第三步已回到基线（+0.0%），KPI 早早达标；", options: { color: INK } },
    { text: "Config2 基线 19cyc 小、内码 t_bch=2 开销占比更重，同样代入 v0/v1/v2 公式，35cyc 起点的 Conv/Direct 差距在小基线下被放大，必须代入到 v2（−2）才压到 20cyc（+5.3%，达标）。", options: { color: INK } },
  ], { x: M + 0.25, y: 5.24, w: W - 2 * M - 0.5, h: 1.4, fontFace: FONT, fontSize: 12,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.22 });

  s.addText("* T_RS-BM、T_BCH-Conv/Direct 公式见「构造原理」「Direct 译码原理」「拉格朗日共享原理」三页；本页只做数值代入，不重复推导。", {
    x: M, y: 6.85, w: W - 2 * M, h: 0.35, fontFace: FONT, fontSize: 10, italic: true, color: GRAY,
    align: "left", margin: 0, valign: "top",
  });

  s.addNotes("时延计算依据表，把前面原理页的公式代入两组config实际参数。Config1: t_rs=15,t_bch=1，基线2×15+3=33cyc；Cascade v0代入BCH-Conv(t=1)=4cyc得37；Direct代入BCH-Direct(t=1)=2cyc得35；Lagrange v1减1得34；v2再减1得33，正好回到基线，+0.0%。Config2: t_rs=8,t_bch=2，基线2×8+3=19cyc；Cascade v0代入BCH-Conv(t=2)=8cyc得27；Direct代入BCH-Direct(t=2)=3cyc得22；v1减1得21；v2再减1得20，+5.3%刚好达标。这张表和第1页KPI表、结果页柱状图的数字是同一套数字，只是这里显式展开了公式代入过程，回应导师说现有PPT时延部分只有柱状图、要加计算依据表格的意见。");
})();

// ============================================================ SLIDE 6 — 结果图
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "仿真结果 · PAM4-AWGN / 突发信道", GREEN);
  title(s, "结果：编码增益、突发抗性", { size: 26 });

  // Top-left: BER-SNR line chart
  s.addText("① BER-SNR：级联 vs 纯 RS（两组码型）", {
    x: M, y: 1.15, w: 6.0, h: 0.3, fontFace: FONT, fontSize: 12.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0 });
  const berAR = 1820 / 700, berW = 6.0, berH = berW / berAR;
  s.addImage({ path: path.join(FIG, "ber_snr_v3.png"), x: M, y: 1.48, w: berW, h: berH });
  s.addText("BER=10⁻⁴：Config1 级联≈8.9dB vs 纯RS≈9.7dB(+0.7dB)；Config2 级联≈8.6dB vs 纯RS≈9.6dB(+1.0dB)", {
    x: M, y: 1.48 + berH + 0.04, w: berW, h: 0.42, fontFace: FONT, fontSize: 9.5, color: INK,
    align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.1 });

  // Top-right: interleave gain bar chart
  const rx = 6.95;
  s.addText("② 矩阵交织抗突发增益（vs existing 符号交织）", {
    x: rx, y: 1.15, w: W - rx - M, h: 0.3, fontFace: FONT, fontSize: 12.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0 });
  const gainAR = 1260 / 700, gainColW = W - rx - M, gainH = 2.3, gainW = gainH * gainAR;
  s.addImage({ path: path.join(FIG, "interleave_gain_v4.png"), x: rx + (gainColW - gainW) / 2, y: 1.48, w: gainW, h: gainH });
  s.addText("AWGN 增益小(+0.09~0.15dB)；突发 p=0.75 增益最大(+0.31~0.34dB)——交织在突发信道下价值最明显", {
    x: rx, y: 1.48 + gainH + 0.04, w: gainColW, h: 0.42, fontFace: FONT, fontSize: 9.5, color: INK,
    align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.1 });

  // Bottom: latency waterfall bar chart (full width, most informative single chart)
  const by = 4.35;
  s.addText("③ 三创新点逐项时延贡献（Cascade → Direct → Lagrange v1/v2）", {
    x: M, y: by, w: W - 2 * M, h: 0.3, fontFace: FONT, fontSize: 12.5, bold: true,
    color: GREEN_DARK, align: "left", margin: 0 });
  const latAR = 1820 / 700, latH = 2.55, latW = latH * latAR;
  s.addImage({ path: path.join(FIG, "latency_innovations_v3.png"), x: (W - latW) / 2, y: by + 0.32, w: latW, h: latH });

  s.addNotes("三张图：①BER-SNR，级联相对纯RS都有0.7~1.0dB编码增益；②矩阵交织在突发信道下的额外增益，p=0.75时最大到+0.34dB；③三创新点(Cascade+Direct+Lagrange)逐项时延贡献柱状图，两个config都能压到10%KPI虚线以下。");
})();

// ============================================================ SLIDE 7 — 时延/复杂度/乘法器/LUT
(function () {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  kicker(s, "时延 · 复杂度 · 硬件资源", GREEN);
  title(s, "时延分析、复杂度分析与硬件资源开销", { size: 26 });

  // Left: latency cycle-cost table
  s.addText("① 时延模型（clock cycles）", { x: M, y: 1.15, w: 6.0, h: 0.3,
    fontFace: FONT, fontSize: 12.5, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  s.addTable([
    [
      { text: "译码器", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
      { text: "时延分解", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
      { text: "cyc", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
    ],
    ["RS-BM(t)", "syndrome+2t·BM+Chien+Forney", "2t+3"],
    ["BCH-Conv(t=1)", "syndrome+ELP+Chien+correct", "4"],
    ["BCH-Conv(t=2)", "1+2t+1+margin", "8"],
    ["BCH-Direct(t=1)", "syndrome+(log查表+修正)", { text: "2", options: { bold: true, color: GREEN_DARK } }],
    ["BCH-Direct(t=2)", "syndrome+precompute+LUT", { text: "3", options: { bold: true, color: GREEN_DARK } }],
  ], {
    x: M, y: 1.48, w: 6.0, h: 1.55, fontFace: FONT, fontSize: 9.5, valign: "middle",
    border: { type: "solid", color: GREEN_75, pt: 0.75 }, autoPage: false,
    colW: [1.3, 3.7, 1.0],
  });
  s.addText([
    { text: "三档共享：", options: { bold: true, color: GREEN_DARK } },
    { text: "T_none→T_v1(−1,共享syndrome)→T_v2(−2,再共享GF乘法器阵列)", options: { color: INK } },
  ], { x: M, y: 3.25, w: 6.0, h: 0.32, fontFace: FONT, fontSize: 9.5, align: "left", margin: 0,
       valign: "top", lineSpacingMultiple: 1.1 });

  // Right: complexity bar chart (image)
  const rx = 6.85;
  s.addText("② 复杂度：每帧平均 𝔉₂ₘ 运算量", { x: rx, y: 1.15, w: W - rx - M, h: 0.3,
    fontFace: FONT, fontSize: 12.5, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  const cxAR = 1260 / 700, cxH = 1.95, cxW = cxH * cxAR;
  s.addImage({ path: path.join(FIG, "complexity_v3.png"), x: rx, y: 1.46, w: cxW, h: cxH });
  s.addText("Direct 相对 Conv 省的是周期数(时延)，总运算量(ops)接近——t 小时 Chien 只贡献 O(t) 个乘法。Direct 价值在关键路径深度，不在总运算量。", {
    x: rx, y: 1.46 + cxH + 0.05, w: W - rx - M, h: 0.55, fontFace: FONT, fontSize: 9.5, color: INK,
    align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.15 });

  // Bottom: hardware resource table — multiplier count / LUT size (native table, synthesized from report.tex/report_v2.tex)
  s.addText("③ 硬件资源开销：GF 乘法器数量 · LUT 大小", { x: M, y: 3.65, w: W - 2 * M, h: 0.3,
    fontFace: FONT, fontSize: 12.5, bold: true, color: GREEN_DARK, align: "left", margin: 0 });
  s.addTable([
    [
      { text: "硬件资源", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
      { text: "Config 1 (KP4, GF(2¹⁰)/GF(2⁸))", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
      { text: "Config 2 (n=255, GF(2⁸))", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
      { text: "说明", options: { fill: { color: GREEN_DARK }, color: WHITE, bold: true } },
    ],
    ["GF 乘法器阵列", "n 路并行（RS: 544 路 / BCH: 144 路，两级不共享物理阵列——不同域）", "n=255 路并行（RS/BCH 同域 GF(2⁸)，v2 可分时复用同一套）",
      "Chien 搜索 + Direct 求值均需并行求值 Λ(α⁻ⁱ)"],
    ["BCH-Direct LUT 大小", "2⁸ = 256 项（t=1 退化为 log 表，非二次查表）", "2⁸ = 256 项 ×(2 根/项) = 512 条目",
      "t=2：k=X²+X 索引，存 (X,X+1)"],
    ["LUT 存储 (ROM/BRAM)", "log 表 256×8bit ≈ 2 kbit（复用 α 幂表）", "256×(2×8bit) ≈ 4 kbit",
      "n 越小 LUT 越小，MUX 深度越浅"],
    ["共享后节省", "Direct 比 Conv 硬件面积降约 62%（Lagendijk 2026）", "同左 + syndrome/乘法器阵列分时复用省 2 cyc",
      "面积节省来自去掉 BM 迭代器 + Chien 专用阵列"],
  ], {
    x: M, y: 3.98, w: W - 2 * M, h: 2.15, fontFace: FONT, fontSize: 9.8, valign: "middle",
    border: { type: "solid", color: GREEN_75, pt: 0.75 }, autoPage: false,
    colW: [1.85, 3.85, 3.85, W - 2 * M - 1.85 - 3.85 - 3.85],
  });

  card(s, M, 6.35, W - 2 * M, 0.85, GREEN_WASH, GREEN_75);
  s.addText([
    { text: "结论：", options: { bold: true, color: GREEN_DARK } },
    { text: "Direct + Lagrange 硬件共享是唯一能同时满足两组 config 10% 时延 KPI 的方案；总运算量(ops)几乎不受 Direct 影响，硬件收益集中在关键路径深度与面积（LUT/乘法器阵列复用），而非算力总量。", options: { color: INK } },
  ], { x: M + 0.25, y: 6.42, w: W - 2 * M - 0.5, h: 0.7, fontFace: FONT, fontSize: 11.5,
       align: "left", margin: 0, valign: "top", lineSpacingMultiple: 1.15 });

  s.addNotes("时延/复杂度/硬件资源三件事。①时延模型：RS-BM 2t+3周期，BCH-Direct比Conv省一半以上周期；②复杂度：每帧F_2m运算量，Direct和Conv总ops接近，价值在时延不在算力；③硬件资源表：GF乘法器阵列路数、Direct LUT大小(256项/512条目)、ROM存储量、共享后面积节省约62%。结论：Direct+Lagrange共享是唯一同时满足两组config KPI的方案。");
})();

pres.writeFile({ fileName: OUT }).then((f) => {
  console.log("WROTE", f);
}).catch((e) => { console.error("ERR", e); process.exit(1); });
