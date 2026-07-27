"""Generate custom technical diagrams for the RS+BCH cascade PPT (Chinese labels).

Produces (into ./assets/):
  A fig_rs_bch_basics.png   —  RS(符号级) vs BCH(比特级) 互补
  B fig_cascade_struct.png  —  级联码构造：外码RS + 内码BCH
  C fig_pipeline.png        —  完整链路：编码→PAM4→AWGN→译码
  D fig_lagrange.png        —  Lagrange 插值取代高斯消元的原理
  E fig_latency.png         —  为什么延迟低：串行关键路径 vs 并行
  F fig_interleave.png      —  交织如何打散突发错误
White-dominant, light-green-bordered palette matching Beijing Institute of
Technology's official emblem colors (VI manual A1-08: PANTONE 347C green,
349C dark green, 1535C brown).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
from matplotlib.lines import Line2D
import numpy as np

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "Songti SC", "STHeiti",
                                   "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

# ---- palette (BIT official VI colors: A1-08, PANTONE 347C/349C/1535C) ----
GREEN      = "#009A44"   # PANTONE 347C 100% — primary school green
GREEN_75   = "#14AE68"   # 75% tint
GREEN_50   = "#89C997"   # 50% tint
GREEN_25   = "#CBE4CE"   # 25% tint — light card wash
GREEN_DARK = "#015C31"   # PANTONE 349C — secondary dark green
BROWN      = "#A23E0A"   # PANTONE 1535C — secondary brown (error/caution accent)
BROWN_TINT = "#FBEEE6"
INK    = "#222222"
GRAY   = "#5C6670"
GRAYL  = "#9AA2AA"
WHITE  = "#FFFFFF"

ASSETS = __import__("os").path.join(__import__("os").path.dirname(__file__), "assets")
__import__("os").makedirs(ASSETS, exist_ok=True)


def _box(ax, x, y, w, h, text, color, tc=None, fs=13, bold=True,
         rad=0.02, lw=2.2, ha="center", va="center", filled=False):
    """White box with colored border + colored text by default.
    filled=True gives the old solid-fill + white-text look (for compact
    cell/heat-strip elements where a filled block reads better)."""
    if filled:
        fc, ec, use_tc = color, color, (tc or WHITE)
    else:
        fc, ec, use_tc = WHITE, color, (tc or color)
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.008,rounding_size={rad}",
                       fc=fc, ec=ec, lw=lw, mutation_aspect=1)
    ax.add_patch(p)
    if text:
        ax.text(x + w / 2 if ha == "center" else x + 0.04,
                y + h / 2, text, ha=ha, va=va, color=use_tc, fontsize=fs,
                fontweight="bold" if bold else "normal", zorder=5, wrap=True)
    return p


def _arrow(ax, x1, y1, x2, y2, color=INK, lw=2.2, style="-|>", ms=18, ls="-"):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, linestyle=ls, shrinkA=2, shrinkB=2,
                        zorder=4)
    ax.add_patch(a)


def _clean(ax, xlim, ylim):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.axis("off"); ax.set_aspect("equal")


# ============================================================ A. RS vs BCH
def fig_rs_bch_basics():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.3))
    fig.patch.set_facecolor(WHITE)

    # ---- LEFT: RS = symbol level ----
    _clean(axL, (0, 10), (0, 8))
    axL.text(5, 7.5, "RS 码：符号级纠错（外码）", ha="center", fontsize=16.5,
             fontweight="bold", color=GREEN_DARK)
    axL.text(5, 6.75, "在 GF($2^m$) 上运算，每 $m$ 比特 = 1 个符号",
             ha="center", fontsize=12, color=GRAY)
    # a codeword = sequence of symbols; some symbols hit by a burst
    n = 10
    x0, w, y = 0.55, 0.88, 4.6
    for i in range(n):
        err = i in (4, 5, 6)          # a burst spanning symbols
        color = BROWN if err else GREEN_DARK
        _box(axL, x0 + i * w, y, w * 0.86, 1.0, f"S{i}", color, fs=11, filled=True)
    axL.text(x0 + 5 * w, y + 1.55, "一段突发错误", ha="center", fontsize=11.5,
             color=BROWN, fontweight="bold")
    _arrow(axL, x0 + 5 * w, y + 1.4, x0 + 5 * w, y + 1.05, color=BROWN, lw=1.8)
    axL.text(5, 3.35, "整段突发只落在少数几个符号里\n→ 只算几个“符号错”，RS 可纠 $t$ 个符号",
             ha="center", fontsize=12, color=INK)
    axL.text(5, 1.6, "★ 抗突发：一个符号内错几比特，对 RS 只算 1 个错",
             ha="center", fontsize=12, color=GREEN, fontweight="bold")

    # ---- RIGHT: BCH = bit level ----
    _clean(axR, (0, 10), (0, 8))
    axR.text(5, 7.5, "BCH 码：比特级纠错（内码）", ha="center", fontsize=16.5,
             fontweight="bold", color=GREEN_DARK)
    axR.text(5, 6.75, "在 GF(2) 上运算，直接对每个比特纠错",
             ha="center", fontsize=12, color=GRAY)
    nb = 16
    x0, w, y = 0.35, 0.57, 4.6
    rng = [3, 11]                     # scattered single-bit errors
    for i in range(nb):
        err = i in rng
        color = BROWN if err else GREEN_75
        _box(axR, x0 + i * w, y, w * 0.82, 1.0, "", color, fs=9, rad=0.015, filled=True)
    axR.text(5, 3.35, "对零散的随机比特错误很高效\n每个码字纠 $t$ 个比特错",
             ha="center", fontsize=12, color=INK)
    axR.text(5, 1.6, "× 怕突发：连续比特错全落进 1 个码字 → 超出 $t$ 失败",
             ha="center", fontsize=12, color=BROWN, fontweight="bold")

    fig.tight_layout()
    fig.savefig(f"{ASSETS}/fig_rs_bch_basics.png", dpi=200,
                bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


# ============================================================ B. cascade structure
def fig_cascade_struct():
    fig, ax = plt.subplots(figsize=(11.5, 5.0))
    fig.patch.set_facecolor(WHITE)
    _clean(ax, (0, 13.3), (0, 6.0))

    ax.text(6.65, 5.6, "级联码构造：外码 RS  ⊕  内码 BCH", ha="center",
            fontsize=17, fontweight="bold", color=GREEN_DARK)

    y = 2.9
    # message
    _box(ax, 0.3, y, 1.7, 1.2, "信息符号\nk 个", GRAY, fs=12)
    _arrow(ax, 2.0, y + 0.6, 2.55, y + 0.6, color=GRAY)
    # RS outer
    _box(ax, 2.55, y, 2.35, 1.2, "外码 RS 编码\n(符号级校验)", GREEN_DARK, fs=12.5)
    _arrow(ax, 4.9, y + 0.6, 5.45, y + 0.6, color=GREEN_DARK)
    # serialize
    _box(ax, 5.45, y, 2.1, 1.2, "符号→比特\n串行化", GRAY, fs=12)
    _arrow(ax, 7.55, y + 0.6, 8.1, y + 0.6, color=GRAY)
    # split to BCH blocks + BCH enc
    _box(ax, 8.1, y, 2.35, 1.2, "分块 → 内码\nBCH 编码", GREEN, fs=12.5)
    _arrow(ax, 10.45, y + 0.6, 11.0, y + 0.6, color=GREEN)
    _box(ax, 11.0, y, 2.0, 1.2, "编码比特\n→ 信道", GRAY, fs=12)

    # annotation band under RS/BCH
    ax.annotate("", xy=(4.9, 2.55), xytext=(2.55, 2.55),
                arrowprops=dict(arrowstyle="-", color=GREEN_DARK, lw=1.2))
    ax.text(3.72, 2.2, "抗突发 · 符号级", ha="center", fontsize=11, color=GREEN_DARK)
    ax.annotate("", xy=(10.45, 2.55), xytext=(8.1, 2.55),
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=1.2))
    ax.text(9.27, 2.2, "纠随机错 · 比特级", ha="center", fontsize=11, color=GREEN)

    # code numbers box
    ax.text(6.65, 1.15,
            "本方案码型（≈0.88 码率）：外码 RS(255,k) / GF($2^8$)，"
            "内码 BCH(n,k,t) / GF($2^7$)～GF($2^8$)",
            ha="center", fontsize=12, color=INK,
            bbox=dict(boxstyle="round,pad=0.4", fc=GREEN_25, ec=GREEN, lw=1.3))
    ax.text(6.65, 0.35, "外码保护内码译码后的残余错误；两级串联把纠错能力叠加",
            ha="center", fontsize=11, color=GRAY, style="italic")

    fig.tight_layout()
    fig.savefig(f"{ASSETS}/fig_cascade_struct.png", dpi=200,
                bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


# ============================================================ C. full pipeline
def fig_pipeline():
    fig, ax = plt.subplots(figsize=(11.5, 5.4))
    fig.patch.set_facecolor(WHITE)
    _clean(ax, (0, 13.3), (0, 6.4))
    ax.text(6.65, 6.0, "完整链路：编码 → PAM4 → AWGN 信道 → 译码", ha="center",
            fontsize=17, fontweight="bold", color=GREEN_DARK)

    w, h = 2.25, 1.2
    # ---- top row: TX (left→right) ----
    yt = 4.1
    tx = [("信息\n比特", GRAY), ("RS 外码\n编码", GREEN_DARK), ("BCH 内码\n编码", GREEN),
          ("PAM4\n调制", GREEN_75)]
    xs = [0.35, 2.95, 5.55, 8.15]
    for i, (t, c) in enumerate(tx):
        _box(ax, xs[i], yt, w, h, t, c, fs=13)
        if i:
            _arrow(ax, xs[i - 1] + w, yt + h / 2, xs[i], yt + h / 2, color=GRAY)
    # channel box
    xch, wch = 10.9, 2.05
    _box(ax, xch, yt, wch, h, "AWGN\n信道", BROWN, fs=13)
    _arrow(ax, xs[-1] + w, yt + h / 2, xch, yt + h / 2, color=GRAY)
    ax.text(xch + wch / 2, yt + h + 0.32, "噪声 / 突发", ha="center", fontsize=11,
            color=BROWN)

    # down connector
    cx = xch + wch / 2
    _arrow(ax, cx, yt, cx, 2.9, color=BROWN, lw=2.4)

    # ---- bottom row: RX (right→left) ----
    yb = 1.7
    rx = [("PAM4 解调\n(软信息 LLR)", GREEN_75), ("BCH 内码译码\nLLOSD", GREEN),
          ("RS 外码译码\nLCC-BR / BM", GREEN_DARK), ("信息\n比特", GRAY)]
    xsb = [10.9, 7.45, 4.0, 0.55]
    for i, (t, c) in enumerate(rx):
        _box(ax, xsb[i], yb, w, h, t, c, fs=12)
        if i:
            _arrow(ax, xsb[i - 1], yb + h / 2, xsb[i] + w, yb + h / 2, color=GRAY)
    # channel down-arrow into PAM4 demod (top-center of first RX box)
    _arrow(ax, cx, 2.9, xsb[0] + w / 2, yb + h, color=BROWN, lw=2.2)

    ax.text(6.65, 0.5, "两级译码是延迟与开销的核心 —— 内码用 LLOSD，外码用 "
            "LCC-BR，二者可共享 Lagrange 结构",
            ha="center", fontsize=11.5, color=INK, style="italic",
            bbox=dict(boxstyle="round,pad=0.35", fc=GREEN_25, ec=GREEN, lw=1.3))

    fig.tight_layout()
    fig.savefig(f"{ASSETS}/fig_pipeline.png", dpi=200,
                bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


# ============================================================ D. Lagrange principle
def fig_lagrange():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.8),
                                   gridspec_kw={"width_ratios": [1, 1.15]})
    fig.patch.set_facecolor(WHITE)

    # LEFT: traditional GE (serial chain) — muted gray = "old / inferior" method
    _clean(axL, (0, 10), (0, 10))
    axL.text(5, 9.4, "传统 OSD：高斯消元", ha="center", fontsize=15.5,
             fontweight="bold", color=GRAY)
    axL.text(5, 8.6, "对 k×k 矩阵逐行消元", ha="center", fontsize=12, color=GRAY)
    # serial chain of pivots
    for i in range(5):
        yy = 7.4 - i * 1.35
        _box(axL, 3.1, yy, 3.8, 0.95, f"消第 {i+1} 行 (依赖上一行)", GRAY, fs=11)
        if i:
            _arrow(axL, 5.0, yy + 1.35, 5.0, yy + 0.95, color=GRAY, lw=1.8)
    axL.text(5, 0.55, "关键路径 ∝ k（串行、需选主元）\n难并行 → 延迟高",
             ha="center", fontsize=12, color=BROWN, fontweight="bold")

    # RIGHT: Lagrange interpolation (parallel) — green = "our method"
    _clean(axR, (0, 11.5), (0, 10))
    axR.text(5.75, 9.4, "本方法 LLOSD：Lagrange 插值", ha="center", fontsize=15.5,
             fontweight="bold", color=GREEN_DARK)
    axR.text(5.75, 8.6, "把重编码看成多项式插值，免高斯消元",
             ha="center", fontsize=12, color=GRAY)

    # k reliable points -> basis (parallel) -> evaluate (parallel)
    _box(axR, 0.5, 5.6, 2.5, 2.6, "k 个最可靠\n位置的值\n(插值节点)", GREEN_DARK, fs=11.5)
    # parallel basis functions
    for j in range(3):
        yy = 7.3 - j * 1.5
        _box(axR, 4.0, yy, 3.0, 1.1, f"基函数 $L_{j}(x)$", GREEN, fs=11)
        _arrow(axR, 3.0, 6.9 - 0.0, 4.0, yy + 0.55, color=GRAYL, lw=1.3)
    axR.text(5.5, 3.15, "…  各基函数彼此独立", ha="center", fontsize=11,
             color=GRAY)
    # evaluate
    _box(axR, 8.1, 5.6, 3.0, 2.6, "并行求值\n得所有位置\n的重编码符号", GREEN_75, fs=11.5)
    for j in range(3):
        yy = 7.3 - j * 1.5
        _arrow(axR, 7.0, yy + 0.55, 8.1, 6.9, color=GRAYL, lw=1.3)

    axR.text(5.75, 1.7, "基函数与求值都可并行、无主元依赖", ha="center",
             fontsize=12, color=INK)
    axR.text(5.75, 0.75, "★ 关键路径大幅缩短 → 延迟低，且可与外码共享中间结果",
             ha="center", fontsize=12, color=GREEN, fontweight="bold")

    fig.tight_layout()
    fig.savefig(f"{ASSETS}/fig_lagrange.png", dpi=200,
                bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


# ============================================================ E. latency
def fig_latency():
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    fig.patch.set_facecolor(WHITE)
    _clean(ax, (0, 13.3), (0, 5.6))
    ax.text(6.65, 5.2, "为什么延迟低：消掉串行关键路径", ha="center",
            fontsize=16.5, fontweight="bold", color=GREEN_DARK)

    # timeline: traditional (long serial)
    ax.text(0.4, 4.3, "传统\n高斯消元", ha="center", fontsize=12.5, color=GRAY,
            fontweight="bold")
    tot = 11.4
    x0 = 1.7
    ncell = 12
    cw = tot / ncell
    for i in range(ncell):
        _box(ax, x0 + i * cw, 3.95, cw * 0.94, 0.7,
             "", GRAY, fs=8, rad=0.01, lw=1.2, filled=True)
    ax.annotate("", xy=(x0 + tot, 3.6), xytext=(x0, 3.6),
                arrowprops=dict(arrowstyle="<->", color=BROWN, lw=1.8))
    ax.text(x0 + tot / 2, 3.25, "关键路径长（每步依赖上一步）", ha="center",
            fontsize=11.5, color=BROWN)

    # timeline: ours (short parallel depth)
    ax.text(0.4, 1.9, "本方法\nLagrange+并行", ha="center", fontsize=12.5,
            color=GREEN_DARK, fontweight="bold")
    ncell2 = 3
    tot2 = tot * ncell2 / ncell
    for i in range(ncell2):
        _box(ax, x0 + i * cw, 1.55, cw * 0.94, 0.7, "", GREEN, fs=8, rad=0.01,
             lw=1.2, filled=True)
    ax.annotate("", xy=(x0 + tot2, 1.2), xytext=(x0, 1.2),
                arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.8))
    ax.text(x0 + tot2 + 0.15, 1.2, "关键路径短", ha="left", fontsize=11.5,
            color=GREEN, va="center", fontweight="bold")
    # the freed parallel work
    ax.text(x0 + tot2 + 1.6, 1.9, "其余计算并行展开（不在关键路径上）",
            ha="left", fontsize=11.5, color=GRAY, va="center")

    ax.text(6.65, 0.4,
            "延迟量化 = 关键路径上的 GF 运算级数 × 单级时钟；并行度 P 提升后，"
            "本方法关键路径显著短于逐行消元",
            ha="center", fontsize=11.3, color=INK, style="italic")

    fig.tight_layout()
    fig.savefig(f"{ASSETS}/fig_latency.png", dpi=200,
                bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


# ============================================================ F. interleaving
def fig_interleave():
    fig, (axT, axB) = plt.subplots(2, 1, figsize=(11.5, 6.4),
                                   gridspec_kw={"height_ratios": [1, 1]})
    fig.patch.set_facecolor(WHITE)

    # ---- TOP: without interleave ----
    _clean(axT, (0, 13.3), (0, 3.4))
    axT.text(6.65, 3.05, "不交织：突发错误全压在一个 BCH 码字里 → 超出 t → 失败",
             ha="center", fontsize=14, fontweight="bold", color=BROWN)
    ncw, cwlen = 6, 10
    bw, bh = 0.34, 0.55
    x0, y0 = 1.2, 1.1
    gap = 1.9
    for c in range(ncw):
        bx = x0 + c * gap
        for j in range(cwlen):
            err = (c == 2 and 2 <= j <= 8)      # whole burst in codeword 2
            fc = BROWN if err else GREEN_50
            axT.add_patch(Rectangle((bx + j * bw, y0), bw * 0.9, bh, fc=fc,
                          ec="white", lw=0.5))
        axT.text(bx + cwlen * bw / 2, y0 - 0.35, f"BCH#{c}", ha="center",
                 fontsize=10.5, color=GRAY)
    axT.text(x0 + 2 * gap + cwlen * bw / 2, y0 + 1.0, "突发", ha="center",
             fontsize=11.5, color=BROWN, fontweight="bold")

    # ---- BOTTOM: with w=1bit interleave ----
    _clean(axB, (0, 13.3), (0, 3.4))
    axB.text(6.65, 3.05, "交织(w=1bit)：矩阵按列写、按行读，突发被打散到多个码字 → 每个≤1错 → 全可纠",
             ha="center", fontsize=13.5, fontweight="bold", color=GREEN_DARK)
    for c in range(ncw):
        bx = x0 + c * gap
        for j in range(cwlen):
            # one error per codeword, staggered — the dispersed burst
            err = (j == 4)
            fc = BROWN if err else GREEN_50
            axB.add_patch(Rectangle((bx + j * bw, y0), bw * 0.9, bh, fc=fc,
                          ec="white", lw=0.5))
        axB.text(bx + cwlen * bw / 2, y0 - 0.35, f"BCH#{c}", ha="center",
                 fontsize=10.5, color=GRAY)
    axB.text(6.65, y0 + 1.15,
             "同一段突发的相邻比特被送进不同 BCH 码字（每字至多 1 个错）",
             ha="center", fontsize=11.5, color=INK)

    fig.tight_layout()
    fig.savefig(f"{ASSETS}/fig_interleave.png", dpi=200,
                bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)


if __name__ == "__main__":
    fig_rs_bch_basics()
    fig_cascade_struct()
    fig_pipeline()
    fig_lagrange()
    fig_latency()
    fig_interleave()
    print("all diagrams written to", ASSETS)
