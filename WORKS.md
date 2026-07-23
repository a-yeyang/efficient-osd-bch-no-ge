# 研究 work 清单与依赖图

## Work 01 — BCH 无高斯消元 OSD 论文复现

目标是复现 Yang et al. 的 BCH 软判译码。核心实现位于 `code/python/src`：

- `gf.py`：GF(2^m) EXP/LOG 表、有限域多项式与 BCH 生成多项式；
- `bch.py`：primitive narrow-sense BCH 编码、BM 硬判、BPSK/AWGN；
- `osd.py`：含 GE 的 OSD 基线；
- `llosd*.py`、`sllosd*.py`、`hsd.py`、`decoders.py`：论文方法与 Numba 加速入口；
- `baselines.py`、`ml.py`、`simulate.py`：基线和 Monte-Carlo 驱动；
- `code/python/experiments`：逐图、逐表脚本和报告生成脚本。

论文、Markdown/DOCX 复现报告和时延分析均位于 `docs/`；原始 JSON、图、表在 `assets/`。

## Work 02 — 软判 RS+BCH 级联与 Lagrange 共享

该研究在 Work 01 的 BCH/LLOSD 上扩展：`cascade_src/upstream.py` 显式导入 Work 01 的 `code/python/src`；本 work 自身实现 PAM4 调制/LLR、RS 编码与 BM/LCC-BR、级联系统、Lagrange 缓存及蒙特卡洛仿真。

方案 A 与 B 的功能输出相同；B 在内层 Lagrange 生成矩阵及外层 RS 代数中实际复用共享缓存（pairwise 差分表、分母与 alpha 查表），并将观测到的复用量动态折算为 F₂ᵐ 运算量节省。报告中也明确记载了该研究的 LCC-BR 是 Chase+BM 简化版，而非完整 Kötter–Vardy/BR 插值器。

## Work 03 — 硬判 RS+BCH 级联与 Direct t=2

该 work 显式复用 Work 01 的 GF/操作计数及 Work 02 的 RS/PAM4 组件。新实现为 BCH t=2 的 syndrome、BM+Chien conventional 解码、LUT direct 根查找、硬判级联和时钟周期模型。

`n=127` 使用 primitive BCH(127,·)，而不是扩展 BCH(128,·)；报告和代码中均保持此约束。

## MATLAB 交付原则

MATLAB 对照代码以算法可验证性为先：同一组 GF/BCH/RS/PAM4、译码器和级联通路共享测试向量。完整 Monte-Carlo 作图的运行时间取决于帧数；交付测试采用固定种子、小样本 smoke/profile，并将理论组合数、无噪声往返、可纠错错误图案、Python/MATLAB 函数级不变量作为通过标准。

## 重要复现边界

仓库原有报告已经声明：Work 01 中 YSVL/CJ 是性能等价的简化基线，HSD 用 BM/Chase 近似完整 LCC-BR；Work 02 的 shared cache 用来量化可重用有限域代数，当前以 F₂ᵐ 运算量统计评估其收益。这些不是目录整理或 MATLAB 翻译能够自动消除的理论/实现差异，测试会把它们标记为“近似路径”而不会误报为原论文的完整算法。
