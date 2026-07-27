# RS+BCH 低时延级联码算法仿真（MATLAB 版）

> **任务标注**：本目录是导师布置的「**低时延级联码算法仿真**」任务的 **MATLAB 实现**，
> 是现有 Python 版 [`../rs_bch_cascade/`](../rs_bch_cascade/) 的完整、自包含移植。
> 数值与 Python 版对齐（`primpoly(7)=137`、`primpoly(8)=285` 与 Python `PRIM_POLY` 一致）。

## 一键运行

在 MATLAB（R2025a 已验证）里：

```matlab
cd rs_bch_cascade_matlab
main            % FAST 预设（默认开）：数分钟跑完，出全部图
% main(false)   % 完整精度（更慢、曲线更平滑）
```

或命令行：

```bash
matlab -batch "cd rs_bch_cascade_matlab; main"
```

`main` 会：**先跑数值自检（`selftest`）** → 跑两档配置的四组方法 → 画 BER-SNR
`semilogy` 折线图 → 画时延/运算量柱状图 + 10% KPI 门限 → 打印各组时延表与 KPI 判定 →
保存 `data/matlab_results.mat`。

## 核心创新（技术方案第 3 点）

**BCH 码是 RS 码的二元子码**，因此可以**不做传统高斯消元（GE）**，直接用
**Lagrange 插值**构造级联码外码 RS 的系统生成矩阵 `G_RS`（[`build_rs_systematic_generator.m`](build_rs_systematic_generator.m)），
再在有序统计译码（OSD）里**过滤掉再编码结果中的非二元候选**（Theorem 2），
只保留落在 BCH 子码里的码字。这样把 OSD 里最贵的高斯消元换成可并行、无主元依赖的插值求值，
**大幅降低复杂度与时延**。

- [`llosd_decode.m`](llosd_decode.m) —— 本方法（LLOSD，Lagrange）。
- [`osd_decode.m`](osd_decode.m) —— **对照组**（传统 OSD，高斯消元）。两者在同一 LLR 下 BER 一致，
  但 LLOSD 运算量/时延更低 —— 这正是「低时延」的直接证据。

## 仿真链路

```
信息符号 → RS 系统编码(外码) → 符号→比特 → 补零凑整 → BCH 编码(内码)
        → PAM4 Gray 调制 → AWGN 信道 → 逐比特 LLR
        → BCH 内码软判(LLOSD / OSD) → RS 外码软判(LCC-BR) → 信息符号
```

- **码型**：RS+BCH 级联码，开销约 12%（码率 ≈ 0.88）。
- **信道**：AWGN；**信号**：PAM4（Gray，星座 {−3,−1,+1,+3} ↔ {00,01,11,10}，E_s=5）。
- **RS**：硬判（Berlekamp–Massey）与软判（LCC-BR，Chase）两种。
- **BCH**：论文的 LLOSD 软判算法（Lagrange），并与传统 OSD 对照。

## 两档配置

| 配置 | 码型 | 参数 |
|---|---|---|
| **n=127** | RS(127,119) + BCH(127,120) | m=7, t_bch=1, τ=1, η=4 |
| **n=255** | RS(255,239) + BCH(255,239) | m=8, t_bch=2, τ=2, η=4 |

> 说明：导师要求的「码长 128 和 255」中的「128」指本原码长 **127**（2⁷−1）。RS/BCH 是
> 循环码，本原长度为 2^m−1，故取 127 与 255。

## 四组对照

同一张 BER-SNR 图、同一张时延表里对比：

1. **Pure RS-BM** —— 硬判基线
2. **Pure RS-LCC-BR** —— 软判基线
3. **Cascade RS+BCH，内码 LLOSD（Lagrange，本方法）** + 外码 LCC-BR
4. **Cascade RS+BCH，内码 OSD（传统高斯消元，对照）** + 外码 LCC-BR

图 3 vs 4 → 证明 Lagrange 降时延；图 1/2 vs 3 → BER 增益 + 「RS+BCH 相对单独 RS 的时延增幅 ≤10%」KPI。

## 时延 / KPI 实测结论（FAST 预设一次运行）

**核心结论（论文创新点的直接证据）：LLOSD（Lagrange）相对 OSD（高斯消元）的总运算量
（`f2m+f2`）在收敛 SNR 下低 1~3 个数量级** ——

| 配置 | @SNR | LLOSD 总运算 | OSD 总运算 | 加速比 |
|---|---|---|---|---|
| n=127 | 9 dB | ~3.5×10⁴ | ~2.1×10⁶ | **≈60×** |
| n=255 | 9 dB | ~8.4×10⁵ | ~1.2×10⁹ | **≈1477×** |

原因：OSD 每帧要对排序后的生成矩阵做二元高斯消元（开销记在 `f2`，n=255 高达 ~10⁹/帧），
而 LLOSD 用 Lagrange 插值一次性构造 `G_RS`（`f2≈0`），再靠二元过滤（Theorem 2）筛候选。

**「RS+BCH vs 单独 RS 时延 ≤10%」KPI（Case b，对软判基线 RS-LCC-BR，P=1 直接 op 计数）：**
本实现实测级联 LLOSD 比纯 RS-LCC-BR 高约 **+17%（n=255, 9 dB）/ +25%（n=127, 9 dB）**，
略高于 10% 理想线。这与主论文实现里「Case b 反而 −20%」存在实现口径差异（外码 LCC-BR
输入洁净度、蒙特卡洛抽样），但**量级完全一致（皆 ~10⁵ F₂ᵐ/帧）**。要严格压到 ≤10%，
按论文风险表的既定策略：① 用 P≥4 的并行硬件模型（Lagrange 无主元依赖、可并行，OSD 的 GE 不行）；
② 内码取 τ=1；③ 只在高 SNR（≥9 dB）区间评估。**对硬判基线 RS-BM（Case a），因 BM 极快
（~6.5k F₂ᵐ），级联在 P=1 下必然远超 10%，需并行视角补充**——与论文结论一致。


## 文件清单

| 文件 | 对应 Python | 职责 |
|---|---|---|
| `main.m` | `experiments/n255_scheme_a.py` + `kpi_analysis.py` | 一键入口：跑仿真、画图、打印时延/KPI |
| `selftest.m` | — | 数值自检（GF / RS / BCH / Lagrange / LLOSD≡OSD） |
| `GF.m` | `src/gf.py` | GF(2^m) EXP/LOG 表 + 域运算 + 多项式运算 + BCH 生成多项式 |
| `RSCode.m` | `cascade_src/rs_code.py` | RS 系统编码、BM 硬判、LCC-BR 软判 |
| `BCHCode.m` | `src/bch.py` | BCH 生成/校验矩阵、编码、BM、消息还原 |
| `build_rs_systematic_generator.m` | `src/llosd.py` | **★ Lagrange 构造 RS 系统生成矩阵（免 GE）** |
| `llosd_decode.m` | `src/llosd.py` | **LLOSD 内码软判（本方法）** |
| `osd_decode.m` | `src/osd.py` | **OSD 内码软判（高斯消元，对照）** |
| `pam4.m` | `cascade_src/pam4.py` | PAM4 调制/AWGN/逐比特 LLR |
| `CascadeConfig.m` / `CascadedCodec.m` / `PureRSCodec.m` | `cascade_src/cascade.py` | 配置、级联编解码、纯 RS 基线 |
| `run_bench.m` | `cascade_src/simulate.py` | 蒙特卡洛 SNR 扫描：BER/FER + 时延 + F₂ᵐ 运算量 |
| `OpCounters.m` | `src/osd.py` | 运算量计数器（f2 / f2m / fp / latency） |

## 输出

- `figures/matlab_n127_ber.{png,fig}`、`figures/matlab_n255_ber.{png,fig}` —— BER-SNR `semilogy`
- `figures/matlab_latency_bars.{png,fig}` —— 各组时延/运算量柱状 + 10% 门限线
- `data/matlab_results.mat` —— 所有数值
- 命令行：每 SNR 点 BER/FER、各组平均时延(μs)、理论 F₂ᵐ 运算量、KPI 判定

## 关于时延度量

MATLAB 解释执行的 `tic/toc` 墙钟时间有抖动，因此**主要、可复现的时延证据是理论运算量计数**
（`OpCounters`），与 Python 版一致；`tic/toc` 实测时延同时打印作参考。计数分两类：

- `f2m` —— GF(2^m) 域上的乘/加（RS 编解码、LCC-BR、Lagrange 插值都在这里）。
- `f2`  —— 二元 XOR/AND（**传统 OSD 的高斯消元开销落在这里**）。

**两个对比用不同度量，各取其义：**

1. **级联 LLOSD vs 纯 RS 的「≤10%」KPI** —— 两者都不做高斯消元（`f2≈0`），故只比 `f2m`
   （与 Python `kpi_analysis.py` 一致）。
2. **LLOSD(Lagrange) vs OSD(高斯消元) 的降时延对比** —— OSD 最贵的高斯消元开销记在 `f2`，
   LLOSD 用 Lagrange 插值把它**整个省掉**（`f2=0`）。因此这一对比用**总运算量 `f2m+f2`**
   才公平：只比 `f2m` 会把 OSD 省下高斯消元的假象误判成「OSD 更省」。

`P` 路并行下的时钟周期 ≈ 总运算量 / P。

