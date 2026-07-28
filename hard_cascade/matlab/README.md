# Matlab port — RS+BCH 硬判决级联码仿真

## 内容

本目录包含论文 Lagendijk 2026 硬判决级联方案的 Matlab 实现，与 `hc_src/`
下的 Python 主实现（v3：`shortened_codes.py` / `cascade_v3.py`）逐行对照移植。

### 有限域 / GF 基础

| 文件 | 说明 |
|---|---|
| `GF_init.m` | GF(2^m) 有限域初始化（EXP/LOG 查表），m=1..10 |
| `gf_add.m` / `gf_mul.m` / `gf_div.m` | GF 标量加/乘/除 |
| `gf_polymul.m` / `gf_polydivmod.m` | GF 多项式乘法 / 除法（余数） |

### BCH（内码，t=1/2，支持缩短码）

| 文件 | 说明 |
|---|---|
| `bch_syndromes.m` | BCH 奇数 syndrome 计算 S_1, S_3（全长）|
| `bch_decode_conventional.m` | BCH t=2 Conventional (BM+Chien 简化版，全长) |
| `bch_decode_direct.m` | BCH t=2 Direct root finding + LUT（全长，Lagendijk §III-A） |
| `build_lut_A.m` | 构建 {}_A LUT（A(Y)=Y²+Y+k 的根表）|
| `bch_init.m` | 构造 BCH(n,k,t) 生成多项式 + 系统形式生成矩阵 G，支持缩短 |
| `bch_encode.m` / `bch_extract_message.m` | 缩短码系统编码 / 取消息 |
| `bch_decode_conventional_shortened.m` | 缩短码 Conventional 译码封装（t=1 闭式 + t=2 复用全长译码器）|
| `bch_decode_direct_shortened.m` | 缩短码 Direct 译码封装（t=1 闭式 + t=2 复用全长译码器）|

### RS（外码，支持缩短码）

| 文件 | 说明 |
|---|---|
| `rs_init.m` | 构造 RS(n,k) 生成多项式（本原码）|
| `rs_encode_systematic.m` / `rs_extract_message.m` | 系统编码 / 取消息 |
| `rs_bm_decode.m` | Berlekamp-Massey + Chien + Forney 硬判决译码（全长）|
| `shortened_rs_init.m` / `shortened_rs_encode.m` / `shortened_rs_extract_message.m` / `shortened_rs_decode.m` | 缩短 RS 码封装（补零/去零）|

### PAM4 调制 + AWGN 信道（硬判决）

| 文件 | 说明 |
|---|---|
| `pam4_bits_to_symbols.m` | Gray 映射：比特 → PAM4 电平 {-3,-1,+1,+3} |
| `pam4_symbols_to_bits_hard.m` | 硬判决：最近电平 → 比特 |
| `pam4_sigma_from_ebn0.m` | Eb/N0 (dB) + 码率 → AWGN σ |
| `awgn_channel.m` | 加性高斯白噪声 |
| `pam4_channel_hard.m` | 调制 + 信道 + 硬判解调一体封装 |

### 级联编解码 + 一键仿真

| 文件 | 说明 |
|---|---|
| `cascade_config.m` | Config 1 (KP4) / Config 2 (255) 两组码型参数 |
| `cascade_init.m` | 构造级联 codec（RS 外码 + BCH 内码，域解耦）|
| `cascade_encode.m` | RS 编码 → 比特串行化（LSB-first）→ 分块 → BCH 编码 |
| `cascade_decode.m` | BCH 逐块译码 → 比特重组 → RS 译码 |
| `symbols_to_bits_lsb.m` | GF(2^m) 符号 → 比特（LSB-first），用于 BER 统计 |
| `smoke_test_v3.m` | v3 正确性自检（RS/BCH/级联往返），`main_hard_cascade` 启动时自动运行 |
| **`main_hard_cascade.m`** | **一键运行主函数**：编码 → PAM4 调制 → AWGN → 硬判解调 → BCH 内码译码 → RS 外码译码，输出 BER-SNR |
| `smoke_test.m` | 旧版单元测试脚本（仅 BCH t=2 全长译码器）|

## 运行方式

```matlab
% 在 Matlab / Octave 中，进入本目录后直接运行：
main_hard_cascade          % 快速模式：少量 Eb/N0 点 + 少量帧数，用于验证流程
main_hard_cascade(false)   % 完整模式：4:1:12 dB 扫描 + 200 帧/点，更接近论文曲线
```

一键运行内部流程：
1. 自动 `addpath` 并运行 `smoke_test_v3()` 正确性自检（RS 缩短码、BCH t=1/t=2
   缩短码、级联零噪声往返全部通过才继续）；
2. 对 Config 1 (KP4) 和 Config 2 (255) 两组码型，在每个 Eb/N0 点做 Monte-Carlo
   仿真：构造信息符号 → 级联编码 → PAM4 调制 → AWGN 信道 → 硬判解调 → BCH 内码
   译码 → RS 外码译码 → 统计 BER/FER；
3. 结果保存到 `matlab/results/main_hard_cascade_results.mat`，并尝试绘制
   BER-SNR 曲线保存为 `matlab/results/ber_snr_hard_cascade.png`（无图形环境时跳过绘图但不报错）。

## Python vs Matlab 对照

Python 完整实现（含 LLOSD/OSD 软判决、复杂度计数等）在同项目的 `hc_src/` 下，
是本项目算法正确性的权威口径。本 Matlab 版本聚焦硬判决主链路的一键仿真，
逐函数对照 `hc_src/shortened_codes.py`、`hc_src/cascade_v3.py`，并复用了
`rs_bch_cascade_matlab/`（姊妹项目）中已验证的 `GF.m` / `RSCode.m` / `pam4.m`
算法逻辑。由于本环境未安装 Matlab/Octave，新增代码未能在本地实际运行验证，
正确性依赖与 Python 参考实现的严格逐行对照；建议使用者在有 Matlab/Octave 的
环境中运行 `main_hard_cascade` 确认输出。
