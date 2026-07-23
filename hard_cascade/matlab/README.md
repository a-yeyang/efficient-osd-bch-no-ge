# Matlab port — BCH t=2 硬判决译码器

## 内容

本目录包含论文 Lagendijk 2026 的**两个 BCH t=2 硬判决译码器**的 Matlab 实现，
作为 Python 主实现的对照移植。

| 文件 | 说明 |
|---|---|
| `GF_init.m` | GF(2^m) 有限域初始化（EXP/LOG 查表）|
| `gf_add.m` | GF 加法（XOR）|
| `gf_mul.m` | GF 乘法（查 EXP/LOG 表） |
| `gf_div.m` | GF 除法 |
| `bch_syndromes.m` | BCH 奇数 syndrome 计算 S_1, S_3 |
| `bch_decode_conventional.m` | BCH t=2 Conventional (BM+Chien 简化版) |
| `bch_decode_direct.m` | BCH t=2 Direct root finding + LUT (Lagendijk §III-A) |
| `build_lut_A.m` | 构建 {}_A LUT（A(Y)=Y²+Y+k 的根表）|
| `smoke_test.m` | 单元测试脚本 |

## 运行方式

```matlab
% 在 Matlab / Octave 中：
addpath('matlab');
smoke_test
```

期望输出：0/1/2-error 全部通过，3-error 意外恢复次数为 0。

## Python vs Matlab 对照

Python 完整实现（含 PAM4 modem, AWGN 信道, RS 编解码, 级联链路, MC 仿真）
在同项目的 `hc_src/` 下。Matlab 版本聚焦于 BCH 译码器的算法正确性验证，
不包含仿真驱动、绘图等外围代码。

用于向导师提交时，`.m` 文件可以独立编译运行；如果需要完整仿真的 Matlab 版本，
可以按需扩展。
