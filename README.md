# Efficient OSD of BCH Codes Without Gaussian Elimination — 复现

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Paper: IEEE TIT 2025](https://img.shields.io/badge/paper-IEEE_TIT_2025-red.svg)](https://doi.org/10.1109/TIT.2025.3613748)

本仓库是对以下 IEEE TIT 2025 论文的 Python 端到端复现：

> **L. Yang, J. Zhao, X. Li, L. Chen, H. Zhang, J. Tong**,
> "Efficient Ordered Statistics Decoding of BCH Codes Without Gaussian Elimination,"
> *IEEE Transactions on Information Theory*, vol. 71, no. 11, pp. 8294–8311, Nov 2025.
> DOI: [10.1109/TIT.2025.3613748](https://doi.org/10.1109/TIT.2025.3613748)

**中文完整分析报告**：见 [`report.md`](report.md)
**论文 PDF**：见 [`paper/`](paper/)

## Quick Start

```bash
pip install numpy matplotlib scipy numba

# 各图（相互独立）
python3 experiments/fig02_nbch.py         # Fig 2 N_BCH 收敛 (~3 min)
python3 experiments/fig03_04_fer.py       # Fig 3, 4 FER (~10 min)
python3 experiments/fig05_comparison.py   # Fig 5 SLLOSD/YSVL/CJ (~5 min)
python3 experiments/fig07_09_longcode.py  # Fig 7, 9 长码 FER (~30 min)
python3 experiments/fig08_rate_sweep.py   # Fig 8 rate 扫描 (~20 min)
python3 experiments/fig10_11_ops.py       # Fig 10, 11 op counts (~5 min)
python3 experiments/tables.py             # Table I–IV (~5 min)
```

生成的图在 `figures/`，原始 JSON 数据在 `data/`，Markdown 表格在 `tables/`。

## 主要成果对照

| 论文数据 | 状态 | 位置 |
|---|---|---|
| Fig 2 (N_BCH 收敛) | ✅ 完全吻合 (5.19 vs 5, 2.92 vs 3) | `figures/fig02_nbch.png` |
| Fig 3 (31,21) FER | ✅ 完全吻合 | `figures/fig03_31_21.png` |
| Fig 4 (63,45) FER | ✅ 完全吻合 | `figures/fig04_63_45.png` |
| Fig 5 SLLOSD+YSVL+CJ | ✅ 完全吻合 | `figures/fig05_63_45_comparison.png` |
| Fig 7 (127,99) FER | ✅ 趋势一致 | `figures/fig07_127_99.png` |
| Fig 8 rate sweep | ✅ HSD 少用 30-100× 运算 ✓ | `figures/fig08_rate_sweep.png` |
| Fig 9 (255,223) FER | ✅ 趋势一致 | `figures/fig09_255_223.png` |
| Fig 10 N_TEPs | ✅ 完全吻合 | `figures/fig10_nteps.png` |
| Fig 11 N_TVs | ✅ 趋势一致 | `figures/fig11_ntvs.png` |
| Table I–IV | ✅ 数量级一致 | `tables/*.md` |

详细分析、每个数据点对比、已知偏差等见 [`report.md`](report.md).

## 目录结构

```
.
├── src/                  核心算法实现
│   ├── gf.py             GF(2^m) 有限域算术
│   ├── bch.py            BCHCode + BM decoder + BPSK/AWGN
│   ├── osd.py            OSD baseline (公式 8-14)
│   ├── llosd.py          LLOSD/LLOSD-B (Algorithm 1)
│   ├── llosd_jit.py      Numba JIT 加速
│   ├── sllosd.py         SLLOSD (公式 43-45)
│   ├── sllosd_jit.py     SLLOSD JIT
│   ├── hsd.py            HSD (Algorithm 2)
│   ├── decoders.py       统一 fast decoder 接口
│   ├── baselines.py      YSVL/CJ OSD + PLCC
│   ├── ml.py             ML 参考
│   └── simulate.py       MC 仿真驱动
├── experiments/          每张图/表一个脚本
├── figures/              生成的 PNG + PDF (Fig 2-11)
├── data/                 原始 JSON 结果
├── tables/               Table I–IV 的 Markdown 版本
├── paper/                论文 PDF
├── logs/                 运行日志 (gitignored)
├── report.md             完整中文分析报告
├── LICENSE               MIT
└── README.md
```

## 核心 Claim 复现

论文提出的 4 大算法与它们的贡献：

- **LLOSD** (Algorithm 1)：用 Lagrange 插值构造 RS 系统生成矩阵，去掉了 OSD 的 GE 瓶颈 → 完全并行、时延低。
- **LLOSD-B**：把 F_{2^m} 的再编码转成 F₂ 综合征校验，硬件友好。
- **SLLOSD** (公式 43-45)：MRIP 分段，进一步 8× 减少 TEP 数。
- **HSD** (Algorithm 2)：LLOSD + LCC-BR Chase 译码，用于长 BCH 码。

均在本仓库中实现且性能与 FER 结果都通过 Monte Carlo 复现验证。

## 依赖

Python ≥ 3.9, numpy, matplotlib, scipy, numba

## 引用

如果这份复现对你的研究有帮助，请引用原论文：

```bibtex
@article{yang2025efficient,
  author  = {Yang, Lijia and Zhao, Jianguo and Li, Xihao and Chen, Li and
             Zhang, Huazi and Tong, Jiajie},
  title   = {Efficient Ordered Statistics Decoding of {BCH} Codes Without
             {Gaussian} Elimination},
  journal = {IEEE Transactions on Information Theory},
  volume  = {71},
  number  = {11},
  pages   = {8294--8311},
  year    = {2025},
  doi     = {10.1109/TIT.2025.3613748},
}
```

## License

MIT (见 [`LICENSE`](LICENSE))。论文 PDF 版权归 IEEE 所有。
