# RS+BCH 级联码低时延译码——问题定义、算法原理、实验方案

## 一、问题定义

### 1.1 系统模型

考虑一个使用 **RS+BCH 级联码**的 FEC 系统，工作在 PAM4-AWGN 信道下。发送端和接收端结构如下：

```
消息 m ∈ F₂^K
     │
     ▼
┌────────────────┐
│ 外码 RS 编码器  │  m → x_RS ∈ F_{2^m}^{n_RS}
│ (n_RS, k_RS)   │
└────────────────┘
     │ (每 m 个二进制位组成一个 F_{2^m} 符号)
     ▼
┌────────────────┐
│ 内码 BCH 编码器 │  x_RS (bit view) → x_BCH ∈ F₂^{n_BCH * (n_RS/k_BCH')} ...
│ (n_BCH, k_BCH) │  (每 k_BCH 位 BCH 编码成 n_BCH 位)
└────────────────┘
     │
     ▼
┌────────────────┐
│ PAM4 modulator │  x_BCH → PAM4 symbols {-3, -1, +1, +3}
└────────────────┘
     │
     ▼
   AWGN Channel
     │
     ▼
┌────────────────┐
│  PAM4 demod +  │  y → per-bit LLR
│  bit-LLR calc  │
└────────────────┘
     │
     ▼
┌────────────────┐
│ 内码 BCH 译码  │  LLR → BCH 硬输出 (LLOSD 软判)
│ (LLOSD)        │
└────────────────┘
     │
     ▼
┌────────────────┐
│ 外码 RS 译码   │  BCH 输出 → RS 译码 (BM 硬判 或 LCC-BR 软判)
│ (BM / LCC-BR)  │
└────────────────┘
     │
     ▼
    m̂ ∈ F₂^K
```

### 1.2 参数选择

三组码率约 0.88 的配置：

| 方案 | 外码 RS | 内码 BCH | 总码率 | GF 域 |
|---|---|---|---|---|
| **n=255-A** | RS(255, 239) | BCH(255, 239) | 0.937 × 0.937 = **0.879** | GF(2^8) |
| **n=255-B** | RS(255, 235) | BCH(255, 247) | 0.921 × 0.969 = **0.891** | GF(2^8) |
| **n=128** | RS(127, 119) | BCH(127, 120) | 0.937 × 0.945 = **0.885** | GF(2^7) |

**注**: n=128 严格是 n=127 (primitive BCH)，若严格要求 n=128 则用 extended BCH。

### 1.3 信道模型

**PAM4 调制**：4 个电平 {-3, -1, +1, +3}，Gray 编码（相邻电平差 1 bit）：

| 电平 | Bit pair (b_1 b_0) |
|---|---|
| −3 | 00 |
| −1 | 01 |
| +1 | 11 |
| +3 | 10 |

**AWGN 信道**：`y = x + n`，其中 `n ~ N(0, σ²)`。

**平均符号能量**：`E_s = E[|x|²] = (9 + 1 + 1 + 9)/4 = 5`。

**每比特能量**：`E_b = E_s / 2` (PAM4 每符号 2 bit)。

**SNR 定义**：`Eb/N0 (dB) = 10·log10(E_b / (σ²·N0))`，其中 `N0 = 2σ²`（双边），或 `Eb/N0 = 10·log10(E_s / (2·2σ²))`。

**bit-LLR 计算**：对每个 PAM4 符号 y，两个 bit 的 LLR 通过 marginalize 另一个 bit 得到：

$$
\text{LLR}(b_i | y) = \log \frac{\sum_{s: b_i(s)=0} \exp(-\|y-s\|^2 / 2\sigma^2)}
                              {\sum_{s: b_i(s)=1} \exp(-\|y-s\|^2 / 2\sigma^2)}
$$

## 二、算法原理

### 2.1 核心创新点：Lagrange 共享

论文 IEEE TIT 2025 的 LLOSD 核心洞察：**BCH ⊂ RS**，所以 BCH 可以用 RS 的代数（Lagrange 插值）来译码。

我们的**扩展**：级联码译码时，**内码 BCH 的 LLOSD 与外码 RS 的 LCC-BR 都需要 Lagrange 插值**。因为两者工作在**同一个 GF(2^m)** 上，可以共享：

1. **α 幂表**：`α^0, α^1, ..., α^{n-1}`（一次计算，两次使用）
2. **Denominator products**：`∏_{j'≠i} (α^i − α^{j'})`（外码 RS 的 code locators 和内码 BCH 的 MRIP 上都可能出现同一位置）
3. **Lagrange 基函数 T_j(α^i)**：至少可以共享**结构**（相同的分母序列）

### 2.2 三级方案对比 (消融实验)

按导师 R1 的要求，做三级对比：

| 方案 | 内码 Lagrange | 外码 Lagrange | 共享代数结构 |
|---|---|---|---|
| **Baseline** | — | 否 (BM 硬判) | 无 (纯 RS) |
| **方案 A** | 是 (LLOSD) | 否 (LCC-BR 独立实现) | 部分（都用 Lagrange 但不共享缓存） |
| **方案 B** | 是 (LLOSD) | 是 (LCC-BR + Lagrange 共享) | 全共享 |

### 2.3 各个译码器的算法

**内码 BCH：LLOSD**（论文 Algorithm 1）
1. 输入：n 个 bit-LLR
2. 按 |LLR| 降序排列位置 → MRIP 集 Θ (前 k' 位)
3. 用 Lagrange 插值构造 RS 系统生成矩阵 G_RS
4. 枚举 TEP（权重 ≤ τ 的错误图案）
5. 对每个 TEP 做**二元过滤**（Theorem 2）
6. 计算通过过滤的候选的相关距离
7. 输出最小距离对应的 BCH 码字

**外码 RS-BM (硬判)**：Berlekamp-Massey + Chien search + Forney
- 输入：n_RS 个 F_{2^m} 符号（BCH 内码输出，每 m bit 一个符号）
- 输出：k_RS 个 F_{2^m} 消息符号
- 复杂度：O(n·t) F_{2^m} ops

**外码 RS-LCC-BR (软判)**：Xing-Chen-Bossert 2020
- 输入：n_RS 个 F_{2^m} 符号 + 每符号的可靠性（从 BCH LLOSD 传递过来）
- η 个最不可靠位置 (LRPs) → 2^η 个测试向量
- 每个测试向量做 Mulders-Storjohann basis reduction + partial root-finding
- 输出：最优 RS 码字候选

### 2.4 Lagrange 共享层的具体实现

**LagrangeCache 类** 提供：

```python
class LagrangeCache:
    def __init__(self, gf: GF, n: int):
        self.gf = gf
        self.n = n
        # α^0, α^1, ..., α^{n-1}
        self.alpha_pow = np.array([gf.pow(2, j) for j in range(n)])
        # ∏_{j'} (α^i - α^{j'}) 的完整表 (对所有 i, j') - 惰性计算
        self._denom_full = None

    def get_alpha_pow(self):
        return self.alpha_pow

    def denominator_product(self, i, exclude_set):
        """计算 ∏_{j' ∉ exclude_set, j' != i} (α^i - α^{j'})"""
        # 复用缓存
        ...

    def lagrange_basis(self, i, support_set):
        """T_j(α^i) for j in support_set - 复用 alpha 表 + denom cache"""
        ...
```

**共享的时机**：内码 LLOSD 完成后，把它的 LagrangeCache 传给外码 LCC-BR，外码复用其中的 alpha 表和 denom cache。

## 三、实验方案

### 3.1 KPI 定义 (导师 R2)

**KPI 1: BER-SNR 性能**
- 4 条曲线：Baseline (纯 RS 硬判) / Baseline (纯 RS 软判) / 方案 A / 方案 B
- 期望：级联码 (方案 A/B) 相比纯 RS 有 1-2 dB 编码增益

**KPI 2: 时延**
- 两个基线：
  - Case (a): `T_{cascade} ≤ 1.10 × T_{RS-BM}` (硬判基线)
  - Case (b): `T_{cascade} ≤ 1.10 × T_{RS-LCC-BR}` (软判基线)
- 两种测量口径：
  - F_{2^m} ops（软件仿真直接数）
  - Clock cycles（用 P ∈ {1, 4, 16, 64} 的抽象并行硬件模型）

### 3.2 仿真参数

| 参数 | 值 |
|---|---|
| SNR 范围 | Eb/N0 ∈ [3.0, 8.0] dB, 步长 0.5 dB (11 点) |
| 每 SNR 帧数 | 10,000 (低 SNR); 100,000 (高 SNR) |
| MC 收敛条件 | ≥60 帧错误或达到最大帧数 |
| 停止条件 | FER < 10⁻⁵ |

### 3.3 交付物

1. **BER-SNR 图**：4 条曲线 × 3 参数方案 = 12 张图（或 1 张图叠 4 条）× 3 = 3 张组合图
2. **时延柱状图**：3 参数方案 × 2 基线 × 3 SNR 点，柱状 + 加速比标注
3. **表格**：每个参数方案一张表，含 F_{2^m} ops、clock cycles (P=1/4/16/64)
4. **消融分析**：方案 A vs 方案 B 的 op 节省量化
5. **最终报告**：LaTeX PDF，含所有图表和讨论

### 3.4 实验矩阵

| # | 参数方案 | 算法 | SNR 范围 | 目标 |
|---|---|---|---|---|
| E1 | n=255-A | 纯 RS 硬判 (BM) | 3-8 dB | 硬判基线 |
| E2 | n=255-A | 纯 RS 软判 (LCC-BR) | 3-8 dB | 软判基线 |
| E3 | n=255-A | 方案 A (LLOSD + LCC-BR) | 3-8 dB | 最小闭环 |
| E4 | n=255-A | 方案 B (Lagrange 共享) | 3-8 dB | 消融对比 |
| E5-E8 | n=255-B | 同上 | 3-8 dB | 横向扩展 |
| E9-E12 | n=128 | 同上 | 3-8 dB | 横向扩展 |

## 四、里程碑

| 阶段 | 时长 | 交付 |
|---|---|---|
| Phase 1 | 2-3 天 | 项目结构 + PAM4 modem + BCH port + RS-BM |
| Phase 2 | 3-5 天 | LCC-BR 完整实现 + 级联链路 + 方案 A 跑通 |
| Phase 3 | 3-5 天 | Lagrange 共享层 + 方案 B + BER-SNR 曲线 |
| Phase 4 | 2-3 天 | 时延测量 + 三个参数方案完整仿真 |
| Phase 5 | 2 天 | LaTeX 报告 + 最终交付 |

**总计约 12-18 天。**

## 五、风险与应对

| 风险 | 应对 |
|---|---|
| LCC-BR 完整实现复杂 | 从 Xing-Chen-Bossert 2020 原文精读, 先做简化版 (partial root-finding 用 exhaustive search) |
| Lagrange 共享收益不够 | 详细做 op profiling, 定位真正瓶颈；退而分析共享的理论上限 |
| Case (a) 达不到 10% | 用 clock-cycle-with-parallelism 视角作补充证明 |
| n=128 extended BCH | 若太复杂, 先做 shortened BCH from n=127 |

## 六、参考

- **主论文**: Yang et al., "Efficient OSD of BCH Codes Without GE," *IEEE TIT* 2025.
- **RS 软判决**: Xing, Chen, Bossert, "Low-complexity Chase decoding of Reed-Solomon codes using module," *IEEE Trans. Commun.*, 2020.
- **BM 算法**: Berlekamp, *Algebraic Coding Theory*, 1968.
- **RS 代数**: MacWilliams & Sloane, *The Theory of Error-Correcting Codes*, 1977.
