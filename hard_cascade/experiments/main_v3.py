"""v3 experiment (advisor's third task, 2026-07-27).

Two cascade configs, decoupled RS/BCH fields, BER (primary) vs SNR + a
three-innovation latency & complexity decomposition (Cascade / Direct / Lagrange).

  Config 1 (NEW):  RS(544,514, t=15)/GF(2^10) + BCH(144,136, t=1)/GF(2^8)   [KP4]
  Config 2 (KEEP): RS(255,239, t=8)/GF(2^8)   + BCH(255,239, t=2)/GF(2^8)

Outputs:
  figures/ber_snr_v3.{png,pdf}         BER-SNR (primary metric)
  figures/latency_innovations_v3.{png,pdf}
  figures/complexity_v3.{png,pdf}
  data/v3_results.json
"""
from __future__ import annotations

import sys
sys.path.insert(0, ".")

import json
import time
import numpy as np
import matplotlib.pyplot as plt

from hc_src.cascade_v3 import (
    CONFIGS_V3, CascadeV3Codec, PureRSV3Codec, run_pam4_channel_hard,
)
from hc_src.latency_model import LatencyModel
from hc_src.upstream import OpCounters


# SNR sweeps chosen from calibration (waterfall regions)
EBN0_SWEEP = {
    "cfg1_kp4": [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5],
    "cfg2_255": [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5],
}


def _bits_of_symbols(syms: np.ndarray, m: int) -> np.ndarray:
    out = np.zeros(syms.size * m, dtype=np.int8)
    for i in range(syms.size):
        si = int(syms[i])
        for b in range(m):
            out[i * m + b] = (si >> b) & 1
    return out


def run_ber(method, codec, cfg, ebn0_list, seed=12345,
            min_frame_errors=40, min_frames=80, max_frames=1500, verbose=True):
    """Monte Carlo BER (primary) + FER + avg F_2m ops per SNR point."""
    out = {"ebn0_db": [], "ber": [], "fer": [], "avg_f2m": [],
           "n_frames": [], "n_frame_err": []}
    m = cfg.m_rs
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
        nf = nfe = nbe = 0
        tot_bits = 0
        sum_f2m = 0.0
        t0 = time.perf_counter()
        while nf < max_frames:
            msg = rng.integers(0, 1 << m, cfg.k_rs)
            coded = codec.encode(msg)
            hard = run_pam4_channel_hard(coded, ebn0, codec.effective_rate, rng)
            ctr = OpCounters()
            mh, _ok, _ = codec.decode(hard, ctr)
            nf += 1
            mb = _bits_of_symbols(msg, m)
            mhb = _bits_of_symbols(mh, m)
            be = int(np.sum(mb != mhb))
            nbe += be
            tot_bits += mb.size
            if be > 0:
                nfe += 1
            sum_f2m += ctr.f2m
            if nfe >= min_frame_errors and nf >= min_frames:
                break
        ber = nbe / max(1, tot_bits)
        fer = nfe / max(1, nf)
        avg_f2m = sum_f2m / max(1, nf)
        dt = time.perf_counter() - t0
        if verbose:
            print(f"  {method:22s} @ {ebn0:4.1f} dB: BER={ber:.3e} "
                  f"FER={fer:.3e} f2m={avg_f2m:.0f} ({nf}f {dt:.1f}s)")
        out["ebn0_db"].append(ebn0)
        out["ber"].append(ber)
        out["fer"].append(fer)
        out["avg_f2m"].append(avg_f2m)
        out["n_frames"].append(nf)
        out["n_frame_err"].append(nfe)
        if ber == 0.0 and nfe == 0 and ebn0 >= ebn0_list[len(ebn0_list) // 2]:
            break
    return out


# =====================================================================
# Latency & complexity decomposition (three innovations)
# =====================================================================
def latency_innovation_breakdown(cfg):
    """Return the cycle count as we switch on each innovation in turn.

    Baseline reference = Pure RS-BM (no inner code).
    Stages:
      S0 Pure RS-BM
      S1 + Cascade only (Conv inner, no share)     -> Cascade innovation cost
      S2 + Direct (Conv->Direct, still no share)    -> Direct innovation saving
      S3 + Lagrange v1 (share syndrome)             -> Lagrange v1 saving
      S4 + Lagrange v2 (share GF mult array)        -> Lagrange v2 saving
    """
    LM = LatencyModel
    t_rs = cfg.t_rs
    t_bch = cfg.t_bch
    m_bch = cfg.m_bch
    conv = LM.bch_conv_cycles(t_bch, m=m_bch)
    direct = LM.bch_direct_cycles(t_bch, m=m_bch)
    pure = LM.rs_bm_cycles(t_rs)
    return {
        "pure_rs": pure,
        "kpi_ceiling": pure * 1.10,
        "bch_conv_cyc": conv,
        "bch_direct_cyc": direct,
        "cascade_conv_none":   LM.cascade_serial(conv, t_rs),
        "cascade_direct_none": LM.cascade_serial(direct, t_rs),
        "cascade_direct_v1":   LM.cascade_lagrange_v1(direct, t_rs),
        "cascade_direct_v2":   LM.cascade_lagrange_v2(direct, t_rs),
        # also conv+v2 for the full table
        "cascade_conv_v1":     LM.cascade_lagrange_v1(conv, t_rs),
        "cascade_conv_v2":     LM.cascade_lagrange_v2(conv, t_rs),
    }


def measure_complexity(cfg, ebn0, n_meas=120, seed=999):
    """Empirical avg F_2m ops per decoded frame for Pure RS vs Cascade
    (Conv vs Direct), at a representative SNR."""
    casc_conv = CascadeV3Codec(cfg, "conv")
    casc_dir = CascadeV3Codec(cfg, "direct")
    pure = PureRSV3Codec(cfg)
    m = cfg.m_rs
    res = {}
    for name, cod in [("pure", pure), ("casc_conv", casc_conv),
                      ("casc_direct", casc_dir)]:
        rng = np.random.default_rng(seed + int(ebn0 * 100))
        s = 0.0
        for _ in range(n_meas):
            msg = rng.integers(0, 1 << m, cfg.k_rs)
            coded = cod.encode(msg)
            hard = run_pam4_channel_hard(coded, ebn0, cod.effective_rate, rng)
            ctr = OpCounters()
            cod.decode(hard, ctr)
            s += ctr.f2m
        res[name] = s / n_meas
    return res


# =====================================================================
# Plots
# =====================================================================
def plot_ber(all_ber, out_base):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    titles = {
        "cfg1_kp4": "Config 1: RS(544,514,t=15)/GF($2^{10}$) + BCH(144,136,t=1)",
        "cfg2_255": "Config 2: RS(255,239,t=8)/GF($2^{8}$) + BCH(255,239,t=2)",
    }
    for ax, key in zip(axes, ["cfg1_kp4", "cfg2_255"]):
        d = all_ber[key]
        ax.semilogy(d["pure"]["ebn0_db"], np.maximum(d["pure"]["ber"], 1e-7),
                    "s--", color="#C0392B", label="Pure RS-BM (baseline)",
                    markersize=6, linewidth=1.6)
        ax.semilogy(d["cascade"]["ebn0_db"], np.maximum(d["cascade"]["ber"], 1e-7),
                    "o-", color="#1E7B34", label="Cascade (Direct + RS-BM)",
                    markersize=6, linewidth=1.8)
        ax.set_xlabel("$E_b/N_0$ (dB)")
        ax.set_ylabel("BER")
        ax.set_title(titles[key], fontsize=9.5)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=9)
        ax.set_ylim(1e-7, 1)
    plt.suptitle("BER-SNR: hard-decision RS+BCH cascade vs pure RS (PAM4/AWGN)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_base + ".png", dpi=140)
    plt.savefig(out_base + ".pdf")
    print(f"Saved {out_base}.png/.pdf")


def plot_latency_innovations(brk, out_base):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    names = {"cfg1_kp4": "Config 1 (KP4)", "cfg2_255": "Config 2 (255)"}
    for ax, key in zip(axes, ["cfg1_kp4", "cfg2_255"]):
        s = brk[key]
        labels = ["Pure\nRS-BM", "+Cascade\n(Conv)", "+Direct\n(no share)",
                  "+Lagrange\nv1", "+Lagrange\nv2"]
        vals = [s["pure_rs"], s["cascade_conv_none"], s["cascade_direct_none"],
                s["cascade_direct_v1"], s["cascade_direct_v2"]]
        colors = ["#C0392B", "#2E6BB0", "#7FB77E", "#3C9A5F", "#1E7B34"]
        ax.bar(range(len(vals)), vals, color=colors, edgecolor="black",
               linewidth=0.6, alpha=0.9)
        kpi = s["kpi_ceiling"]
        ax.axhline(kpi, color="black", ls="--", lw=1,
                   label=f"+10% KPI = {kpi:.1f} cyc")
        for i, v in enumerate(vals):
            r = (v - s["pure_rs"]) / s["pure_rs"] * 100
            txt = "baseline" if i == 0 else f"{r:+.1f}%"
            bold = (i > 0 and r <= 10)
            ax.text(i, v + max(vals) * 0.015, f"{v}\n({txt})", ha="center",
                    fontsize=8, fontweight="bold" if bold else "normal")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Clock cycles")
        ax.set_title(f"{names[key]}: 3-innovation latency\n"
                     f"(RS t={s['bch_direct_cyc']*0+CONFIGS_V3[key].t_rs}, "
                     f"BCH-Direct={s['bch_direct_cyc']} cyc)", fontsize=9)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(True, alpha=0.3, axis="y")
        ax.set_ylim(0, max(vals) * 1.28)
    plt.suptitle("Latency contribution of each innovation: Cascade + Direct + Lagrange",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_base + ".png", dpi=140)
    plt.savefig(out_base + ".pdf")
    print(f"Saved {out_base}.png/.pdf")


def plot_complexity(cx, out_base):
    fig, ax = plt.subplots(figsize=(9, 5))
    keys = ["cfg1_kp4", "cfg2_255"]
    names = ["Config 1 (KP4)", "Config 2 (255)"]
    x = np.arange(len(keys))
    w = 0.26
    pure = [cx[k]["pure"] for k in keys]
    conv = [cx[k]["casc_conv"] for k in keys]
    direct = [cx[k]["casc_direct"] for k in keys]
    ax.bar(x - w, pure, w, label="Pure RS-BM", color="#C0392B", alpha=0.9,
           edgecolor="black", linewidth=0.5)
    ax.bar(x, conv, w, label="Cascade (Conv inner)", color="#2E6BB0", alpha=0.9,
           edgecolor="black", linewidth=0.5)
    ax.bar(x + w, direct, w, label="Cascade (Direct inner)", color="#1E7B34",
           alpha=0.9, edgecolor="black", linewidth=0.5)
    for i in range(len(keys)):
        for off, v in [(-w, pure[i]), (0, conv[i]), (w, direct[i])]:
            ax.text(i + off, v * 1.01, f"{v:.0f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Avg $\\mathbb{F}_{2^m}$ ops per frame")
    ax.set_title("Decoding complexity: pure RS vs cascade (Conv/Direct inner)",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(out_base + ".png", dpi=140)
    plt.savefig(out_base + ".pdf")
    print(f"Saved {out_base}.png/.pdf")


# =====================================================================
def main():
    all_ber = {}
    for key, cfg in CONFIGS_V3.items():
        print("=" * 78)
        print(cfg.describe())
        casc = CascadeV3Codec(cfg, "direct")
        pure = PureRSV3Codec(cfg)
        sweep = EBN0_SWEEP[key]
        r_casc = run_ber("Cascade(Direct)", casc, cfg, sweep)
        r_pure = run_ber("Pure RS-BM", pure, cfg, sweep)
        all_ber[key] = {"cascade": r_casc, "pure": r_pure}

    # Latency breakdown (analytic)
    brk = {key: latency_innovation_breakdown(cfg)
           for key, cfg in CONFIGS_V3.items()}

    # Complexity (empirical) at a mid-waterfall SNR
    cx_snr = {"cfg1_kp4": 8.0, "cfg2_255": 7.0}
    cx = {key: measure_complexity(cfg, cx_snr[key])
          for key, cfg in CONFIGS_V3.items()}

    print("\n" + "=" * 78)
    print("LATENCY (3 innovations) — clock cycles")
    for key, s in brk.items():
        base = s["pure_rs"]
        print(f"\n### {CONFIGS_V3[key].name}  (Pure RS-BM = {base} cyc, "
              f"KPI ceiling = {s['kpi_ceiling']:.1f})")
        for lbl, k in [("Cascade Conv (no share)", "cascade_conv_none"),
                       ("Cascade Direct (no share)", "cascade_direct_none"),
                       ("Cascade Direct + Lagrange v1", "cascade_direct_v1"),
                       ("Cascade Direct + Lagrange v2", "cascade_direct_v2")]:
            v = s[k]
            r = (v - base) / base * 100
            print(f"   {lbl:34s}: {v:3d} cyc  {r:+6.1f}%  "
                  f"{'PASS' if r <= 10 else 'FAIL'}")

    plot_ber(all_ber, "figures/ber_snr_v3")
    plot_latency_innovations(brk, "figures/latency_innovations_v3")
    plot_complexity(cx, "figures/complexity_v3")

    with open("data/v3_results.json", "w") as f:
        json.dump({"ber": all_ber, "latency": brk, "complexity": cx,
                   "cx_snr": cx_snr}, f, indent=2)
    print("\nSaved data/v3_results.json")


if __name__ == "__main__":
    main()
