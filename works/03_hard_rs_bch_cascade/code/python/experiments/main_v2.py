"""v2 experiment: verify KPI with aggressive Lagrange sharing (v2) +
n=127 Direct compressed to 2 cycles.

BER results reuse the v1 run (data/all_results.json). This script only
recomputes the latency numbers under the v2 optimizations and produces
new figures + KPI tables.
"""
from work_paths import setup
setup()

import json
import numpy as np
import matplotlib.pyplot as plt

from hc_src.cascade import HardCascadeConfig
from hc_src.latency_model import LatencyModel


CONFIGS = {
    "n255": HardCascadeConfig(m=8, k_rs=239),  # RS(255,239) t=8 + BCH(255,239) t=2
    "n127": HardCascadeConfig(m=7, k_rs=113),  # RS(127,113) t=7 + BCH(127,113) t=2
}


def summarize(cfg: HardCascadeConfig) -> dict:
    """All latency numbers for one config."""
    LM = LatencyModel
    m = cfg.m
    t = cfg.t_rs
    rs = LM.rs_bm_cycles(t)
    conv_cyc = LM.bch_conv_cycles(2, m=m)
    dir_cyc = LM.bch_direct_cycles(2, m=m)
    return {
        "pure_rs": rs,
        "kpi_ceiling": rs * 1.10,
        "bch_conv_cyc": conv_cyc,
        "bch_direct_cyc": dir_cyc,
        # Conv variants
        "conv_none":  LM.cascade_serial(conv_cyc, t),
        "conv_v1":    LM.cascade_lagrange_v1(conv_cyc, t),
        "conv_v2":    LM.cascade_lagrange_v2(conv_cyc, t),
        # Direct variants
        "direct_none": LM.cascade_serial(dir_cyc, t),
        "direct_v1":   LM.cascade_lagrange_v1(dir_cyc, t),
        "direct_v2":   LM.cascade_lagrange_v2(dir_cyc, t),
    }


def print_kpi_table(all_lat: dict):
    print("=" * 90)
    print("KPI Table (Clock Cycles) — v1 vs v2 optimization")
    print("=" * 90)
    for cfg_name, s in all_lat.items():
        rs = s["pure_rs"]
        kpi = s["kpi_ceiling"]
        print(f"\n### {cfg_name} (RS-BM baseline = {rs} cyc, 10% KPI ceiling = {kpi:.1f} cyc)")
        print(f"BCH-Direct cycles: {s['bch_direct_cyc']}  (BCH-Conv: {s['bch_conv_cyc']})")
        print(f"{'Config':<40} {'Cycles':>8} {'vs baseline':>12} {'KPI':>8}")
        print("-" * 72)
        for label, key in [
            ("Cascade Conv (no share)",       "conv_none"),
            ("Cascade Conv + v1 share",       "conv_v1"),
            ("Cascade Conv + v2 share",       "conv_v2"),
            ("Cascade Direct (no share)",     "direct_none"),
            ("Cascade Direct + v1 share",     "direct_v1"),
            ("Cascade Direct + v2 share",     "direct_v2"),
        ]:
            cyc = s[key]
            ratio = (cyc - rs) / rs * 100
            passes = "PASS" if ratio <= 10.0 else "FAIL"
            print(f"{label:<40} {cyc:>8} {ratio:>+11.1f}% {passes:>8}")


def make_v2_latency_plot(all_lat: dict, out_path_base: str):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (cfg_name, s) in zip(axes, all_lat.items()):
        labels = [
            "Pure RS-BM\n(baseline)",
            "Conv\nno share",
            "Direct\nno share",
            "Direct\nv1 share",
            "Direct\nv2 share",
        ]
        values = [
            s["pure_rs"],
            s["conv_none"],
            s["direct_none"],
            s["direct_v1"],
            s["direct_v2"],
        ]
        # Color code: red = baseline, blue = Conv, greens graduated for Direct variants
        colors = ["red", "blue", "#90EE90", "#3CB371", "#228B22"]
        bars = ax.bar(range(len(labels)), values, color=colors, alpha=0.85,
                       edgecolor="black", linewidth=0.6)

        kpi = s["kpi_ceiling"]
        ax.axhline(kpi, color="black", linestyle="dashed", linewidth=1,
                    label=f"+10% KPI = {kpi:.1f} cyc")
        for i, v in enumerate(values):
            ratio = (v - s["pure_rs"]) / s["pure_rs"] * 100
            label_str = "baseline" if i == 0 else f"{ratio:+.1f}%"
            ax.text(i, v + 0.3, f"{v}\n({label_str})", ha="center", fontsize=8,
                    fontweight="bold" if ratio <= 10 and i > 0 else "normal")

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Clock cycles")
        ax.set_title(f"{cfg_name}: BCH t=2, BCH-Direct = {s['bch_direct_cyc']} cyc")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3, axis="y")
        ax.set_ylim(0, max(values) * 1.25)
    plt.suptitle("Cascade Latency v2: aggressive Lagrange share + smaller-GF Direct",
                  fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path_base + ".png", dpi=140)
    plt.savefig(out_path_base + ".pdf")
    print(f"\nSaved {out_path_base}.png/.pdf")


def make_evolution_plot(all_lat: dict, out_path: str):
    """Show KPI evolution: v0 (Conv no share) → v1 → v2, with 10% KPI line."""
    fig, ax = plt.subplots(figsize=(10, 5))
    stages = ["Conv\n(no share)", "Direct\n(no share)",
              "Direct\n+ v1 share", "Direct\n+ v2 share"]
    for cfg_name, s in all_lat.items():
        cycles = [s["conv_none"], s["direct_none"], s["direct_v1"], s["direct_v2"]]
        ratios = [(c - s["pure_rs"]) / s["pure_rs"] * 100 for c in cycles]
        ax.plot(stages, ratios, "o-", label=f"{cfg_name} (RS baseline={s['pure_rs']} cyc)",
                 linewidth=2, markersize=10)
        for i, (stage, r) in enumerate(zip(stages, ratios)):
            ax.text(i, r + 1.5, f"{r:+.1f}%", ha="center", fontsize=9)

    ax.axhline(10, color="red", linestyle="dashed", linewidth=1.5, label="10% KPI ceiling")
    ax.axhline(0, color="gray", linestyle="dotted", linewidth=1)
    ax.set_ylabel("Latency increase vs Pure RS-BM (%)")
    ax.set_title("Cascade latency evolution: optimizations from Conv to Direct v2 share",
                  fontsize=11)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_ylim(-5, 55)
    plt.tight_layout()
    plt.savefig(out_path + ".png", dpi=140)
    plt.savefig(out_path + ".pdf")
    print(f"Saved {out_path}.png/.pdf")


def main():
    all_lat = {name: summarize(cfg) for name, cfg in CONFIGS.items()}
    print_kpi_table(all_lat)

    make_v2_latency_plot(all_lat, "figures/latency_bars_v2")
    make_evolution_plot(all_lat, "figures/latency_evolution")

    with open("data/v2_latency.json", "w") as f:
        json.dump(all_lat, f, indent=2)
    print("\nSaved data/v2_latency.json")


if __name__ == "__main__":
    main()
