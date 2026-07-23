"""Main experiment: n=255 and n=127 hard-decision cascade."""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

from hc_src.cascade import HardCascadeConfig, HardCascadedCodec, PureRSHardCodec
from hc_src.simulate import run_bench, result_to_dict
from hc_src.latency_model import LatencyModel


CONFIGS = {
    "n255": HardCascadeConfig(m=8, k_rs=239),  # RS(255,239) t=8 + BCH(255,239) t=2
    "n127": HardCascadeConfig(m=7, k_rs=113),  # RS(127,113) t=7 + BCH(127,113) t=2 (R=0.792)
}


def build_codecs(cfg: HardCascadeConfig) -> dict:
    return {
        "pure_rs":       PureRSHardCodec(cfg),
        "cascade_conv":  HardCascadedCodec(cfg, bch_decoder="conv"),
        "cascade_direct": HardCascadedCodec(cfg, bch_decoder="direct"),
    }


def run_experiment(cfg_name: str, cfg: HardCascadeConfig, ebn0_list: np.ndarray) -> dict:
    print(f"\n\n########## {cfg_name}: {cfg.describe()} ##########")
    codecs = build_codecs(cfg)
    results = {}
    for name, codec in codecs.items():
        print(f"\n--- {cfg_name} / {name} ---")
        max_frames = 800 if "pure" in name else 400
        res = run_bench(
            method_name=f"{cfg_name}/{name}",
            codec=codec,
            ebn0_list=ebn0_list,
            n_info_symbols=cfg.k_rs,
            m_bits_per_symbol=cfg.m,
            min_frame_errors=20,
            max_frames=max_frames,
        )
        results[name] = result_to_dict(res)
    return results


def latency_summary(cfg: HardCascadeConfig) -> dict:
    """Compute all latency numbers for a given config."""
    LM = LatencyModel
    rs_only = LM.rs_bm_cycles(cfg.t_rs)
    conv_cyc = LM.bch_conv_cycles(2)
    dir_cyc = LM.bch_direct_cycles(2)
    return {
        "pure_rs":            rs_only,
        "conv_no_share":      LM.cascade_serial(conv_cyc, cfg.t_rs),
        "direct_no_share":    LM.cascade_serial(dir_cyc, cfg.t_rs),
        "conv_lagrange":      LM.cascade_lagrange_shared(conv_cyc, cfg.t_rs),
        "direct_lagrange":    LM.cascade_lagrange_shared(dir_cyc, cfg.t_rs),
        "kpi_target":         int(rs_only * 1.10),  # 110% of baseline
    }


def make_fer_plot(all_results: dict, out_path_base: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    styles = {
        "pure_rs":         ("D-", "red",   "Pure RS-BM"),
        "cascade_conv":    ("o-", "blue",  "Cascade (BCH-Conv + RS-BM)"),
        "cascade_direct":  ("^-", "green", "Cascade (BCH-Direct + RS-BM)"),
    }
    for ax, (cfg_name, results) in zip(axes, all_results.items()):
        for k, (mk, c, label) in styles.items():
            r = results[k]
            ebn = r["ebn0_db"][:len(r["fer"])]
            fer = [max(v, 1e-6) for v in r["fer"]]
            ax.semilogy(ebn, fer, mk, color=c, label=label, markersize=6)
        cfg = CONFIGS[cfg_name]
        ax.set_xlabel("Eb/N0 (dB)")
        ax.set_ylabel("FER")
        ax.set_title(f"{cfg_name}: {cfg.describe()}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path_base + ".png", dpi=140)
    plt.savefig(out_path_base + ".pdf")
    print(f"Saved {out_path_base}.png/.pdf")


def make_latency_bar(all_lat: dict, out_path_base: str):
    fig, ax = plt.subplots(figsize=(10, 5))
    cfgs = list(all_lat.keys())
    x = np.arange(len(cfgs))
    w = 0.16
    metrics = [
        ("pure_rs",         "Pure RS-BM (baseline)",           "red"),
        ("conv_no_share",   "Cascade Conv (no share)",          "blue"),
        ("direct_no_share", "Cascade Direct (no share)",        "green"),
        ("conv_lagrange",   "Cascade Conv (Lagrange shared)",   "darkblue"),
        ("direct_lagrange", "Cascade Direct (Lagrange shared)", "darkgreen"),
    ]
    for i, (key, label, color) in enumerate(metrics):
        values = [all_lat[cfg][key] for cfg in cfgs]
        ax.bar(x + (i - 2) * w, values, w, label=label, color=color, alpha=0.85)
    for i, cfg in enumerate(cfgs):
        kpi = all_lat[cfg]["kpi_target"]
        ax.hlines(kpi, x[i] - 2.5 * w, x[i] + 2.5 * w,
                  colors="black", linestyles="dashed", linewidth=1)
        ax.annotate(f"+10% KPI={kpi}", xy=(x[i], kpi), xytext=(x[i] + 0.05, kpi + 1),
                    fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(cfgs)
    ax.set_ylabel("Clock cycles")
    ax.set_title("Cascade latency: clock cycles vs Pure RS-BM baseline")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(out_path_base + ".png", dpi=140)
    plt.savefig(out_path_base + ".pdf")
    print(f"Saved {out_path_base}.png/.pdf")


def print_kpi_table(all_lat: dict):
    print("\n\n=== Latency KPI Table ===")
    print(f"{'config':<10}{'pure_rs':>10}{'Conv':>12}{'Direct':>12}{'Conv+share':>14}{'Direct+share':>14}{'KPI':>10}")
    for cfg_name, lat in all_lat.items():
        rs = lat["pure_rs"]
        kpi = lat["kpi_target"]
        conv = lat["conv_no_share"]
        direct = lat["direct_no_share"]
        conv_s = lat["conv_lagrange"]
        direct_s = lat["direct_lagrange"]
        ratio_direct_s = (direct_s - rs) / rs * 100
        print(f"{cfg_name:<10}{rs:>10}{conv:>12}{direct:>12}{conv_s:>14}{direct_s:>14}{kpi:>10}")
        print(f"{'  →ratio':<10}{'':>10}"
              f"{(conv-rs)/rs*100:>+11.1f}%"
              f"{(direct-rs)/rs*100:>+11.1f}%"
              f"{(conv_s-rs)/rs*100:>+13.1f}%"
              f"{(direct_s-rs)/rs*100:>+13.1f}%"
              f"{' (10%)':>10}")


def main():
    all_results = {}
    all_lat = {}
    for cfg_name, cfg in CONFIGS.items():
        if cfg.m == 8:
            ebn0_list = np.arange(6.0, 12.5, 0.5)
        else:
            ebn0_list = np.arange(5.0, 12.0, 0.5)
        results = run_experiment(cfg_name, cfg, ebn0_list)
        all_results[cfg_name] = results
        all_lat[cfg_name] = latency_summary(cfg)

        with open(f"data/{cfg_name}_results.json", "w") as f:
            json.dump({"results": results, "latency": all_lat[cfg_name],
                       "config": {"m": cfg.m, "k_rs": cfg.k_rs, "n_rs": cfg.n_rs,
                                  "t_rs": cfg.t_rs}}, f, indent=2, default=str)

    with open("data/all_results.json", "w") as f:
        json.dump({"results": all_results, "latency": all_lat}, f, indent=2, default=str)

    make_fer_plot(all_results, "figures/fer_all")
    make_latency_bar(all_lat, "figures/latency_bars")
    print_kpi_table(all_lat)


if __name__ == "__main__":
    main()
