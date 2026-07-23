"""Comprehensive experiment: all 3 configurations x all 4 methods.

Configurations:
    n=255-A: RS(255,239) + BCH(255,239)   [t_bch=2, R=0.879]
    n=255-B: RS(255,235) + BCH(255,247)   [t_bch=1, R=0.891]
    n=127:   RS(127,119) + BCH(127,120)   [t_bch=1, R=0.885]

Methods:
    pure_rs_bm (hard baseline)
    pure_rs_lccbr (soft baseline)
    scheme_a (LLOSD inner + LCC-BR outer, no cache sharing)
    scheme_b (LLOSD + LCC-BR + Lagrange shared cache)
"""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

from cascade_src.cascade import CascadeConfig
from cascade_src.simulate import run_bench, get_codec_pack, result_to_dict


CONFIGS = {
    "n255_A": CascadeConfig(m=8, k_rs=239, t_bch=2, llosd_tau=2, lcc_eta=4),
    "n255_B": CascadeConfig(m=8, k_rs=235, t_bch=1, llosd_tau=1, lcc_eta=4),
    "n127":   CascadeConfig(m=7, k_rs=119, t_bch=1, llosd_tau=1, lcc_eta=4),
}


def run_all():
    all_results = {}
    for cfg_name, cfg in CONFIGS.items():
        print(f"\n########## {cfg_name}: {cfg.describe()} ##########\n")
        results_by_method = {}
        # Choose SNR range based on config size
        if cfg.m == 8:
            ebn0_list = np.arange(6.0, 12.5, 0.5)
        else:
            ebn0_list = np.arange(5.0, 11.5, 0.5)

        for method_key in ["pure_rs_bm", "pure_rs_lccbr", "scheme_a", "scheme_b"]:
            print(f"\n--- {cfg_name} / {method_key} ---")
            pack = get_codec_pack(cfg, method_key)
            max_frames = 200 if method_key.startswith("pure") else 150
            res = run_bench(
                method_name=f"{cfg_name}/{method_key}",
                codec_pack=pack,
                ebn0_list=ebn0_list,
                n_info_symbols=cfg.k_rs,
                m_bits_per_symbol=cfg.m,
                min_frame_errors=15,
                max_frames=max_frames,
            )
            results_by_method[method_key] = result_to_dict(res)
        all_results[cfg_name] = results_by_method

        # Save intermediate
        with open(f"data/{cfg_name}_all_methods.json", "w") as f:
            json.dump(results_by_method, f, indent=2, default=str)

    # Combined plot: one panel per config
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    styles = {
        "pure_rs_bm": ("D-", "red", "Pure RS-BM"),
        "pure_rs_lccbr": ("s-", "orange", "Pure RS-LCC-BR"),
        "scheme_a": ("^-", "blue", "Cascade A"),
        "scheme_b": ("v--", "green", "Cascade B"),
    }
    for ax, (cfg_name, cfg) in zip(axes, CONFIGS.items()):
        rr = all_results[cfg_name]
        for mk_key, (mk, c, label) in styles.items():
            if mk_key not in rr: continue
            r = rr[mk_key]
            ebn = r["ebn0_db"][:len(r["fer"])]
            fer = [max(v, 1e-5) for v in r["fer"]]
            ax.semilogy(ebn, fer, mk, color=c, label=label, markersize=5)
        ax.set_xlabel("Eb/N0 (dB)")
        ax.set_ylabel("FER")
        ax.set_title(f"{cfg_name}: {cfg.describe()}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("figures/all_configs_fer.png", dpi=140)
    plt.savefig("figures/all_configs_fer.pdf")
    print("\nSaved figures/all_configs_fer.{png,pdf}")

    with open("data/all_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)


if __name__ == "__main__":
    run_all()
