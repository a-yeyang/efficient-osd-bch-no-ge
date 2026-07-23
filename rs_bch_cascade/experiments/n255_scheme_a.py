"""Main experiment: n=255 Scheme A minimum closed loop.

This is the "smallest closed loop" recommended by the advisor's answers:
  - n=255 config: RS(255, 239) + BCH(255, 239)
  - Baselines: pure RS-BM, pure RS-LCC-BR
  - Method under test: Scheme A (LLOSD inner + LCC-BR outer, no cache sharing)

Later experiments (fig_scheme_ab_compare.py) will add Scheme B.
"""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

from cascade_src.cascade import CascadeConfig
from cascade_src.simulate import run_bench, get_codec_pack, result_to_dict


def main():
    # n=255 Scheme A config
    cfg = CascadeConfig(m=8, k_rs=239, t_bch=2, llosd_tau=2, lcc_eta=4)
    print(f"[Main experiment: n=255 Scheme A] {cfg.describe()}")

    ebn0_list = np.arange(6.0, 12.5, 0.5)
    m_bits = cfg.m

    all_results = {}

    for method_key, label in [
        ("pure_rs_bm", "Pure RS-BM"),
        ("pure_rs_lccbr", "Pure RS-LCC-BR"),
        ("scheme_a", "Cascade A (LLOSD+LCC-BR)"),
    ]:
        print(f"\n--- {label} ---")
        pack = get_codec_pack(cfg, method_key)
        # For n=255, BM/LLOSD is much slower — use fewer frames
        max_frames = 400 if method_key.startswith("pure") else 200
        res = run_bench(
            method_name=label,
            codec_pack=pack,
            ebn0_list=ebn0_list,
            n_info_symbols=cfg.k_rs,
            m_bits_per_symbol=m_bits,
            min_frame_errors=15,
            max_frames=max_frames,
        )
        all_results[method_key] = result_to_dict(res)

    # Save
    with open("data/n255_scheme_a.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    # Plot BER-SNR
    plt.figure(figsize=(7, 5))
    styles = {
        "pure_rs_bm": ("D-", "red", "Pure RS-BM"),
        "pure_rs_lccbr": ("s-", "orange", "Pure RS-LCC-BR"),
        "scheme_a": ("^-", "blue", "Cascade A"),
    }
    for k, (mk, c, label) in styles.items():
        r = all_results[k]
        ebn = r["ebn0_db"]
        fer = [max(v, 1e-5) for v in r["fer"]]
        plt.semilogy(ebn[:len(fer)], fer, mk, color=c, label=label, markersize=6)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.title(f"n=255 minimum closed loop: {cfg.describe()}")
    plt.tight_layout()
    plt.savefig("figures/n255_scheme_a_fer.png", dpi=140)
    plt.savefig("figures/n255_scheme_a_fer.pdf")

    # BER plot
    plt.figure(figsize=(7, 5))
    for k, (mk, c, label) in styles.items():
        r = all_results[k]
        ebn = r["ebn0_db"]
        ber = [max(v, 1e-6) for v in r["ber"]]
        plt.semilogy(ebn[:len(ber)], ber, mk, color=c, label=label, markersize=6)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("BER (info bits)")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.title(f"n=255 minimum closed loop: {cfg.describe()}")
    plt.tight_layout()
    plt.savefig("figures/n255_scheme_a_ber.png", dpi=140)
    plt.savefig("figures/n255_scheme_a_ber.pdf")

    # Op comparison table
    print("\n\n=== Per-frame F_2^m ops (avg) ===")
    print(f"{'Eb/N0 (dB)':<12} {'RS-BM':>12} {'RS-LCC-BR':>12} {'Cascade A':>12}")
    for i, ebn in enumerate(all_results["pure_rs_bm"]["ebn0_db"]):
        if i >= len(all_results["scheme_a"]["ebn0_db"]):
            break
        print(f"{ebn:<12.2f} {all_results['pure_rs_bm']['avg_f2m_ops'][i]:>12.0f} "
              f"{all_results['pure_rs_lccbr']['avg_f2m_ops'][i]:>12.0f} "
              f"{all_results['scheme_a']['avg_f2m_ops'][i]:>12.0f}")

    print("\nSaved figures/n255_scheme_a_*.{png,pdf}")


if __name__ == "__main__":
    main()
