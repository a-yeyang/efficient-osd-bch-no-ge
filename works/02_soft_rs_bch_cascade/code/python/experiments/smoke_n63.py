"""Smoke test on n=63 (small) config to verify end-to-end simulation pipeline."""
from work_paths import setup
setup()

import json
import numpy as np
import matplotlib.pyplot as plt

from cascade_src.cascade import CascadeConfig
from cascade_src.simulate import run_bench, get_codec_pack, result_to_dict


def main():
    # Small n=63 config for smoke test
    cfg = CascadeConfig(m=6, k_rs=57, t_bch=1, llosd_tau=2, lcc_eta=3)
    print(f"[SMOKE TEST] {cfg.describe()}")

    ebn0_list = np.arange(4.0, 8.5, 0.5)
    m_bits = cfg.m

    all_results = {}

    for method_key, label in [
        ("pure_rs_bm", "Pure RS-BM (hard)"),
        ("pure_rs_lccbr", "Pure RS-LCC-BR (soft)"),
        ("scheme_a", "Cascade A (LLOSD+LCC-BR)"),
        ("scheme_b", "Cascade B (Lagrange shared)"),
    ]:
        print(f"\n--- {label} ---")
        pack = get_codec_pack(cfg, method_key)
        res = run_bench(
            method_name=label,
            codec_pack=pack,
            ebn0_list=ebn0_list,
            n_info_symbols=cfg.k_rs,
            m_bits_per_symbol=m_bits,
            min_frame_errors=30,
            max_frames=1500,
        )
        all_results[method_key] = result_to_dict(res)

    with open("data/smoke_n63.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    # Plot
    plt.figure(figsize=(7, 5))
    styles = {
        "pure_rs_bm": ("D-", "red", "Pure RS-BM"),
        "pure_rs_lccbr": ("s-", "orange", "Pure RS-LCC-BR"),
        "scheme_a": ("^-", "blue", "Cascade A"),
        "scheme_b": ("v--", "green", "Cascade B"),
    }
    for k, (mk, c, label) in styles.items():
        r = all_results[k]
        ebn = r["ebn0_db"][:len(r["fer"])]
        fer = [max(v, 1e-6) for v in r["fer"]]
        plt.semilogy(ebn, fer, mk, color=c, label=label, markersize=6)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.title(f"Smoke test: {cfg.describe()}")
    plt.tight_layout()
    plt.savefig("figures/smoke_n63.png", dpi=140)
    plt.savefig("figures/smoke_n63.pdf")
    print(f"\nSaved figures/smoke_n63.{{png,pdf}}")


if __name__ == "__main__":
    main()
