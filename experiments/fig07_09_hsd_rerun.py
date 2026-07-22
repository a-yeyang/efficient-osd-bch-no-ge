"""Re-run only the HSD portion of Fig 7 & 9 with the fixed HSD (no LLOSD early
terminate).
"""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.decoders import hsd_fast


def run_fer(code, decoder_fn, ebn0_list, min_errors, max_frames, seed=0, label=""):
    zero = np.zeros(code.n, dtype=np.int8)
    x0 = bpsk_modulate(zero)
    rate = code.k / code.n
    fers = []
    import time
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
        sigma = sigma_from_ebn0(ebn0, rate)
        n_frames = 0
        n_errors = 0
        t0 = time.perf_counter()
        while n_frames < max_frames:
            y = x0 + sigma * rng.standard_normal(code.n)
            L = 2.0 * y / (sigma * sigma)
            c_hat, _ = decoder_fn(code, L)
            n_frames += 1
            if not np.array_equal(c_hat, zero):
                n_errors += 1
            if n_errors >= min_errors and n_frames >= 200:
                break
        fer = n_errors / n_frames if n_frames > 0 else 0.0
        dt = time.perf_counter() - t0
        print(f"  {label} @ {ebn0:.1f} dB: FER = {fer:.2e}  ({n_errors}/{n_frames})  {dt:.1f}s", flush=True)
        fers.append(fer)
    return fers


def fig7_hsd():
    """Re-run HSD variants for Fig 7."""
    code = BCHCode(m=7, t=4)
    print(f"=== Fig 7 HSD re-run: ({code.n}, {code.k}) ===")
    ebn0_list = np.arange(3.0, 7.5, 1.0)

    with open("data/fig07_127_99.json") as f:
        old = json.load(f)

    results = {k: v for k, v in old.items() if k != "ebn0_db"}
    results["HSD(1,4)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=4),
                                    ebn0_list, 40, 3000, label="HSD(1,4)")
    results["HSD(1,6)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=6),
                                    ebn0_list, 40, 2000, label="HSD(1,6)")
    results["HSD(1,8)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=8),
                                    ebn0_list, 40, 1500, label="HSD(1,8)")

    plt.figure(figsize=(6, 4.5))
    styles = {
        "OSD(1)":   ("o-", "red"),
        "OSD(2)":   ("s-", "red"),
        "LLOSD(2)": ("D-", "royalblue"),
        "LLOSD(3)": ("D--", "royalblue"),
        "HSD(1,4)": ("^-", "black"),
        "HSD(1,6)": ("v-", "black"),
        "HSD(1,8)": ("x-", "black"),
    }
    for name in ["OSD(1)", "OSD(2)", "LLOSD(2)", "LLOSD(3)", "HSD(1,4)", "HSD(1,6)", "HSD(1,8)"]:
        if name not in results:
            continue
        fers = results[name]
        st, c = styles[name]
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-6), st, color=c, label=name, markersize=4)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8)
    plt.title("Fig. 7: (127, 99) BCH (fixed HSD)")
    plt.tight_layout()
    plt.savefig("figures/fig07_127_99.png", dpi=140)
    plt.savefig("figures/fig07_127_99.pdf")
    with open("data/fig07_127_99.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **results}, f, indent=2)


def fig9_hsd():
    code = BCHCode(m=8, t=4)
    print(f"\n=== Fig 9 HSD re-run: ({code.n}, {code.k}) ===")
    ebn0_list = np.arange(3.5, 7.5, 0.5)

    with open("data/fig09_255_223.json") as f:
        old = json.load(f)

    results = {k: v for k, v in old.items() if k != "ebn0_db"}
    results["HSD(1,4)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=4),
                                    ebn0_list, 30, 500, label="HSD(1,4)")
    results["HSD(1,6)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=6),
                                    ebn0_list, 30, 400, label="HSD(1,6)")
    results["HSD(1,8)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=8),
                                    ebn0_list, 30, 300, label="HSD(1,8)")

    plt.figure(figsize=(6, 4.5))
    styles = {
        "BM":       ("D-", "red"),
        "OSD(1)":   ("o-", "red"),
        "LLOSD(2)": ("x-", "royalblue"),
        "PLCC(6)":  ("s-", "black"),
        "PLCC(8)":  ("s--", "black"),
        "HSD(1,4)": ("^-", "black"),
        "HSD(1,6)": ("v-", "black"),
        "HSD(1,8)": ("<-", "black"),
    }
    for name in ["BM", "OSD(1)", "LLOSD(2)", "PLCC(6)", "PLCC(8)",
                 "HSD(1,4)", "HSD(1,6)", "HSD(1,8)"]:
        if name not in results:
            continue
        fers = results[name]
        st, c = styles[name]
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-5), st, color=c, label=name, markersize=4)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8, loc="lower left")
    plt.title("Fig. 9: (255, 223) BCH (fixed HSD)")
    plt.tight_layout()
    plt.savefig("figures/fig09_255_223.png", dpi=140)
    plt.savefig("figures/fig09_255_223.pdf")
    with open("data/fig09_255_223.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **results}, f, indent=2)


if __name__ == "__main__":
    fig7_hsd()
    fig9_hsd()
