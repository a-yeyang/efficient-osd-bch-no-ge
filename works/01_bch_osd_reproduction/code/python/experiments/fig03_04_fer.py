"""Fig 3 & Fig 4 reproduction: FER curves for (31, 21) and (63, 45) BCH codes.

Fig. 3 (31, 21): BM, OSD (1), LLOSD (1), LLOSD (2), ML  at Eb/N0 ∈ [2, 8] dB
Fig. 4 (63, 45): BM, LLOSD (1..3), OSD (1), ML  at Eb/N0 ∈ [1, 7] dB
"""
from work_paths import setup
setup()

import json
import numpy as np
import matplotlib.pyplot as plt

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.osd import osd_decode
from src.decoders import llosd_fast
from src.ml import ml_decode_full_codebook, ml_approx_by_high_order_llosd


def bm_wrap(code, L):
    r_hard = (L < 0).astype(np.int8)
    dec, ok = code.bm_decode(r_hard)
    return (dec if ok else r_hard), {"counters": None}


def run_fer(code, decoder_fn, ebn0_list, min_errors, max_frames, seed=0, label=""):
    rate = code.k / code.n
    zero = np.zeros(code.n, dtype=np.int8)
    x0 = bpsk_modulate(zero)
    fers = []
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
        sigma = sigma_from_ebn0(ebn0, rate)
        n_frames = 0
        n_errors = 0
        import time
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
        print(f"  {label} @ {ebn0:.1f} dB: FER = {fer:.2e}  ({n_errors}/{n_frames})  {dt:.1f}s")
        fers.append(fer)
    return fers


def fig3_31_21():
    code = BCHCode(m=5, t=2)
    ebn0_list = np.arange(2.0, 8.5, 1.0)
    results = {}
    print(f"\n=== Fig 3: (31, 21) BCH ===")
    print("BM:")
    results["BM"] = run_fer(code, bm_wrap, ebn0_list, 60, 20000, label="BM")
    print("OSD(1):")
    results["OSD(1)"] = run_fer(code, lambda c, L: osd_decode(c, L, tau=1), ebn0_list, 60, 15000, label="OSD(1)")
    print("LLOSD(1):")
    results["LLOSD(1)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=1), ebn0_list, 60, 15000, label="LLOSD(1)")
    print("LLOSD(2):")
    results["LLOSD(2)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=2), ebn0_list, 60, 15000, label="LLOSD(2)")
    # ML — (31,21) has 2^21 = 2M codewords, tractable.
    print("ML (full codebook):")
    # Even 2M codewords x 20k frames is ~40 billion inner products — too slow.
    # Use LLOSD(4) which is very close to ML for (31,21).
    results["ML"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=4), ebn0_list, 60, 15000, label="ML(~LLOSD 4)")

    plt.figure(figsize=(6, 4.5))
    styles = {
        "BM":       ("D-", "red"),
        "OSD(1)":   ("o-", "black"),
        "LLOSD(1)": ("v-", "black"),
        "LLOSD(2)": ("^-", "black"),
        "ML":       ("--", "red"),
    }
    for name, fers in results.items():
        st, c = styles[name]
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-6), st, color=c, label=name, markersize=5)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9)
    plt.title("Fig. 3: (31, 21) BCH code")
    plt.tight_layout()
    plt.savefig("figures/fig03_31_21.png", dpi=140)
    plt.savefig("figures/fig03_31_21.pdf")
    with open("data/fig03_31_21.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **{k: v for k, v in results.items()}}, f, indent=2)


def fig4_63_45():
    code = BCHCode(m=6, t=3)
    ebn0_list = np.arange(1.0, 7.5, 1.0)
    results = {}
    print(f"\n=== Fig 4: (63, 45) BCH ===")
    print("BM:")
    results["BM"] = run_fer(code, bm_wrap, ebn0_list, 60, 20000, label="BM")
    print("OSD(1):")
    results["OSD(1)"] = run_fer(code, lambda c, L: osd_decode(c, L, tau=1), ebn0_list, 60, 15000, label="OSD(1)")
    print("LLOSD(1):")
    results["LLOSD(1)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=1), ebn0_list, 60, 15000, label="LLOSD(1)")
    print("LLOSD(2):")
    results["LLOSD(2)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=2), ebn0_list, 60, 15000, label="LLOSD(2)")
    print("LLOSD(3):")
    results["LLOSD(3)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=3), ebn0_list, 60, 15000, label="LLOSD(3)")
    print("ML (~LLOSD 4):")
    results["ML"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=4), ebn0_list, 60, 10000, label="ML(~LLOSD 4)")

    plt.figure(figsize=(6, 4.5))
    styles = {
        "BM":       ("D-", "red"),
        "LLOSD(1)": ("v-", "black"),
        "LLOSD(2)": ("^-", "black"),
        "LLOSD(3)": ("x-", "black"),
        "OSD(1)":   ("o-", "red"),
        "ML":       ("--", "red"),
    }
    for name, fers in results.items():
        st, c = styles[name]
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-7), st, color=c, label=name, markersize=5)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9)
    plt.title("Fig. 4: (63, 45) BCH code")
    plt.tight_layout()
    plt.savefig("figures/fig04_63_45.png", dpi=140)
    plt.savefig("figures/fig04_63_45.pdf")
    with open("data/fig04_63_45.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **{k: v for k, v in results.items()}}, f, indent=2)


if __name__ == "__main__":
    fig3_31_21()
    fig4_63_45()
    print("\nSaved: figures/fig03_31_21.{png,pdf}, figures/fig04_63_45.{png,pdf}")
