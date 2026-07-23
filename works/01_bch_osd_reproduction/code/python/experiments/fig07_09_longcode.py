"""Fig 7 & Fig 9 reproduction: long BCH codes (127, 99) and (255, 223).

Fig. 7 (127, 99): OSD (1), OSD (2), LLOSD (2), LLOSD (3), HSD (1, 4), HSD (1, 6), HSD (1, 8)
Fig. 9 (255, 223): BM, OSD (1), LLOSD (2), PLCC (6), PLCC (8), HSD (1, 4..8), ML

Compute is much heavier for these codes; we shrink Monte Carlo per point.
"""
from work_paths import setup
setup()

import json
import numpy as np
import matplotlib.pyplot as plt

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.osd import osd_decode
from src.decoders import llosd_fast, hsd_fast
from src.baselines import plcc_decode


def bm_wrap(code, L):
    r_hard = (L < 0).astype(np.int8)
    dec, ok = code.bm_decode(r_hard)
    return (dec if ok else r_hard), {"counters": None}


def run_fer(code, decoder_fn, ebn0_list, min_errors, max_frames, seed=0, label=""):
    rate = code.k / code.n
    zero = np.zeros(code.n, dtype=np.int8)
    x0 = bpsk_modulate(zero)
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
        print(f"  {label} @ {ebn0:.1f} dB: FER = {fer:.2e}  ({n_errors}/{n_frames})  {dt:.1f}s")
        fers.append(fer)
    return fers


def fig7_127_99():
    code = BCHCode(m=7, t=4)
    print(f"\n=== Fig 7: ({code.n}, {code.k}) ===")
    ebn0_list = np.arange(3.0, 7.5, 1.0)
    results = {}
    results["OSD(1)"] = run_fer(code, lambda c, L: osd_decode(c, L, tau=1),
                                 ebn0_list, 40, 3000, label="OSD(1)")
    results["OSD(2)"] = run_fer(code, lambda c, L: osd_decode(c, L, tau=2),
                                 ebn0_list, 40, 3000, label="OSD(2)")
    results["LLOSD(2)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=2),
                                   ebn0_list, 40, 3000, label="LLOSD(2)")
    results["LLOSD(3)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=3),
                                   ebn0_list, 40, 2000, label="LLOSD(3)")
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
    for name, fers in results.items():
        st, c = styles[name]
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-6), st, color=c, label=name, markersize=4)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8)
    plt.title("Fig. 7: (127, 99) BCH")
    plt.tight_layout()
    plt.savefig("figures/fig07_127_99.png", dpi=140)
    plt.savefig("figures/fig07_127_99.pdf")
    with open("data/fig07_127_99.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **results}, f, indent=2)


def fig9_255_223():
    code = BCHCode(m=8, t=4)
    print(f"\n=== Fig 9: ({code.n}, {code.k}) ===")
    ebn0_list = np.arange(3.5, 7.5, 0.5)
    results = {}
    # For the 255-length code, per-frame decoding is quite a bit slower.
    # Cap iterations tighter.
    results["BM"] = run_fer(code, bm_wrap, ebn0_list, 30, 3000, label="BM")
    results["OSD(1)"] = run_fer(code, lambda c, L: osd_decode(c, L, tau=1),
                                 ebn0_list, 30, 1000, label="OSD(1)")
    results["LLOSD(2)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=2),
                                   ebn0_list, 30, 800, label="LLOSD(2)")
    results["PLCC(6)"] = run_fer(code, lambda c, L: plcc_decode(c, L, eta=6),
                                  ebn0_list, 30, 500, label="PLCC(6)")
    results["PLCC(8)"] = run_fer(code, lambda c, L: plcc_decode(c, L, eta=8),
                                  ebn0_list, 30, 400, label="PLCC(8)")
    results["HSD(1,4)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=4),
                                   ebn0_list, 30, 500, label="HSD(1,4)")
    results["HSD(1,6)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=6),
                                   ebn0_list, 30, 400, label="HSD(1,6)")
    results["HSD(1,8)"] = run_fer(code, lambda c, L: hsd_fast(c, L, tau=1, eta=8),
                                   ebn0_list, 30, 300, label="HSD(1,8)")
    # No true ML for 255; use HSD(1,8) as ML proxy.
    results["ML"] = results["HSD(1,8)"]

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
        "ML":       ("r--", None),
    }
    for name, fers in results.items():
        st, c = styles.get(name, ("-", None))
        kwargs = dict(label=name, markersize=4)
        if c: kwargs["color"] = c
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-5), st, **kwargs)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8, loc="lower left")
    plt.title("Fig. 9: (255, 223) BCH")
    plt.tight_layout()
    plt.savefig("figures/fig09_255_223.png", dpi=140)
    plt.savefig("figures/fig09_255_223.pdf")
    with open("data/fig09_255_223.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **results}, f, indent=2)


if __name__ == "__main__":
    fig7_127_99()
    fig9_255_223()
