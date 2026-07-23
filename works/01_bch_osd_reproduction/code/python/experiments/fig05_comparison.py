"""Fig 5 reproduction: (63, 45) BCH with LLOSD/SLLOSD + YSVL OSD + CJ OSD."""
from work_paths import setup
setup()

import json
import numpy as np
import matplotlib.pyplot as plt

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.osd import osd_decode
from src.decoders import llosd_fast, sllosd_fast
from src.baselines import ysvl_osd_decode, cj_osd_decode


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
        print(f"  {label} @ {ebn0:.1f} dB: FER = {fer:.2e}  ({n_errors}/{n_frames})")
        fers.append(fer)
    return fers


def main():
    code = BCHCode(m=6, t=3)
    ebn0_list = np.arange(1.0, 7.5, 1.0)
    results = {}
    print("=== Fig 5: (63, 45) — LLOSD/SLLOSD/YSVL/CJ/OSD/ML ===")
    results["LLOSD(3)"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=3),
                                   ebn0_list, 60, 15000, label="LLOSD(3)")
    results["SLLOSD(3,2)"] = run_fer(code, lambda c, L: sllosd_fast(c, L, theta_tuple=(3,2)),
                                      ebn0_list, 60, 15000, label="SLLOSD(3,2)")
    results["YSVL OSD(1)"] = run_fer(code, lambda c, L: ysvl_osd_decode(c, L, tau=1),
                                      ebn0_list, 60, 15000, label="YSVL(1)")
    results["CJ OSD(1)"] = run_fer(code, lambda c, L: cj_osd_decode(c, L, tau=1),
                                    ebn0_list, 60, 15000, label="CJ(1)")
    results["OSD(1)"] = run_fer(code, lambda c, L: osd_decode(c, L, tau=1),
                                 ebn0_list, 60, 15000, label="OSD(1)")
    results["ML"] = run_fer(code, lambda c, L: llosd_fast(c, L, tau=4),
                             ebn0_list, 60, 10000, label="ML(~LLOSD 4)")

    plt.figure(figsize=(7, 5))
    styles = {
        "LLOSD(3)":   ("k^-",     None),
        "SLLOSD(3,2)":("s-",  "royalblue"),
        "YSVL OSD(1)":("D-",  "gold"),
        "CJ OSD(1)":  ("x-",  "orange"),
        "OSD(1)":     ("o-",  "red"),
        "ML":         ("r--",     None),
    }
    for name, fers in results.items():
        st, c = styles[name]
        kwargs = dict(label=name, markersize=5)
        if c: kwargs["color"] = c
        plt.semilogy(ebn0_list, np.maximum(fers, 1e-7), st, **kwargs)
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel("FER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9, loc="lower left")
    plt.title("Fig. 5: (63, 45) BCH — LLOSD/SLLOSD/YSVL/CJ vs OSD/ML")
    plt.tight_layout()
    plt.savefig("figures/fig05_63_45_comparison.png", dpi=140)
    plt.savefig("figures/fig05_63_45_comparison.pdf")
    with open("data/fig05_63_45.json", "w") as f:
        json.dump({"ebn0_db": ebn0_list.tolist(), **{k: v for k, v in results.items()}}, f, indent=2)


if __name__ == "__main__":
    main()
