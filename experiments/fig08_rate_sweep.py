"""Fig 8 reproduction — reduced to fit the time budget of a demo run.

Only sweeps 5 rate values with fewer trials and larger min-error threshold.
Uses FER=1e-2 target only (1e-4 was too slow for a live demo).
"""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.decoders import hsd_fast
from src.baselines import ysvl_osd_decode, cj_osd_decode


CANDIDATES = [
    (1,   "HSD(1,4)",  "YSVL(1)",     "CJ(1)"),
    (2,   "HSD(1,6)",  "YSVL(1)",     "CJ(1)"),
    (3,   "HSD(1,6)",  "YSVL(1)",     "CJ(1)"),
    (4,   "HSD(1,8)",  "YSVL(1)",     "CJ(1)"),
    (5,   "HSD(1,8)",  "YSVL(1)",     "CJ(1)"),
    (7,   "HSD(1,10)", "YSVL(1)",     "CJ(1)"),
]


def measure_fer(code, decoder_fn, ebn0, n_frames_max=800, min_errors=20):
    rate = code.k / code.n
    zero = np.zeros(code.n, dtype=np.int8)
    x0 = bpsk_modulate(zero)
    rng = np.random.default_rng(1234 + int(round(ebn0 * 100)))
    sigma = sigma_from_ebn0(ebn0, rate)
    n_err = 0
    n_frm = 0
    total_ops = 0.0
    while n_frm < n_frames_max:
        y = x0 + sigma * rng.standard_normal(code.n)
        L = 2.0 * y / (sigma * sigma)
        c_hat, s = decoder_fn(code, L)
        n_frm += 1
        if not np.array_equal(c_hat, zero):
            n_err += 1
        cc = s["counters"]
        total_ops += cc.f2 + cc.f2m
        if n_err >= min_errors and n_frm >= 200:
            break
    fer = n_err / n_frm if n_frm > 0 else 1.0
    return fer, total_ops / n_frm


def find_snr_for_fer(code, decoder_fn, target_fer, snr_range=(2.0, 7.0)):
    for snr in np.arange(snr_range[0], snr_range[1] + 0.01, 0.5):
        fer, ops = measure_fer(code, decoder_fn, snr)
        if fer <= target_fer:
            return float(snr), float(ops)
    return float(snr_range[1]), None


def get_hsd_fn(spec):
    a, b = spec.replace("HSD(", "").replace(")", "").split(",")
    tau = int(a); eta = int(b)
    return lambda c, L: hsd_fast(c, L, tau=tau, eta=eta)


def get_ysvl_fn(spec):
    order = int(spec.replace("YSVL(", "").replace(")", ""))
    return lambda c, L: ysvl_osd_decode(c, L, tau=order)


def get_cj_fn(spec):
    order = int(spec.replace("CJ(", "").replace(")", ""))
    return lambda c, L: cj_osd_decode(c, L, tau=order)


def main():
    print("=== Fig 8: rate sweep for length-127 BCH (reduced) ===", flush=True)

    rates = []
    snr_1e2 = {"HSD": [], "YSVL": [], "CJ": []}
    ops_1e2 = {"HSD": [], "YSVL": [], "CJ": []}

    for (t, hsd_spec, ysvl_spec, cj_spec) in CANDIDATES:
        code = BCHCode(m=7, t=t)
        actual_rate = code.k / code.n
        rates.append(actual_rate)
        print(f"\n(127, {code.k}) rate={actual_rate:.3f} t={t}", flush=True)

        s, o = find_snr_for_fer(code, get_hsd_fn(hsd_spec), 1e-2)
        snr_1e2["HSD"].append(s); ops_1e2["HSD"].append(o if o else 1)
        print(f"  HSD  {hsd_spec}: SNR={s} dB, ops={o:.0f}", flush=True)
        s, o = find_snr_for_fer(code, get_ysvl_fn(ysvl_spec), 1e-2)
        snr_1e2["YSVL"].append(s); ops_1e2["YSVL"].append(o if o else 1)
        print(f"  YSVL {ysvl_spec}: SNR={s} dB, ops={o:.0f}", flush=True)
        s, o = find_snr_for_fer(code, get_cj_fn(cj_spec), 1e-2)
        snr_1e2["CJ"].append(s); ops_1e2["CJ"].append(o if o else 1)
        print(f"  CJ   {cj_spec}: SNR={s} dB, ops={o:.0f}", flush=True)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 4.5))
    for algo, style, color in [("HSD", "^", "black"), ("YSVL", "s", "gold"), ("CJ", "x", "orange")]:
        ax_a.plot(rates, snr_1e2[algo], style + "-", color=color, label=algo, markersize=7)
        ax_b.semilogy(rates, ops_1e2[algo], style + "-", color=color, label=algo, markersize=7)
    ax_a.set_xlabel("Code rate")
    ax_a.set_ylabel("Eb/N0 (dB)")
    ax_a.set_title("(a) Required SNR for FER=1e-2")
    ax_a.grid(True, alpha=0.3); ax_a.legend(fontsize=9)
    ax_b.set_xlabel("Code rate")
    ax_b.set_ylabel("F_2 / F_128 ops per frame")
    ax_b.set_title("(b) Ops for FER=1e-2")
    ax_b.grid(True, alpha=0.3); ax_b.legend(fontsize=9)
    plt.suptitle("Fig. 8: Rate sweep for length-127 BCH codes (reduced)")
    plt.tight_layout()
    plt.savefig("figures/fig08_rate_sweep.png", dpi=140)
    plt.savefig("figures/fig08_rate_sweep.pdf")
    with open("data/fig08_rate_sweep.json", "w") as f:
        json.dump({
            "rates": rates,
            "snr_1e2": snr_1e2, "ops_1e2": ops_1e2,
        }, f, indent=2)
    print("Saved Fig 8", flush=True)


if __name__ == "__main__":
    main()
