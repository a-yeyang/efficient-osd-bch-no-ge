"""Fig 10, 11 reproduction:
  Fig 10: average N_TEPs actually processed by LLOSD with tau ∈ {1, 2}, for
          (31, 21) and (63, 45), across Eb/N0.
  Fig 11: average N_TVs actually processed by LCC-BR portion of HSD with
          tau = 1 and η ∈ {4, 6}, for (63, 39) and (255, 223).
"""
from work_paths import setup
setup()

import json
import numpy as np
import matplotlib.pyplot as plt

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.decoders import llosd_fast, hsd_fast


def avg_stat(code, decoder_fn, ebn0_list, n_trials, stat_key="n_teps", seed=0, label=""):
    rate = code.k / code.n
    zero = np.zeros(code.n, dtype=np.int8)
    x0 = bpsk_modulate(zero)
    avgs = []
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
        sigma = sigma_from_ebn0(ebn0, rate)
        vals = []
        for _ in range(n_trials):
            y = x0 + sigma * rng.standard_normal(code.n)
            L = 2.0 * y / (sigma * sigma)
            _, s = decoder_fn(code, L)
            vals.append(s.get(stat_key, 0))
        m = float(np.mean(vals))
        avgs.append(m)
        print(f"  {label} @ {ebn0:.1f} dB: avg {stat_key} = {m:.3f}")
    return avgs


def fig10():
    print("=== Fig 10: avg N_TEPs in LLOSD ===")
    ebn0_list = np.arange(3.0, 7.5, 0.5)

    code31 = BCHCode(m=5, t=2)
    r31_t1 = avg_stat(code31, lambda c, L: llosd_fast(c, L, tau=1),
                      ebn0_list, 400, "n_teps", label="(31,21) tau=1")
    r31_t2 = avg_stat(code31, lambda c, L: llosd_fast(c, L, tau=2),
                      ebn0_list, 400, "n_teps", label="(31,21) tau=2")

    code63 = BCHCode(m=6, t=3)
    r63_t1 = avg_stat(code63, lambda c, L: llosd_fast(c, L, tau=1),
                      ebn0_list, 200, "n_teps", label="(63,45) tau=1")
    r63_t2 = avg_stat(code63, lambda c, L: llosd_fast(c, L, tau=2),
                      ebn0_list, 200, "n_teps", label="(63,45) tau=2")

    plt.figure(figsize=(6, 4))
    plt.semilogy(ebn0_list, r63_t2, "s-", color="k", label="(63,45) BCH, τ=2")
    plt.semilogy(ebn0_list, r63_t1, "o-", color="k", label="(63,45) BCH, τ=1")
    plt.semilogy(ebn0_list, r31_t2, "s--", color="royalblue", label="(31,21) BCH, τ=2")
    plt.semilogy(ebn0_list, r31_t1, "o--", color="royalblue", label="(31,21) BCH, τ=1")
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel(r"$\overline{N}_{\mathrm{TEPs}}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8)
    plt.title("Fig. 10: avg N_TEPs processed in LLOSD")
    plt.tight_layout()
    plt.savefig("figures/fig10_nteps.png", dpi=140)
    plt.savefig("figures/fig10_nteps.pdf")
    with open("data/fig10_nteps.json", "w") as f:
        json.dump({
            "ebn0_db": ebn0_list.tolist(),
            "n_teps_63_45_tau1": r63_t1,
            "n_teps_63_45_tau2": r63_t2,
            "n_teps_31_21_tau1": r31_t1,
            "n_teps_31_21_tau2": r31_t2,
        }, f, indent=2)


def fig11():
    print("\n=== Fig 11: avg N_TVs in HSD ===")
    ebn0_list = np.arange(2.0, 6.5, 0.5)

    # (63, 39) has t = 4 (design distance 9).
    code63_39 = BCHCode(m=6, t=4)  # this gives (63, k) — need to check
    # Actually we want a (63, 39) code; that requires t=4 for narrow-sense
    # primitive BCH. Verify:
    print(f"  BCHCode(m=6, t=4) has k = {code63_39.k}")
    if code63_39.k != 39:
        # (63,39) is a different design; try t=5? Skip if wrong.
        print(f"  Falling back to k={code63_39.k}")
    r63_t1_eta4 = avg_stat(code63_39, lambda c, L: hsd_fast(c, L, tau=1, eta=4),
                            ebn0_list, 100, "n_tvs", label="(63,39) η=4")
    r63_t1_eta6 = avg_stat(code63_39, lambda c, L: hsd_fast(c, L, tau=1, eta=6),
                            ebn0_list, 100, "n_tvs", label="(63,39) η=6")

    code255 = BCHCode(m=8, t=4)
    # (255,223) is much slower; use fewer trials.
    r255_t1_eta4 = avg_stat(code255, lambda c, L: hsd_fast(c, L, tau=1, eta=4),
                             ebn0_list, 30, "n_tvs", label="(255,223) η=4")
    r255_t1_eta6 = avg_stat(code255, lambda c, L: hsd_fast(c, L, tau=1, eta=6),
                             ebn0_list, 30, "n_tvs", label="(255,223) η=6")

    plt.figure(figsize=(6, 4))
    plt.plot(ebn0_list, r63_t1_eta6, "s-", color="k", label="(63,39) BCH, η=6")
    plt.plot(ebn0_list, r63_t1_eta4, "o-", color="k", label="(63,39) BCH, η=4")
    plt.plot(ebn0_list, r255_t1_eta6, "s--", color="royalblue", label="(255,223) BCH, η=6")
    plt.plot(ebn0_list, r255_t1_eta4, "o--", color="royalblue", label="(255,223) BCH, η=4")
    plt.xlabel("Eb/N0 (dB)")
    plt.ylabel(r"$\overline{N}_{\mathrm{TVs}}$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=8)
    plt.title("Fig. 11: avg N_TVs processed in HSD (τ=1)")
    plt.tight_layout()
    plt.savefig("figures/fig11_ntvs.png", dpi=140)
    plt.savefig("figures/fig11_ntvs.pdf")
    with open("data/fig11_ntvs.json", "w") as f:
        json.dump({
            "ebn0_db": ebn0_list.tolist(),
            "n_tvs_63_39_eta4": r63_t1_eta4,
            "n_tvs_63_39_eta6": r63_t1_eta6,
            "n_tvs_255_223_eta4": r255_t1_eta4,
            "n_tvs_255_223_eta6": r255_t1_eta6,
        }, f, indent=2)


if __name__ == "__main__":
    fig10()
    fig11()
