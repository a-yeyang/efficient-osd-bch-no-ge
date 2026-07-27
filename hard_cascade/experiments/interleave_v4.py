"""v4 experiment: interleaving before/after on burst & AWGN channels.

Advisor's 4th task (2026-07-27): add the Liu/Song/Wang APCC-2022 matrix
interleaving to our RS+BCH cascade (Config 1 / KP4 = the paper's exact code) and
measure how much it improves performance, before vs after interleaving.

Compares 5 schemes  x  3 channels:
  schemes  : none (no interleave), existing (symbol interleave),
             proposed w=1bit, w=2bits, w=1symbol
  channels : AWGN (p=0), burst p=0.3, burst p=0.75  (one-tap DFE propagation)

Outputs:
  figures/interleave_ber_v4.{png,pdf}   3-panel BER-SNR (paper Fig-5 style)
  figures/interleave_gain_v4.{png,pdf}  coding-gain summary bar chart
  data/v4_results.json
"""
from __future__ import annotations

import sys
sys.path.insert(0, ".")

import json
import time
import numpy as np
import matplotlib.pyplot as plt

from hc_src.cascade_v4 import InterleavedCascadeCodec, run_channel
from hc_src.upstream import OpCounters


MODES = ["none", "existing", "w1sym", "w2bits", "w1bit"]
MODE_LABEL = {
    "none": "No interleave",
    "existing": "Existing (symbol)",
    "w1sym": "Proposed w=1 symbol",
    "w2bits": "Proposed w=2 bits",
    "w1bit": "Proposed w=1 bit",
}
MODE_STYLE = {
    "none":     ("#7F8C8D", "x:", 5),
    "existing": ("#000000", "s-", 5),
    "w1sym":    ("#2E6BB0", "^-", 5),
    "w2bits":   ("#1E7B34", "d-", 5),
    "w1bit":    ("#C0392B", "o-", 5),
}

CHANNELS = [
    ("awgn",   0.0,  [7.5, 8.0, 8.5, 9.0, 9.5, 10.0]),
    ("p0.3",   0.3,  [8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0]),
    ("p0.75",  0.75, [9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0]),
]
CHAN_TITLE = {
    "awgn": "AWGN channel ($p=0$)",
    "p0.3": "Burst channel ($p=0.3$)",
    "p0.75": "Burst channel ($p=0.75$)",
}


def bits_of(sym, m):
    return ((sym[:, None] >> np.arange(m)) & 1).astype(np.int8).reshape(-1)


def run_ber(mode, p, ebn0_list, seed=20260727,
            min_frame_errors=60, min_frames=100, max_frames=3000, verbose=True):
    cod = InterleavedCascadeCodec(mode=mode, bch_decoder="direct")
    m = cod.m_rs
    mode_off = MODES.index(mode) * 101      # deterministic per-mode seed offset
    out = {"ebn0_db": [], "ber": [], "fer": [], "n_frames": [], "n_frame_err": []}
    for ebn0 in ebn0_list:
        rng = np.random.default_rng(seed + int(round(ebn0 * 100)) + mode_off)
        nf = nfe = nbe = tot = 0
        t0 = time.perf_counter()
        while nf < max_frames:
            msg = rng.integers(0, 1 << m, cod.n_info_symbols)
            coded = cod.encode(msg)
            rx = run_channel(coded, ebn0, cod.effective_rate, p, rng)
            mh, ok, info = cod.decode(rx)
            nf += 1
            e = int((bits_of(msg, m) != bits_of(mh, m)).sum())
            nbe += e
            tot += msg.size * m
            if e > 0:
                nfe += 1
            if nfe >= min_frame_errors and nf >= min_frames:
                break
        ber = nbe / max(1, tot)
        fer = nfe / max(1, nf)
        dt = time.perf_counter() - t0
        if verbose:
            print(f"  [{mode:9s} {('AWGN' if p==0 else f'p={p}'):6s}] "
                  f"@ {ebn0:4.1f} dB: BER={ber:.3e} FER={fer:.3e} ({nf}f {dt:.1f}s)")
        out["ebn0_db"].append(ebn0)
        out["ber"].append(ber)
        out["fer"].append(fer)
        out["n_frames"].append(nf)
        out["n_frame_err"].append(nfe)
        if ber == 0.0 and nfe == 0 and ebn0 >= ebn0_list[len(ebn0_list) // 2]:
            break
    return out


def crossing(ebn0, ber, target):
    """Eb/N0 at which BER == target via log-linear interpolation."""
    ebn0 = np.asarray(ebn0, float)
    ber = np.asarray(ber, float)
    for i in range(len(ber) - 1):
        b0, b1 = ber[i], ber[i + 1]
        if b0 >= target and (b1 <= target or b1 == 0.0):
            e0, e1 = ebn0[i], ebn0[i + 1]
            b1e = (min(b0 / 10.0, target / 10.0)) if b1 <= 0 else b1
            l0, l1 = np.log10(max(b0, 1e-12)), np.log10(max(b1e, 1e-12))
            lt = np.log10(target)
            frac = (lt - l0) / (l1 - l0) if l1 != l0 else 0.0
            return float(e0 + frac * (e1 - e0))
    return None


# =====================================================================
def plot_ber_panels(results, out_base):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    for ax, (chan, p, _) in zip(axes, CHANNELS):
        for mode in MODES:
            d = results[chan][mode]
            color, mk, ms = MODE_STYLE[mode]
            ax.semilogy(d["ebn0_db"], np.maximum(d["ber"], 1e-7), mk,
                        color=color, label=MODE_LABEL[mode], markersize=ms,
                        linewidth=1.5)
        ax.set_xlabel("$E_b/N_0$ (dB)")
        ax.set_ylabel("BER")
        ax.set_title(CHAN_TITLE[chan], fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_ylim(1e-7, 1)
        ax.legend(fontsize=8, loc="lower left")
    plt.suptitle("RS(544,514)+BCH(144,136) cascade: interleaving on AWGN vs burst channels "
                 "(hard-decision, PAM4)", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_base + ".png", dpi=140)
    plt.savefig(out_base + ".pdf")
    print(f"Saved {out_base}.png/.pdf")


def plot_gain(results, out_base, target=1e-4):
    """Coding gain vs the 'existing' scheme at a reachable target BER."""
    fig, ax = plt.subplots(figsize=(9, 5))
    chans = [c[0] for c in CHANNELS]
    proposed = ["w1sym", "w2bits", "w1bit"]
    x = np.arange(len(chans))
    w = 0.24
    gains = {mode: [] for mode in proposed}
    for chan in chans:
        e_exist = crossing(results[chan]["existing"]["ebn0_db"],
                            results[chan]["existing"]["ber"], target)
        for mode in proposed:
            e_m = crossing(results[chan][mode]["ebn0_db"],
                           results[chan][mode]["ber"], target)
            g = (e_exist - e_m) if (e_exist is not None and e_m is not None) else 0.0
            gains[mode].append(g)
    colors = {"w1sym": "#2E6BB0", "w2bits": "#1E7B34", "w1bit": "#C0392B"}
    for i, mode in enumerate(proposed):
        ax.bar(x + (i - 1) * w, gains[mode], w, label=MODE_LABEL[mode],
               color=colors[mode], alpha=0.9, edgecolor="black", linewidth=0.5)
        for j, g in enumerate(gains[mode]):
            ax.text(x[j] + (i - 1) * w, g + 0.01, f"{g:+.2f}", ha="center",
                    fontsize=7)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([CHAN_TITLE[c] for c in chans], fontsize=9)
    ax.set_ylabel(f"Coding gain vs Existing @ BER={target:.0e} (dB)")
    ax.set_title("Extra coding gain of proposed interleaving over the existing scheme",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(out_base + ".png", dpi=140)
    plt.savefig(out_base + ".pdf")
    print(f"Saved {out_base}.png/.pdf")
    return gains


# =====================================================================
def main():
    results = {}
    t_start = time.perf_counter()
    for chan, p, sweep in CHANNELS:
        print("=" * 78)
        print(f"CHANNEL {chan}  (p={p})")
        results[chan] = {}
        for mode in MODES:
            results[chan][mode] = run_ber(mode, p, sweep)
    print("=" * 78)
    print(f"Total sim time: {(time.perf_counter()-t_start)/60:.1f} min")

    plot_ber_panels(results, "figures/interleave_ber_v4")
    gains = plot_gain(results, "figures/interleave_gain_v4", target=1e-4)

    print("\nCoding gain vs existing @ BER=1e-4 (dB):")
    for mode in ("w1sym", "w2bits", "w1bit"):
        print(f"  {MODE_LABEL[mode]:22s}: " +
              "  ".join(f"{c[0]}={g:+.2f}" for c, g in zip(CHANNELS, gains[mode])))

    with open("data/v4_results.json", "w") as f:
        json.dump({"results": results, "gains": gains,
                   "channels": [(c[0], c[1]) for c in CHANNELS],
                   "modes": MODES}, f, indent=2)
    print("\nSaved data/v4_results.json")


if __name__ == "__main__":
    main()
