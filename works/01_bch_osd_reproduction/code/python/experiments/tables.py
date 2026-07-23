"""Tables I–IV reproduction: decoding complexity and latency.

Table I  ((63, 45) at Eb/N0 ∈ {4, 5, 6} dB) — OSD(1), LLOSD(3), LLOSD-B(3)
Table II ((63, 45) at Eb/N0 ∈ {4, 5, 6} dB) — LLOSD-B(3), SLLOSD-B(3,2), YSVL OSD(1), CJ OSD(1)
Table III ((127, 99) at Eb/N0 ∈ {4, 5, 6} dB) — OSD(1), LLOSD(3), HSD(1,8), HSD(1,6)
Table IV  ((255, 223) at Eb/N0 ∈ {4, 5, 6} dB) — OSD(1), LLOSD(2), PLCC(8), HSD(1,8)

Each table has (A) Decoding Complexity (F_2, F_2^m, floating ops) and (B)
Decoding Latency (μs).  We instrument each decoder and average over MC frames.
"""
from work_paths import setup
setup()

import json
import numpy as np

from src.bch import BCHCode, bpsk_modulate, sigma_from_ebn0
from src.osd import osd_decode
from src.decoders import llosd_fast, sllosd_fast, hsd_fast
from src.baselines import ysvl_osd_decode, cj_osd_decode, plcc_decode


def instrument(code, decoder_fn, ebn0, n_trials=500, seed=0):
    """Return avg (f2, f2m, fp, latency_us) over n_trials frames."""
    rate = code.k / code.n
    zero = np.zeros(code.n, dtype=np.int8)
    x0 = bpsk_modulate(zero)
    rng = np.random.default_rng(seed + int(round(ebn0 * 100)))
    sigma = sigma_from_ebn0(ebn0, rate)
    f2, f2m, fp, lat = 0.0, 0.0, 0.0, 0.0
    for _ in range(n_trials):
        y = x0 + sigma * rng.standard_normal(code.n)
        L = 2.0 * y / (sigma * sigma)
        _, s = decoder_fn(code, L)
        c = s["counters"]
        f2  += c.f2
        f2m += c.f2m
        fp  += c.fp
        lat += c.latency_us
    return (f2 / n_trials, f2m / n_trials, fp / n_trials, lat / n_trials)


def build_table(code, decoders, ebn0s, n_trials=500):
    """decoders: list of (name, decoder_fn). Returns dict of dicts."""
    tbl = {}
    for name, fn in decoders:
        row = {}
        for ebn0 in ebn0s:
            f2, f2m, fp, lat = instrument(code, fn, ebn0, n_trials=n_trials)
            row[f"{ebn0}dB"] = {"F2": f2, "F2m": f2m, "FP": fp, "latency_us": lat}
            print(f"  {name} @ {ebn0} dB: F2={f2:.2e}  F2m={f2m:.2e}  FP={fp:.2f}  lat={lat:.1f}μs")
        tbl[name] = row
    return tbl


def format_scientific(x):
    if x < 1:
        return f"{x:.2f}"
    return f"{x:.2e}"


def dump_markdown(name, tbl, ebn0s, path):
    lines = [f"# {name}\n"]
    lines.append("| Algorithm | Eb/N0 (dB) | F_2 ops | F_2^m ops | Floating ops | Latency (μs) |")
    lines.append("|---|---|---|---|---|---|")
    for algo, row in tbl.items():
        first = True
        for ebn0 in ebn0s:
            r = row[f"{ebn0}dB"]
            name_col = algo if first else ""
            lines.append(f"| {name_col} | {ebn0} | {format_scientific(r['F2'])} | "
                          f"{format_scientific(r['F2m'])} | {r['FP']:.1f} | {r['latency_us']:.1f} |")
            first = False
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ebn0s = [4, 5, 6]

    # Table I
    print("=== Table I: (63, 45) ===")
    code = BCHCode(m=6, t=3)
    decoders = [
        ("OSD(1)", lambda c, L: osd_decode(c, L, tau=1)),
        ("LLOSD(3)", lambda c, L: llosd_fast(c, L, tau=3, use_binary_reencoding=False)),
        ("LLOSD-B(3)", lambda c, L: llosd_fast(c, L, tau=3, use_binary_reencoding=True)),
    ]
    tbl1 = build_table(code, decoders, ebn0s, n_trials=500)
    with open("tables/table_I.json", "w") as f:
        json.dump(tbl1, f, indent=2)
    dump_markdown("Table I. Numerical Results in Decoding the (63, 45) BCH Code",
                  tbl1, ebn0s, "tables/table_I.md")

    # Table II
    print("\n=== Table II: (63, 45) — with YSVL / CJ ===")
    decoders2 = [
        ("LLOSD-B(3)", lambda c, L: llosd_fast(c, L, tau=3, use_binary_reencoding=True)),
        ("SLLOSD-B(3,2)", lambda c, L: sllosd_fast(c, L, theta_tuple=(3, 2), use_binary_reencoding=True)),
        ("YSVL OSD(1)", lambda c, L: ysvl_osd_decode(c, L, tau=1)),
        ("CJ OSD(1)", lambda c, L: cj_osd_decode(c, L, tau=1)),
    ]
    tbl2 = build_table(code, decoders2, ebn0s, n_trials=500)
    with open("tables/table_II.json", "w") as f:
        json.dump(tbl2, f, indent=2)
    dump_markdown("Table II. Numerical Results in Decoding the (63, 45) BCH Code",
                  tbl2, ebn0s, "tables/table_II.md")

    # Table III
    print("\n=== Table III: (127, 99) ===")
    code = BCHCode(m=7, t=4)
    decoders3 = [
        ("OSD(1)", lambda c, L: osd_decode(c, L, tau=1)),
        ("LLOSD(3)", lambda c, L: llosd_fast(c, L, tau=3, use_binary_reencoding=True)),
        ("HSD(1,8)", lambda c, L: hsd_fast(c, L, tau=1, eta=8)),
        ("HSD(1,6)", lambda c, L: hsd_fast(c, L, tau=1, eta=6)),
    ]
    tbl3 = build_table(code, decoders3, ebn0s, n_trials=200)
    with open("tables/table_III.json", "w") as f:
        json.dump(tbl3, f, indent=2)
    dump_markdown("Table III. Numerical Results in Decoding the (127, 99) BCH Code",
                  tbl3, ebn0s, "tables/table_III.md")

    # Table IV
    print("\n=== Table IV: (255, 223) ===")
    code = BCHCode(m=8, t=4)
    decoders4 = [
        ("OSD(1)", lambda c, L: osd_decode(c, L, tau=1)),
        ("LLOSD(2)", lambda c, L: llosd_fast(c, L, tau=2, use_binary_reencoding=True)),
        ("PLCC(8)", lambda c, L: plcc_decode(c, L, eta=8)),
        ("HSD(1,8)", lambda c, L: hsd_fast(c, L, tau=1, eta=8)),
    ]
    tbl4 = build_table(code, decoders4, ebn0s, n_trials=100)
    with open("tables/table_IV.json", "w") as f:
        json.dump(tbl4, f, indent=2)
    dump_markdown("Table IV. Numerical Results in Decoding the (255, 223) BCH Code",
                  tbl4, ebn0s, "tables/table_IV.md")


if __name__ == "__main__":
    main()
