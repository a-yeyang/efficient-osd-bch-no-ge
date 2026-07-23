"""KPI verification: does Cascade A / B satisfy the 10% latency budget?

Reports:
- Per-method F_2^m ops at SNR points where FER converges (e.g., 9-10 dB)
- Cascade A / B vs Pure RS-BM ratio (Case a KPI)
- Cascade A / B vs Pure RS-LCC-BR ratio (Case b KPI)
- Clock cycle estimation under P ∈ {1, 4, 16, 64} parallelism
"""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

# Load existing n=255 scheme A data
with open("data/n255_scheme_a.json") as f:
    data = json.load(f)


def clock_cycles(f2m_ops: float, parallelism: int) -> float:
    """Estimate clock cycles under P-way parallelism.

    Simple model: each F_2^m mult = 1 cycle, but we can do P in parallel.
    Total cycles = ceil(f2m_ops / P).
    """
    return f2m_ops / parallelism


def kpi_check(baseline_ops, method_ops, threshold=0.10):
    """Check if method_ops <= (1 + threshold) * baseline_ops."""
    ratio = (method_ops - baseline_ops) / baseline_ops if baseline_ops > 0 else float("inf")
    return ratio, ratio <= threshold


def main():
    print("=" * 70)
    print("KPI Verification: n=255 Scheme A vs baselines")
    print("=" * 70)

    ebn0_list = data["pure_rs_bm"]["ebn0_db"]

    for i, ebn0 in enumerate(ebn0_list):
        if i >= len(data["scheme_a"]["ebn0_db"]):
            break
        # skip low SNR where FER is dominated by uncoded errors
        if data["scheme_a"]["fer"][i] > 0.5:
            continue

        bm_ops = data["pure_rs_bm"]["avg_f2m_ops"][i]
        lcc_ops = data["pure_rs_lccbr"]["avg_f2m_ops"][i]
        cascade_ops = data["scheme_a"]["avg_f2m_ops"][i]

        r_a, ok_a = kpi_check(bm_ops, cascade_ops)
        r_b, ok_b = kpi_check(lcc_ops, cascade_ops)

        print(f"\n@ {ebn0} dB (Cascade FER={data['scheme_a']['fer'][i]:.2e}):")
        print(f"  RS-BM:       {bm_ops:>10.0f} F₂ᵐ ops")
        print(f"  RS-LCC-BR:   {lcc_ops:>10.0f} F₂ᵐ ops")
        print(f"  Cascade A:   {cascade_ops:>10.0f} F₂ᵐ ops")
        print(f"  Case (a) vs BM:      +{r_a*100:+6.1f}%  KPI ≤+10%: {'✓' if ok_a else '✗'}")
        print(f"  Case (b) vs LCC-BR:  +{r_b*100:+6.1f}%  KPI ≤+10%: {'✓' if ok_b else '✗'}")

        # Clock cycles at various parallelism
        print(f"  Clock cycles:")
        for P in [1, 4, 16, 64]:
            c_bm = clock_cycles(bm_ops, P)
            c_lcc = clock_cycles(lcc_ops, P)
            c_cas = clock_cycles(cascade_ops, P)
            print(f"    P={P:<3d}: RS-BM={c_bm:>10.0f}  RS-LCC-BR={c_lcc:>10.0f}  Cascade={c_cas:>10.0f}")

    # Save KPI summary as JSON
    summary = {"ebn0_db": [], "kpi_case_a": [], "kpi_case_b": [],
               "rs_bm_ops": [], "rs_lccbr_ops": [], "cascade_a_ops": []}
    for i, ebn0 in enumerate(ebn0_list):
        if i >= len(data["scheme_a"]["ebn0_db"]):
            break
        bm_ops = data["pure_rs_bm"]["avg_f2m_ops"][i]
        lcc_ops = data["pure_rs_lccbr"]["avg_f2m_ops"][i]
        cascade_ops = data["scheme_a"]["avg_f2m_ops"][i]
        r_a, _ = kpi_check(bm_ops, cascade_ops)
        r_b, _ = kpi_check(lcc_ops, cascade_ops)
        summary["ebn0_db"].append(ebn0)
        summary["kpi_case_a"].append(r_a)
        summary["kpi_case_b"].append(r_b)
        summary["rs_bm_ops"].append(bm_ops)
        summary["rs_lccbr_ops"].append(lcc_ops)
        summary["cascade_a_ops"].append(cascade_ops)
    with open("data/kpi_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Plot latency bar chart at key SNR
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    key_snr_indices = [i for i, e in enumerate(ebn0_list)
                       if e in [9.0, 9.5, 10.0] and i < len(data["scheme_a"]["ebn0_db"])]

    for ax, ratio_key, title, threshold_label in [
        (ax1, "kpi_case_a", "Case (a): Cascade vs Pure RS-BM (hard baseline)", "+10%"),
        (ax2, "kpi_case_b", "Case (b): Cascade vs Pure RS-LCC-BR (soft baseline)", "+10%"),
    ]:
        snrs = [ebn0_list[i] for i in key_snr_indices]
        ratios = [summary[ratio_key][i] * 100 for i in key_snr_indices]
        colors = ["green" if r <= 10 else "red" for r in ratios]
        ax.bar(range(len(snrs)), ratios, color=colors, alpha=0.7)
        ax.axhline(10, color="black", linestyle="--", label=threshold_label)
        ax.set_xticks(range(len(snrs)))
        ax.set_xticklabels([f"{s} dB" for s in snrs])
        ax.set_ylabel("Latency increase (%)")
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis="y")
        ax.legend()
        for i, r in enumerate(ratios):
            ax.text(i, r + 1, f"{r:+.1f}%", ha="center", fontsize=10)

    plt.tight_layout()
    plt.savefig("figures/n255_kpi_bars.png", dpi=140)
    plt.savefig("figures/n255_kpi_bars.pdf")
    print("\nSaved figures/n255_kpi_bars.{png,pdf}")


if __name__ == "__main__":
    main()
