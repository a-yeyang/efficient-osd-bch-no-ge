"""Scheme A vs B latency comparison: quantifies the Lagrange sharing benefit."""
import sys
sys.path.insert(0, '.')

import json
import numpy as np
import matplotlib.pyplot as plt

with open("data/all_results.json") as f:
    data = json.load(f)


def main():
    print("=" * 70)
    print("Scheme A vs Scheme B ops comparison (Lagrange sharing benefit)")
    print("=" * 70)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax, cfg_name in zip(axes, ["n255_A", "n255_B", "n127"]):
        rr = data[cfg_name]
        ebn = rr["scheme_a"]["ebn0_db"]
        ops_a = rr["scheme_a"]["avg_f2m_ops"]
        ops_b = rr["scheme_b"]["avg_f2m_ops"]

        # Compute savings
        savings = [(a - b) / a * 100 if a > 0 else 0
                   for a, b in zip(ops_a, ops_b)]

        ax.bar(range(len(ebn)), savings, color="green", alpha=0.7)
        ax.set_xticks(range(len(ebn)))
        ax.set_xticklabels([f"{s:.1f}" for s in ebn], rotation=45)
        ax.set_xlabel("Eb/N0 (dB)")
        ax.set_ylabel("Ops savings (%)")
        ax.set_title(f"{cfg_name}: Scheme B vs A (Lagrange sharing)")
        ax.grid(True, alpha=0.3, axis="y")
        ax.axhline(0, color="black", linewidth=0.5)

        for i, s in enumerate(savings):
            if abs(s) > 0.1:
                ax.text(i, s + 0.1 if s > 0 else s - 0.5,
                        f"{s:.1f}%", ha="center", fontsize=7)

        print(f"\n{cfg_name}:")
        for e, a, b, s in zip(ebn, ops_a, ops_b, savings):
            print(f"  @ {e:.1f} dB: A={a:.0f}, B={b:.0f}, saved={s:.2f}%")

    plt.tight_layout()
    plt.savefig("figures/scheme_ab_savings.png", dpi=140)
    plt.savefig("figures/scheme_ab_savings.pdf")
    print("\nSaved figures/scheme_ab_savings.{png,pdf}")


if __name__ == "__main__":
    main()
