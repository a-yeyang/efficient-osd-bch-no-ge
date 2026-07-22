"""Generate latency comparison figures for the delay analysis section.

Two figures:
  1. Bar chart: OSD(1) vs LLOSD/LLOSD-B/HSD latency at 4/5/6 dB, for both
     (63, 45) and (127, 99).
  2. Latency-vs-SNR line plot for (63, 45).
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path("/Users/chenshiyang.10/workspace/llosd_reproduction")

t1 = json.loads((PROJECT / "tables/table_I.json").read_text())
t3 = json.loads((PROJECT / "tables/table_III.json").read_text())
t4 = json.loads((PROJECT / "tables/table_IV.json").read_text())

snrs = ["4dB", "5dB", "6dB"]

# --- Figure 1: Bar chart -----------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# (63, 45)
algos1 = ["OSD(1)", "LLOSD(3)", "LLOSD-B(3)"]
data1 = np.array([[t1[a][s]["latency_us"] for s in snrs] for a in algos1])
x = np.arange(len(snrs))
width = 0.25
colors = ["#c62828", "#1976d2", "#2e7d32"]
for i, algo in enumerate(algos1):
    offset = (i - 1) * width
    ax1.bar(x + offset, data1[i], width, label=algo, color=colors[i])
ax1.set_xticks(x)
ax1.set_xticklabels(snrs)
ax1.set_ylabel("Latency (μs)")
ax1.set_title("(63, 45) BCH — per-frame decoding latency")
ax1.set_yscale("log")
ax1.legend()
ax1.grid(True, alpha=0.3, axis="y")

# Annotate speedup
for i in range(len(snrs)):
    osd = data1[0, i]
    llosd_b = data1[2, i]
    ax1.annotate(f"{osd/llosd_b:.0f}x",
                 xy=(x[i] + width, llosd_b),
                 xytext=(x[i] + width, llosd_b * 1.4),
                 ha="center", fontsize=9, color="#2e7d32",
                 fontweight="bold")

# (127, 99)
algos2 = ["OSD(1)", "LLOSD(3)", "HSD(1,6)"]
data2 = np.array([[t3[a][s]["latency_us"] for s in snrs] for a in algos2])
for i, algo in enumerate(algos2):
    offset = (i - 1) * width
    ax2.bar(x + offset, data2[i], width, label=algo, color=colors[i])
ax2.set_xticks(x)
ax2.set_xticklabels(snrs)
ax2.set_ylabel("Latency (μs)")
ax2.set_title("(127, 99) BCH — per-frame decoding latency")
ax2.set_yscale("log")
ax2.legend()
ax2.grid(True, alpha=0.3, axis="y")

# Annotate speedup
for i in range(len(snrs)):
    osd = data2[0, i]
    llosd = data2[1, i]
    ax2.annotate(f"{osd/llosd:.0f}x",
                 xy=(x[i], llosd),
                 xytext=(x[i], llosd * 1.4),
                 ha="center", fontsize=9, color="#1976d2",
                 fontweight="bold")

plt.suptitle("Decoding Latency Comparison — LLOSD/LLOSD-B/HSD vs OSD", fontsize=12)
plt.tight_layout()
plt.savefig(PROJECT / "figures" / "latency_bars.png", dpi=140)
plt.savefig(PROJECT / "figures" / "latency_bars.pdf")

# --- Figure 2: Latency vs SNR ------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
markers = {"OSD(1)": "o-", "LLOSD(3)": "s-", "LLOSD-B(3)": "^-"}
color_map = {"OSD(1)": "#c62828", "LLOSD(3)": "#1976d2", "LLOSD-B(3)": "#2e7d32"}
snr_vals = [4, 5, 6]
for algo in ["OSD(1)", "LLOSD(3)", "LLOSD-B(3)"]:
    ys = [t1[algo][s]["latency_us"] for s in snrs]
    ax.semilogy(snr_vals, ys, markers[algo], color=color_map[algo],
                label=algo, markersize=8, linewidth=2)
ax.set_xlabel("Eb/N0 (dB)")
ax.set_ylabel("Per-frame decoding latency (μs)")
ax.set_title("(63, 45) BCH: latency decreases with SNR (early-terminate effect)")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
ax.set_xticks(snr_vals)
plt.tight_layout()
plt.savefig(PROJECT / "figures" / "latency_vs_snr.png", dpi=140)
plt.savefig(PROJECT / "figures" / "latency_vs_snr.pdf")

# --- Speedup summary text -----------------------------------------------
print("=== Latency Speedup Summary ===")
print(f"\n(63, 45) @ 5 dB:")
osd = t1["OSD(1)"]["5dB"]["latency_us"]
llosd = t1["LLOSD(3)"]["5dB"]["latency_us"]
llosd_b = t1["LLOSD-B(3)"]["5dB"]["latency_us"]
print(f"  OSD(1)   : {osd:.1f} μs")
print(f"  LLOSD(3) : {llosd:.1f} μs  ({osd/llosd:.1f}x faster)")
print(f"  LLOSD-B(3): {llosd_b:.1f} μs  ({osd/llosd_b:.1f}x faster)")

print(f"\n(127, 99) @ 5 dB:")
osd3 = t3["OSD(1)"]["5dB"]["latency_us"]
llosd3 = t3["LLOSD(3)"]["5dB"]["latency_us"]
hsd = t3["HSD(1,6)"]["5dB"]["latency_us"]
print(f"  OSD(1)   : {osd3:.1f} μs")
print(f"  LLOSD(3) : {llosd3:.1f} μs  ({osd3/llosd3:.1f}x faster)")
print(f"  HSD(1,6) : {hsd:.1f} μs  ({osd3/hsd:.1f}x faster)")

print(f"\n(255, 223) @ 5 dB:")
osd4 = t4["OSD(1)"]["5dB"]["latency_us"]
llosd4 = t4["LLOSD(2)"]["5dB"]["latency_us"]
hsd4 = t4["HSD(1,8)"]["5dB"]["latency_us"]
print(f"  OSD(1)   : {osd4:.1f} μs")
print(f"  LLOSD(2) : {llosd4:.1f} μs  ({osd4/llosd4:.1f}x faster)")
print(f"  HSD(1,8) : {hsd4:.1f} μs  ({osd4/hsd4:.1f}x faster)")

print("\nSaved: figures/latency_bars.png, figures/latency_vs_snr.png")
