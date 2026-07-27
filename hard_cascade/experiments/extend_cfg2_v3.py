"""Extend cfg2 pure BER sweep by [10.0, 10.5] dB and regenerate figures.

Config 1 (slow) is untouched; only the fast cfg2 pure baseline is extended so
its BER clearly crosses 1e-4. Merges into data/v3_results.json and rebuilds
all three v3 figures from the merged data.
"""
import sys
sys.path.insert(0, ".")
import json
import numpy as np

from hc_src.cascade_v3 import CONFIGS_V3, PureRSV3Codec
from experiments.main_v3 import (
    run_ber, plot_ber, plot_latency_innovations, plot_complexity,
)

D = json.load(open("data/v3_results.json"))

cfg = CONFIGS_V3["cfg2_255"]
pure = PureRSV3Codec(cfg)
extra = run_ber("Pure RS-BM (ext)", pure, cfg, [10.0, 10.5])

p = D["ber"]["cfg2_255"]["pure"]
# append only new SNR points (avoid duplicates)
existing = set(p["ebn0_db"])
for i, e in enumerate(extra["ebn0_db"]):
    if e in existing:
        continue
    for kk in ("ebn0_db", "ber", "fer", "avg_f2m", "n_frames", "n_frame_err"):
        p[kk].append(extra[kk][i])

json.dump(D, open("data/v3_results.json", "w"), indent=2)
print("Merged cfg2 pure extension. New cfg2 pure BER:")
print("  ", " ".join(f"{e}:{b:.1e}" for e, b in zip(p["ebn0_db"], p["ber"])))

# regenerate figures from merged data
plot_ber(D["ber"], "figures/ber_snr_v3")
plot_latency_innovations(D["latency"], "figures/latency_innovations_v3")
plot_complexity(D["complexity"], "figures/complexity_v3")
