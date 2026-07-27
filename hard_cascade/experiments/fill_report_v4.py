"""Fill gain placeholders in report_v4.tex from data/v4_results.json.

Computes coding gain of each proposed interleaving scheme relative to the
'existing' scheme, at BER=1e-4 (log-linear crossing interpolation), for each of
the three channels (AWGN, burst p=0.3, burst p=0.75).
"""
import json
import re
import numpy as np

with open("data/v4_results.json") as f:
    D = json.load(f)

results = D["results"]


def crossing(ebn0, ber, target=1e-4):
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


def gain(chan, mode, target=1e-4):
    e_exist = crossing(results[chan]["existing"]["ebn0_db"],
                       results[chan]["existing"]["ber"], target)
    e_m = crossing(results[chan][mode]["ebn0_db"],
                   results[chan][mode]["ber"], target)
    if e_exist is None or e_m is None:
        return None
    return e_exist - e_m


def fmt(x):
    return "N/A" if x is None else f"{x:+.2f}"


chans = ["awgn", "p0.3", "p0.75"]
subs = {
    "GAINsymAWGN": fmt(gain("awgn", "w1sym")),
    "GAINsymPthree": fmt(gain("p0.3", "w1sym")),
    "GAINsymPseven": fmt(gain("p0.75", "w1sym")),
    "GAINtwoAWGN": fmt(gain("awgn", "w2bits")),
    "GAINtwoPthree": fmt(gain("p0.3", "w2bits")),
    "GAINtwoPseven": fmt(gain("p0.75", "w2bits")),
    "GAINbitAWGN": fmt(gain("awgn", "w1bit")),
    "GAINbitPthree": fmt(gain("p0.3", "w1bit")),
    "GAINbitPseven": fmt(gain("p0.75", "w1bit")),
}

# best gain at p=0.75 across proposed schemes (for the headline number)
best = max(
    (g for g in (gain("p0.75", m) for m in ("w1sym", "w2bits", "w1bit")) if g is not None),
    default=None)
subs["GAINbest"] = "N/A" if best is None else f"{best:.2f}"

print("Computed substitutions:")
for k, v in subs.items():
    print(f"  \\{k} -> {v}")

path = "docs/report_v4.tex"
with open(path) as f:
    tex = f.read()
for macro, val in subs.items():
    tex = re.sub(r"\\" + macro + r"(?![a-zA-Z])", val, tex)
with open(path, "w") as f:
    f.write(tex)

leftover = re.findall(r"\\GAIN[a-zA-Z]+", tex)
print(f"\nFilled {path}")
print("Leftover placeholders:", set(leftover) if leftover else "none")
