"""Fill LaTeX placeholders in report_v3.tex from data/v3_results.json.

Computes the Eb/N0 at which BER crosses 1e-4 (log-linear interpolation) for
cascade vs pure, per config, and the coding gain; fills complexity numbers.
"""
import json
import re
import numpy as np

with open("data/v3_results.json") as f:
    D = json.load(f)


def crossing(ebn0, ber, target=1e-4):
    """Eb/N0 where BER == target, via log-linear interpolation. Robust to zeros."""
    ebn0 = np.asarray(ebn0, float)
    ber = np.asarray(ber, float)
    # replace zeros with a small floor so log works; but crossing detection
    # uses the bracket around target
    for i in range(len(ber) - 1):
        b0, b1 = ber[i], ber[i + 1]
        if b0 >= target and (b1 <= target or b1 == 0.0):
            e0, e1 = ebn0[i], ebn0[i + 1]
            if b1 <= 0:
                # extrapolate a little past e0 using previous slope if possible
                b1_eff = min(b0 / 10.0, target / 10.0)
            else:
                b1_eff = b1
            # log-linear
            l0, l1 = np.log10(max(b0, 1e-12)), np.log10(max(b1_eff, 1e-12))
            lt = np.log10(target)
            frac = (lt - l0) / (l1 - l0) if l1 != l0 else 0.0
            return e0 + frac * (e1 - e0)
    return None


def fmt(x):
    return "N/A" if x is None else f"{x:.1f}"


subs = {}
gains = {}
for key, macros in [("cfg1_kp4", ("BERoneCASC", "BERonePURE", "BERoneGAIN")),
                    ("cfg2_255", ("BERtwoCASC", "BERtwoPURE", "BERtwoGAIN"))]:
    c = D["ber"][key]["cascade"]
    p = D["ber"][key]["pure"]
    ec = crossing(c["ebn0_db"], c["ber"])
    ep = crossing(p["ebn0_db"], p["ber"])
    g = (ep - ec) if (ec is not None and ep is not None) else None
    subs[macros[0]] = fmt(ec)
    subs[macros[1]] = fmt(ep)
    subs[macros[2]] = fmt(g)
    gains[key] = (ec, ep, g)

# complexity
cx = D["complexity"]
subs["CXoneP"] = f"{cx['cfg1_kp4']['pure']:.0f}"
subs["CXoneC"] = f"{cx['cfg1_kp4']['casc_conv']:.0f}"
subs["CXoneD"] = f"{cx['cfg1_kp4']['casc_direct']:.0f}"
subs["CXtwoP"] = f"{cx['cfg2_255']['pure']:.0f}"
subs["CXtwoC"] = f"{cx['cfg2_255']['casc_conv']:.0f}"
subs["CXtwoD"] = f"{cx['cfg2_255']['casc_direct']:.0f}"

print("Computed substitutions:")
for k, v in subs.items():
    print(f"  \\{k} -> {v}")

# Apply to report
path = "docs/report_v3.tex"
with open(path) as f:
    tex = f.read()
for macro, val in subs.items():
    tex = re.sub(r"\\" + macro + r"(?![a-zA-Z])", val, tex)
with open(path, "w") as f:
    f.write(tex)
print(f"\nFilled {path}")

# sanity: any leftover placeholders?
leftover = re.findall(r"\\BER[a-zA-Z]+|\\CX[a-zA-Z]+", tex)
print("Leftover placeholders:", set(leftover) if leftover else "none")
