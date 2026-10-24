"""Session 19, Step 4: device numbers -> quantum figures of merit.

From steps 1-3 (silica-dispersion baseline): photon coherence time, maximum
heralding efficiency (eta_esc), and the frequency-bin entangled state
    |psi> = sum_k c_k |p+k>_s |p-k>_i ,  c_k ~ sqrt(brightness_k)
for D = 2, 3, 4 (ASSUMED: the D bin pairs nearest the pump, k = 1..D).
QuTiP fidelity to the maximally entangled state, checked against the analytic
F = (sum_k sqrt(p_k))^2 / D, for D scalings x0.5, x1, x2, x5.
Compared with the measured fidelities of Kues et al. 2017.
Results -> results/s19_quantum.json
"""
import os, json
import numpy as np
import qutip as qt
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")   # project results folder
OUT = os.path.join(RESULTS, "s19_quantum.json")


def load(name):
    with open(os.path.join(RESULTS, name)) as fh:
        return json.load(fh)


S1, S2, S3 = load("s19_waveguide.json"), load("s19_coupler.json"), load("s19_dispersion.json")

# ---------- baseline device numbers ----------
CH = S2["chosen"]
GAP = CH["gap_m"]
QL, QI, QC = CH["QL"], CH["Qi"], CH["Qc"]
LINEWIDTH = 823e6                       # Hz, step 2 (nu/Q_L = 823.0 MHz)
ETA_ESC = CH["eta_esc"]                 # Q_L / Q_c
D_BASE = S3["pump"]["D_lumerical_ps_nm_km"]
FSR = S3["pump"]["FSR_Hz"]
SCALES = ("0.5", "1", "2", "5")
DIMS = (2, 3, 4)
F_MEASURED = {2: 0.885, 3: 0.809, 4: 0.766}     # Kues et al., Nature 546, 622 (2017)
LINEWIDTH_PAPER = 800e6                          # "about 800 MHz"
TAU_PAPER = 0.6e-9                               # "about 0.6 ns"
USABLE = 0.5                                     # ASSUMED: a bin is usable if brightness >= 0.5

results = []


def check(label, sim, theory, tol=1e-9):
    """Compare a computed number with its independent value (course-notebook style)."""
    ok = bool(np.isclose(sim, theory, atol=tol, rtol=0))
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {label}: simulation = {sim:.10f}, theory = {theory:.10f}")
    return ok


# ---------- 1. coherence time and heralding efficiency ----------
# Lorentzian line of FWHM Gamma (Hz): energy decays at rate 2 pi Gamma, the field
# correlation |g1(tau)| = exp(-pi Gamma |tau|).
tau_inv = 1 / LINEWIDTH                          # 1 / Gamma
tau_2pi = 1 / (2 * np.pi * LINEWIDTH)            # photon (energy) lifetime in the cavity
tau_pi = 1 / (np.pi * LINEWIDTH)                 # 1/e time of |g1| (Lorentzian coherence time)
tau_half = 1 / (2 * LINEWIDTH)                   # the paper gives no formula; this one reproduces ~0.6 ns
print("--- coherence time (Gamma = 823 MHz) ---")
for lab, t in (("1/Gamma", tau_inv), ("1/(2 pi Gamma)", tau_2pi), ("1/(pi Gamma)", tau_pi),
               ("1/(2 Gamma)", tau_half)):
    t800 = t * LINEWIDTH / LINEWIDTH_PAPER        # same definition at the paper's 800 MHz
    print(f"  {lab:15s} = {t * 1e9:.3f} ns   ({t / TAU_PAPER:.2f} x paper's ~0.6 ns; "
          f"{t800 * 1e9:.3f} ns at 800 MHz)")
print(f"  linewidth {LINEWIDTH / 1e6:.0f} MHz vs paper ~{LINEWIDTH_PAPER / 1e6:.0f} MHz "
      f"({(LINEWIDTH / LINEWIDTH_PAPER - 1) * 100:+.1f} %)")

print("\n--- maximum heralding efficiency ---")
# eta_esc = Q_L/Q_c = 1 - Q_L/Q_i (same thing, since 1/Q_L = 1/Q_i + 1/Q_c)
check("eta_esc = Q_L/Q_c vs 1 - Q_L/Q_i", ETA_ESC, 1 - QL / QI, tol=1e-12)
print(f"  eta_esc = {ETA_ESC:.4f}  (Q_L/Q_i = {QL / QI:.4f})")


# ---------- 2. frequency-bin states ----------
def bin_state(c):
    """sum_k c_k |k>_s |k>_i in a D x D two-photon space (k = 0 .. D-1 <-> bins 1..D)."""
    D = len(c)
    return sum(ck * qt.tensor(qt.basis(D, k), qt.basis(D, k)) for k, ck in enumerate(c))


def fid_analytic(p):
    """|<Phi_D|psi>|^2 for |psi> = sum sqrt(p_k)|kk>, |Phi_D> = sum |kk>/sqrt(D)."""
    return np.sum(np.sqrt(p)) ** 2 / len(p)


def brightness(scale):
    return np.array([b["brightness_numeric"] for b in S3["bins_vs_D_scale"][scale]["bins"]])


print("\n--- bin-imbalance fidelity (QuTiP vs analytic) ---")
fid = {}
for s in SCALES:
    B = brightness(s)
    fid[s] = {}
    for D in DIMS:
        p = B[:D] / B[:D].sum()                       # ASSUMED: bins k = 1..D
        c = np.sqrt(p)
        psi = bin_state(c)
        phi = bin_state(np.ones(D) / np.sqrt(D))      # maximally entangled
        F_expect = float(np.real(qt.expect(qt.ket2dm(phi), psi)))   # <Phi|rho|Phi>, course convention
        F_uhl = float(qt.fidelity(psi, phi)) ** 2     # qt.fidelity returns |<a|b>|, so square it
        F_an = fid_analytic(p)
        check(f"D x {s:>3s}, dim {D}: <Phi|rho|Phi> vs (sum sqrt p)^2/D", F_expect, F_an)
        check(f"D x {s:>3s}, dim {D}: qt.fidelity^2 vs analytic", F_uhl, F_an)
        fid[s][D] = {"p_k": p.tolist(), "c_k": c.tolist(), "F_qutip": F_expect,
                     "F_qutip_uhlmann_sq": F_uhl, "F_analytic": F_an,
                     "infidelity": 1 - F_expect}

print("\n  fidelity table (1 - F in brackets)")
print("  D-scale " + "".join(f"     dim {D}        " for D in DIMS))
for s in SCALES:
    print(f"  x{s:6s} " + "".join(f"  {fid[s][D]['F_qutip']:.7f} ({fid[s][D]['infidelity']:.1e})" for D in DIMS))

# ---------- 3. comparison with the measurement ----------
print("\n--- how much of the measured infidelity does bin imbalance explain? ---")
comp = {}
for D in DIMS:
    lost = 1 - F_MEASURED[D]
    comp[D] = {"F_measured": F_MEASURED[D], "measured_infidelity": lost,
               "explained_fraction": {s: fid[s][D]["infidelity"] / lost for s in SCALES}}
    print(f"  dim {D}: measured 1-F = {lost:.3f}; explained by imbalance: "
          + ", ".join(f"x{s}: {comp[D]['explained_fraction'][s] * 100:.1e} %" for s in SCALES))

# Illustration: how strong would the dispersion have to be for bin imbalance alone
# to reach the measured fidelity?  Lorentzian brightness with Delta_k = s * D2 k^2 / 2pi.
D2_HZ = S3["bins_vs_D_scale"]["1"]["D2_rad_s"] / (2 * np.pi)


def F_model(scale, D):
    k = np.arange(1, D + 1)
    B = 1 / (1 + (2 * scale * D2_HZ * k**2 / LINEWIDTH) ** 2)
    return fid_analytic(B / B.sum())


print("\n  illustration: D scaling needed for imbalance alone to give the measured F")
for D in DIMS:
    s_need = brentq(lambda s: F_model(s, D) - F_MEASURED[D], 1, 1e5)
    p_need = 1 / (1 + (2 * s_need * D2_HZ * np.arange(1, D + 1) ** 2 / LINEWIDTH) ** 2)
    comp[D]["D_scale_needed"] = s_need
    comp[D]["brightness_needed"] = p_need.tolist()
    print(f"  dim {D}: D x {s_need:.0f} (D = {s_need * D_BASE:,.0f} ps/(nm km)); brightness k=1..{D}: "
          f"{np.round(p_need, 3).tolist()}")

# ---------- 4. usable bins ----------
print(f"\n--- usable bins (brightness >= {USABLE}, ASSUMED threshold; computed up to k = 10) ---")
usable = {}
for s in ("0.5", "1", "2", "5", "10"):
    B = np.array([b["brightness_numeric"] for b in S3["bins_vs_D_scale"][s]["bins"]])
    n = int(np.sum(B >= USABLE))
    usable[s] = {"k_max": n, "at_edge": bool(n == len(B)), "brightness_k10": float(B[-1])}
    print(f"  D x {s:>3s} ({S3['bins_vs_D_scale'][s]['D_ps_nm_km']:.0f} ps/(nm km)): k = 1..{n}"
          f"{' (all computed bins, so >= 10)' if n == len(B) else ''}, i.e. {2 * n} resonances "
          f"around the pump; brightness at k=10 = {B[-1]:.3f}")

print(f"\n{sum(results)}/{len(results)} checks passed")

out = {
    "step": "s19_quantum",
    "baseline": {"gap_m": GAP, "Q_L": QL, "Q_i": QI, "Q_c": QC, "linewidth_Hz": LINEWIDTH,
                 "eta_esc": ETA_ESC, "eta_esc_check_1_minus_QL_over_Qi": 1 - QL / QI,
                 "D_ps_nm_km": D_BASE, "FSR_Hz": FSR},
    "coherence_time_s": {"1/Gamma": tau_inv, "1/(2 pi Gamma)": tau_2pi, "1/(pi Gamma)": tau_pi,
                         "1/(2 Gamma)": tau_half, "paper_approx": TAU_PAPER,
                         "paper_wording": "narrow spectral linewidth (~800 MHz) and the related long "
                                          "coherence time (~0.6 ns); no formula given"},
    "assumptions": {"bins": "D-level state uses bin pairs k = 1..D (nearest the pump)",
                    "amplitudes": "c_k = sqrt(brightness_k / sum brightness), real and in phase",
                    "fidelity": "F = <Phi_D|rho|Phi_D> (squared convention, as course notebooks)",
                    "usable_bin_threshold": USABLE},
    "fidelity": {s: {str(D): v for D, v in d.items()} for s, d in fid.items()},
    "comparison_with_measurement": {str(D): v for D, v in comp.items()},
    "usable_bins": usable,
    "checks_passed": f"{sum(results)}/{len(results)}",
}
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=2)
print(f"saved {OUT}")
