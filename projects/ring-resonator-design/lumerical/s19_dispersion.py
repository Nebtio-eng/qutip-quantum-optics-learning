"""Session 19, Step 3: material dispersion, resonance alignment and bin brightness.

1. Same waveguide as step 1, but the silica cladding now follows the Malitson
   Sellmeier formula (Lumerical native Sellmeier material).  Core stays n = 1.70
   (Hydex dispersion unknown).  FDE sweep 1500-1600 nm -> n_eff, n_g, D.
   Checks: n_g, D Lumerical vs our derivative (as step 1); material index vs
   Malitson in Python; n_eff vs constant-index solves at n_clad = Malitson(lam).
2. Resonances from n_eff(lam) L = m lam with L fixed from step 1.  For bins
   k = 1..10: Delta_k = nu_(p+k) + nu_(p-k) - 2 nu_p, brightness
   1/(1 + (2 Delta_k / linewidth)^2), for D scaled x0.5, 1, 2, 5, 10.
   Check: Delta_k vs analytic D2 k^2 / 2pi, D2 = lam^2 D (2 pi FSR)^2 / (2 pi n_g).
3. TM0 and the guided higher-order modes: their resonances in the same ring,
   listed when within one linewidth of a TE0 resonance (bins -10..+10).
   Sensitivity check: bent-waveguide FDE at R (n_eff shift, bend loss).
4. Plot brightness vs k -> results/s19_dispersion_brightness.png
Results -> results/s19_dispersion.json
"""
import os, json, time, importlib.util
import numpy as np
from scipy.constants import c
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")   # project results folder
OUT = os.path.join(RESULTS, "s19_dispersion.json")
PNG = os.path.join(RESULTS, "s19_dispersion_brightness.png")
LUMAPI = r"C:\Program Files\Lumerical\v261\api\python\lumapi.py"


def load_lumapi():
    spec = importlib.util.spec_from_file_location("lumapi", LUMAPI)
    lumapi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lumapi)
    return lumapi


# ---------- inputs from steps 1 and 2 ----------
with open(os.path.join(RESULTS, "s19_waveguide.json")) as fh:
    S1 = json.load(fh)
with open(os.path.join(RESULTS, "s19_coupler.json")) as fh:
    S2 = json.load(fh)
W, H = S1["design"]["W_m"], S1["design"]["H_m"]
N_CORE = S1["design"]["n_core"]                   # ASSUMED constant: Hydex dispersion unknown
L_RING = S1["at_1550nm"]["L_m"]                   # accepted geometry, kept fixed
R = S1["at_1550nm"]["R_m"]
LINEWIDTH = 823e6                                 # Hz, step 2 (Q_L = 235,000)
LAM0 = 1.55e-6
LAM_START, LAM_STOP, NPTS = 1.50e-6, 1.60e-6, 21
SIM_SPAN, CELLS = 6e-6, 300                       # step 1 fine mesh (20 nm)
K_BINS = np.arange(1, 11)
D_SCALES = (0.5, 1, 2, 5, 10)
PS_NM_KM = 1e6                                    # 1 s/m^2 = 1e6 ps/(nm km)

# Malitson (1965) fused silica, lambda in um: n^2 = 1 + sum B lam^2 / (lam^2 - C^2)
SELL_B = (0.6961663, 0.4079426, 0.8974794)
SELL_C = (0.0684043, 0.1162414, 9.896161)        # um


def malitson(lam_m):
    l2 = (np.asarray(lam_m) * 1e6) ** 2
    return np.sqrt(1 + sum(b * l2 / (l2 - ci**2) for b, ci in zip(SELL_B, SELL_C)))


def poly_fit(lam, neff, deg=6, s=50e-9):
    """Smooth n_eff(lam): degree-6 polynomial in x = (lam - 1550 nm)/50 nm."""
    x = (lam - LAM0) / s
    p = np.polyfit(x, neff, deg)
    resid = float(np.max(np.abs(np.polyval(p, x) - neff)))
    n = lambda l: np.polyval(p, (np.asarray(l) - LAM0) / s)
    dn = lambda l: np.polyval(np.polyder(p, 1), (np.asarray(l) - LAM0) / s) / s
    d2n = lambda l: np.polyval(np.polyder(p, 2), (np.asarray(l) - LAM0) / s) / s**2
    return {"p": p, "resid": resid, "n": n, "ng": lambda l: n(l) - l * dn(l),
            "D": lambda l: -(l / c) * d2n(l)}


# ======================= Lumerical =======================
def mode_list(mode, ntrial):
    n = int(np.real(mode.findmodes()))
    out = []
    for k in range(1, n + 1):
        name = f"FDE::data::mode{k}"
        te = float(np.real(mode.getdata(name, "TE polarization fraction")))
        ne = float(np.real(mode.getdata(name, "neff")).flatten()[0])
        loss = float(np.real(mode.getdata(name, "loss")))       # dB/m
        y = np.real(mode.getdata(name, "y")).flatten()
        z = np.real(mode.getdata(name, "z")).flatten()
        Ey = np.squeeze(mode.getdata(name, "Ey"))
        Ez = np.squeeze(mode.getdata(name, "Ez"))
        # z-mirror class: "A" = same symmetry as TE0 (Ey even in z, Ez odd).
        # In an up-down symmetric ring only class-A modes can couple to TE0.
        def pz(E):
            k0 = np.unravel_index(np.argmax(np.abs(E)), E.shape)
            Er = np.real(E * np.conj(E[k0]) / np.abs(E[k0]))
            Em = np.array([np.interp(-z, z, Er[i, :]) for i in range(Er.shape[0])])
            return float(np.sum(Er * Em) / np.sum(Er * Er))
        dom_y = np.sum(np.abs(Ey) ** 2) >= np.sum(np.abs(Ez) ** 2)
        cls = ("A" if pz(Ey) > 0 else "B") if dom_y else ("A" if pz(Ez) < 0 else "B")
        I = np.abs(Ey) ** 2 + np.abs(Ez) ** 2
        yc = float(np.sum(I.sum(axis=1) * y) / I.sum())
        zc = float(np.sum(I.sum(axis=0) * z) / I.sum())
        out.append({"mode": k, "neff": ne, "TE_fraction": te, "z_class": cls,
                    "loss_dB_per_cm": loss / 100, "y_centroid_m": yc, "z_centroid_m": zc})
    return out


def run_lumerical():
    lumapi = load_lumapi()
    res = {}
    t0 = time.time()
    with lumapi.MODE(hide=True) as mode:
        # ---- materials ----
        m = mode.addmaterial("(n,k) Material")
        mode.setmaterial(m, "name", "core_glass")
        mode.setmaterial("core_glass", "Refractive Index", N_CORE)
        m = mode.addmaterial("Sellmeier")                 # eps = A0 + sum B lam^2/(lam^2 - C)
        mode.setmaterial(m, "name", "silica_malitson")
        mode.setmaterial("silica_malitson", "A0", 1.0)
        for i, (b, ci) in enumerate(zip(SELL_B, SELL_C), 1):
            mode.setmaterial("silica_malitson", f"B{i}", b)
            mode.setmaterial("silica_malitson", f"C{i}", ci**2)   # v261: C in um^2 (verified)
        m = mode.addmaterial("(n,k) Material")             # for the constant-index check
        mode.setmaterial(m, "name", "clad_const")
        mode.setmaterial("clad_const", "Refractive Index", 1.444)

        lam_chk = np.linspace(LAM_START, LAM_STOP, 41)
        n_lum = np.real(np.asarray(mode.getindex("silica_malitson", c / lam_chk))).flatten()
        res["material_check_max_diff"] = float(np.max(np.abs(n_lum - malitson(lam_chk))))
        print(f"silica index Lumerical vs Malitson (1500-1600 nm): max diff "
              f"{res['material_check_max_diff']:.1e}")

        # ---- geometry (as step 1) ----
        mode.addrect(name="clad", x=0, x_span=2e-6, y=0, y_span=10e-6, z=0, z_span=10e-6,
                     material="silica_malitson")
        mode.addrect(name="core", x=0, x_span=2e-6, y=0, y_span=W, z=0, z_span=H,
                     material="core_glass")
        mode.setnamed("core", "override mesh order from material database", 1)
        mode.setnamed("core", "mesh order", 1)
        mode.setnamed("clad", "override mesh order from material database", 1)
        mode.setnamed("clad", "mesh order", 3)
        mode.addfde(solver_type="2D X normal", x=0, y=0, y_span=SIM_SPAN, z=0, z_span=SIM_SPAN)
        mode.set("define y mesh by", "number of mesh cells")
        mode.set("define z mesh by", "number of mesh cells")
        mode.set("mesh cells y", CELLS)
        mode.set("mesh cells z", CELLS)

        def find_at(lam, ntrial=8):
            mode.setanalysis("wavelength", lam)
            mode.setanalysis("number of trial modes", ntrial)
            mode.setanalysis("search", "near n")          # v261 spelling (step 1)
            mode.setanalysis("use max index", 1)
            return mode_list(mode, ntrial)

        # ---- which modes are guided at 1500 nm ----
        mode.setanalysis("bent waveguide", 0)
        ms = find_at(LAM_START)
        n_clad = float(malitson(LAM_START))
        print(f"\nmodes at 1500 nm (n_clad = {n_clad:.5f}):")
        for q in ms:
            print(f"  mode{q['mode']}: n_eff {q['neff']:.5f}  TE {q['TE_fraction']:.3f}  "
                  f"class {q['z_class']}  {'guided' if q['neff'] > n_clad else 'not guided'}")
        guided = [q for q in ms if q["neff"] > n_clad]
        te0 = max((q for q in guided if q["TE_fraction"] > 0.5), key=lambda q: q["neff"])
        n_floor = 0.5 * (te0["neff"] + n_clad)
        tm0 = max((q for q in guided if q["TE_fraction"] < 0.5 and q["neff"] > n_floor),
                  key=lambda q: q["neff"])
        hom = [q for q in guided if q["neff"] < n_floor]
        families = [("TE0", te0), ("TM0", tm0)] + [(f"HOM{i + 1}", q) for i, q in enumerate(hom)]
        res["modes_1500nm"] = ms

        # ---- track each family over the sweep ----
        res["sweeps"] = {}
        for label, q in families:
            t = time.time()
            find_at(LAM_START)                               # same ordering as above
            mode.selectmode(q["mode"])
            mode.setanalysis("track selected mode", 1)
            mode.setanalysis("stop wavelength", LAM_STOP)  # sweep starts at "wavelength"
            mode.setanalysis("number of points", NPTS)
            mode.setanalysis("number of test modes", 8)
            mode.setanalysis("detailed dispersion calculation", 1 if label == "TE0" else 0)
            mode.frequencysweep()
            f = np.real(mode.getdata("frequencysweep", "f")).flatten()
            ne = np.real(mode.getdata("frequencysweep", "neff")).flatten()
            lam = c / f
            o = np.argsort(lam)
            sw = {"lam_m": lam[o].tolist(), "neff": ne[o].tolist(),
                  "TE_fraction_1500": q["TE_fraction"], "z_class": q["z_class"]}
            if label == "TE0":
                vg = np.real(mode.getdata("frequencysweep", "vg")).flatten()
                D = np.real(mode.getdata("frequencysweep", "D")).flatten()
                sw["ng_lumerical"] = (c / vg[o]).tolist()
                sw["D_lumerical_s_m2"] = D[o].tolist()
            res["sweeps"][label] = sw
            print(f"  {label}: swept {lam.min() * 1e9:.0f}-{lam.max() * 1e9:.0f} nm in {time.time() - t:.0f} s")

        # ---- check: constant-index solves with n_clad = Malitson(lam) ----
        mode.switchtolayout()
        mode.setnamed("clad", "material", "clad_const")
        res["const_index_check"] = []
        for lam in (1.50e-6, 1.55e-6, 1.60e-6):
            mode.switchtolayout()
            mode.setmaterial("clad_const", "Refractive Index", float(malitson(lam)))
            qq = find_at(lam, ntrial=4)
            n_te = max((x["neff"] for x in qq if x["TE_fraction"] > 0.5))
            res["const_index_check"].append({"lam_m": lam, "neff_const": n_te})
        mode.switchtolayout()
        mode.setnamed("clad", "material", "silica_malitson")

        # ---- check: bent waveguide at R (Malitson cladding, 1550 nm) ----
        mode.setanalysis("wavelength", LAM0)
        sweep_straight = find_at(LAM0, ntrial=10)
        bent = None
        for orient in (0, 90):
            mode.setanalysis("bent waveguide", 1)
            mode.setanalysis("bend radius", R)
            mode.setanalysis("bend orientation", orient)
            bent = find_at(LAM0, ntrial=10)
            b0 = max((x for x in bent if x["TE_fraction"] > 0.5), key=lambda x: x["neff"])
            # in-plane bend pushes the mode outward along y, not z
            if abs(b0["y_centroid_m"]) > abs(b0["z_centroid_m"]):
                break
        res["bend_check"] = {"orientation_deg": orient, "R_m": R,
                             "straight_1550": sweep_straight, "bent_1550": bent}
        mode.setanalysis("bent waveguide", 0)
    res["lumerical_time_s"] = time.time() - t0
    print(f"Lumerical total {res['lumerical_time_s']:.0f} s")
    return res


# ======================= ring resonances =======================
def resonance(psi, m, nu_guess, fsr):
    """nu solving psi(nu) L / c = m, with psi(nu) = nu * n_eff(c/nu) (beta = 2 pi psi / c)."""
    return brentq(lambda v: psi(v) * L_RING / c - m, nu_guess - 0.45 * fsr, nu_guess + 0.45 * fsr,
                  xtol=1e-3, rtol=1e-15)


def te0_combs(fit):
    """Pump resonance near 1550 nm and bins -10..10, for each D scaling.

    D scaling (ASSUMED definition): keep psi0 and the group index at the pump,
    scale everything beyond first order,
        psi_s(nu) = psi0 + ng_p (nu - nu_p) + s [psi(nu) - psi0 - ng_p (nu - nu_p)]
    so the pump resonance and FSR stay put and beta2 (and higher orders) scale by s.
    """
    psi = lambda v: v * fit["n"](c / v)
    m_p = int(np.round(fit["n"](LAM0) * L_RING / LAM0))
    fsr0 = c / (fit["ng"](LAM0) * L_RING)
    nu_p = resonance(psi, m_p, c / LAM0, fsr0)
    lam_p = c / nu_p
    ng_p = float(fit["ng"](lam_p))
    psi0 = psi(nu_p)
    combs = {}
    for s in D_SCALES:
        ps = lambda v, s=s: psi0 + ng_p * (v - nu_p) + s * (psi(v) - psi0 - ng_p * (v - nu_p))
        nus = {mu: resonance(ps, m_p + mu, nu_p + mu * fsr0, fsr0) for mu in range(-11, 12)}
        combs[s] = nus
    return {"m_p": m_p, "nu_p": nu_p, "lam_p": lam_p, "ng_p": ng_p, "combs": combs}


def bins_table(tc, D_p_lum):
    """Delta_k from the resonances vs analytic D2 k^2, and the bin brightness."""
    nus1 = tc["combs"][1]
    fsr = 0.5 * (nus1[1] - nus1[-1])                     # FSR at the pump
    rows = {}
    for s in D_SCALES:
        nus = tc["combs"][s]
        D2 = tc["lam_p"] ** 2 * (s * D_p_lum) * (2 * np.pi * fsr) ** 2 / (2 * np.pi * tc["ng_p"])
        r = []
        for k in K_BINS:
            dk = nus[k] + nus[-k] - 2 * nus[0]             # Hz
            da = D2 * k**2 / (2 * np.pi)                   # Hz
            r.append({"k": int(k), "Delta_numeric_Hz": dk, "Delta_analytic_Hz": da,
                      "brightness_numeric": 1 / (1 + (2 * dk / LINEWIDTH) ** 2),
                      "brightness_analytic": 1 / (1 + (2 * da / LINEWIDTH) ** 2)})
        rows[s] = {"D_ps_nm_km": s * D_p_lum * PS_NM_KM, "D2_rad_s": D2, "bins": r}
    return fsr, rows


def crossings(fits, tc):
    """TE0 bins -10..10 vs resonances of the other mode families in the same ring."""
    nus = tc["combs"][1]
    te0 = {mu: nus[mu] for mu in range(-10, 11)}
    nu_lo, nu_hi = te0[-10] - 250e9, te0[10] + 250e9
    out = {}
    for label, fit in fits.items():
        if label == "TE0":
            continue
        psi = lambda v, fit=fit: v * fit["n"](c / v)
        ng = float(fit["ng"](LAM0))
        fsr = c / (ng * L_RING)
        m_lo = int(np.ceil(psi(nu_lo) * L_RING / c))
        m_hi = int(np.floor(psi(nu_hi) * L_RING / c))
        res_h = []
        for m in range(m_lo, m_hi + 1):
            # first guess from the local linear model, then the exact root
            g = nu_lo + (m - psi(nu_lo) * L_RING / c) * fsr
            res_h.append(resonance(psi, m, g, fsr))
        res_h = np.array(res_h)
        offs = []
        for mu, v in te0.items():
            j = int(np.argmin(np.abs(res_h - v)))
            offs.append({"bin": mu, "offset_Hz": float(res_h[j] - v)})
        near = [o for o in offs if abs(o["offset_Hz"]) < LINEWIDTH]
        best = min(offs, key=lambda o: abs(o["offset_Hz"]))
        out[label] = {"neff_1550": float(fit["n"](LAM0)), "ng_1550": ng, "FSR_Hz": fsr,
                      "FSR_minus_TE0_Hz": fsr - (te0[1] - te0[-1]) / 2,
                      "offsets": offs, "within_linewidth": near, "closest": best,
                      "offset_shift_per_1e-5_neff_Hz": (c / LAM0) * 1e-5 / ng}
    return out


# ======================= plot =======================
def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ramp = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]   # blue 250..700, validated ordinal
    ink, ink2, muted, grid, surf = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
    fig, ax = plt.subplots(figsize=(7, 4.3), dpi=200)
    fig.patch.set_facecolor(surf)
    ax.set_facecolor(surf)
    for col, s in zip(ramp, D_SCALES):
        b = [r["brightness_numeric"] for r in rows[s]["bins"]]
        ax.plot(K_BINS, b, color=col, lw=2, marker="o", ms=6, mec=surf, mew=1.5,
                label=f"D × {s:g}  ({rows[s]['D_ps_nm_km']:.0f} ps/(nm·km))")
    ax.set_xlim(0.5, 10.5)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(K_BINS)
    ax.set_xlabel("bin pair k (resonances p ± k)", color=ink2)
    ax.set_ylabel("relative pair brightness", color=ink2)
    ax.set_title("Bin brightness vs k for scaled dispersion\n"
                 f"1 / (1 + (2Δ$_k$/Γ)$^2$), Γ = {LINEWIDTH / 1e6:.0f} MHz, Malitson cladding, core n = 1.70",
                 color=ink, fontsize=10, loc="left")
    ax.grid(True, color=grid, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(muted)
    ax.tick_params(colors=muted, labelcolor=ink2)
    leg = ax.legend(loc="lower left", frameon=False, fontsize=8.5)
    for t in leg.get_texts():
        t.set_color(ink)
    fig.tight_layout()
    fig.savefig(PNG, facecolor=surf)
    plt.close(fig)


# ======================= main =======================
def _js(o):
    """numpy -> plain JSON types."""
    if isinstance(o, dict):
        return {str(k): _js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_js(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


RAW = os.path.join(RESULTS, "s19_dispersion_raw.json")


def main(stage="all"):
    """stage "all": run Lumerical, save raw data, analyse; "analyze": reuse raw data."""
    if stage == "analyze":
        with open(RAW) as fh:
            res = json.load(fh)
    else:
        res = run_lumerical()
        with open(RAW, "w") as fh:
            json.dump(_js(res), fh, indent=2)
    fits = {k: poly_fit(np.array(v["lam_m"]), np.array(v["neff"])) for k, v in res["sweeps"].items()}
    te = res["sweeps"]["TE0"]
    lam = np.array(te["lam_m"])
    f0 = fits["TE0"]

    # ---------- TE0 dispersion: Malitson vs step 1 ----------
    ng_lum = float(np.interp(LAM0, lam, te["ng_lumerical"]))
    D_lum = float(np.interp(LAM0, lam, te["D_lumerical_s_m2"]))
    s1 = S1["at_1550nm"]
    print("\n--- TE0 at 1550 nm: constant indices (step 1) | Malitson cladding (this step) ---")
    print(f"n_eff                  {s1['neff']:.5f} | {float(f0['n'](LAM0)):.5f}")
    print(f"n_g  Lumerical / ours  {s1['ng_lumerical']:.5f} | {ng_lum:.5f} / {float(f0['ng'](LAM0)):.5f}")
    print(f"D    Lumerical / ours  {s1['D_lumerical_ps_nm_km']:.2f} | {D_lum * PS_NM_KM:.2f} / "
          f"{float(f0['D'](LAM0)) * PS_NM_KM:.2f}  ps/(nm km)")
    print(f"n_eff poly-fit residual {f0['resid']:.1e}")
    print(f"silica alone at 1550: D_mat = {-(LAM0 / c) * np.gradient(np.gradient(malitson(lam), lam), lam)[10] * PS_NM_KM:.1f} ps/(nm km) (bulk)")
    cc = [{"lam_m": r["lam_m"], "neff_const": r["neff_const"],
           "neff_sweep": float(f0["n"](r["lam_m"])),
           "diff": r["neff_const"] - float(f0["n"](r["lam_m"]))} for r in res["const_index_check"]]
    for r in cc:
        print(f"check n_eff at {r['lam_m'] * 1e9:.0f} nm: constant n_clad = Malitson(lam) "
              f"{r['neff_const']:.6f} | dispersive sweep {r['neff_sweep']:.6f}  diff {r['diff']:.1e}")
    # independent D: second difference of the three constant-index solves (50 nm steps)
    n3 = [r["neff_const"] for r in cc]
    D_3pt = -(LAM0 / c) * (n3[0] + n3[2] - 2 * n3[1]) / (50e-9) ** 2
    print(f"check D at 1550 nm: 3-point second difference of constant-index solves "
          f"{D_3pt * PS_NM_KM:.2f} ps/(nm km) (truncation error O(50 nm^2))")

    # ---------- resonances, Delta_k, brightness ----------
    tc = te0_combs(f0)
    D_p_lum = float(np.interp(tc["lam_p"], lam, te["D_lumerical_s_m2"]))
    fsr, rows = bins_table(tc, D_p_lum)
    print(f"\npump resonance m = {tc['m_p']}, lam_p = {tc['lam_p'] * 1e9:.4f} nm, "
          f"n_g = {tc['ng_p']:.5f}, FSR = {fsr / 1e9:.3f} GHz (L fixed at {L_RING * 1e3:.4f} mm; "
          f"200 GHz would need R = {c / (tc['ng_p'] * 200e9) / (2 * np.pi) * 1e6:.2f} um)")
    print(f"D at lam_p: Lumerical {D_p_lum * PS_NM_KM:.2f} ps/(nm km)")
    for s in D_SCALES:
        r = rows[s]
        print(f"\nD x {s:g} ({r['D_ps_nm_km']:.1f} ps/(nm km), D2/2pi = {r['D2_rad_s'] / (2 * np.pi) / 1e6:.4f} MHz)")
        print("  k   Delta numeric [MHz]  Delta analytic [MHz]  rel.diff   brightness num / analytic")
        for b in r["bins"]:
            rd = (b["Delta_numeric_Hz"] - b["Delta_analytic_Hz"]) / b["Delta_analytic_Hz"]
            print(f"  {b['k']:2d}  {b['Delta_numeric_Hz'] / 1e6:12.4f}  {b['Delta_analytic_Hz'] / 1e6:16.4f}"
                  f"      {rd:+.2e}    {b['brightness_numeric']:.4f} / {b['brightness_analytic']:.4f}")

    # ---------- avoided crossings ----------
    cr = crossings(fits, tc)
    print(f"\n--- other mode families vs TE0 bins -10..10 (window +-{LINEWIDTH / 1e6:.0f} MHz) ---")
    for label, x in cr.items():
        sw = res["sweeps"][label]
        print(f"{label}: n_eff {x['neff_1550']:.5f}, n_g {x['ng_1550']:.4f}, class {sw['z_class']}, "
              f"FSR {x['FSR_Hz'] / 1e9:.3f} GHz (TE0 - {-x['FSR_minus_TE0_Hz'] / 1e9:.3f} GHz), "
              f"fit resid {fits[label]['resid']:.1e}")
        print(f"   closest: bin {x['closest']['bin']:+d}, offset {x['closest']['offset_Hz'] / 1e9:+.3f} GHz; "
              f"within linewidth: {[(o['bin'], round(o['offset_Hz'] / 1e6)) for o in x['within_linewidth']] or 'none'}"
              f"; offset moves {x['offset_shift_per_1e-5_neff_Hz'] / 1e9:.2f} GHz per 1e-5 in n_eff")

    # ---------- bend check ----------
    bc = res["bend_check"]
    print(f"\n--- bent waveguide at R = {R * 1e6:.2f} um (orientation {bc['orientation_deg']} deg), 1550 nm ---")
    pairs = []
    for q in bc["straight_1550"]:
        if q["neff"] <= malitson(LAM0):
            continue
        cand = [b for b in bc["bent_1550"] if b["z_class"] == q["z_class"]
                and abs(b["TE_fraction"] - q["TE_fraction"]) < 0.35]
        if not cand:
            continue
        b = min(cand, key=lambda b: abs(b["neff"] - q["neff"]))
        dn = b["neff"] - q["neff"]
        pairs.append({"straight": q, "bent": b, "dneff": dn,
                      "offset_shift_Hz": -(c / LAM0) * dn / (fits["TE0"]["ng"](LAM0))})
        print(f"  straight n_eff {q['neff']:.5f} (TE {q['TE_fraction']:.2f}, class {q['z_class']}) -> bent "
              f"{b['neff']:.5f}  dn {dn:+.1e}  bend loss {b['loss_dB_per_cm']:.2e} dB/cm  "
              f"y-shift {b['y_centroid_m'] * 1e9:+.0f} nm")

    # sensitivity: shift each family's comb by its bend-induced n_eff change (taken
    # constant over the band), relative to TE0, and redo the crossing scan.
    # A resonance moves by -nu * dn / n_g when n_eff changes by dn.
    te0_pair = min(pairs, key=lambda p: abs(p["straight"]["neff"] - float(f0["n"](LAM0))))
    dnu = {}
    for label, x in cr.items():
        p = min(pairs, key=lambda p: abs(p["straight"]["neff"] - x["neff_1550"]))
        rel = -(c / LAM0) * (p["dneff"] / x["ng_1550"] - te0_pair["dneff"] / tc["ng_p"])
        fsr_h = x["FSR_Hz"]
        new = [{"bin": o["bin"],
                "offset_Hz": (o["offset_Hz"] + rel + 0.5 * fsr_h) % fsr_h - 0.5 * fsr_h}
               for o in x["offsets"]]
        dnu[label] = {"relative_shift_Hz": rel, "offsets": new,
                      "within_linewidth": [o for o in new if abs(o["offset_Hz"]) < LINEWIDTH],
                      "closest": min(new, key=lambda o: abs(o["offset_Hz"]))}
    print("\n--- crossing scan with bend-induced n_eff shifts applied (sensitivity) ---")
    for label, x in dnu.items():
        print(f"{label}: comb moves {x['relative_shift_Hz'] / 1e9:+.1f} GHz vs TE0; closest bin "
              f"{x['closest']['bin']:+d} at {x['closest']['offset_Hz'] / 1e9:+.3f} GHz; within linewidth: "
              f"{[(o['bin'], round(o['offset_Hz'] / 1e6)) for o in x['within_linewidth']] or 'none'}")

    plot(rows)
    print(f"\nsaved {PNG}")

    out = {
        "step": "s19_dispersion",
        "inputs": {"L_m": L_RING, "R_m": R, "linewidth_Hz": LINEWIDTH, "n_core": N_CORE,
                   "cladding": "Malitson Sellmeier (Lumerical native Sellmeier material)",
                   "D_scaling_definition": "psi_s = psi0 + ng_p dnu + s (psi - psi0 - ng_p dnu), psi = nu n_eff",
                   "brightness_formula": "1 / (1 + (2 Delta_k / linewidth)^2)"},
        "material_check_max_diff": res["material_check_max_diff"],
        "modes_1500nm": res["modes_1500nm"],
        "te0_1550": {"neff": float(f0["n"](LAM0)), "ng_lumerical": ng_lum,
                     "ng_ours": float(f0["ng"](LAM0)), "D_lumerical_ps_nm_km": D_lum * PS_NM_KM,
                     "D_ours_ps_nm_km": float(f0["D"](LAM0)) * PS_NM_KM, "fit_resid": f0["resid"],
                     "step1_constant_index": {k: s1[k] for k in ("neff", "ng_lumerical", "D_lumerical_ps_nm_km")}},
        "const_index_check": cc,
        "D_3pt_const_index_ps_nm_km": D_3pt * PS_NM_KM,
        "pump": {"m": tc["m_p"], "nu_Hz": tc["nu_p"], "lam_m": tc["lam_p"], "ng": tc["ng_p"],
                 "FSR_Hz": fsr, "D_lumerical_ps_nm_km": D_p_lum * PS_NM_KM},
        "bins_vs_D_scale": {str(s): rows[s] for s in D_SCALES},
        "te0_resonances_Hz": {str(s): {str(k): v for k, v in tc["combs"][s].items()} for s in D_SCALES},
        "crossings": cr,
        "crossings_bend_shifted": dnu,
        "bend_check": {"orientation_deg": bc["orientation_deg"], "pairs": pairs,
                       "bent_modes": bc["bent_1550"]},
        "sweeps": res["sweeps"],
        "neff_polys": {k: {"coeffs_highest_first": v["p"].tolist(), "resid": v["resid"],
                           "x_definition": "x = (lam - 1.55e-6) / 50e-9"} for k, v in fits.items()},
        "lumerical_time_s": res["lumerical_time_s"],
    }
    with open(OUT, "w") as fh:
        json.dump(_js(out), fh, indent=2)
    print(f"saved {OUT}")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
