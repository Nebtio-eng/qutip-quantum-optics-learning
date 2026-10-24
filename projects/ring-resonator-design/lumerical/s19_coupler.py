"""Session 19, Step 2: bus-ring coupler and Q budget.

Stage "fde" (default):
  1. MODE FDE on the bus + ring cross-section (two identical parallel cores) for
     gaps 0.4-1.2 um: even/odd TE supermodes -> kappa = pi*dn/lam per unit length.
     Run at two meshes (40 nm, 20 nm) as a convergence check.
  2. Straight-coupler approximation for a ring of radius R: integrate kappa(g(x))
     along the curved gap -> coupling phase Phi, power coupling kappa^2 = sin^2(Phi).
  3. Q budget with ASSUMED loss 0.06 dB/cm: Q_i, Q_c(g), Q_L(g), linewidth,
     eta_esc; solve for the gap giving Q_L = 235,000.
     Checks: Q_i by hand vs exact ring formula; Q_L (1/Q sum) vs exact all-pass
     ring formula vs Lorentzian fit of the computed transmission spectrum.
  4. Supermode list at the chosen gap: does bus TE0 mix with higher-order modes?
Stage "eme" (python s19_coupler.py eme):
  EME of the real curved coupler at the chosen gap (independent check of kappa^2,
  and direct |S|^2 from bus TE0 into each ring mode).
Results -> results/s19_coupler.json
"""
import os, sys, json, time, tempfile, importlib.util
import numpy as np
from scipy.constants import c
from scipy.optimize import brentq, curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")   # project results folder
OUT = os.path.join(RESULTS, "s19_coupler.json")
LUMAPI = r"C:\Program Files\Lumerical\v261\api\python\lumapi.py"


def load_lumapi():
    spec = importlib.util.spec_from_file_location("lumapi", LUMAPI)
    lumapi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lumapi)
    return lumapi


# ---------- inputs from step 1 ----------
with open(os.path.join(RESULTS, "s19_waveguide.json")) as fh:
    S1 = json.load(fh)
W, H = S1["design"]["W_m"], S1["design"]["H_m"]
N_CORE, N_CLAD = S1["design"]["n_core"], S1["design"]["n_clad"]
NG = S1["at_1550nm"]["ng_lumerical"]
L_RING = S1["at_1550nm"]["L_m"]
R = S1["at_1550nm"]["R_m"]                      # 136.93 um, ring centre-line radius
P_NEFF = np.array(S1["neff_poly"]["coeffs_highest_first"])
LAM0 = 1.55e-6
NU0 = c / LAM0

# ---------- step-2 parameters ----------
# ASSUMED: bus waveguide identical to the ring waveguide (same W, H, layer);
# gap = edge-to-edge distance in the plane of the chip.
GAPS = np.round(np.arange(0.4, 1.2001, 0.1), 2) * 1e-6
LOSS_DB_CM = 0.06                               # ASSUMED propagation loss
Q_TARGET = 235_000                              # Kues et al. 2017, loaded Q
MARGIN = 2.25e-6                                # cladding margin outside the cores (as step 1)
Z_SPAN = 6e-6
MESHES = (40e-9, 20e-9)                         # coarse, fine (fine = reported)


def neff_1(lam):
    """Single-waveguide n_eff(lam) from step 1 (polynomial fit of the FDE sweep)."""
    return np.polyval(P_NEFF, (lam - LAM0) / 50e-9)


# ======================= FDE: coupler cross-section =======================
def build_xsec(mode):
    for name, n in [("core_glass", N_CORE), ("clad_silica", N_CLAD)]:
        m = mode.addmaterial("(n,k) Material")
        mode.setmaterial(m, "name", name)
        mode.setmaterial(name, "Refractive Index", n)
    mode.addrect(name="clad", x=0, x_span=2e-6, y=0, y_span=30e-6, z=0, z_span=12e-6,
                 material="clad_silica")
    for nm in ("bus", "ringwg"):
        mode.addrect(name=nm, x=0, x_span=2e-6, y=0, y_span=W, z=0, z_span=H,
                     material="core_glass")
    for nm, order in (("bus", 1), ("ringwg", 1), ("clad", 3)):
        mode.setnamed(nm, "override mesh order from material database", 1)
        mode.setnamed(nm, "mesh order", order)
    mode.addfde(solver_type="2D X normal", x=0, y=0, z=0, z_span=Z_SPAN)
    mode.set("define y mesh by", "number of mesh cells")
    mode.set("define z mesh by", "number of mesh cells")


def parity(mode, name, comp):
    """Mirror symmetry of the dominant field component about y = 0 (the gap centre).

    +1 = even supermode, -1 = odd.  Computed as sum E(y)E(-y) / sum E(y)^2 after
    removing the global phase, so it does not depend on the sign convention.
    """
    E = np.squeeze(mode.getdata(name, comp))                   # (ny, nz)
    y = np.real(mode.getdata(name, "y")).flatten()
    k = np.unravel_index(np.argmax(np.abs(E)), E.shape)
    Er = np.real(E * np.conj(E[k]) / np.abs(E[k]))
    Em = np.array([np.interp(-y, y, Er[:, j]) for j in range(Er.shape[1])]).T
    return float(np.sum(Er * Em) / np.sum(Er * Er))


def solve_gap(mode, gap, step, ntrial=8):
    """All supermodes of the two-core cross-section at one gap (lam = 1550 nm)."""
    mode.switchtolayout()
    yc = (W + gap) / 2                                          # core centres at +-yc
    mode.setnamed("bus", "y", -yc)
    mode.setnamed("ringwg", "y", +yc)
    y_span = 2 * W + gap + 2 * MARGIN
    mode.select("FDE")
    mode.set("y span", y_span)
    mode.set("mesh cells y", int(round(y_span / step)))
    mode.set("mesh cells z", int(round(Z_SPAN / step)))
    mode.setanalysis("wavelength", LAM0)
    mode.setanalysis("number of trial modes", ntrial)
    mode.setanalysis("search", "near n")                      # v261 spelling (see step 1)
    mode.setanalysis("use max index", 1)
    n = int(np.real(mode.findmodes()))
    modes = []
    for k in range(1, n + 1):
        name = f"FDE::data::mode{k}"
        te = float(np.real(mode.getdata(name, "TE polarization fraction")))
        ne = float(np.real(mode.getdata(name, "neff")).flatten()[0])
        p = parity(mode, name, "Ey" if te > 0.5 else "Ez")
        modes.append({"mode": k, "neff": ne, "TE_fraction": te, "parity": p})
    return modes


def te_pair(modes, n_floor):
    """Even and odd fundamental TE supermodes (TE fraction > 0.5, n_eff above the
    higher-order-mode group, parity sign decides even/odd)."""
    fund = [m for m in modes if m["TE_fraction"] > 0.5 and m["neff"] > n_floor]
    ev = [m for m in fund if m["parity"] > 0.5]
    od = [m for m in fund if m["parity"] < -0.5]
    if not ev or not od:
        raise RuntimeError(f"TE even/odd pair not found: {modes}")
    return max(ev, key=lambda m: m["neff"]), max(od, key=lambda m: m["neff"])


# ======================= coupler and ring theory =======================
def kappa_model(gaps, kap):
    """kappa(g) for any gap: log-linear interpolation between FDE points, and the
    exponential fit kappa ~ exp(-gamma g) beyond the sweep (only reached in the
    tails of the curved coupler, where kappa is tiny)."""
    lk = np.log(kap)
    gamma, lnk0 = np.polyfit(gaps, lk, 1)
    gamma = -gamma
    resid = float(np.max(np.abs(lk - (lnk0 - gamma * gaps))))

    def kfun(g):
        g = np.asarray(g, float)
        inside = np.exp(np.interp(g, gaps, lk))
        above = kap[-1] * np.exp(-gamma * (g - gaps[-1]))
        below = kap[0] * np.exp(-gamma * (g - gaps[0]))
        return np.where(g > gaps[-1], above, np.where(g < gaps[0], below, inside))
    return kfun, gamma, resid


def coupling_phase(kfun, g0, R=R, xmax=50e-6, n=20001):
    """Straight-coupler approximation for a ring of radius R next to a straight bus.

    Centre-to-centre separation along the bus: d(x) = d0 + R - sqrt(R^2 - x^2),
    so the local gap is g(x) = g0 + R - sqrt(R^2 - x^2).  Each slice dx is treated
    as a straight parallel coupler with the local kappa(g(x)) (adiabatic, phase
    matched because bus and ring waveguides are identical; bend neglected):
        Phi = integral kappa(g(x)) dx,   kappa^2_power = sin^2(Phi)
    Effective interaction length: L_eff = Phi / kappa(g0).
    """
    x = np.linspace(-xmax, xmax, n)
    g = g0 + R - np.sqrt(R**2 - x**2)
    phi = float(np.trapezoid(kfun(g), x))
    return phi, phi / float(kfun(g0))


def q_budget(phi):
    """Ring Q values for coupling phase phi (power coupling sin^2 phi)."""
    k2 = np.sin(phi) ** 2
    alpha = LOSS_DB_CM * 100 * np.log(10) / 10                   # power loss, 1/m
    Qi = 2 * np.pi * NG / (LAM0 * alpha)                         # = omega * n_g / (alpha c)
    # coupling Q: fraction k2 of the energy leaves per round trip T = n_g L / c
    # (small-k2 form; exact all-pass formula printed as a check below)
    Qc = 2 * np.pi * NG * L_RING / (LAM0 * k2)
    QL = 1 / (1 / Qi + 1 / Qc)
    return {"kappa2": k2, "alpha_per_m": alpha, "Qi": Qi, "Qc": Qc, "QL": QL,
            "linewidth_Hz": NU0 / QL, "eta_esc": QL / Qc}


def q_exact(a, t):
    """Exact all-pass ring FWHM -> Q = pi n_g L sqrt(a t) / (lam (1 - a t)).
    a = round-trip field amplitude, t = self-coupling (field)."""
    return np.pi * NG * L_RING * np.sqrt(a * t) / (LAM0 * (1 - a * t))


def spectrum_fit(phi, alpha):
    """Transmission of the all-pass ring near 1550 nm with the full n_eff(lam)
    phase, then a Lorentzian fit of the dip -> Q_L,fit = nu0 / FWHM."""
    t, a = np.cos(phi), np.exp(-alpha * L_RING / 2)
    m = int(np.round(neff_1(LAM0) * L_RING / LAM0))
    lam_res = brentq(lambda l: neff_1(l) * L_RING - m * l, 1.545e-6, 1.555e-6)
    nu_res = c / lam_res
    fw = nu_res / Q_TARGET
    nu = np.linspace(nu_res - 15 * fw, nu_res + 15 * fw, 4001)
    th = 2 * np.pi * neff_1(c / nu) * L_RING * nu / c
    T = (a**2 - 2 * a * t * np.cos(th) + t**2) / (1 - 2 * a * t * np.cos(th) + (a * t) ** 2)

    def lor(v, v0, gam, A, B):
        return B - A / (1 + ((v - v0) / (gam / 2)) ** 2)
    p, _ = curve_fit(lor, nu - nu_res, T, p0=[0, fw, 1 - T.min(), 1])
    return {"lam_res_m": lam_res, "Q_fit": nu_res / abs(p[1]), "FWHM_fit_Hz": abs(p[1]),
            "T_min": float(T.min()), "T_min_analytic": float(((a - t) / (1 - a * t)) ** 2),
            "max_fit_error": float(np.max(np.abs(lor(nu - nu_res, *p) - T)))}


# ======================= stage: FDE + analysis =======================
def stage_fde():
    lumapi = load_lumapi()
    n_floor = 0.5 * (neff_1(LAM0) + N_CLAD)   # separates fundamental from higher-order modes
    sweep = {}
    t0 = time.time()
    with lumapi.MODE(hide=True) as mode:
        build_xsec(mode)
        for step in MESHES:
            tag = f"{step * 1e9:.0f}nm"
            sweep[tag] = []
            print(f"\n=== FDE gap sweep, mesh {tag} ===")
            for g in GAPS:
                ms = solve_gap(mode, g, step)
                ev, od = te_pair(ms, n_floor)
                dn = ev["neff"] - od["neff"]
                sweep[tag].append({"gap_m": g, "n_even": ev["neff"], "n_odd": od["neff"],
                                   "dn": dn, "kappa_per_m": np.pi * dn / LAM0,
                                   "parity_even": ev["parity"], "parity_odd": od["parity"],
                                   "TE_even": ev["TE_fraction"], "TE_odd": od["TE_fraction"]})
                print(f"  gap {g * 1e6:.1f} um: n_even {ev['neff']:.6f}  n_odd {od['neff']:.6f}  "
                      f"dn {dn:.3e}  kappa {np.pi * dn / LAM0 * 1e-6:.4e} /um  "
                      f"(parity {ev['parity']:+.2f}/{od['parity']:+.2f})")
        print(f"FDE sweep total {time.time() - t0:.0f} s")

        # ---------- coupler model + Q budget per mesh ----------
        res = {}
        for tag in sweep:
            kap = np.array([r["kappa_per_m"] for r in sweep[tag]])
            kfun, gamma, resid = kappa_model(GAPS, kap)
            rows = []
            for g in GAPS:
                phi, leff = coupling_phase(kfun, g)
                rows.append({"gap_m": g, "phi": phi, "L_eff_m": leff, **q_budget(phi)})
            g_star = brentq(lambda g: q_budget(coupling_phase(kfun, g)[0])["QL"] - Q_TARGET,
                            GAPS[0], GAPS[-1])
            res[tag] = {"kfun": kfun, "gamma": gamma, "fit_resid": resid, "rows": rows,
                        "g_star": g_star}

        fine, coarse = res["20nm"], res["40nm"]
        g_star = fine["g_star"]
        phi_s, leff_s = coupling_phase(fine["kfun"], g_star)
        qb = q_budget(phi_s)
        a = np.exp(-qb["alpha_per_m"] * L_RING / 2)
        t = np.cos(phi_s)
        Qi_exact = q_exact(a, 1.0)
        QL_exact = q_exact(a, t)
        Qc_exact = q_exact(1.0, t)
        spec = spectrum_fit(phi_s, qb["alpha_per_m"])
        leff_analytic = np.sqrt(2 * np.pi * R / fine["gamma"])
        # same target with the exact all-pass Q instead of the small-kappa 1/Q sum
        g_exact = brentq(lambda g: q_exact(a, np.cos(coupling_phase(fine["kfun"], g)[0]))
                         - Q_TARGET, GAPS[0], GAPS[-1])

        # ---------- supermodes at the chosen gap (multimode check) ----------
        print(f"\n=== supermodes at chosen gap {g_star * 1e6:.3f} um (20 nm mesh, 16 trial modes) ===")
        sm = solve_gap(mode, g_star, MESHES[1], ntrial=16)
        for m in sm:
            grp = ("fundamental" if m["neff"] > n_floor else
                   "higher-order" if m["neff"] > N_CLAD else "box mode (n_eff < n_clad)")
            print(f"  mode{m['mode']:2d}: n_eff {m['neff']:.6f}  TE {m['TE_fraction']:.3f}  "
                  f"parity {m['parity']:+.2f}  ({grp})")
        ev, od = te_pair(sm, n_floor)
        kap_s = np.pi * (ev["neff"] - od["neff"]) / LAM0

    # ---------- phase-mismatch estimate for bus TE0 -> ring higher-order mode ----------
    # Coupled-mode theory along the curved gap: amplitude transferred to a mode with
    # propagation-constant mismatch dbeta is
    #   kappa_x(g0) * integral exp(-gamma x^2 / 2R) exp(i dbeta x) dx
    #   = kappa_x L_eff * exp(-dbeta^2 R / (2 gamma))
    # ASSUMED cross-coupling kappa_x ~ sqrt(kappa_00 * kappa_hh), with kappa_hh from
    # the even/odd splitting of each higher-order pair (order-of-magnitude proxy).
    # n_eff < n_clad: not guided, these are modes of the metal-walled simulation box
    hom = [m for m in sm if N_CLAD < m["neff"] < n_floor]
    box = [m["mode"] for m in sm if m["neff"] <= N_CLAD]
    pairs = []
    used = set()
    for m in sorted(hom, key=lambda m: -m["neff"]):
        if m["mode"] in used:
            continue
        partner = min((q for q in hom if q["mode"] != m["mode"] and q["mode"] not in used
                       and np.sign(q["parity"]) != np.sign(m["parity"])
                       and abs(q["TE_fraction"] - m["TE_fraction"]) < 0.3),
                      key=lambda q: abs(q["neff"] - m["neff"]), default=None)
        if partner is None:
            continue
        used |= {m["mode"], partner["mode"]}
        n_h = 0.5 * (m["neff"] + partner["neff"])
        k_hh = np.pi * abs(m["neff"] - partner["neff"]) / LAM0
        k_x = np.sqrt(kap_s * k_hh)
        dbeta = 2 * np.pi * (neff_1(LAM0) - n_h) / LAM0
        amp = k_x * leff_s * np.exp(-dbeta**2 * R / (2 * fine["gamma"]))
        pairs.append({"modes": [m["mode"], partner["mode"]], "n_eff_pair": n_h,
                      "TE_fraction": 0.5 * (m["TE_fraction"] + partner["TE_fraction"]),
                      "kappa_hh_per_m": k_hh, "kappa_x_per_m": k_x,
                      "dbeta_per_m": dbeta, "power_transfer_est": float(min(amp, 1) ** 2)})

    # ---------- report ----------
    print("\n--- coupling and Q budget vs gap (20 nm mesh) ---")
    print(" gap[um]   dn         kappa[1/um]  L_eff[um]  kappa^2     Q_c        Q_L      dnu[MHz]  eta_esc")
    for r, s in zip(fine["rows"], sweep["20nm"]):
        print(f"  {r['gap_m'] * 1e6:.1f}   {s['dn']:.3e}   {s['kappa_per_m'] * 1e-6:.4e}   "
              f"{r['L_eff_m'] * 1e6:6.2f}   {r['kappa2']:.3e}  {r['Qc']:9.0f}  {r['QL']:9.0f}  "
              f"{r['linewidth_Hz'] * 1e-6:7.1f}   {r['eta_esc']:.4f}")
    print(f"\nexp fit: kappa ~ exp(-gamma g), gamma = {fine['gamma'] * 1e-6:.3f} /um "
          f"(max |ln residual| {fine['fit_resid']:.3f})")
    print(f"L_eff at g*: numeric integral {leff_s * 1e6:.2f} um | analytic sqrt(2 pi R / gamma) "
          f"{leff_analytic * 1e6:.2f} um")

    print(f"\n--- chosen gap for Q_L = {Q_TARGET:,} ---")
    print(f"g* = {g_star * 1e6:.3f} um (20 nm mesh) | {coarse['g_star'] * 1e6:.3f} um (40 nm mesh) | "
          f"{g_exact * 1e6:.3f} um if exact ring Q is used")
    print(f"kappa(g*) interpolated {fine['kfun'](g_star) * 1e-6:.4e} /um | direct FDE at g* "
          f"{kap_s * 1e-6:.4e} /um")
    print(f"Phi = {phi_s:.4f} rad, kappa^2 = {qb['kappa2']:.4f}, t = {t:.5f}, a = {a:.6f}")
    print(f"Q_i  hand 2 pi n_g/(lam alpha) = {qb['Qi']:,.0f} | exact ring formula = {Qi_exact:,.0f}")
    print(f"Q_c  small-kappa = {qb['Qc']:,.0f} | exact ring formula = {Qc_exact:,.0f}")
    print(f"Q_L  1/Q sum = {qb['QL']:,.0f} | exact ring formula = {QL_exact:,.0f} | "
          f"Lorentzian fit of spectrum = {spec['Q_fit']:,.0f}")
    print(f"linewidth = {qb['linewidth_Hz'] * 1e-6:.1f} MHz, eta_esc = {qb['eta_esc']:.4f}, "
          f"T_min = {spec['T_min']:.4f} (analytic {spec['T_min_analytic']:.4f})")

    print("\n--- mesh convergence of kappa (40 -> 20 nm) ---")
    rel = []
    for s40, s20 in zip(sweep["40nm"], sweep["20nm"]):
        rc = (s20["kappa_per_m"] - s40["kappa_per_m"]) / s20["kappa_per_m"]
        rel.append(rc)
        print(f"  gap {s20['gap_m'] * 1e6:.1f} um: kappa change {rc * 100:+.2f} %")

    print("\n--- multimode: bus TE0 -> ring higher-order modes (phase-mismatch estimate) ---")
    for p in pairs:
        print(f"  pair {p['modes']}: n_eff {p['n_eff_pair']:.4f} (TE {p['TE_fraction']:.2f}), "
              f"dbeta {p['dbeta_per_m'] * 1e-6:.3f} /um, est. power transfer {p['power_transfer_est']:.1e}")
    print(f"  excluded box modes (n_eff < n_clad): {box}")

    out = {
        "step": "s19_coupler",
        "inputs": {"from_step1": {"ng": NG, "L_m": L_RING, "R_m": R, "neff_1550": float(neff_1(LAM0))},
                   "loss_dB_per_cm": LOSS_DB_CM, "Q_target": Q_TARGET,
                   "bus_equals_ring_waveguide": True, "gaps_m": GAPS.tolist(),
                   "meshes_m": list(MESHES)},
        "fde_sweep": sweep,
        "coupler_model": {
            "method": "straight-coupler approximation: Phi = integral kappa(g0 + R - sqrt(R^2-x^2)) dx "
                      "over |x| <= 50 um; kappa(g) log-linear interpolated, exponential beyond sweep; "
                      "L_eff = Phi / kappa(g0); kappa^2 = sin^2(Phi)",
            "gamma_per_m": fine["gamma"], "fit_max_ln_residual": fine["fit_resid"],
            "L_eff_analytic_m": leff_analytic},
        "q_vs_gap": [{k: v for k, v in r.items()} for r in fine["rows"]],
        "chosen": {"gap_m": g_star, "gap_m_coarse_mesh": coarse["g_star"],
                   "gap_m_exact_ring_Q": g_exact, "phi": phi_s,
                   "L_eff_m": leff_s, "kappa_per_m_interp": float(fine["kfun"](g_star)),
                   "kappa_per_m_direct": kap_s, "t": t, "a": a, **qb,
                   "Qi_exact": Qi_exact, "Qc_exact": Qc_exact, "QL_exact": QL_exact,
                   "spectrum": spec},
        "mesh_convergence_kappa_rel": rel,
        "supermodes_at_gap": sm,
        "higher_order_estimate": pairs,
        "box_modes_excluded": box,
    }
    os.makedirs(RESULTS, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nsaved {OUT}")


# ======================= stage: EME check =======================
def port_modes(mode, port, yb, yr):
    """n_eff, TE fraction and location (bus / ring) of each selected port mode.

    Location = power-weighted centroid of |E|^2 in y, compared with the bus centre
    (yb) and the ring centre line at the port (yr)."""
    ne = mode.getresult(port, "neff")
    ne = np.real(np.asarray(ne["neff"] if isinstance(ne, dict) else ne)).flatten()
    prof = mode.getresult(port, "mode profiles")
    y = np.real(np.asarray(prof["y"])).flatten()
    out = []
    for j in range(ne.size):
        E = np.asarray(prof[f"E{j + 1}"])               # (x, y, z, f, 3 components)
        I = np.abs(E) ** 2
        Iy, Iz = np.sum(I[..., 1]), np.sum(I[..., 2])
        w = np.sum(I, axis=tuple(k for k in range(I.ndim) if k != 1)).flatten()
        yc = float(np.sum(w * y) / np.sum(w))
        where = "bus" if abs(yc - yb) < abs(yc - yr) else "ring"
        out.append({"mode": j + 1, "neff": float(ne[j]), "TE_fraction": float(Iy / (Iy + Iz)),
                    "y_centroid_m": yc, "location": where})
    return out


def stage_eme(x_half=20e-6, cells=40, nmodes=12, step=40e-9, n_port=6, save=True):
    """EME of the actual curved coupler at the chosen gap.

    Bus: straight core along x at y = 0.  Ring: ring object, centre-line radius R,
    centre at y = W + g + R, only the arc near the bus.  EME spans |x| <= 20 um
    (the gap has opened by 1.5 um there; the model above puts 99.6 % of Phi inside).
    Full-span ports: at x = +-20 um the cross-section modes are already localised
    in the bus or in the ring (the tilted ring's n_eff differs from the bus by ~8e-4,
    far more than the residual coupling), so each port mode is labelled bus/ring by
    its field centroid.  Inject bus TE0 at port 1:
      |S(ring TE0)|^2 = kappa^2,  |S(bus TE0)|^2 = t^2,  and their sum checks energy.
    (Restricted-span ports were tried first: in v261 their modes came out of the
    x = 0 cross-section, not the port plane, so they are not used.)
    """
    with open(OUT) as fh:
        out = json.load(fh)
    g = out["chosen"]["gap_m"]
    yc = W + g + R
    y_ring = lambda x: yc - np.sqrt(R**2 - x**2)                 # ring centre line
    y_min = -W / 2 - MARGIN
    y_max = y_ring(x_half) + W / 2 + MARGIN

    lumapi = load_lumapi()
    t0 = time.time()
    with lumapi.MODE(hide=True) as mode:
        for name, n in [("core_glass", N_CORE), ("clad_silica", N_CLAD)]:
            m = mode.addmaterial("(n,k) Material")
            mode.setmaterial(m, "name", name)
            mode.setmaterial(name, "Refractive Index", n)
        mode.addrect(name="clad", x=0, x_span=4 * x_half, y=0.5 * (y_min + y_max),
                     y_span=(y_max - y_min) + 10e-6, z=0, z_span=12e-6, material="clad_silica")
        mode.addrect(name="bus", x=0, x_span=4 * x_half, y=0, y_span=W, z=0, z_span=H,
                     material="core_glass")
        mode.addring(name="ring", x=0, y=yc, z=0, z_span=H, inner_radius=R - W / 2,
                     outer_radius=R + W / 2, theta_start=240, theta_stop=300,
                     material="core_glass")
        for nm, order in (("bus", 1), ("ring", 1), ("clad", 3)):
            mode.setnamed(nm, "override mesh order from material database", 1)
            mode.setnamed(nm, "mesh order", order)

        mode.addeme()
        mode.set("x min", -x_half)        # v261: x max is inactive, set by x min + group spans
        mode.set("y min", y_min)
        mode.set("y max", y_max)
        mode.set("z", 0)
        mode.set("z span", Z_SPAN)
        mode.set("wavelength", LAM0)
        mode.set("define y mesh by", "maximum mesh step")
        mode.set("define z mesh by", "maximum mesh step")
        mode.set("dy", step)
        mode.set("dz", step)
        mode.set("number of cell groups", 1)
        mode.set("group spans", np.array([[2 * x_half]]))
        mode.set("cells", np.array([[cells]]))
        mode.set("subcell method", np.array([[1]]))              # 1 = CVCS
        mode.set("number of modes for all cell groups", nmodes)

        p1, p2 = "EME::Ports::port_1", "EME::Ports::port_2"
        for p in (p1, p2):
            mode.setnamed(p, "mode selection", "user select")
            mode.select(p)
            mode.updateportmodes(np.arange(1, n_port + 1))      # the property alone is not enough

        # v261: EME run() fails ("no simulation data") unless the project is saved first.
        # The file holds all cell modes (~0.8 GB), so it goes to the temp folder, not results/.
        mode.save(os.path.join(tempfile.gettempdir(), "s19_coupler_eme.lms"))
        mode.run()
        mode.emepropagate()
        runtime = time.time() - t0
        print(f"EME ({cells} cells, {nmodes} modes/cell, {step * 1e9:.0f} nm mesh) done in {runtime:.0f} s")
        S = np.asarray(mode.getresult("EME", "user s matrix"))
        m1 = port_modes(mode, p1, 0.0, y_ring(-x_half))
        m2 = port_modes(mode, p2, 0.0, y_ring(+x_half))

    for p, ms in (("port 1", m1), ("port 2", m2)):
        for m in ms:
            print(f"  {p} mode {m['mode']}: n_eff {m['neff']:.5f}  TE {m['TE_fraction']:.2f}  "
                  f"centroid y {m['y_centroid_m'] * 1e6:6.2f} um -> {m['location']}")

    def pick(ms, where):
        c = [m for m in ms if m["location"] == where and m["TE_fraction"] > 0.5]
        return max(c, key=lambda m: m["neff"]) if c else None
    src, bus_out, ring_out = pick(m1, "bus"), pick(m2, "bus"), pick(m2, "ring")
    if not (src and bus_out and ring_out):
        raise RuntimeError("could not identify bus/ring TE0 port modes")
    n1 = len(m1)
    col = src["mode"] - 1                                        # input: bus TE0 at port 1
    P = np.abs(S[n1:, col]) ** 2                                 # power into each port-2 mode
    for m, pw in zip(m2, P):
        m["power_from_bus_TE0"] = float(pw)
        print(f"  bus TE0 -> port 2 mode {m['mode']} ({m['location']}, TE {m['TE_fraction']:.2f}): "
              f"{pw:.4e}")
    k2 = float(P[ring_out["mode"] - 1])
    t2 = float(P[bus_out["mode"] - 1])
    refl = float(np.sum(np.abs(S[:n1, col]) ** 2))
    k2_model = out["chosen"]["kappa2"]
    print(f"kappa^2: FDE + straight-coupler model {k2_model:.4f} | EME {k2:.4f} "
          f"(ratio {k2 / k2_model:.3f})")
    print(f"energy: kappa^2 + t^2 = {k2 + t2:.5f}, all port-2 modes {P.sum():.5f}, reflection {refl:.1e}")
    res = {"gap_m": g, "x_span_m": 2 * x_half, "cells": cells, "subcell": "CVCS",
           "modes_per_cell": nmodes, "mesh_step_m": step, "runtime_s": runtime,
           "kappa2_eme": k2, "t2_eme": t2, "sum_port2": float(P.sum()), "reflection": refl,
           "kappa2_model": k2_model, "port1_modes": m1, "port2_modes": m2}
    if save:
        out.setdefault("eme_check", []).append(res)
        with open(OUT, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"saved {OUT}")
    return res


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "fde"
    {"fde": stage_fde, "eme": stage_eme}[stage]()
