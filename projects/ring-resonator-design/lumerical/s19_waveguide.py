"""Session 19, Step 1: waveguide cross-section of the Kues-type microring.

MODE FDE on the straight-waveguide cross-section, fundamental TE mode,
wavelength sweep 1500-1600 nm.  Outputs n_eff, n_g and D at 1550 nm and the
ring length / radius that give FSR = 200 GHz.

Checks (printed side by side):
  n_g : Lumerical c/v_g  vs  n_eff - lam*dn_eff/dlam  (our derivative)
  D   : Lumerical D      vs  -(lam/c)*d2n_eff/dlam2   (our derivative)
  mesh: rerun at a 2x finer mesh, report change in n_eff, n_g, D, R
  FSR : resonance condition n_eff(lam)*L = m*lam solved directly -> FSR
Results -> results/s19_waveguide.json
"""
import os, json, time, importlib.util
import numpy as np
from scipy.constants import c
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")   # project results folder
os.makedirs(RESULTS, exist_ok=True)

LUMAPI = r"C:\Program Files\Lumerical\v261\api\python\lumapi.py"


def load_lumapi():
    spec = importlib.util.spec_from_file_location("lumapi", LUMAPI)
    lumapi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lumapi)
    return lumapi


# ---------- design ----------
# ASSUMED geometry and indices (earlier papers on this Hydex chip family, not Kues 2017)
W, H = 1.50e-6, 1.45e-6          # ASSUMED core width, height
N_CORE, N_CLAD = 1.70, 1.444     # ASSUMED constant indices: material dispersion ignored,
                                 # so D below is the waveguide (geometric) dispersion only
LAM0 = 1.55e-6                   # pump wavelength
LAM_START, LAM_STOP, NPTS = 1.50e-6, 1.60e-6, 21
FSR_TARGET = 200e9               # Kues et al. 2017

# numerical settings (not physics)
SIM_SPAN = 6e-6                  # ASSUMED FDE window 6x6 um, default metal boundaries
MESH_COARSE, MESH_FINE = 150, 300   # cells per side: 40 nm and 20 nm steps


def build(mode):
    """Rectangular core buried in silica; x is the propagation direction."""
    for name, n in [("core_glass", N_CORE), ("clad_silica", N_CLAD)]:
        m = mode.addmaterial("(n,k) Material")
        mode.setmaterial(m, "name", name)
        mode.setmaterial(name, "Refractive Index", n)

    mode.addrect(name="clad", x=0, x_span=2e-6, y=0, y_span=10e-6, z=0, z_span=10e-6,
                 material="clad_silica")
    mode.addrect(name="core", x=0, x_span=2e-6, y=0, y_span=W, z=0, z_span=H,
                 material="core_glass")
    # core must win where it overlaps the cladding block
    mode.setnamed("core", "override mesh order from material database", 1)
    mode.setnamed("core", "mesh order", 1)
    mode.setnamed("clad", "override mesh order from material database", 1)
    mode.setnamed("clad", "mesh order", 3)

    mode.addfde(solver_type="2D X normal", x=0, y=0, y_span=SIM_SPAN, z=0, z_span=SIM_SPAN)
    mode.set("define y mesh by", "number of mesh cells")
    mode.set("define z mesh by", "number of mesh cells")


def solve(mode, cells):
    """Find the fundamental TE mode at 1500 nm, then track it across the sweep."""
    mode.switchtolayout()
    mode.select("FDE")
    mode.set("mesh cells y", cells)
    mode.set("mesh cells z", cells)

    mode.setanalysis("wavelength", LAM_START)
    mode.setanalysis("number of trial modes", 4)
    # v261 spelling: "search" takes 'near n' / 'in range'; searching near the
    # structure's max index is the separate flag "use max index"
    mode.setanalysis("search", "near n")
    mode.setanalysis("use max index", 1)
    nmodes = int(np.real(mode.findmodes()))

    # fundamental TE: TE fraction > 0.5 and largest n_eff.  The core is almost
    # square (1.50 x 1.45 um), so TE0 and TM0 are nearly degenerate: select by
    # polarisation, not by n_eff alone.
    best, best_n, best_te, modes = None, 0.0, 0.0, []
    for k in range(1, nmodes + 1):
        name = f"FDE::data::mode{k}"
        te = float(np.real(mode.getdata(name, "TE polarization fraction")))
        ne = float(np.real(mode.getdata(name, "neff")).flatten()[0])
        modes.append({"mode": k, "neff": ne, "TE_fraction": te})
        print(f"  mode{k}: n_eff = {ne:.5f}, TE fraction = {te:.3f}")
        if te > 0.5 and ne > best_n:
            best, best_n, best_te = f"mode{k}", ne, te
    if best is None:
        raise RuntimeError("no TE-like mode found")
    print(f"  -> using {best} (TE fraction {best_te:.3f})")

    mode.selectmode(best)
    mode.setanalysis("track selected mode", 1)
    # v261: "start wavelength" is inactive; the sweep starts at "wavelength" (set above)
    mode.setanalysis("stop wavelength", LAM_STOP)
    mode.setanalysis("number of points", NPTS)
    mode.setanalysis("number of test modes", 4)
    mode.setanalysis("detailed dispersion calculation", 1)
    mode.frequencysweep()

    f = np.real(mode.getdata("frequencysweep", "f")).flatten()
    neff = np.real(mode.getdata("frequencysweep", "neff")).flatten()
    vg = np.real(mode.getdata("frequencysweep", "vg")).flatten()
    D = np.real(mode.getdata("frequencysweep", "D")).flatten()   # SI: s/m^2

    lam = c / f
    o = np.argsort(lam)
    print(f"  sweep covered {lam.min() * 1e9:.1f}-{lam.max() * 1e9:.1f} nm, {lam.size} points")
    return {"lam": lam[o], "neff": neff[o], "ng_lum": c / vg[o], "D_lum": D[o],
            "TE_fraction": best_te, "modes_at_1500nm": modes}


def our_derivatives(lam, neff):
    """Independent n_g and D from the n_eff(lam) curve alone.

    n_g = n_eff - lam dn/dlam      (group index)
    D   = -(lam/c) d2n/dlam2       (GVD parameter, s/m^2)
    Two numerical routes: finite differences on the solver grid, and a smooth
    polynomial fit (degree 6, centred and scaled for conditioning).
    """
    dn_fd = np.gradient(neff, lam)
    d2n_fd = np.gradient(dn_fd, lam)

    s = 50e-9
    x = (lam - LAM0) / s
    p = np.polyfit(x, neff, 6)
    fit_resid = float(np.max(np.abs(np.polyval(p, x) - neff)))  # mode hop would show here
    dn_p = np.polyval(np.polyder(p, 1), x) / s
    d2n_p = np.polyval(np.polyder(p, 2), x) / s**2

    return {"ng_fd": neff - lam * dn_fd, "D_fd": -(lam / c) * d2n_fd,
            "ng_fit": neff - lam * dn_p, "D_fit": -(lam / c) * d2n_p,
            "poly": p, "scale": s, "fit_resid": fit_resid}


def fsr_from_resonances(p, s, L):
    """FSR from the resonance condition n_eff(lam) L = m lam (no n_g used).

    Solve for the two resonances bracketing 1550 nm and take the frequency gap.
    """
    neff = lambda l: np.polyval(p, (l - LAM0) / s)
    m0 = int(np.floor(neff(LAM0) * L / LAM0))          # resonance just above 1550 nm
    lams = [brentq(lambda l: neff(l) * L - m * l, LAM_START, LAM_STOP) for m in (m0, m0 + 1)]
    return c / lams[1] - c / lams[0], lams


def at(lam, arr):
    return float(np.interp(LAM0, lam, arr))


PS_NM_KM = 1e6   # 1 s/m^2 = 1e6 ps/(nm km)


def summarise(r):
    lam = r["lam"]
    d = our_derivatives(lam, r["neff"])
    ng = at(lam, r["ng_lum"])
    # Ring length from FSR = c / (n_g L).  ASSUMED: straight-waveguide n_g holds in
    # the ring (bend correction at R ~ 135 um neglected).
    L = c / (ng * FSR_TARGET)
    fsr_chk, lam_res = fsr_from_resonances(d["poly"], d["scale"], L)
    monotonic = bool(np.all(np.diff(r["neff"]) < 0))   # n_eff must fall with lam
    return {
        "neff": at(lam, r["neff"]),
        "ng_lumerical": ng,
        "ng_ours_fit": at(lam, d["ng_fit"]),
        "ng_ours_fd": at(lam, d["ng_fd"]),
        "D_lumerical_ps_nm_km": at(lam, r["D_lum"]) * PS_NM_KM,
        "D_ours_fit_ps_nm_km": at(lam, d["D_fit"]) * PS_NM_KM,
        "D_ours_fd_ps_nm_km": at(lam, d["D_fd"]) * PS_NM_KM,
        "L_m": L, "R_m": L / (2 * np.pi),
        "FSR_check_Hz": fsr_chk, "resonances_check_m": lam_res,
        "TE_fraction": r["TE_fraction"], "neff_fit_max_residual": d["fit_resid"],
        "neff_monotonic": monotonic,
    }, d


def main():
    lumapi = load_lumapi()
    runs = {}
    t0 = time.time()
    with lumapi.MODE(hide=True) as mode:
        build(mode)
        for cells in (MESH_COARSE, MESH_FINE):
            print(f"\n=== mesh {cells} x {cells} ({SIM_SPAN / cells * 1e9:.0f} nm step) ===")
            t = time.time()
            runs[cells] = solve(mode, cells)
            print(f"  sweep done in {time.time() - t:.0f} s")
    print(f"\nLumerical total {time.time() - t0:.0f} s")

    s_c, _ = summarise(runs[MESH_COARSE])
    s_f, d_f = summarise(runs[MESH_FINE])   # fine mesh = reported values

    # ---------- side-by-side report ----------
    print(f"\n--- results at 1550 nm (fine mesh {MESH_FINE}, reported) ---")
    print(f"n_eff                         = {s_f['neff']:.5f}")
    print(f"n_g   Lumerical / fit / FD    = {s_f['ng_lumerical']:.5f} / "
          f"{s_f['ng_ours_fit']:.5f} / {s_f['ng_ours_fd']:.5f}")
    print(f"D     Lumerical / fit / FD    = {s_f['D_lumerical_ps_nm_km']:.2f} / "
          f"{s_f['D_ours_fit_ps_nm_km']:.2f} / {s_f['D_ours_fd_ps_nm_km']:.2f}  ps/(nm km)")
    print(f"ring length L                 = {s_f['L_m'] * 1e3:.4f} mm")
    print(f"ring radius R                 = {s_f['R_m'] * 1e6:.2f} um")
    print(f"FSR target / from resonances  = {FSR_TARGET / 1e9:.2f} / {s_f['FSR_check_Hz'] / 1e9:.3f} GHz")
    print(f"TE fraction {s_f['TE_fraction']:.3f}, n_eff monotonic {s_f['neff_monotonic']}, "
          f"poly-fit residual {s_f['neff_fit_max_residual']:.1e}")

    print(f"\n--- mesh convergence ({MESH_COARSE} -> {MESH_FINE} cells) ---")
    conv = {}
    for key, label, scale, fmt in [("neff", "n_eff", 1, ".2e"),
                                   ("ng_lumerical", "n_g", 1, ".2e"),
                                   ("D_lumerical_ps_nm_km", "D [ps/(nm km)]", 1, ".3f"),
                                   ("R_m", "R [um]", 1e6, ".3f")]:
        a, b = s_c[key] * scale, s_f[key] * scale
        conv[key] = {"coarse": s_c[key], "fine": s_f[key], "abs_change": b - a,
                     "rel_change": (b - a) / b}
        print(f"{label:16s} {a:.6f} -> {b:.6f}   change {b - a:{fmt}} ({(b - a) / b:.1e} rel)")

    rf = runs[MESH_FINE]
    out = {
        "step": "s19_waveguide",
        "design": {"W_m": W, "H_m": H, "n_core": N_CORE, "n_clad": N_CLAD,
                   "material_dispersion": False, "FSR_target_Hz": FSR_TARGET, "lam0_m": LAM0},
        "numerics": {"solver": "MODE FDE 2D X normal, metal boundaries",
                     "window_m": SIM_SPAN, "mesh_cells": [MESH_COARSE, MESH_FINE],
                     "reported_mesh": MESH_FINE, "sweep_points": NPTS},
        "at_1550nm": s_f,
        "coarse_mesh_at_1550nm": s_c,
        "mesh_convergence": conv,
        # full curves for step 3 (resonance positions need n_eff(lam), not just D)
        "sweep_fine": {"lam_m": rf["lam"].tolist(), "neff": rf["neff"].tolist(),
                       "ng_lumerical": rf["ng_lum"].tolist(),
                       "D_lumerical_s_m2": rf["D_lum"].tolist(),
                       "ng_ours_fit": d_f["ng_fit"].tolist(),
                       "D_ours_fit_s_m2": d_f["D_fit"].tolist()},
        "neff_poly": {"coeffs_highest_first": d_f["poly"].tolist(),
                      "x_definition": "x = (lam - 1.55e-6) / 50e-9"},
        "modes_at_1500nm_fine": rf["modes_at_1500nm"],
    }
    path = os.path.join(RESULTS, "s19_waveguide.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nsaved {path}")


if __name__ == "__main__":
    main()
