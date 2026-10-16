import numpy as np
import qutip as qt

def f(v, d=4):
    return round(float(np.real(v)), d)

N = 10                       # photon-number cutoff per mode
REP = 100e6                  # ASSUMED pump repetition rate: 100 MHz
ETA_H = 0.5                  # ASSUMED herald-arm efficiency (filters + coupling + detector)
ETA_S = 0.5                  # ASSUMED signal-arm efficiency (delivered to the user)

a = qt.destroy(N)
n = qt.num(N)

# ---------- building blocks ----------
def pair_state(mu):
    """Two-mode squeezed vacuum (Session 12), mean mu photons per arm."""
    x = mu / (1 + mu)                                   # = tanh^2 r
    psi = sum(np.sqrt((1 - x) * x**k) * qt.tensor(qt.basis(N, k), qt.basis(N, k))
              for k in range(N))
    return qt.ket2dm(psi)

def herald_weights(eta, kind):
    """P(detector fires | k photons in idler)."""
    if kind == "click":                                 # doorbell: 1 or more
        return np.array([1 - (1 - eta)**k for k in range(N)])
    return np.array([k * eta * (1 - eta)**(k - 1) if k > 0 else 0.0   # counter: exactly 1
                     for k in range(N)])

def herald(mu, eta_h, kind):
    """Signal state given that the idler detector fired, and the herald probability."""
    w = herald_weights(eta_h, kind)
    M = qt.tensor(qt.qeye(N), qt.Qobj(np.diag(np.sqrt(w))))
    unnorm = (M * pair_state(mu) * M).ptrace(0)         # Session 6 partial trace
    p = float(np.real(unnorm.tr()))
    return unnorm / p, p

A = qt.tensor(qt.destroy(N), qt.qeye(N))
B = qt.tensor(qt.qeye(N), qt.destroy(N))
def lossy(rho, eta):
    """Session 10: loss = beam splitter with hidden vacuum mode, then partial trace."""
    th = np.arccos(np.sqrt(eta))
    U = (th * (A.dag() * B - A * B.dag())).expm()
    return (U * qt.tensor(rho, qt.ket2dm(qt.basis(N, 0))) * U.dag()).ptrace(0)

def stats(rho):
    P = np.real(rho.diag())
    m = qt.expect(n, rho)
    g2 = qt.expect(a.dag() * a.dag() * a * a, rho) / m**2
    return P[0], P[1], 1 - P[0] - P[1], g2

def g2_series(mu, eta, kind):
    """Theory: exact g2 of the heralded state, summed by hand (no QuTiP)."""
    x = mu / (1 + mu)
    k = np.arange(200)
    wk = x**k * np.where(k > 0, 1 - (1 - eta)**k, 0) if kind == "click" else \
         x**k * k * eta * (1 - eta)**np.maximum(k - 1, 0)
    return np.sum(k * (k - 1) * wk) * np.sum(wk) / np.sum(k * wk)**2

# ---------- A. Validate the heralding model ----------
for mu in [0.01, 0.1]:
    rho, p = herald(mu, ETA_H, "click")
    g2 = stats(rho)[3]
    print(mu, f(p), f(ETA_H * mu / (1 + ETA_H * mu)),            # herald prob: sim, formula
          f(g2), f(g2_series(mu, ETA_H, "click")),                # g2: sim, exact series
          f(2 * mu * (2 - ETA_H)))                                # g2: small-mu approximation
# expect: 0.01 0.005 0.005 0.0294 0.0294 0.03
#         0.1 0.0476 0.0476 0.25 0.25 0.3   (approximation drifts at higher mu)

rho, _ = herald(0.1, ETA_H, "click")
print(f(stats(rho)[3]), f(stats(lossy(rho, ETA_S))[3]))           # expect 0.25 0.25 (g2 survives loss)

# ---------- B. Compare four designs ----------
print("design | herald MHz | P1 | P2+ | g2 | delivered single photons MHz")

# D1: attenuated laser, mean 0.1 at the chip, then the same signal-arm loss
laser = lossy(qt.ket2dm(qt.coherent(N, np.sqrt(0.1))), ETA_S)
P0, P1, P2, g2 = stats(laser)
print("D1 laser mu=0.1      | n/a  ", f(P1), f(P2), f(g2, 3), f(P1 * REP / 1e6, 2))

for name, mu, kind in [("D2 ring mu=0.01 click", 0.01, "click"),
                       ("D3 ring mu=0.1  click", 0.1, "click"),
                       ("D4 ring mu=0.1  PNR  ", 0.1, "pnr")]:
    rho, p = herald(mu, ETA_H, kind)
    P0, P1, P2, g2 = stats(lossy(rho, ETA_S))
    print(name, "|", f(p * REP / 1e6, 2), f(P1), f(P2), f(g2, 3), f(p * P1 * REP / 1e6, 2))
# expect:
# D1 laser mu=0.1      | n/a   0.0476 0.0012 1.0 4.76
# D2 ring mu=0.01 click | 0.5 0.5 0.0037 0.029 0.25
# D3 ring mu=0.1  click | 4.76 0.4981 0.0368 0.25 2.37
# D4 ring mu=0.1  PNR   | 4.54 0.4992 0.0238 0.17 2.26

# ---------- C. Rate vs purity trade-off ----------
print("mu | herald MHz | g2 click | g2 PNR")
for mu in [0.001, 0.003, 0.01, 0.03, 0.1]:
    _, p = herald(mu, ETA_H, "click")
    print(mu, f(p * REP / 1e6, 3),
          f(stats(herald(mu, ETA_H, "click")[0])[3]),
          f(stats(herald(mu, ETA_H, "pnr")[0])[3]))
# expect: 0.001 0.05 0.003 0.002 | 0.003 0.15 0.0089 0.006 | 0.01 0.498 0.0294 0.0197
#         0.03 1.478 0.0849 0.057 | 0.1 4.762 0.25 0.1701

# ---------- D. Spectral purity sets HOM between two independent sources ----------
M3 = 3
def op(k):
    ops = [qt.qeye(M3)] * 4
    ops[k] = qt.destroy(M3)
    return qt.tensor(ops)

a1, a2, b1, b2 = [op(k) for k in range(4)]       # port a / port b, spectral mode 1 / 2
G = (a1.dag() * b1 - a1 * b1.dag()) + (a2.dag() * b2 - a2 * b2.dag())
U = (np.pi / 4 * G).expm()
Na = a1.dag() * a1 + a2.dag() * a2
Nb = b1.dag() * b1 + b2.dag() * b2
vac4 = qt.tensor([qt.basis(M3, 0)] * 4)
ina = [a1.dag(), a2.dag()]
inb = [b1.dag(), b2.dag()]

def hom_visibility(lam):
    """Each heralded photon is in spectral mode 1 with prob lam, mode 2 with 1-lam."""
    w = [lam, 1 - lam]
    C = sum(w[j] * w[k] * qt.expect(Na * Nb, U * ina[j] * inb[k] * vac4)
            for j in range(2) for k in range(2))
    return 1 - float(np.real(C)) / 0.5

for lam in [1.0, 0.95, 0.9, 0.5]:
    purity = lam**2 + (1 - lam)**2
    print(lam, f(hom_visibility(lam)), f(purity), f(1 / purity, 3))   # V, purity, Schmidt number K
# expect: 1.0 1.0 1.0 1.0 | 0.95 0.905 0.905 1.105 | 0.9 0.82 0.82 1.22 | 0.5 0.5 0.5 2.0
