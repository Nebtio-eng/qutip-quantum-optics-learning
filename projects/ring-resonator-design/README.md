# Ring resonator design

Lumerical MODE simulation of the microring of Kues et al., Nature 546, 622 (2017), from the waveguide cross-section to the quantum figures of merit of the photon-pair source. The analysis and plots are in [notebook 13](../../notebooks/13_ring_resonator_design.ipynb), which reads the saved results and runs without Lumerical.

## Scripts (`lumerical/`, run in order with Ansys Lumerical MODE)

| Script | What it computes | Runtime |
|---|---|---|
| `s19_waveguide.py` | FDE mode solve of the core cross-section: n_eff, n_g and dispersion vs wavelength, mesh convergence, ring radius for 200 GHz FSR | ~3 min |
| `s19_coupler.py` | Bus–ring supermodes vs gap, coupling along the curved gap, Q budget, gap for Q_L = 235,000, EME check, coupling to higher-order modes | ~10 min |
| `s19_dispersion.py` | Silica material dispersion, resonance alignment and bin brightness vs dispersion, avoided-crossing search; `python s19_dispersion.py analyze` reruns only the analysis | ~8 min |
| `s19_quantum.py` | Bin-state fidelity in QuTiP, coherence time conventions, design → figure-of-merit table (no Lumerical needed) | seconds |

Each script writes `results/s19_<step>.json`, which the next step reads.

## Assumptions

- Core n = 1.70, 1.50 µm × 1.45 µm, silica cladding (from earlier papers on this platform).
- Propagation loss 0.06 dB/cm.
- Bus waveguide identical to the ring waveguide.
- Frequency-bin states use the bin pairs nearest the pump.

## Key results

| Quantity | Value |
|---|---|
| n_eff, n_g at 1550 nm | 1.608, 1.745 |
| Ring radius (200 GHz FSR) | 136.9 µm |
| Dispersion D (with silica dispersion) | 44.3 ps/(nm·km) |
| Gap for Q_L = 235,000 | 0.70 µm (0.63–0.70 µm including the EME check) |
| Linewidth | 823 MHz |
| Escape efficiency (max heralding efficiency) | 0.954 |
| Brightness of bin pair 10 | 0.74 of bin pair 1 |
