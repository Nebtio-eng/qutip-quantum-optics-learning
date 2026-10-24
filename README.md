# Quantum Optics with QuTiP: from a single photon to an on-chip light source

A self-built course in quantum optics for integrated photonics, written while moving from classical photonics (Zemax, Lumerical) into quantum photonics. It runs from the first quantum state to a microring photon-pair source designed in Lumerical, modelled in QuTiP and checked against published experiments.

Every result is **validated**: each notebook compares its simulations with analytic formulas or an independent method through a `check()` helper and ends with a pass count. Notebooks 09, 10 and 12 compare the models with published experiments.

**213 automated theory-vs-simulation checks, all passing.**

## Headline results

**1. [Ring resonator design → quantum performance](notebooks/13_ring_resonator_design.ipynb)** (Lumerical MODE + QuTiP). The microring of Kues et al. (Nature 2017) rebuilt from its waveguide cross-section: the model reproduces the 200 GHz FSR with a 136.9 µm radius and the 235,000 loaded Q with an ~800 MHz linewidth. It shows that 95 % of the generated photons can leave the ring (the ceiling on heralding efficiency), and that dispersion keeps at least 10 frequency-bin pairs usable.

**2. [Source model validated against a real experiment](notebooks/10_validation_heralded_source.ipynb).** Fitted with one parameter to the heralded g⁽²⁾ of a silicon pair source (Paesani et al., Nat. Commun. 2020), the model predicts the measured unheralded g⁽²⁾ (1.99 vs 1.97 ± 0.03) and the two-source HOM visibility (0.944 vs 0.96 ± 0.02).

**3. [Why fidelity falls with dimension in Kues et al. 2017](notebooks/12_kues2017_fidelity_analysis.ipynb).** Fixed white noise and unequal bin brightness are ruled out, the latter by both the measured spectrum and the ring's own dispersion. The remaining explanations point to the measurement setup rather than the chip, and the notebook proposes the measurement that separates them.

**4. [Quantum light source design](projects/quantum-light-source-design/design_report.md).** Five heralded-source designs compared: a gently pumped ring with a photon-number-resolving herald gives g⁽²⁾(0) = 0.020, against 1.0 for an attenuated laser.

![Heralded g2 against herald rate](projects/quantum-light-source-design/figures/g2_vs_rate.png)

## Learning path

| # | Notebook | What it covers | Photonics link |
|---|---|---|---|
| 01 | [Quantum foundations](notebooks/01_quantum_foundations.ipynb) | States, Born rule, operators, eigenvalues, expectation values, Bloch sphere | Polarization: a quantum state is a Jones vector |
| 02 | [Single-photon and nonclassical states](notebooks/02_single_photon_states.ipynb) | Fock, coherent and squeezed states, g⁽²⁾(0), Mandel Q, Wigner function | Laser vs single-photon source vs squeezed light |
| 03 | [Gates and operations](notebooks/03_quantum_gates_operations.ipynb) | Unitarity, Hadamard, phase shifter, Mach–Zehnder, CNOT | Beam splitters, phase shifters, directional couplers |
| 04 | [Entanglement](notebooks/04_entanglement_multiparticle.ipynb) | Tensor products, Bell states, entanglement entropy | Entangled photon-pair sources |
| 05 | [Measurements and density matrices](notebooks/05_measurements_density_matrices.ipynb) | Mixed states, partial trace, loss, detector efficiency, tomography | Real chips lose light, real detectors miss photons |
| 06 | [Hong–Ou–Mandel effect](notebooks/06_hong_ou_mandel_effect.ipynb) | Two-photon interference, distinguishability, dip vs delay | The standard test of photon indistinguishability |
| 07 | [Photonic quantum circuits](notebooks/07_photonic_quantum_circuits.ipynb) | Chips as unitaries, permanents, boson sampling, teleportation | Programmable MZI meshes |
| 08 | [Quantum frequency combs](notebooks/08_quantum_frequency_combs.ipynb) | Four-wave mixing, two-mode squeezing, frequency-bin entanglement | Microring sources (Kues et al., Nature 2017) |
| 09 | [Validation against experiment](notebooks/09_validation_and_comparison.ipynb) | Theory vs simulation vs published data | Hong–Ou–Mandel 1987; Kues et al. 2017 |
| 10 | [Fitting the source model to data](notebooks/10_validation_heralded_source.ipynb) | Multimode pair statistics, purity from g⁽²⁾, fitting, multi-photon HOM | Paesani et al. 2020 silicon source |
| 11 | [Time evolution and g⁽²⁾(τ)](notebooks/11_time_evolution_g2_tau.ipynb) | Schrödinger and master equations, photon lifetime, bunching and antibunching | Directional coupler, ring Q, HBT measurement |
| 12 | [Fidelity vs dimension](notebooks/12_kues2017_fidelity_analysis.ipynb) | Noise, imbalance and dephasing models fitted to tomography data | Kues et al. 2017 frequency-bin qudits |
| 13 | [Ring resonator design](notebooks/13_ring_resonator_design.ipynb) | Mode solving, coupling, Q budget, dispersion, bin brightness | Lumerical design to quantum figures of merit |

**Projects**
- [Ring resonator design](projects/ring-resonator-design/): Lumerical MODE scripts (FDE and EME) and their saved results, read by notebook 13.
- [Quantum light source design](projects/quantum-light-source-design/): [report](projects/quantum-light-source-design/design_report.md), [notebook](projects/quantum-light-source-design/simulations.ipynb), [script](projects/quantum-light-source-design/analysis.py).

Every concept is explained the same way: the idea in plain words with an optics analogy, the honest gap in that explanation, a one-sentence summary, and then code that checks the simulation against the analytic answer.

## Validation against published experiments

| Quantity | Model | Experiment | Verdict |
|---|---|---|---|
| HOM visibility, weak coherent light | 0.5 | limit quoted in the literature | agrees |
| Unheralded g⁽²⁾, silicon pair source (Paesani et al.) | 1.990 | 1.97 ± 0.03 | agrees |
| Two-source HOM visibility (Paesani et al.) | 0.944 | 0.96 ± 0.02 | agrees within error |
| Ring photon lifetime / linewidth, Q = 92,000 (Ma et al.) | 75.7 ps / 2.1 GHz | ~76 ps / ~2.1 GHz | agrees |
| Ring linewidth, Q = 235,000 (Kues et al.) | 823 MHz | ~800 MHz | agrees |
| Microring effective dimension (Kues et al.) | g⁽²⁾ = 1 + 1/K | 10.45 ± 0.53 | ideal K = 10 within 0.85 error bars |
| Qudit fidelity vs dimension (Kues et al.) | white noise: almost flat | falls 0.885 → 0.766 | white noise ruled out; analysed in notebook 12 |

## Run it

```bash
git clone https://github.com/<your-username>/qutip-quantum-optics-learning.git
cd qutip-quantum-optics-learning
pip install -r requirements.txt
jupyter lab
```

Every notebook runs top to bottom on a laptop (notebook 09 takes under a minute) and ends with `N of N checks passed`. Notebook 13 reads saved Lumerical results, so it needs no Lumerical license; the scripts in `projects/ring-resonator-design/lumerical/` regenerate them with Ansys Lumerical MODE. Built and tested with Python 3, QuTiP 5.3.1.

## Limitations and next steps

- The ring geometry and core index are taken from earlier papers on the same platform; the core material's dispersion is unknown and is treated as a range.
- The bus–ring coupling model and an EME check differ by up to 35 %; this moves the design gap, not the escape efficiency.
- Avoided crossings with higher-order ring modes cannot be located at the model's precision.
- Dark counts and Raman noise are not yet in the source model.

## References

- C. K. Hong, Z. Y. Ou, L. Mandel, Phys. Rev. Lett. 59, 2044 (1987).
- M. Kues et al., "On-chip generation of high-dimensional entangled quantum states and their coherent control," Nature 546, 622 (2017).
- C. Reimer et al., Frontiers in Optics FTh3E.2 (2017).
- S. Paesani et al., "Near-ideal spontaneous photon sources in silicon quantum photonics," Nat. Commun. 11, 2505 (2020). [arXiv:2005.09579](https://arxiv.org/abs/2005.09579)
- C. Ma et al., "Silicon photonic entangled photon-pair and heralded single photon generation with CAR > 12,000 and g⁽²⁾(0) < 0.006," [arXiv:1710.01001](https://arxiv.org/abs/1710.01001).
- J. R. Johansson, P. D. Nation, F. Nori, QuTiP, [qutip.org](https://qutip.org).
