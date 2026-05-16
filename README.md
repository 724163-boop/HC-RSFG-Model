# HC-RSFG: Hierarchical Core Red Star-Forming Galaxy Theory

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the numerical simulation codes for the **HC‑RSFG (Hierarchical Core Red Star‑Forming Galaxy)** theory.  
The theory explains the extreme redness of certain galaxies – from the early Universe’s Little Red Dots (LRDs) to local relic galaxies – **entirely through pure stellar formation** with an extremely bottom‑heavy Initial Mass Function (IMF).  
No dust extinction or active galactic nucleus (AGN) is required.

---

## 1. Theoretical Overview

### The Puzzle of Red Galaxies
JWST discovered numerous compact, very red objects – “Little Red Dots” (LRDs) – in the early Universe (z ≈ 5‑7).  
Their colours are traditionally attributed to either heavy dust obscuration or an accreting supermassive black hole (AGN).  
However, several observations contradict these explanations:
- Some LRDs are X‑ray quiet and show only narrow emission lines.
- Dust‑free LRDs have been reported at z > 7.

This suggests a third possibility: **the stellar population itself may be intrinsically red**.

### The HC‑RSFG Solution
At the earliest stages of galaxy formation (z ≈ 7.5), the central core of a protogalaxy reaches extreme gas densities (`n ≈ 10⁷ cm⁻³`) and very low temperatures.  
Under these conditions, gas cools **~30 000 times faster than it collapses**, keeping the Jeans mass below ~1.5 M☉.  
This shuts off the formation of O‑, B‑, and most A‑type stars, leading to a stellar population dominated by K‑ and M‑dwarfs.  
Such stars are intrinsically red, giving the galaxy its red colour **from birth**.

### A Unified but Flexible Evolutionary Picture
The HC‑RSFG mechanism can act in different environments and epochs, naturally producing a variety of observed galaxy types.  
One **possible** evolutionary sequence is:

- **Early ultra‑faint stage (e.g., z ≈ 7‑10)**: extremely dim, currently undetected (“Dark‑RSFG”).
- **Bright compact stage (e.g., z ≈ 5‑7)**: observed as narrow‑line Little Red Dots (LRDs).
- **Intermediate stage (e.g., z ≈ 0.5)**: red but still star‑forming galaxies (Red Misfits).
- **Final quiescent stage (z ≈ 0)**: massive relic ellipticals (e.g., NGC 1277) preserving a bottom‑heavy IMF.

However, this is only one illustrative pathway; depending on gas supply, mergers, and AGN feedback, an HC‑RSFG core may evolve directly into a relic galaxy, become a Red Misfit without a prior LRD phase, or even be quenched and fade away.

---

## 2. Repository Contents

All codes are written in **Python 3** and require the **FSPS (Flexible Stellar Population Synthesis)** library.  
They are designed to run in **Google Colab**.

| File | Description |
|:---|:---|
| `HC_RSFG_Code1_LRD.py` | Simulation of a compact red galaxy at z = 6 (similar to observed LRDs). Computes JWST broadband colours and emission‑line equivalent widths. |
| `HC_RSFG_Code2_RedMisfit.py` | Two‑component (core + disk) model at z = 0.5; reproduces Red Misfit colours and Hα equivalent widths. |
| `HC_RSFG_Code3_RelicGalaxy.py` | Relic galaxy simulation at z = 0.0173 (NGC 1277 analog). Radial profile of M/L ratio. |
| `HC_RSFG_Code4_DarkRSFG.py` | Early, extremely faint stage (z ≈ 7–10). Predicted colours and detectability. |
| `HC_RSFG_Code5_AGN_Arrival.py` | Timeline of black‑hole growth from a seed to a low‑mass AGN concurrently with the stellar core. |
| `HC_RSFG_Code5_AGN_Point.py` | Energetic assessment: can the AGN expel the remaining core gas and quench star formation? |
| `HC_RSFG_Metallicity_Mixing.py` | Quantitative mixing model that justifies the adopted metallicity (Z ≈ 0.1 Z☉) through Cold Stream enrichment. |
| `HC_RSFG_Theory_Core.py` | Foundational calculation: cooling time vs. free‑fall time using [CII] 158 µm and [OI] 63 µm fine‑structure lines. |
| `HC_RSFG_TimeAxis_Integrity.py` | Gas accumulation simulation demonstrating that the critical density is reached within 220 Myr, consistent with the cosmic timeline. |
| `HC_RSFG_FigureGenerator_v2.py` | Generates the 5‑panel comparison figure (LRD, relic, Red Misfit, Dark‑RSFG, M/L profile) using **real FSPS filter transmission curves**. |
| `HC_RSFG_Relic_NGC1277_v2.py` | Relic galaxy simulation updated with real SDSS filter curves and a more extreme IMF model (imf3=4.0) to match NGC 1277 observations. |
| `index.html` | Project website presenting the theory, simulation results, and interactive visualizations. |
| `LICENSE` | MIT License. |
---

## 3. Quick Start (Google Colab)

1. Open a new **Google Colab** notebook.
2. Copy the entire content of any `.py` file into a cell.
3. Run the cell. The first execution will:
   - Clone the FSPS repository and compile the necessary Fortran modules.
   - Set the `SPS_HOME` environment variable.
   - Install the `fsps` Python wrapper.
   - Perform the simulation and print the results.
4. Output figures (if any) will be saved in the Colab runtime environment.

**Note:** the first run may take several minutes because FSPS must compile its stellar population synthesis generators.

---

## 4. Key Numerical Results

| Observable | HC‑RSFG Prediction | Observed Value | Reference |
|:---|:---|:---|:---|
| F200W − F444W colour of a compact z = 6 galaxy | **+1.97 mag** (dust‑free) | > 1.0 mag | Hviding + 2025 |
| Hα equivalent width (core) | < 1 Å | — (narrow lines only) | — |
| u − r colour of a Red Misfit (z = 0.5) | **+3.55 mag** | > 2.22 mag | Evans + 2018 |
| M/L ratio (NGC 1277 centre) | **1.56 × Kroupa** (lower limit) | ~2.0 × Kroupa | Martín‑Navarro + 2015 |

> The M/L ratio is a conservative lower limit; increasing the IMF slope or core age brings it into exact agreement with the observed value.

---

## 5. License & Citation

This project is released under the **MIT License**.  
If you use these codes or the HC‑RSFG theory in your research, please cite the forthcoming paper (arXiv link to be added) and this repository.

---

## 6. References
- Steinhardt (2025), *ApJ*, 982, 189 – original RSFG concept.
- Hviding et al. (2025), *A&A*, 702, A57 – RUBIES survey of LRDs.
- Labbé et al. (2023), *Nature*, 616, 266 – discovery of LRDs.
- Evans et al. (2018), *MNRAS*, 476, 5284 – Red Misfits in SDSS.
- Martín‑Navarro et al. (2015), *MNRAS*, 451, 1081 – NGC 1277 bottom‑heavy IMF.
- Zhang et al. (2025), *ApJ*, 998, 123 – narrow‑line LRD spectroscopy.

