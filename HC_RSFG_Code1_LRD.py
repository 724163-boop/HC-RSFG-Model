#### Code 1: LRD State Simulation
```python
#!/usr/bin/env python3
"""
HC-RSFG Theory Code 1: Little Red Dot (LRD) State Simulation
Phase 4: HC-RSFG Full Activation Phase (z=6.0, tage=0.9 Gyr)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

import fsps

# 1. Parameters
Z = 6.0
TAGE = 0.9 
LOGZSOL = -1.0 
DUST2 = 0.08 
TAU = 0.10 
GAS_LOGU = -3.0 

# CRITICAL FIX: imf_type must be 0 to use custom imf1, imf2, imf3
IMF_TYPE = 0 
IMF1 = 1.8
IMF2 = 1.8
IMF3 = 3.0 
IMF_UMAX = 1.5 
LRD_THRESH = 1.0 

BANDS = ["jwst_f115w", "jwst_f150w", "jwst_f200w", "jwst_f277w", "jwst_f356w", "jwst_f444w"]
SEED = 42

# 2. FSPS Instance Creation
def make_sps(neb=False, imf_type=IMF_TYPE, logzsol=LOGZSOL, dust2=DUST2, zred=Z, tage=TAGE, imf_umax=IMF_UMAX):
    kw = dict(compute_vega_mags=False, zcontinuous=1)
    if neb:
        kw["add_neb_emission"] = True
    sps = fsps.StellarPopulation(**kw)
    sps.params["imf_type"] = imf_type
    if imf_type == 0:
        sps.params["imf1"] = IMF1
        sps.params["imf2"] = IMF2
        sps.params["imf3"] = IMF3
    sps.params["logzsol"] = logzsol
    sps.params["sfh"] = 1
    sps.params["tau"] = TAU
    sps.params["tage"] = tage
    sps.params["dust2"] = dust2
    sps.params["zred"] = zred
    if imf_type == IMF_TYPE:
        sps.params["imf_upper_limit"] = imf_umax
    if neb:
        sps.params["gas_logz"] = logzsol
        sps.params["gas_logu"] = GAS_LOGU
    return sps

print("=================================================================")
print(" HC-RSFG Theory - Code 1: LRD State Simulation")
print(f" z={Z}, tage={TAGE} Gyr, logzsol={LOGZSOL}, dust2={DUST2}")
print(f" IMF: Custom Broken Power Law (imf3={IMF3})")
print("=================================================================")

sps_hc = make_sps(neb=False)
sps_neb = make_sps(neb=True)
sps_std = make_sps(neb=False, imf_type=1) # 1 = Chabrier

wave_obs, spec_hc = sps_hc.get_spectrum(tage=TAGE)
wave_neb, spec_neb = sps_neb.get_spectrum(tage=TAGE)
wave_std, spec_std = sps_std.get_spectrum(tage=TAGE)

wave_rest = wave_obs / (1 + Z)
wave_rest_neb = wave_neb / (1 + Z)

# 3. Band Colors
mags_hc = sps_hc.get_mags(tage=TAGE, bands=BANDS)
mags_std = sps_std.get_mags(tage=TAGE, bands=BANDS)

def get_color(mags, b1, b2):
    return mags[BANDS.index(b1)] - mags[BANDS.index(b2)]

c200_444_hc = get_color(mags_hc, "jwst_f200w", "jwst_f444w")
c277_444_hc = get_color(mags_hc, "jwst_f277w", "jwst_f444w")
c115_200_hc = get_color(mags_hc, "jwst_f115w", "jwst_f200w")
c200_444_std = get_color(mags_std, "jwst_f200w", "jwst_f444w")
c277_444_std = get_color(mags_std, "jwst_f277w", "jwst_f444w")

print("\n[JWST Band Colors (HC-RSFG vs Chabrier)]")
print(f" {'Color':<15} {'HC-RSFG':>10} {'Chabrier':>10} {'Diff':>8}")
print("-" * 50)
for label, hc, std in [
    ("F115W-F200W", c115_200_hc, None),
    ("F200W-F444W", c200_444_hc, c200_444_std),
    ("F277W-F444W", c277_444_hc, c277_444_std),
]:
    if std is not None:
        diff = hc - std
        lrd_flag = " LRD!" if (label == "F200W-F444W" and hc > LRD_THRESH) else ""
        print(f" {label:<15} {hc:>+10.3f} {std:>+10.3f} {diff:>+8.3f} {lrd_flag}")
    else:
        print(f" {label:<15} {hc:>+10.3f}")

# 4. Emission Line EW
EMISSION_LINES = {
    "[OII] 3727": 3727.1, "[NeIII] 3869": 3869.9, "Hd 4102": 4102.9,
    "Hg 4341": 4341.7, "Hb 4863": 4862.7, "[OIII] 5007": 5007.0,
    "Ha 6565": 6564.6, "[NII] 6583": 6583.4,
}

def compute_ew(wave_r, spec, w0, dw_line=20, dw_cont=80):
    m_line = (wave_r > w0 - dw_line) & (wave_r < w0 + dw_line)
    m_cont = ((wave_r > w0 - dw_cont) & (wave_r < w0 - dw_line * 2)) | \
             ((wave_r > w0 + dw_line * 2) & (wave_r < w0 + dw_cont))
    if m_line.sum() == 0 or m_cont.sum() == 0: return 0.0
    peak = np.max(spec[m_line])
    cont = np.median(spec[m_cont])
    if cont <= 0: return 0.0
    return (peak - cont) / cont * (2 * dw_line)

print("\n[Emission Line EW (logzsol=-1.0, gas_logu=-3.0)]")
print(f" {'Line':<15} {'EW [A]':>8} {'Classification'}")
print("-" * 45)
ew_results = {}
for name, w0 in EMISSION_LINES.items():
    ew = compute_ew(wave_rest_neb, spec_neb, w0)
    ew_results[name] = ew
    if ew > 2.0: tag = "Strong (AGN-like)"
    elif ew > 0.5: tag = "Moderate"
    elif ew > 0.05: tag = "Weak (Narrow LRD)"
    else: tag = "Undetected"
    print(f" {name:<15} {ew:>8.3f} {tag}")

# 5 & 6. Grids for Plots
logzsol_grid = np.array([-3.0, -2.5, -2.0, -1.5, -1.3, -1.0, -0.7, -0.5, -0.3, 0.0])
c200_444_grid, c277_444_grid = [], []
for lz in logzsol_grid:
    sp = make_sps(neb=False, logzsol=lz)
    m = sp.get_mags(tage=TAGE, bands=BANDS)
    c200_444_grid.append(get_color(m, "jwst_f200w", "jwst_f444w"))
    c277_444_grid.append(get_color(m, "jwst_f277w", "jwst_f444w"))
c200_444_grid, c277_444_grid = np.array(c200_444_grid), np.array(c277_444_grid)

dust2_grid = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0])
c200_444_dust = []
for d2 in dust2_grid:
    sp = make_sps(neb=False, logzsol=LOGZSOL, dust2=d2)
    m = sp.get_mags(tage=TAGE, bands=BANDS)
    c200_444_dust.append(get_color(m, "jwst_f200w", "jwst_f444w"))
c200_444_dust = np.array(c200_444_dust)

# 7. Plotting
DARK_BG, PANEL_BG = "#090912", "#0d0d20"
C_HC, C_STD, C_NEB, C_OBS, C_TXT = "#ff7043", "#4dd0e1", "#ef5350", "#ffd54f", "#e8eaf6"
plt.rcParams.update({"figure.facecolor": DARK_BG, "text.color": C_TXT, "axes.labelcolor": C_TXT, 
                     "xtick.color": C_TXT, "ytick.color": C_TXT, "axes.edgecolor": "#3d3d6b", "font.family": "sans-serif"})

fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.3)

# Panel A
axA = fig.add_subplot(gs[0, :], facecolor=PANEL_BG)
mask_jwst = (wave_obs > 9000) & (wave_obs < 55000)
norm_A = np.percentile(spec_hc[mask_jwst], 95)
axA.plot(wave_obs[mask_jwst]/1e4, spec_hc[mask_jwst]/norm_A, color=C_HC, lw=1.8, label=f"HC-RSFG (imf3={IMF3})")
axA.plot(wave_std[mask_jwst]/1e4, spec_std[mask_jwst]/norm_A, color=C_STD, lw=1.2, ls="--", label="Chabrier")
axA.plot(wave_neb[mask_jwst]/1e4, spec_neb[mask_jwst]/norm_A, color=C_NEB, lw=0.9, alpha=0.5, label="HC-RSFG+Nebular")
axA.set(xlabel="Observed Wavelength [μm]", ylabel="Normalized Flux", title=f"Phase 4 LRD Spectrum (z={Z})")
axA.legend(framealpha=0.2)

# Panel B
axB = fig.add_subplot(gs[1, 0], facecolor=PANEL_BG)
mask_opt = (wave_rest > 2200) & (wave_rest < 7200)
norm_B = np.percentile(spec_hc[mask_opt], 98)
axB.plot(wave_rest[mask_opt], spec_hc[mask_opt]/norm_B, color=C_HC, label="HC-RSFG")
axB.plot(wave_rest[mask_opt], spec_std[mask_opt]/norm_B, color=C_STD, ls="--", label="Chabrier")
axB.set(xlabel="Rest Wavelength [A]", ylabel="Normalized Flux", title="Rest-frame UV-Optical")
axB.legend(framealpha=0.2)

# Panel C
axC = fig.add_subplot(gs[1, 1], facecolor=PANEL_BG)
ew_names = list(ew_results.keys())
ew_vals = [ew_results[k] for k in ew_names]
axC.barh(ew_names, ew_vals, color=C_NEB, alpha=0.8)
axC.axvline(1.0, color=C_OBS, ls="--", label="EW=1A Threshold")
axC.set(xlabel="Equivalent Width [A]", title="Emission Line EW (< 1A)")
axC.legend(framealpha=0.2)

# Panel D
axD = fig.add_subplot(gs[1, 2], facecolor=PANEL_BG)
band_waves = np.array([1.15, 1.50, 2.00, 2.77, 3.56, 4.44])
flux_hc = 10**(-0.4*(np.array([mags_hc[BANDS.index(b)] for b in BANDS]) - mags_hc[4]))
flux_std = 10**(-0.4*(np.array([mags_std[BANDS.index(b)] for b in BANDS]) - mags_std[4]))
axD.plot(band_waves, flux_hc, "o-", color=C_HC, label="HC-RSFG")
axD.plot(band_waves, flux_std, "s--", color=C_STD, label="Chabrier")
axD.set(xlabel="Observed Wavelength [μm]", ylabel="Relative Flux (F356W=1)", title="Mock SED Fitting")
axD.legend(framealpha=0.2)

# Panel E
axE = fig.add_subplot(gs[2, 0], facecolor=PANEL_BG)
axE.plot(logzsol_grid, c200_444_grid, "o-", color=C_HC)
axE.axhline(LRD_THRESH, color=C_OBS, ls="--", label=f"LRD Threshold")
axE.set(xlabel="Metallicity [log Z/Zsun]", ylabel="F200W - F444W [mag]", title="Metallicity vs Color")
axE.legend(framealpha=0.2)

# Panel F
axF = fig.add_subplot(gs[2, 1], facecolor=PANEL_BG)
av_grid = dust2_grid * 3.1
axF.plot(av_grid, c200_444_dust, "o-", color=C_HC)
axF.axhline(LRD_THRESH, color=C_OBS, ls="--")
axF.set(xlabel="Dust Av [mag]", ylabel="F200W - F444W [mag]", title="Dust vs Color (logZ=-1.0)")

# Panel G
axG = fig.add_subplot(gs[2, 2], facecolor=PANEL_BG)
sc = axG.scatter(c200_444_grid, c277_444_grid, c=logzsol_grid, cmap="plasma")
plt.colorbar(sc, ax=axG).set_label("log Z/Zsun")
axG.scatter([c200_444_hc], [c277_444_hc], s=150, color=C_OBS, marker="*", label="This Model")
axG.set(xlabel="F200W - F444W [mag]", ylabel="F277W - F444W [mag]", title="Color-Color Diagram")
axG.legend(framealpha=0.2)

fig.suptitle("HC-RSFG Theory - Code 1: LRD State Simulation", fontsize=14, y=0.95)
plt.savefig("HC_RSFG_Code1_LRD.png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()
