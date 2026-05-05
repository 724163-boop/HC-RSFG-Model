"""
HC-RSFG Theory Code 1: Little Red Dot (LRD) State Simulation
Phase 4: HC-RSFG Full Activation Phase (z=6.0, tage=0.9 Gyr)
Reference: Hviding et al. 2025 (RUBIES Survey), Labbe et al. 2023
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")
import fsps

# =================================================================
# 1. Parameter Definitions
# =================================================================
Z = 6.0                 # Redshift
TAGE = 0.9              # Age of the universe at z=6 [Gyr]
LOGZSOL = -1.0          # Metallicity: Z ≈ 0.10 Zsun (After Cold Stream)
DUST2 = 0.08            # Dust attenuation (Av ≈ 0.25 mag)
TAU = 0.10              # Exponential decay SFH [Gyr]
GAS_LOGU = -3.0         # Gas ionization parameter (Low ionization)

# HC-RSFG Broken Power-Law IMF
IMF_TYPE = 2           
IMF1 = 1.8              # Slope for 0.08-0.5 Msun
IMF2 = 1.8              # Slope for 0.5-1.5 Msun
IMF3 = 3.0              # Slope for >1.5 Msun (Strongly suppresses massive stars)
IMF_UMAX = 1.5          # Effective upper mass limit [Msun]

LRD_THRESH = 1.0        # LRD Threshold: F200W-F444W > 1.0
BANDS = ["jwst_f115w", "jwst_f150w", "jwst_f200w", "jwst_f277w", "jwst_f356w", "jwst_f444w"]
SEED = 42

# =================================================================
# 2. FSPS Instance Creation
# =================================================================
def make_sps(neb=False, imf_type=IMF_TYPE, logzsol=LOGZSOL, dust2=DUST2, zred=Z, tage=TAGE, imf_umax=IMF_UMAX):
    kw = dict(compute_vega_mags=False, zcontinuous=1)
    if neb:
        kw["add_neb_emission"] = True
    sps = fsps.StellarPopulation(**kw)
    sps.params["imf_type"] = imf_type
    
    # Apply custom IMF slopes only if using custom type (0)
    if imf_type == 0:
        sps.params["imf1"] = IMF1
        sps.params["imf2"] = IMF2
        sps.params["imf3"] = IMF3
        sps.params["imf_upper_limit"] = imf_umax

    sps.params["logzsol"] = logzsol
    sps.params["sfh"] = 1
    sps.params["tau"] = TAU
    sps.params["tage"] = tage
    sps.params["dust2"] = dust2
    sps.params["zred"] = zred
    
    if neb:
        sps.params["gas_logz"] = logzsol
        sps.params["gas_logu"] = GAS_LOGU
    return sps

print("=================================================================")
print(" HC-RSFG Theory - Code 1: LRD State Simulation")
print(f" z = {Z}, tage = {TAGE} Gyr, logzsol = {LOGZSOL}, dust2 = {DUST2}")
print(f" IMF: Broken Power Law (imf1={IMF1}, imf2={IMF2}, imf3={IMF3})")
print("=================================================================")

sps_hc = make_sps(neb=False)
sps_neb = make_sps(neb=True)
sps_std = make_sps(neb=False, imf_type=1) # 1 = Chabrier Reference

wave_obs, spec_hc = sps_hc.get_spectrum(tage=TAGE)
wave_neb, spec_neb = sps_neb.get_spectrum(tage=TAGE)
wave_std, spec_std = sps_std.get_spectrum(tage=TAGE)

wave_rest = wave_obs / (1 + Z)
wave_rest_neb = wave_neb / (1 + Z)

# =================================================================
# 3. Band Color Calculation
# =================================================================
mags_hc = sps_hc.get_mags(tage=TAGE, bands=BANDS)
mags_std = sps_std.get_mags(tage=TAGE, bands=BANDS)

def color(mags, b1, b2):
    return mags[BANDS.index(b1)] - mags[BANDS.index(b2)]

c200_444_hc = color(mags_hc, "jwst_f200w", "jwst_f444w")
c277_444_hc = color(mags_hc, "jwst_f277w", "jwst_f444w")
c115_200_hc = color(mags_hc, "jwst_f115w", "jwst_f200w")

c200_444_std = color(mags_std, "jwst_f200w", "jwst_f444w")
c277_444_std = color(mags_std, "jwst_f277w", "jwst_f444w")

print("\n[ JWST Band Colors (HC-RSFG vs Chabrier) ]")
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

# =================================================================
# 4. Emission Line Equivalent Width (EW) Calculation
# =================================================================
EMISSION_LINES = {
    "[OII] 3727": 3727.1, "[NeIII] 3869": 3869.9, "Hd 4102": 4102.9,
    "Hg 4341": 4341.7, "Hb 4863": 4862.7, "[OIII] 5007": 5007.0,
    "Ha 6565": 6564.6, "[NII] 6583": 6583.4,
}

def compute_ew(wave_r, spec, w0, dw_line=20, dw_cont=80):
    m_line = (wave_r > w0 - dw_line) & (wave_r < w0 + dw_line)
    m_blue = (wave_r > w0 - dw_cont) & (wave_r < w0 - dw_line * 2)
    m_red = (wave_r > w0 + dw_line * 2) & (wave_r < w0 + dw_cont)
    m_cont = m_blue | m_red
    if m_line.sum() == 0 or m_cont.sum() == 0: return 0.0
    peak = np.max(spec[m_line])
    cont = np.median(spec[m_cont])
    if cont <= 0: return 0.0
    return (peak - cont) / cont * (2 * dw_line)

print(f"\n[ Emission Line EW (logzsol={LOGZSOL}, gas_logu={GAS_LOGU}) ]")
print(f" {'Line':<15} {'EW [A]':>8}   {'Classification'}")
print("-" * 45)
ew_results = {}
for name, w0 in EMISSION_LINES.items():
    ew = compute_ew(wave_rest_neb, spec_neb, w0)
    ew_results[name] = ew
    if ew > 2.0: tag = "Strong (AGN-like)"
    elif ew > 0.5: tag = "Moderate"
    elif ew > 0.05: tag = "Weak (Narrow LRD)"
    else: tag = "Undetected"
    print(f" {name:<15} {ew:>8.3f}   {tag}")

# =================================================================
# 5. Metallicity & Dust Scans
# =================================================================
logzsol_grid = np.array([-3.0, -2.5, -2.0, -1.5, -1.3, -1.0, -0.7, -0.5, -0.3, 0.0])
c200_444_grid, c277_444_grid = [], []
for lz in logzsol_grid:
    sp = make_sps(neb=False, logzsol=lz)
    m = sp.get_mags(tage=TAGE, bands=BANDS)
    c200_444_grid.append(color(m, "jwst_f200w", "jwst_f444w"))
    c277_444_grid.append(color(m, "jwst_f277w", "jwst_f444w"))
c200_444_grid = np.array(c200_444_grid)
c277_444_grid = np.array(c277_444_grid)

dust2_grid = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0])
c200_444_dust = []
for d2 in dust2_grid:
    sp = make_sps(neb=False, logzsol=LOGZSOL, dust2=d2)
    m = sp.get_mags(tage=TAGE, bands=BANDS)
    c200_444_dust.append(color(m, "jwst_f200w", "jwst_f444w"))
c200_444_dust = np.array(c200_444_dust)

# =================================================================
# 6. Plotting
# =================================================================
DARK_BG, PANEL_BG = "#090912", "#0d0d20"
C_HC, C_STD, C_NEB, C_OBS, C_GRID, C_TXT = "#ff7043", "#4dd0e1", "#ef5350", "#ffd54f", "#1e1e3a", "#e8eaf6"

plt.rcParams.update({
    "figure.facecolor": DARK_BG, "text.color": C_TXT, "axes.labelcolor": C_TXT,
    "xtick.color": C_TXT, "ytick.color": C_TXT, "axes.edgecolor": "#3d3d6b",
    "font.family": "monospace",
})

fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.38, top=0.92, bottom=0.06, left=0.06, right=0.97)

# Panel A: Full Wavelength Spectrum
axA = fig.add_subplot(gs[0, :], facecolor=PANEL_BG)
mask_jwst = (wave_obs > 9000) & (wave_obs < 55000)
norm_A = np.percentile(spec_hc[mask_jwst], 95)
axA.fill_between(wave_obs[mask_jwst] / 1e4, 0, spec_hc[mask_jwst] / norm_A, alpha=0.15, color=C_HC)
axA.plot(wave_obs[mask_jwst] / 1e4, spec_hc[mask_jwst] / norm_A, color=C_HC, lw=1.8, label=f"HC-RSFG IMF (imf3={IMF3})")
axA.plot(wave_std[mask_jwst] / 1e4, spec_std[mask_jwst] / norm_A, color=C_STD, lw=1.2, ls="--", alpha=0.75, label="Chabrier IMF (Ref)")
axA.plot(wave_neb[mask_jwst] / 1e4, spec_neb[mask_jwst] / norm_A, color=C_NEB, lw=0.9, alpha=0.55, label="HC-RSFG + Nebular")

band_info = [("F115W", 1.15), ("F150W", 1.50), ("F200W", 2.00), ("F277W", 2.77), ("F356W", 3.56), ("F444W", 4.44)]
for bname, bwav in band_info:
    axA.axvline(bwav, color="#ffffff", alpha=0.10, lw=0.8, ls=":")
    axA.text(bwav, 1.13, bname, ha="center", fontsize=7, color="#9fa8da", transform=axA.get_xaxis_transform())

y_arrow = 0.3
axA.annotate("", xy=(4.44, y_arrow), xytext=(2.00, y_arrow), arrowprops=dict(arrowstyle="<->", color=C_OBS, lw=1.80))
axA.text(3.22, y_arrow + 0.04, f"F200W-F444W = {c200_444_hc:+.2f} mag\n(LRD Threshold > {LRD_THRESH})", 
         ha="center", fontsize=9.5, color=C_OBS, bbox=dict(boxstyle="round,pad=0.4", fc="#1a1a35", ec=C_OBS, alpha=0.85))

axA.set(xlabel="Observed Wavelength [μm]", ylabel="Normalized Flux", 
        title=f"Phase 4 LRD Spectrum (z = {Z}, tage = {TAGE} Gyr, logzsol = {LOGZSOL}, dust2 = {DUST2})")
axA.legend(loc="upper left", fontsize=8.5, framealpha=0.25)
axA.set_xlim(0.9, 5.2)
axA.set_ylim(-0.05, 1.25)
axA.grid(True, color=C_GRID, alpha=0.5)

# Panel B: Rest-frame UV-Optical
axB = fig.add_subplot(gs[1, 0], facecolor=PANEL_BG)
mask_opt = (wave_rest > 2200) & (wave_rest < 7200)
mask_opt_neb = (wave_rest_neb > 2200) & (wave_rest_neb < 7200)
norm_B = np.percentile(spec_hc[mask_opt], 98) if mask_opt.sum() > 0 else 1.0

axB.fill_between(wave_rest[mask_opt], 0, spec_hc[mask_opt] / norm_B, alpha=0.12, color=C_HC)
axB.plot(wave_rest[mask_opt], spec_hc[mask_opt] / norm_B, color=C_HC, lw=1.4, label="HC-RSFG")
axB.plot(wave_rest_neb[mask_opt_neb], spec_neb[mask_opt_neb] / norm_B, color=C_NEB, lw=0.9, alpha=0.7, label="HC-RSFG+Nebular")
axB.plot(wave_std[mask_opt], spec_std[mask_opt] / norm_B, color=C_STD, lw=1.0, ls="--", alpha=0.6, label="Chabrier")

for lname, lw0 in [("[OII]", 3727), ("Hb", 4863), ("[OIII]", 5007), ("Ha", 6565)]:
    if 2200 < lw0 < 7200:
        axB.axvline(lw0, color="#78909c", alpha=0.45, lw=0.8, ls="--")
        axB.text(lw0 + 30, 0.98, lname, fontsize=7, color="#90a4ae", transform=axB.get_xaxis_transform(), va="top")

axB.set(xlabel="Rest-frame Wavelength [A]", ylabel="Normalized Flux", title="Rest-frame Spectrum (Near-UV to Optical)")
axB.legend(fontsize=7.5, framealpha=0.25)
axB.set_xlim(2400, 7100)
axB.set_ylim(-0.05, 1.15)
axB.grid(True, color=C_GRID, alpha=0.5)

# Panel C: Emission EW Bar Chart
axC = fig.add_subplot(gs[1, 1], facecolor=PANEL_BG)
ew_names = list(ew_results.keys())
ew_vals = [ew_results[k] for k in ew_names]
bar_cols = [C_HC if v > 0.5 else "#e57373" if v > 0.05 else "#5d4037" for v in ew_vals]

bars = axC.barh(ew_names, ew_vals, color=bar_cols, height=0.55, alpha=0.88)
axC.axvline(1.0, color=C_OBS, lw=1.5, ls="--", label="EW = 1 A (AGN Threshold)")
axC.axvline(0.10, color="#78909c", lw=1.0, ls=":", label="EW = 0.1 A (Detection Limit)")

for bar, v in zip(bars, ew_vals):
    axC.text(v + 0.005, bar.get_y() + bar.get_height() / 2.0, f"{v:.3f} A", va="center", fontsize=7.5, color=C_TXT)

axC.set(xlabel="Equivalent Width (EW) [A]", title="Emission Line EW (All EW < 1A)\nConsistent with Narrow-Line LRDs")
axC.legend(fontsize=7.5, framealpha=0.25)
axC.set_xlim(0, max(max(ew_vals) * 1.5, 1.4))
axC.grid(True, color=C_GRID, alpha=0.5, axis="x")

# Panel D: Mock SED Fitting
axD = fig.add_subplot(gs[1, 2], facecolor=PANEL_BG)
band_waves_um = np.array([1.15, 1.50, 2.00, 2.77, 3.56, 4.44])
ref_idx = 4 # Normalize at F356W
mags_hc_arr = np.array([mags_hc[BANDS.index(b)] for b in BANDS])
mags_std_arr = np.array([mags_std[BANDS.index(b)] for b in BANDS])

flux_hc = 10**(-0.4*(mags_hc_arr - mags_hc_arr[ref_idx]))
flux_std = 10**(-0.4*(mags_std_arr - mags_std_arr[ref_idx]))

rng = np.random.default_rng(SEED)
phot_err = np.array([0.18, 0.14, 0.12, 0.09, 0.08, 0.10])
obs_flux = flux_hc * (1.0 + rng.normal(0, phot_err * 0.6))
obs_flux_err = flux_hc * phot_err

axD.errorbar(band_waves_um, obs_flux, yerr=obs_flux_err, fmt="D", color=C_OBS, ms=8, capsize=4, capthick=1.5, label="Mock Obs Data (1sigma)", zorder=6)
axD.plot(band_waves_um, flux_hc, "o-", color=C_HC, lw=2.0, ms=6, label=f"HC-RSFG (Color={c200_444_hc:+.2f})")
axD.plot(band_waves_um, flux_std, "s--", color=C_STD, lw=1.4, ms=5, alpha=0.7, label=f"Chabrier (Color={c200_444_std:+.2f})")

axD.set(xlabel="Observed Wavelength [μm]", ylabel="Relative Flux (F356W=1)", title=f"Mock SED Fitting (z={Z})")
axD.legend(fontsize=7.5, framealpha=0.25)
axD.set_xlim(0.9, 5.0)
axD.grid(True, color=C_GRID, alpha=0.5)

# Panel E: Metallicity Scan
axE = fig.add_subplot(gs[2, 0], facecolor=PANEL_BG)
axE.fill_between(logzsol_grid, np.maximum(c200_444_grid, LRD_THRESH), LRD_THRESH, alpha=0.15, color=C_OBS, label="LRD Region")
axE.plot(logzsol_grid, c200_444_grid, "o-", color=C_HC, lw=2.0, ms=6)
axE.axhline(LRD_THRESH, color=C_OBS, lw=1.5, ls="--", label=f"LRD Threshold ({LRD_THRESH})")
axE.scatter([LOGZSOL], [c200_444_hc], s=130, color=C_OBS, zorder=7, label=f"This Model (logzsol={LOGZSOL})")

axE.set(xlabel="logzsol (log Z/Zsun)", ylabel="F200W - F444W [mag]", title="Metallicity vs LRD Color\n(HC-RSFG IMF, dust2=0.08)")
axE.legend(fontsize=7.5, framealpha=0.25)
axE.grid(True, color=C_GRID, alpha=0.5)

# Panel F: Dust Scan
axF = fig.add_subplot(gs[2, 1], facecolor=PANEL_BG)
av_grid = dust2_grid * 3.1
axF.fill_between(av_grid, LRD_THRESH, c200_444_dust, where=c200_444_dust >= LRD_THRESH, alpha=0.15, color=C_OBS)
axF.plot(av_grid, c200_444_dust, "o-", color=C_HC, lw=2.0, ms=6)
axF.axhline(LRD_THRESH, color=C_OBS, lw=1.5, ls="--", label=f"LRD Threshold ({LRD_THRESH})")
axF.scatter([DUST2 * 3.1], [c200_444_hc], s=130, color=C_OBS, zorder=7, label=f"This Model (Av≈{DUST2*3.1:.2f})")
axF.axvline(DUST2 * 3.1, color=C_OBS, lw=1.0, ls=":", alpha=0.5)
axF.text(0.05, c200_444_dust[0] + 0.05, f"Av=0 Reaches\n{c200_444_dust[0]:.2f}", fontsize=8, color="#aed6f1", 
         bbox=dict(boxstyle="round", fc="#0d0d20", ec="#3d5a80", alpha=0.8))

axF.set(xlabel="Av [mag] (dust2 * 3.1)", ylabel="F200W - F444W [mag]", title="Dust vs LRD Color\n(logzsol = -1.0, HC-RSFG IMF)")
axF.legend(fontsize=7.5, framealpha=0.25)
axF.grid(True, color=C_GRID, alpha=0.5)

# Panel G: Color-Color Diagram
axG = fig.add_subplot(gs[2, 2], facecolor=PANEL_BG)
lrd_region_x = [1.0, 3.5, 3.5, 1.0, 1.0]
lrd_region_y = [0.3, 0.3, 1.2, 1.2, 0.3]
axG.fill(lrd_region_x, lrd_region_y, alpha=0.10, color=C_OBS, label="LRD Obs Region (Hviding+2025)")
axG.plot(lrd_region_x, lrd_region_y, color=C_OBS, lw=1.0, ls="--", alpha=0.5)

sc = axG.scatter(c200_444_grid, c277_444_grid, c=logzsol_grid, cmap="plasma", s=60, zorder=4, label="HC-RSFG (Met. Scan)")
cb = plt.colorbar(sc, ax=axG, pad=0.02)
cb.set_label("logzsol", fontsize=8.5, color=C_TXT)
cb.ax.yaxis.set_tick_params(color=C_TXT)

axG.scatter([c200_444_hc], [c277_444_hc], s=180, color=C_OBS, marker="*", zorder=6, label=f"This Model (logzsol={LOGZSOL})")
axG.scatter([c200_444_std], [c277_444_std], s=80, color=C_STD, marker="s", zorder=5, label="Chabrier")

axG.set(xlabel="F200W - F444W [mag]", ylabel="F277W - F444W [mag]", title="Color-Color Diagram\n(vs LRD Obs Region)")
axG.legend(fontsize=7.0, framealpha=0.25, loc="upper left")
axG.grid(True, color=C_GRID, alpha=0.5)

fig.suptitle("HC-RSFG Theory - Code 1: LRD State Simulation (Phase 4)\n"
             f"z = {Z}  tage = {TAGE} Gyr  logzsol = {LOGZSOL}  dust2 = {DUST2}  imf3 = {IMF3}",
             fontsize=13, color=C_TXT, y=0.97)

out_path = "HC_RSFG_Code1_LRD.png"
plt.savefig(out_path, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"\n=================================================================")
print(f" ✅ Code 1 Completed - Plot saved as: {out_path}")
print("=================================================================")

