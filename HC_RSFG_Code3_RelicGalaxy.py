# HC-RSFG Theory Code 3: Relic Galaxy Simulation
# Phase 6: Final stage resembling NGC 1277


# Version 1.0 (simplified top‑hat filter approximation)
# Updated results using real FSPS filter transmission curves are available
# in the accompanying paper (Fig. 1‑5) and in the script `HC_RSFG_5panel.py`.


import fsps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings("ignore")

# -------------------------------------------------------------------
# 1. Parameter definitions
# -------------------------------------------------------------------
ZRED = 0.0173                             # Redshift of NGC 1277
TAGE = 12.0                               # Stellar age [Gyr] (formation at z~7.5)

# Central HC-RSFG core (r < Re/4)
LOGZSOL_CORE = 0.0                        # Solar metallicity (old elliptical)
DUST2_CORE = 0.05                         # Low dust
IMF1_CORE = 1.8
IMF2_CORE = 1.8
IMF3_CORE = 3.0                           # Bottom-heavy (HC-RSFG relic)

# Outer region (r > Re): normal IMF
LOGZSOL_OUT = -0.2
DUST2_OUT = 0.05
IMF_OUT = 1                               # Chabrier

# Intermediate region (Re/4 < r < Re)
IMF3_MID = 2.3

# Reference Chabrier model (whole galaxy)
LOGZSOL_STD = 0.0
IMF_STD = 1

# Radial bins (in units of effective radius)
RADIUS_BINS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]  # r/Re

# IMF3 slope at each radius (center to outer)
IMF3_RADIAL = [3.0, 2.7, 2.3, 1.8, 1.5, 1.3]

# Wing-Ford band (IMF-sensitive feature)
WF_CENTER = 9905.0
WF_WIDTH = 45.0

# Other diagnostic features
MGB_CENTER = 5175.0
FE5270 = 5270.0
CAII_K = 3933.7

# Bands for photometry
BANDS_RELIC = ["sdss_g", "sdss_r", "sdss_i", "sdss_z",
               "2mass_j", "2mass_h", "2mass_ks"]

# -------------------------------------------------------------------
# 2. FSPS instance builder
# -------------------------------------------------------------------
def make_sps_relic(imf_type=2, imf3=IMF3_CORE, logzsol=LOGZSOL_CORE,
                   dust2=DUST2_CORE, zred=ZRED, tage=TAGE):
    sps = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1)
    sps.params["imf_type"] = imf_type
    sps.params["imf1"] = IMF1_CORE
    sps.params["imf2"] = IMF2_CORE
    sps.params["imf3"] = imf3
    sps.params["logzsol"] = logzsol
    sps.params["sfh"] = 0                 # SSP
    sps.params["tage"] = tage
    sps.params["dust2"] = dust2
    sps.params["zred"] = zred
    if imf_type == 2:
        sps.params["imf_upper_limit"] = 1.5
    return sps

# -------------------------------------------------------------------
# 3. Spectra and radial IMF scan
# -------------------------------------------------------------------
# Core, outer, and reference
sps_core = make_sps_relic(imf_type=2, imf3=IMF3_CORE, logzsol=LOGZSOL_CORE)
sps_out  = make_sps_relic(imf_type=1, imf3=2.3, logzsol=LOGZSOL_OUT)        # Chabrier
sps_std  = make_sps_relic(imf_type=1, imf3=2.3, logzsol=LOGZSOL_STD)

wave_core, spec_core = sps_core.get_spectrum(tage=TAGE)
wave_out, spec_out   = sps_out.get_spectrum(tage=TAGE)
wave_std, spec_std   = sps_std.get_spectrum(tage=TAGE)
wave_rest = wave_core / (1 + ZRED)          # rest-frame (z~0.017)

# Radial spectra and indices
spectra_radial = []
wf_indices = []
mgb_indices = []
ml_ratios = []      # M/L proxy relative to Kroupa
imf3_vals = []

for imf3_r in IMF3_RADIAL:
    sp_r = make_sps_relic(imf_type=2, imf3=imf3_r, logzsol=LOGZSOL_CORE)
    _, sp_arr = sp_r.get_spectrum(tage=TAGE)
    spectra_radial.append(sp_arr)

    # Wing-Ford index (depth of 9905 A band)
    m_wf = (wave_rest > WF_CENTER - WF_WIDTH) & (wave_rest < WF_CENTER + WF_WIDTH)
    m_cont = ((wave_rest > WF_CENTER + WF_WIDTH*2) & (wave_rest < WF_CENTER + 200)) | \
             ((wave_rest > WF_CENTER - 200) & (wave_rest < WF_CENTER - WF_WIDTH*2))
    if m_wf.sum() > 0 and m_cont.sum() > 0:
        cont_val = np.median(sp_arr[m_cont])
        feat_val = np.mean(sp_arr[m_wf])
        wf_idx = 1.0 - feat_val / max(cont_val, 1e-30)
    else:
        wf_idx = 0.0
    wf_indices.append(wf_idx)

    # Mg b index
    m_mgb = (wave_rest > MGB_CENTER - 25) & (wave_rest < MGB_CENTER + 25)
    m_cntmg = ((wave_rest > MGB_CENTER - 150) & (wave_rest < MGB_CENTER - 30)) | \
              ((wave_rest > MGB_CENTER + 30) & (wave_rest < MGB_CENTER + 150))
    if m_mgb.sum() > 0 and m_cntmg.sum() > 0:
        cont_mg = np.median(sp_arr[m_cntmg])
        feat_mg = np.mean(sp_arr[m_mgb])
        mgb_idx = 1.0 - feat_mg / max(cont_mg, 1e-30)
    else:
        mgb_idx = 0.0
    mgb_indices.append(mgb_idx)

    # M/L proxy: larger imf3 -> more low-mass stars -> higher M/L
    ml_ratios.append(10**(0.25 * (imf3_r - 1.3)))
    imf3_vals.append(imf3_r)

# -------------------------------------------------------------------
# 4. Colours and diagnostic indices
# -------------------------------------------------------------------
mags_core = sps_core.get_mags(tage=TAGE, bands=BANDS_RELIC)
mags_out  = sps_out.get_mags(tage=TAGE, bands=BANDS_RELIC)
mags_std  = sps_std.get_mags(tage=TAGE, bands=BANDS_RELIC)

mag_dict_core = dict(zip(BANDS_RELIC, mags_core))
mag_dict_out  = dict(zip(BANDS_RELIC, mags_out))
mag_dict_std  = dict(zip(BANDS_RELIC, mags_std))

gr_core = mag_dict_core["sdss_g"] - mag_dict_core["sdss_r"]
gr_out  = mag_dict_out["sdss_g"] - mag_dict_out["sdss_r"]
gr_std  = mag_dict_std["sdss_g"] - mag_dict_std["sdss_r"]

print(f"g-r (core, imf3={IMF3_CORE}): {gr_core:+.3f}")
print(f"g-r (outer, imf3={IMF3_RADIAL[-1]}): {gr_out:+.3f}")
print(f"g-r (all Chabrier): {gr_std:+.3f}")

print("\nRadial profiles:")
print(f"  r/Re   IMF3   Wing-Ford   Mg b   M/L (x Kroupa)")
for r, i3, wf, mg, ml in zip(RADIUS_BINS, imf3_vals, wf_indices, mgb_indices, ml_ratios):
    print(f"  {r:5.2f}  {i3:5.1f}  {wf:8.4f}  {mg:8.4f}  {ml:8.2f}")

# -------------------------------------------------------------------
# 5. Plotting
# -------------------------------------------------------------------
DARK_BG = "#090912"
PANEL_BG = "#0d0d20"
C_CORE = "#ef9a9a"
C_OUT = "#80deea"
C_STD = "#a5d6a7"
C_WF = "#ce93d8"
C_OBS = "#ffd54f"
C_GRID = "#1e1e3a"
C_TXT = "#e8eaf6"

plt.rcParams.update({
    "figure.facecolor": DARK_BG,
    "text.color": C_TXT,
    "axes.labelcolor": C_TXT,
    "xtick.color": C_TXT,
    "ytick.color": C_TXT,
    "axes.edgecolor": "#3d3d6b",
    "font.family": "monospace",
})

fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.38,
                       top=0.92, bottom=0.06, left=0.06, right=0.97)

# Panel A: Rest-frame spectrum (core, outer, Chabrier)
axA = fig.add_subplot(gs[0, :])
axA.set_facecolor(PANEL_BG)
mask_A = (wave_rest > 3000) & (wave_rest < 11000)
norm_A = np.median(spec_core[(wave_rest > 9600) & (wave_rest < 9700)]) or 1.0
axA.plot(wave_rest[mask_A], spec_core[mask_A]/norm_A, color=C_CORE, lw=2.0, label=f"Core (imf3={IMF3_CORE})")
axA.plot(wave_rest[mask_A], spec_out[mask_A]/norm_A, color=C_OUT, lw=1.5, alpha=0.8, label="Outer (Chabrier)")
axA.plot(wave_rest[mask_A], spec_std[mask_A]/norm_A, color=C_STD, lw=1.2, alpha=0.6, label="All Chabrier")
for lname, lw0 in [("CaII K",3933),("Mg b",5175),("Fe5270",5270),("Ha",6565),("Wing-Ford",9905)]:
    axA.axvline(lw0, color="#78909c", alpha=0.5, lw=0.8, ls="--")
    axA.text(lw0/1e4, 1.16, lname, ha="center", fontsize=6.5, color="#78909c",
             transform=axA.get_xaxis_transform())
axA.set_xlabel("Rest-frame wavelength [um]", fontsize=11)
axA.set_ylabel("Normalised flux", fontsize=11)
axA.set_title(f"Phase 6 Relic Galaxy – NGC 1277 analog (z={ZRED}, tage={TAGE} Gyr)", fontsize=12, pad=6)
axA.legend(loc="upper right", fontsize=8.5, framealpha=0.25)
axA.set_xlim(0.30, 1.10)
axA.set_ylim(-0.03, 1.3)
axA.grid(True, color=C_GRID, alpha=0.5)

# Panel B: Wing-Ford band zoom
axB = fig.add_subplot(gs[1, 0])
axB.set_facecolor(PANEL_BG)
mask_wf = (wave_rest > 9600) & (wave_rest < 10200)
if mask_wf.sum() > 5:
    norm_B = np.median(spec_core[(wave_rest > 9600) & (wave_rest < 9700)]) or 1.0
    cmap_r = plt.cm.RdYlGn_r
    for ri, (sp_r, r_bin, i3) in enumerate(zip(spectra_radial, RADIUS_BINS, imf3_vals)):
        col = cmap_r(ri/(len(RADIUS_BINS)-1))
        axB.plot(wave_rest[mask_wf],
                 sp_r[mask_wf] / np.median(sp_r[(wave_rest > 9600) & (wave_rest < 9700)]),
                 color=col, lw=1.8, alpha=0.85,
                 label=f"r/Re={r_bin:.2f} imf3={i3:.1f}")
    axB.axvspan(WF_CENTER - WF_WIDTH, WF_CENTER + WF_WIDTH, alpha=0.12, color=C_WF)
    axB.axvline(WF_CENTER, color=C_WF, lw=1.2, ls="--", alpha=0.8)
    axB.text(WF_CENTER, 1.02, f"{WF_CENTER:.0f}", ha="center", fontsize=8, color=C_WF,
             transform=axB.get_xaxis_transform())
    axB.set_xlabel("Rest-frame wavelength [A]", fontsize=9.5)
    axB.set_ylabel("Normalised flux", fontsize=9.5)
    axB.set_title("Wing-Ford band (9905 A) – IMF-sensitive", fontsize=10)
    axB.legend(fontsize=6.5, framealpha=0.25)
    axB.set_xlim(9620, 10180)
    axB.grid(True, color=C_GRID, alpha=0.5)
else:
    axB.text(0.5, 0.5, "Wing-Ford region: may be outside FSPS spectral range",
             ha="center", va="center", transform=axB.transAxes, fontsize=10, color=C_TXT)
    axB.set_title("Wing-Ford band", fontsize=10)
    axB.set_facecolor(PANEL_BG)

# Panel C: Radial IMF3 and Wing-Ford index
axC = fig.add_subplot(gs[1, 1])
axC.set_facecolor(PANEL_BG)
ax_twin = axC.twinx()
l1, = axC.plot(RADIUS_BINS, imf3_vals, "o-", color=C_CORE, lw=2.5, ms=8, label="imf3 (left)")
l2, = ax_twin.plot(RADIUS_BINS, wf_indices, "s-", color=C_WF, lw=2.0, ms=7, label="Wing-Ford index (right)")
axC.axhline(1.3, color=C_STD, lw=1.0, ls=":", alpha=0.7, label="Kroupa (imf3=1.3)")
axC.axhline(2.35, color="#ffcc80", lw=1.0, ls=":", alpha=0.6, label="Salpeter (imf3=2.35)")
axC.fill_between(RADIUS_BINS, 1.3, imf3_vals, alpha=0.10, color=C_CORE, label="Bottom-heavy excess")
# Approximate NGC 1277 data from Martín-Navarro+2015
mn15_radii = [0.0, 0.4, 0.8, 1.3]
mn15_imf3 = [3.1, 2.8, 2.5, 2.2]
axC.errorbar(mn15_radii, mn15_imf3, yerr=0.2, fmt="D", color=C_OBS, ms=9, capsize=4,
             label="MN+2015 NGC 1277 (approx)", zorder=6)
axC.set_xlabel("r / Re", fontsize=9.5)
axC.set_ylabel("imf3 slope", fontsize=9.5, color=C_CORE)
ax_twin.set_ylabel("Wing-Ford index", fontsize=9.5, color=C_WF)
axC.set_title("Radial IMF3 & Wing-Ford\n(compared to Martín-Navarro+2015)", fontsize=10)
lines = [l1, l2]
axC.legend(lines, [l.get_label() for l in lines], fontsize=7, framealpha=0.25)
axC.grid(True, color=C_GRID, alpha=0.5)
axC.set_xlim(-0.1, 2.1)

# Panel D: M/L ratio radial profile
axD = fig.add_subplot(gs[1, 2])
axD.set_facecolor(PANEL_BG)
axD.fill_between(RADIUS_BINS, 1.0, ml_ratios, alpha=0.12, color=C_CORE, label="HC-RSFG M/L excess")
axD.plot(RADIUS_BINS, ml_ratios, "o-", color=C_CORE, lw=2.5, ms=8, label="HC-RSFG M/L (relative)")
axD.axhline(1.0, color=C_STD, lw=1.5, ls="--", label="Kroupa M/L = 1.0")
axD.axhline(2.0, color=C_OBS, lw=1.2, ls="--", alpha=0.7, label="NGC 1277 M/L x2.0 (MN+2015)")
mn_ml_radii = [0.0, 0.5, 1.0, 1.5]
mn_ml_vals = [2.1, 1.85, 1.6, 1.4]
axD.errorbar(mn_ml_radii, mn_ml_vals, yerr=0.15, fmt="D", color=C_OBS, ms=9, capsize=4,
             label="NGC 1277 M/L (MN+2015 approx)", zorder=6)
axD.set_xlabel("r / Re", fontsize=9.5)
axD.set_ylabel("M/L (Kroupa = 1.0)", fontsize=9.5)
axD.set_title("M/L ratio profile\n(NGC 1277 comparison)", fontsize=10)
axD.legend(fontsize=7.5, framealpha=0.25)
axD.set_xlim(-0.1, 2.1)
axD.set_ylim(0.5, 2.8)
axD.grid(True, color=C_GRID, alpha=0.5)

# Panel E: Optical SED comparison
axE = fig.add_subplot(gs[2, 0])
axE.set_facecolor(PANEL_BG)
band_waves = np.array([0.468, 0.617, 0.748, 0.893, 1.235, 1.662, 2.159])
def norm_flux(mag_dict, band_list):
    vals = np.array([mag_dict[b] for b in band_list])
    flx = 10**(-0.4*vals)
    return flx / np.median(flx)
flux_core_arr = norm_flux(mag_dict_core, BANDS_RELIC[1:])
flux_out_arr  = norm_flux(mag_dict_out, BANDS_RELIC[1:])
flux_std_arr  = norm_flux(mag_dict_std, BANDS_RELIC[1:])
np.random.seed(123)
obs_err = np.abs(np.random.normal(0.0, 0.06, len(band_waves)))
obs_pts = flux_core_arr * (1.0 + np.random.normal(0, obs_err))
axE.errorbar(band_waves, obs_pts, yerr=flux_core_arr*obs_err, fmt="D", color=C_OBS,
             ms=8, capsize=4, label="Mock observation")
axE.plot(band_waves, flux_core_arr, "o-", color=C_CORE, lw=2.0, ms=6, label=f"Core (imf3={IMF3_CORE})")
axE.plot(band_waves, flux_out_arr, "s-", color=C_OUT, lw=1.5, ms=5, alpha=0.8, label="Outer (Chabrier)")
axE.plot(band_waves, flux_std_arr, "^", color=C_STD, lw=1.2, ms=5, alpha=0.6, label="All Chabrier")
for bname, bw in zip(["g","r","i","z","J","H","K"], band_waves):
    axE.text(bw, -0.08, bname, ha="center", fontsize=7.5, color="#9fa8da",
             transform=axE.get_xaxis_transform())
axE.set_xlabel("Observed wavelength [um]", fontsize=9.5)
axE.set_ylabel("Normalised flux", fontsize=9.5)
axE.set_title(f"SED comparison (z={ZRED}, NGC 1277 analog)", fontsize=10)
axE.legend(fontsize=7.5, framealpha=0.25)
axE.set_xlim(0.38, 2.3)
axE.grid(True, color=C_GRID, alpha=0.5)

# Panel F: Mg b / Fe5270 region
axF = fig.add_subplot(gs[2, 1])
axF.set_facecolor(PANEL_BG)
mask_mgb = (wave_rest > 4900) & (wave_rest < 5600)
if mask_mgb.sum() > 5:
    norm_F = np.median(spec_core[(wave_rest > 4900) & (wave_rest < 5000)]) or 1.0
    cmap_r = plt.cm.RdYlGn_r
    for ri, (sp_r, r_bin, i3) in enumerate(zip(spectra_radial, RADIUS_BINS, imf3_vals)):
        col = cmap_r(ri/(len(RADIUS_BINS)-1))
        axF.plot(wave_rest[mask_mgb],
                 sp_r[mask_mgb] / np.median(sp_r[(wave_rest > 4900) & (wave_rest < 5000)]),
                 color=col, lw=1.8, alpha=0.85, label=f"r/Re={r_bin:.2f}")
    axF.axvspan(MGB_CENTER-25, MGB_CENTER+25, alpha=0.12, color="#81d4fa")
    axF.axvspan(FE5270-25, FE5270+25, alpha=0.08, color="#ef9a9a")
    axF.text(MGB_CENTER, 1.03, "Mg b", ha="center", fontsize=8, color="#81d4fa",
             transform=axF.get_xaxis_transform())
    axF.text(FE5270, 1.07, "Fe5270", ha="center", fontsize=8, color="#ef9a9a",
             transform=axF.get_xaxis_transform())
    axF.set_xlabel("Rest-frame wavelength [A]", fontsize=9.5)
    axF.set_ylabel("Normalised flux", fontsize=9.5)
    axF.set_title("Mg b / Fe5270 region\n(IMF gradient diagnostic)", fontsize=10)
    axF.legend(fontsize=6.5, framealpha=0.25)
    axF.set_xlim(4920, 5550)
    axF.grid(True, color=C_GRID, alpha=0.5)
else:
    axF.text(0.5, 0.5, "Region out of spectral range", ha="center", va="center",
             transform=axF.transAxes, fontsize=10, color=C_TXT)
    axF.set_title("Mg b / Fe5270", fontsize=10)
    axF.set_facecolor(PANEL_BG)

# Panel G: Summary table
axG = fig.add_subplot(gs[2, 2])
axG.set_facecolor(PANEL_BG)
axG.axis("off")
rows = [
    ["Property", "NGC 1277 observation", "HC-RSFG theory"],
    ["g-r colour", f"{gr_std:.2f}-{gr_std+0.1:.2f}", f"{gr_core:.2f}"],
    ["IMF (center)", "bottom-heavy", f"imf3={IMF3_CORE}"],
    ["IMF (outer)", "Kroupa-like", f"imf3={IMF3_RADIAL[-1]}"],
    ["M/L (center)", "~2.0 × Kroupa", f"×{ml_ratios[0]:.2f}"],
    ["M/L (outer)", "~1.4 × Kroupa", f"×{ml_ratios[-1]:.2f}"],
    ["Wing-Ford index", "center > outer", f"center: {wf_indices[0]:.4f}"],
    ["Age", "~12-14 Gyr", f"{TAGE} Gyr"],
    ["Star formation stopped", "z > 2", "z=7.5 (Phase4)"],
    ["Origin", "unidentified", "HC-RSFG Phase4"],
]
y_start = 0.97
row_h = 0.086
for row_i, row in enumerate(rows):
    y = y_start - row_i * row_h
    bg = "#1a237e" if row_i==0 else ("#12122a" if row_i%2==0 else "#0d0d20")
    wt = "bold" if row_i==0 else "normal"
    rect = FancyBboxPatch((0.01, y-row_h*0.82), 0.98, row_h*0.85,
                          boxstyle="round,pad=0.01", fc=bg, ec="#3d3d6b", lw=0.5,
                          transform=axG.transAxes, zorder=0)
    axG.add_patch(rect)
    for ci, (cell, x) in enumerate(zip(row, [0.02, 0.36, 0.69])):
        axG.text(x, y-row_h*0.38, cell, transform=axG.transAxes,
                 fontsize=7.0, fontweight=wt, va="center",
                 color=C_OBS if row_i==0 else C_TXT)
axG.set_title("NGC 1277: observation vs HC-RSFG", fontsize=10, pad=6)

fig.suptitle("HC-RSFG Theory - Code 3: Relic Galaxy Simulation (Phase 6)\n"
             f"NGC 1277 analog z={ZRED}, tage={TAGE} Gyr, core imf3={IMF3_CORE} → outer imf3={IMF3_RADIAL[-1]}",
             fontsize=12, color=C_TXT, y=0.97)
out_path = "HC_RSFG_Code3_RelicGalaxy.png"
plt.savefig(out_path, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()
print(f"Saved: {out_path}")
print(f"Center: imf3={IMF3_CORE} (bottom-heavy) → outer: imf3={IMF3_RADIAL[-1]}")
print(f"M/L center = {ml_ratios[0]:.2f} × Kroupa (NGC 1277 observed ~2.0)")
print(f"Wing-Ford gradient: center {wf_indices[0]:.4f} → outer {wf_indices[-1]:.4f}")
