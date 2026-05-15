# HC-RSFG Theory Code 4: Dark-RSFG Simulation
# Phase 2-3: Early stage (z ~ 7.5-10, tage < 0.5 Gyr)
# Predicts extremely faint objects – currently undetectable but may be found with future deep JWST surveys.


# Version 1.0 (simplified top‑hat filter approximation)
# Updated results using real FSPS filter transmission curves are available
# in the accompanying paper (Fig. 1‑5) and in the script `HC_RSFG_5panel.py`.


import fsps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Circle
import warnings
warnings.filterwarnings("ignore")

# -------------------------------------------------------------------
# 1. Parameter definitions
# -------------------------------------------------------------------
# Dark-RSFG phase: age grid and corresponding redshifts
TAGE_STAGES = [0.05, 0.1, 0.2, 0.3, 0.5]          # Gyr
ZRED_STAGES = [10.4, 9.5, 8.3, 7.9, 7.3]          # Approximate redshifts

# Common parameters for Dark-RSFG
LOGZSOL_DARK = -2.0              # Metallicity Z ~ 0.01 Zsun (early low metallicity)
DUST2_DARK = 0.00                # No dust in Phase 2
TAU_DARK = 0.05                  # Very short SFH (rapid core formation)
IMF1_DARK = 1.8
IMF2_DARK = 1.8
IMF3_DARK = 3.5                  # Extremely bottom-heavy (M_J ~ 1.0 Msun)

# Phase 4 LRD for comparison (already discovered)
LOGZSOL_LRD = -1.0
DUST2_LRD = 0.08
TAGE_LRD = 0.9
ZRED_LRD = 6.0
IMF3_LRD = 3.0

# JWST sensitivity (typical deep survey limits)
JWST_LIMIT_10HR = 5.0            # nJy (5σ) for 10-hour exposure
JWST_LIMIT_100HR = 1.5           # nJy (5σ) for 100-hour exposure

# Assumed core mass for flux scaling
LOG_MASS_CORE = 8.0              # log(M*/Msun)

# -------------------------------------------------------------------
# 2. FSPS instances for Dark-RSFG stages and LRD
# -------------------------------------------------------------------
def make_dark_sps(tage, zred, logzsol=LOGZSOL_DARK, dust2=DUST2_DARK, imf3=IMF3_DARK):
    sps = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1,
                                 add_neb_emission=True)
    sps.params["imf_type"] = 2          # Kroupa broken power law
    sps.params["imf1"] = IMF1_DARK
    sps.params["imf2"] = IMF2_DARK
    sps.params["imf3"] = imf3
    sps.params["imf_upper_limit"] = 1.5 # Effective upper mass cutoff
    sps.params["logzsol"] = logzsol
    sps.params["sfh"] = 0               # SSP
    sps.params["tage"] = tage
    sps.params["dust2"] = dust2
    sps.params["zred"] = zred
    sps.params["gas_logz"] = logzsol
    sps.params["gas_logu"] = -4.0
    return sps

# Create sps for each Dark-RSFG stage
sps_stages = []
for tage, zred in zip(TAGE_STAGES, ZRED_STAGES):
    print(f"Computing: tage={tage} Gyr, z={zred}...", flush=True)
    sps_stages.append(make_dark_sps(tage, zred))

# Phase 4 LRD for comparison
sps_lrd = make_dark_sps(tage=TAGE_LRD, zred=ZRED_LRD,
                        logzsol=LOGZSOL_LRD, dust2=DUST2_LRD, imf3=IMF3_LRD)

# -------------------------------------------------------------------
# 3. Spectra and broadband colours
# -------------------------------------------------------------------
BANDS_NIR = ["jwst_f115w", "jwst_f150w", "jwst_f200w", "jwst_f277w", "jwst_f356w", "jwst_f444w"]

# Extract spectra and store
wave_stages = []
spec_stages = []
for sp in sps_stages:
    w, s = sp.get_spectrum(tage=sp.params["tage"])
    wave_stages.append(w)
    spec_stages.append(s)

w_lrd, s_lrd = sps_lrd.get_spectrum(tage=TAGE_LRD)

# Colours over time
def colour(mag_dict, b1, b2):
    return mag_dict[b1] - mag_dict[b2]

colors_200_444 = []
colors_150_277 = []
mab_f277w = []          # apparent F277W mag (flux proxy)

for i, (sp, tage, zred) in enumerate(zip(sps_stages, TAGE_STAGES, ZRED_STAGES)):
    m = sp.get_mags(tage=tage, bands=BANDS_NIR)
    md = dict(zip(BANDS_NIR, m))
    colors_200_444.append(colour(md, "jwst_f200w", "jwst_f444w"))
    colors_150_277.append(colour(md, "jwst_f150w", "jwst_f277w"))
    mab_f277w.append(md["jwst_f277w"])

m_lrd = sps_lrd.get_mags(tage=TAGE_LRD, bands=BANDS_NIR)
md_lrd = dict(zip(BANDS_NIR, m_lrd))
c200_lrd = colour(md_lrd, "jwst_f200w", "jwst_f444w")

# -------------------------------------------------------------------
# 4. Emission line EW evolution
# -------------------------------------------------------------------
HA_REST = 6564.6
OIII_REST = 5007.0
OII_REST = 3727.1

def compute_ew_from_rest(wave_obs, spec, w0_rest, zred, dw=30, dw_cont=120):
    wave_rest = wave_obs / (1 + zred)
    m_line = (wave_rest > w0_rest - dw) & (wave_rest < w0_rest + dw)
    m_cont = (((wave_rest > w0_rest - dw_cont) & (wave_rest < w0_rest - dw*2)) |
              ((wave_rest > w0_rest + dw*2) & (wave_rest < w0_rest + dw_cont)))
    if m_line.sum() == 0 or m_cont.sum() == 0:
        return 0.0
    peak = np.max(spec[m_line])
    cont = np.median(spec[m_cont])
    return max((peak - cont) / cont * (2 * dw), 0.0) if cont > 0 else 0.0

ew_ha_stages = []
ew_oiii_stages = []
ew_oii_stages = []

for wave, spec, tage, zred in zip(wave_stages, spec_stages, TAGE_STAGES, ZRED_STAGES):
    ew_ha_stages.append(compute_ew_from_rest(wave, spec, HA_REST, zred))
    ew_oiii_stages.append(compute_ew_from_rest(wave, spec, OIII_REST, zred))
    ew_oii_stages.append(compute_ew_from_rest(wave, spec, OII_REST, zred))

ew_ha_lrd = compute_ew_from_rest(w_lrd, s_lrd, HA_REST, ZRED_LRD)
ew_oiii_lrd = compute_ew_from_rest(w_lrd, s_lrd, OIII_REST, ZRED_LRD)

print("\nEmission line EW evolution:")
print(f"  tage [Gyr]   z     EW(Ha)    EW([OIII])   EW([OII])")
for tage, zred, eha, eoiii, eoii in zip(TAGE_STAGES, ZRED_STAGES,
                                        ew_ha_stages, ew_oiii_stages, ew_oii_stages):
    print(f"  {tage:7.2f}  {zred:5.1f}  {eha:10.4f}  {eoiii:12.4f}  {eoii:10.4f}")
print(f"  LRD comp   {ZRED_LRD:5.1f}  {ew_ha_lrd:10.4f}  {ew_oiii_lrd:12.4f}")

# -------------------------------------------------------------------
# 5. Flux estimation and JWST detectability
# -------------------------------------------------------------------
def luminosity_distance_mpc(z, H0=67.4, Om=0.315):
    c_kms = 2.998e5
    dz = 0.001
    z_arr = np.arange(0, z+dz, dz)
    integrand = 1.0 / np.sqrt(Om * (1+z_arr)**3 + (1-Om))
    DC = c_kms / H0 * np.trapz(integrand, z_arr)
    return (1+z) * DC

def flux_njy(mag_ab, DL_mpc, log_mass):
    mag_scaled = mag_ab - 2.5 * log_mass
    flux_jy = 10**(-0.4 * (mag_scaled + 48.60))
    return flux_jy * 1e9

fluxes_f277w = []
for tage, zred, mab in zip(TAGE_STAGES, ZRED_STAGES, mab_f277w):
    DL = luminosity_distance_mpc(zred)
    fl = flux_njy(mab, DL, LOG_MASS_CORE)
    fluxes_f277w.append(fl)
    det_10h = "Y" if fl > JWST_LIMIT_10HR else "N"
    det_100h = "Y" if fl > JWST_LIMIT_100HR else "N"
    print(f" z={zred}, tage={tage}: F277W={fl:.2f} nJy  10h:{det_10h}  100h:{det_100h}", flush=True)

DL_LRD = luminosity_distance_mpc(ZRED_LRD)
fl_lrd = flux_njy(md_lrd["jwst_f277w"], DL_LRD, LOG_MASS_CORE)
print(f" z={ZRED_LRD} LRD: F277W={fl_lrd:.2f} nJy")

# -------------------------------------------------------------------
# 6. Plotting
# -------------------------------------------------------------------
DARK_BG = "#090912"
PANEL_BG = "#0d0d20"
C_OBS = "#ffd54f"
C_GRID = "#1e1e3a"
C_TXT = "#e8eaf6"
C_LRD = "#ff8a65"
C_LIMIT = "#ef5350"

STAGE_COLORS = ["#ce93d8", "#9575cd", "#7986cb", "#64b5f6", "#4dd0e1"]

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

# Panel A: Observed-frame spectral time evolution
axA = fig.add_subplot(gs[0, :])
axA.set_facecolor(PANEL_BG)
for i, (wave, spec, tage, zred, col) in enumerate(zip(wave_stages, spec_stages,
                                                      TAGE_STAGES, ZRED_STAGES, STAGE_COLORS)):
    mask = (wave > 8000) & (wave < 55000)
    norm = np.percentile(spec[mask], 95) if mask.sum()>0 else 1.0
    axA.plot(wave[mask]/1e4, spec[mask]/norm, color=col, lw=1.5, alpha=0.85,
             label=f"Dark-RSFG z={zred} tage={tage} Gyr")
mask_lrd = (w_lrd > 8000) & (w_lrd < 55000)
norm_lrd = np.percentile(s_lrd[mask_lrd], 95) or 1.0
axA.plot(w_lrd[mask_lrd]/1e4, s_lrd[mask_lrd]/norm_lrd, color=C_LRD, lw=2.2, ls="--",
         alpha=0.85, label=f"Phase4 LRD z={ZRED_LRD} tage={TAGE_LRD} Gyr (ref)")
for bname, bw in [("F150W",1.50),("F200W",2.00),("F277W",2.77),("F356W",3.56),("F444W",4.44)]:
    axA.axvline(bw, color="#ffffff", alpha=0.08, lw=0.7, ls=".")
    axA.text(bw, 1.12, bname, ha="center", fontsize=7, color="#9fa8da",
             transform=axA.get_xaxis_transform())
axA.set_xlabel("Observed wavelength [um]", fontsize=11)
axA.set_ylabel("Normalised flux", fontsize=11)
axA.set_title(f"Phase 2-3 Dark-RSFG – Spectral evolution\n"
              f"(imf3={IMF3_DARK}, logzsol={LOGZSOL_DARK}, dust2={DUST2_DARK})", fontsize=12, pad=6)
axA.legend(loc="upper right", fontsize=7.5, framealpha=0.25, ncol=2)
axA.set_xlim(0.9, 5.5)
axA.set_ylim(-0.03, 1.3)
axA.grid(True, color=C_GRID, alpha=0.5)

# Panel B: F200W-F444W colour evolution
axB = fig.add_subplot(gs[1, 0])
axB.set_facecolor(PANEL_BG)
axB.plot(TAGE_STAGES, colors_200_444, "o-", color=STAGE_COLORS[2], lw=2.5, ms=8, label="Dark-RSFG")
axB.scatter([TAGE_LRD], [c200_lrd], s=180, color=C_LRD, marker="*", zorder=6,
            label=f"Phase4 LRD (z={ZRED_LRD})")
axB.axhline(1.0, color=C_OBS, lw=1.5, ls="--", label="LRD threshold (F200W-F444W=1.0)")
for i, (tage, c200, zred) in enumerate(zip(TAGE_STAGES, colors_200_444, ZRED_STAGES)):
    axB.annotate(f"z={zred}", (tage, c200), textcoords="offset points", xytext=(5,4),
                 fontsize=7.5, color=STAGE_COLORS[i])
axB.set_xlabel("Stellar age [Gyr]", fontsize=9.5)
axB.set_ylabel("F200W - F444W [mag]", fontsize=9.5)
axB.set_title("F200W-F444W colour evolution\n(Dark-RSFG -> LRD)", fontsize=10)
axB.legend(fontsize=7.5, framealpha=0.25)
axB.grid(True, color=C_GRID, alpha=0.5)

# Panel C: Emission line EW time evolution
axC = fig.add_subplot(gs[1, 1])
axC.set_facecolor(PANEL_BG)
axC.semilogy(TAGE_STAGES, [max(e, 1e-5) for e in ew_ha_stages], "o-", color="#ef9a9a",
             lw=2.0, ms=7, label="EW(Ha)")
axC.semilogy(TAGE_STAGES, [max(e, 1e-5) for e in ew_oiii_stages], "s-", color="#80deea",
             lw=2.0, ms=7, label="EW([OIII])")
axC.semilogy(TAGE_STAGES, [max(e, 1e-5) for e in ew_oii_stages], "^", color="#a5d6a7",
             lw=1.8, ms=7, label="EW([OII])")
axC.scatter([TAGE_LRD], [max(ew_ha_lrd, 1e-5)], s=180, color=C_LRD, marker="*",
            zorder=6, label="Phase4 LRD Ha")
axC.axhline(1.0, color=C_OBS, lw=1.2, ls="--", alpha=0.7, label="EW = 1 A")
axC.set_xlabel("Stellar age [Gyr]", fontsize=9.5)
axC.set_ylabel("EW [A] (log scale)", fontsize=9.5)
axC.set_title("Emission line EW evolution\n(consistently tiny – no O/B stars)", fontsize=10)
axC.legend(fontsize=7.5, framealpha=0.25)
axC.grid(True, color=C_GRID, alpha=0.5)

# Panel D: JWST detectability
axD = fig.add_subplot(gs[1, 2])
axD.set_facecolor(PANEL_BG)
bar_x = np.array(TAGE_STAGES + [TAGE_LRD])
bar_y = np.array(fluxes_f277w + [fl_lrd])
bar_cols = STAGE_COLORS + [C_LRD]
bar_labels = [f"z={z}" for z in ZRED_STAGES] + [f"LRD z={ZRED_LRD}"]
bars = axD.bar(bar_x, bar_y, width=0.07, color=bar_cols, alpha=0.85, align="center")
for bar, val, lbl in zip(bars, bar_y, bar_labels):
    axD.text(bar.get_x()+bar.get_width()/2, val*1.05, f"{val:.1f}", ha="center",
             fontsize=7.5, color=C_TXT)
    axD.text(bar.get_x()+bar.get_width()/2, -0.6, lbl, ha="center", fontsize=6.5,
             color=C_TXT, rotation=30)
axD.axhline(JWST_LIMIT_10HR, color=C_LIMIT, lw=2.0, ls="--",
            label=f"JWST 5σ 10h ({JWST_LIMIT_10HR} nJy)")
axD.axhline(JWST_LIMIT_100HR, color="#ff8a65", lw=1.5, ls="--",
            label=f"JWST 5σ 100h ({JWST_LIMIT_100HR} nJy)")
axD.fill_between([-0.1, bar_x[-1]+0.2], JWST_LIMIT_10HR, 0, alpha=0.07, color=C_LIMIT,
                 label="Below current detection limit")
axD.set_xlabel("Stellar age [Gyr]", fontsize=9.5)
axD.set_ylabel("F277W flux [nJy]", fontsize=9.5)
axD.set_title(f"JWST detectability\n(assuming M={10**LOG_MASS_CORE:.0e} Msun)", fontsize=10)
axD.legend(fontsize=7.5, framealpha=0.25, loc="upper left")
axD.set_xlim(-0.05, TAGE_LRD+0.15)
axD.set_ylim(-0.5, max(bar_y)*1.4)
axD.grid(True, color=C_GRID, alpha=0.5, axis="y")

# Panel E: Darkest stage vs Phase4 LRD (rest-frame)
axE = fig.add_subplot(gs[2, 0])
axE.set_facecolor(PANEL_BG)
wave_r_dark = wave_stages[0] / (1 + ZRED_STAGES[0])
wave_r_lrd = w_lrd / (1 + ZRED_LRD)
mask_d = (wave_r_dark > 500) & (wave_r_dark < 8000)
mask_l = (wave_r_lrd > 500) & (wave_r_lrd < 8000)
norm_d = np.percentile(spec_stages[0][mask_d], 95) or 1.0
norm_l = np.percentile(s_lrd[mask_l], 95) or 1.0
axE.fill_between(wave_r_dark[mask_d], 0, spec_stages[0][mask_d]/norm_d,
                 alpha=0.10, color=STAGE_COLORS[0])
axE.plot(wave_r_dark[mask_d], spec_stages[0][mask_d]/norm_d, color=STAGE_COLORS[0], lw=1.8,
         label=f"Dark-RSFG z={ZRED_STAGES[0]} tage={TAGE_STAGES[0]} Gyr")
axE.plot(wave_r_lrd[mask_l], s_lrd[mask_l]/norm_l, color=C_LRD, lw=1.8, ls="--",
         label=f"Phase4 LRD z={ZRED_LRD} tage={TAGE_LRD} Gyr")
for lname, lw0 in [("Lya",1216),("CIV",1549),("CIII",1909),("MgII",2798),
                   ("OII",3727),("Hβ",4863),("OIII",5007),("Ha",6565)]:
    axE.axvline(lw0, color="#546e7a", alpha=0.35, lw=0.7, ls=":")
    if lw0 < 8000:
        axE.text(lw0+30, 0.97, lname, fontsize=6, color="#78909c",
                 transform=axE.get_xaxis_transform(), va="top")
axE.set_xlabel("Rest-frame wavelength [A]", fontsize=9.5)
axE.set_ylabel("Normalised flux", fontsize=9.5)
axE.set_title("Darkest stage vs Phase4 LRD\nRest-frame spectrum", fontsize=10)
axE.legend(fontsize=7.5, framealpha=0.25)
axE.set_xlim(900, 7500)
axE.set_ylim(-0.02, 1.1)
axE.grid(True, color=C_GRID, alpha=0.5)

# Panel F: Colour-colour evolution track
axF = fig.add_subplot(gs[2, 1])
axF.set_facecolor(PANEL_BG)
all_c200 = colors_200_444 + [c200_lrd]
all_c150 = colors_150_277 + [md_lrd["jwst_f150w"] - md_lrd["jwst_f277w"]]
for i in range(len(all_c200)-1):
    axF.annotate("", xy=(all_c200[i+1], all_c150[i+1]),
                 xytext=(all_c200[i], all_c150[i]),
                 arrowprops=dict(arrowstyle="->", color="#ffd54f", lw=1.2, alpha=0.7))
sc = axF.scatter(all_c200[:-1], all_c150[:-1], c=TAGE_STAGES, cmap="plasma_r",
                 s=80, zorder=5, label="Dark-RSFG track")
cb = plt.colorbar(sc, ax=axF, pad=0.02)
cb.set_label("tage [Gyr]", fontsize=8.5, color=C_TXT)
cb.ax.yaxis.set_tick_params(color=C_TXT)
axF.scatter([c200_lrd], [all_c150[-1]], s=200, color=C_LRD, marker="*", zorder=7, label="Phase4 LRD")
axF.axhline(0, color=C_OBS, lw=0.8, ls=":", alpha=0.5)
axF.axhline(1.0, color=C_OBS, lw=1.2, ls="--", alpha=0.6, label="F200W-F444W > 1 (LRD threshold)")
for i, (c200, c150, tage, zred) in enumerate(zip(all_c200[:-1], all_c150[:-1],
                                                TAGE_STAGES, ZRED_STAGES)):
    axF.annotate(f"z={zred}\n{tage} Gyr", (c200, c150), textcoords="offset points",
                 xytext=(5,3), fontsize=6.5, color=STAGE_COLORS[i])
axF.set_xlabel("F200W - F444W [mag]", fontsize=9.5)
axF.set_ylabel("F150W - F277W [mag]", fontsize=9.5)
axF.set_title("Colour-colour evolution\nDark-RSFG -> LRD", fontsize=10)
axF.legend(fontsize=7.5, framealpha=0.25)
axF.grid(True, color=C_GRID, alpha=0.5)

# Panel G: Theory prediction vs observation strategy
axG = fig.add_subplot(gs[2, 2])
axG.set_facecolor(PANEL_BG)
axG.axis("off")
rows = [
    ["Property", "Dark-RSFG prediction", "Phase4 LRD (ref)"],
    ["Redshift", "z = 7.5-10", f"z = 5-7"],
    ["Stellar age", "tage < 0.5 Gyr", f"tage = {TAGE_LRD} Gyr"],
    ["Metallicity", f"Z = {10**LOGZSOL_DARK:.0e} Zsun", f"Z = 0.1 Zsun"],
    ["imf3", f"{IMF3_DARK} (extreme bottom-heavy)", f"{IMF3_LRD} (bottom-heavy)"],
    ["F200W-F444W", f"{colors_200_444[0]:.2f}", f"{c200_lrd:.2f}"],
    ["EW(Ha)", f"< {ew_ha_stages[0]:.3f} A", "< 1 A"],
    ["F277W flux", f"{fluxes_f277w[0]:.1f} nJy", f"{fl_lrd:.1f} nJy"],
    ["JWST 10h", "Difficult", "Detectable"],
    ["JWST 100h", "To be confirmed", "Easy"],
    ["Strategy", "Ultra-deep lensing field", "Standard deep survey"],
]
y_start = 0.97
row_h = 0.082
for row_i, row in enumerate(rows):
    y = y_start - row_i * row_h
    bg = "#1a237e" if row_i==0 else ("#12122a" if row_i%2==0 else "#0d0d20")
    wt = "bold" if row_i==0 else "normal"
    rect = FancyBboxPatch((0.01, y-row_h*0.80), 0.98, row_h*0.83,
                          boxstyle="round,pad=0.01", fc=bg, ec="#3d3d6b", lw=0.5,
                          transform=axG.transAxes, zorder=0)
    axG.add_patch(rect)
    for ci, (cell, x) in enumerate(zip(row, [0.02, 0.36, 0.69])):
        axG.text(x, y-row_h*0.36, cell, transform=axG.transAxes,
                 fontsize=6.8, fontweight=wt, va="center",
                 color=C_OBS if row_i==0 else C_TXT)
axG.set_title("Dark-RSFG: theory prediction & strategy", fontsize=10, pad=6)

fig.suptitle("HC-RSFG Theory - Code 4: Dark-RSFG Simulation (Phase 2-3, undetected)\n"
             f"z ~ 7.5-10, tage = {TAGE_STAGES[0]}-{TAGE_STAGES[-1]} Gyr, "
             f"logzsol = {LOGZSOL_DARK}, imf3 = {IMF3_DARK}",
             fontsize=12, color=C_TXT, y=0.97)
out_path = "HC_RSFG_Code4_DarkRSFG.png"
plt.savefig(out_path, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()
print(f"Saved: {out_path}")
print(f"Darkest stage (z={ZRED_STAGES[0]}, tage={TAGE_STAGES[0]}):")
print(f"  F277W = {fluxes_f277w[0]:.1f} nJy (detection difficult)")
print(f"  EW(Ha) = {ew_ha_stages[0]:.4f} A (nearly zero)")
print(f"Phase4 LRD comparison (z={ZRED_LRD}): F277W = {fl_lrd:.1f} nJy")
print("Strategy: Use gravitational lensing (mu x10-100) in ultra-deep surveys.")
