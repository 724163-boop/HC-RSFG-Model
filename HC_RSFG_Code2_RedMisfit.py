"""
HC-RSFG Theory Code 2: Red Misfit Intermediate State Simulation
Phase 5 Transition (z ≈ 1-3, Universe Age 2-6 Gyr)
Reference: Evans et al. (2018) - Discovery of Red Misfits
"""


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

# =================================================================
# 1. Parameter Definitions
# =================================================================
# -- HC-RSFG Core (Old stars formed at z=7.5) --
ZRED_OBS = 0.5          # Observation Redshift
TAGE_CORE = 4.0         # Core Age [Gyr]
LOGZSOL_CORE = -0.5     # Core Metallicity
DUST2_CORE = 0.10       # Core Dust
IMF_TYPE_CORE = 2       
IMF1_HC = 1.8
IMF2_HC = 1.8
IMF3_HC = 3.0           # Bottom-Heavy IMF

# -- Initial Disk (Normal star formation starting in Phase 5) --
TAGE_DISK = 0.8         # Disk Age [Gyr]
LOGZSOL_DISK = -0.3     # Disk Metallicity
DUST2_DISK = 0.25       # Disk Dust
TAU_DISK = 0.5          # Exponential decay tau [Gyr]
IMF_DISK = 1            # Chabrier (Standard)

FRAC_CORE = 0.75        # Core mass fraction
FRAC_DISK = 0.25        # Disk mass fraction
UR_RED_SEQ = 2.22       # Red sequence threshold (Evans+2018)

BANDS_OPT = ["sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z"]
BANDS_NIR = ["2mass_j", "2mass_h", "2mass_ks"]

# =================================================================
# 2. FSPS Instances
# =================================================================
def make_core(neb=False):
    kw = dict(compute_vega_mags=False, zcontinuous=1)
    if neb: kw["add_neb_emission"] = True
    sps = fsps.StellarPopulation(**kw)
    sps.params["imf_type"] = IMF_TYPE_CORE
    sps.params["imf1"] = IMF1_HC
    sps.params["imf2"] = IMF2_HC
    sps.params["imf3"] = IMF3_HC
    sps.params["imf_upper_limit"] = 1.5 
    sps.params["logzsol"] = LOGZSOL_CORE
    sps.params["sfh"] = 0
    sps.params["tage"] = TAGE_CORE
    sps.params["dust2"] = DUST2_CORE
    sps.params["zred"] = ZRED_OBS
    if neb:
        sps.params["gas_logz"] = LOGZSOL_CORE
        sps.params["gas_logu"] = -4.0
    return sps

def make_disk(neb=False):
    kw = dict(compute_vega_mags=False, zcontinuous=1)
    if neb: kw["add_neb_emission"] = True
    sps = fsps.StellarPopulation(**kw)
    sps.params["imf_type"] = IMF_DISK
    sps.params["logzsol"] = LOGZSOL_DISK
    sps.params["sfh"] = 1 
    sps.params["tau"] = TAU_DISK
    sps.params["tage"] = TAGE_DISK
    sps.params["dust2"] = DUST2_DISK
    sps.params["zred"] = ZRED_OBS
    if neb:
        sps.params["gas_logz"] = LOGZSOL_DISK
        sps.params["gas_logu"] = -2.5
    return sps

print("=================================================================")
print(" HC-RSFG Theory - Code 2: Red Misfit Simulation")
print(f" Obs Redshift: z = {ZRED_OBS}")
print(f" Core: SSP, tage={TAGE_CORE} Gyr, IMF3={IMF3_HC} (Bottom-Heavy)")
print(f" Disk: tau={TAU_DISK} Gyr, tage={TAGE_DISK} Gyr, Chabrier")
print(f" Mass Ratio -> Core:Disk = {FRAC_CORE:.0%} : {FRAC_DISK:.0%}")
print("=================================================================")

sps_core = make_core(neb=False)
sps_core_neb = make_core(neb=True)
sps_disk = make_disk(neb=False)
sps_disk_neb = make_disk(neb=True)

sps_all_std = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1)
sps_all_std.params["imf_type"] = 1
sps_all_std.params["logzsol"] = LOGZSOL_CORE
sps_all_std.params["sfh"] = 0
sps_all_std.params["tage"] = TAGE_CORE
sps_all_std.params["dust2"] = DUST2_CORE
sps_all_std.params["zred"] = ZRED_OBS

# =================================================================
# 3. Spectrum Synthesis
# =================================================================
wave_c, spec_c = sps_core.get_spectrum(tage=TAGE_CORE)
wave_cn, spec_cn = sps_core_neb.get_spectrum(tage=TAGE_CORE)
wave_d, spec_d = sps_disk.get_spectrum(tage=TAGE_DISK)
wave_dn, spec_dn = sps_disk_neb.get_spectrum(tage=TAGE_DISK)
wave_s, spec_s = sps_all_std.get_spectrum(tage=TAGE_CORE)

norm_c = np.trapz(spec_c, wave_c) or 1.0
norm_d = np.trapz(spec_d, wave_d) or 1.0
norm_cn = np.trapz(spec_cn, wave_cn) or 1.0
norm_dn = np.trapz(spec_dn, wave_dn) or 1.0

spec_mix = FRAC_CORE * (spec_c / norm_c) + FRAC_DISK * (spec_d / norm_d)
spec_mix_neb = FRAC_CORE * (spec_cn / norm_cn) + FRAC_DISK * (spec_dn / norm_dn)
spec_std_n = spec_s / (np.trapz(spec_s, wave_s) or 1.0)
wave_rest = wave_c / (1 + ZRED_OBS)

# =================================================================
# 4. Band Color Calculation
# =================================================================
def get_color_composite(sps_a, tage_a, sps_b, tage_b, frac_a, frac_b, bands):
    mags_a = np.array(sps_a.get_mags(tage=tage_a, bands=bands))
    mags_b = np.array(sps_b.get_mags(tage=tage_b, bands=bands))
    flux_a = 10 ** (-0.4 * mags_a) * frac_a
    flux_b = 10 ** (-0.4 * mags_b) * frac_b
    flux_tot = flux_a + flux_b
    return dict(zip(bands, -2.5 * np.log10(np.clip(flux_tot, 1e-30, None))))

bands_all = BANDS_OPT + BANDS_NIR
mag_mix = get_color_composite(sps_core, TAGE_CORE, sps_disk, TAGE_DISK, FRAC_CORE, FRAC_DISK, bands_all)
mag_core = dict(zip(bands_all, sps_core.get_mags(tage=TAGE_CORE, bands=bands_all)))
mag_std = dict(zip(bands_all, sps_all_std.get_mags(tage=TAGE_CORE, bands=bands_all)))

ur_mix = mag_mix["sdss_u"] - mag_mix["sdss_r"]
ur_core = mag_core["sdss_u"] - mag_core["sdss_r"]
ur_std = mag_std["sdss_u"] - mag_std["sdss_r"]

print(f"\n[ Optical Band Colors @ z={ZRED_OBS} ]")
print(f" u-r (Mixed Model) = {ur_mix:+.3f} ({'Red Sequence' if ur_mix > UR_RED_SEQ else 'Blue Cloud'})")

# =================================================================
# 5. Emission Line EW Calculation
# =================================================================
HA_REST = 6564.6
def compute_ew_simple(wave_r, spec, w0, dw=30, dw_cont=120):
    m_line = (wave_r > w0 - dw) & (wave_r < w0 + dw)
    m_cont = (((wave_r > w0 - dw_cont) & (wave_r < w0 - dw * 2)) |
              ((wave_r > w0 + dw * 2) & (wave_r < w0 + dw_cont)))
    if m_line.sum() == 0 or m_cont.sum() == 0: return 0.0
    peak = np.max(spec[m_line])
    cont = np.median(spec[m_cont])
    return max((peak - cont) / cont * (2 * dw), 0.0) if cont > 0 else 0.0

ew_ha_mix = compute_ew_simple(wave_rest, spec_mix_neb, HA_REST)
print(f" Mixed Model Ha EW = {ew_ha_mix:.2f} A")

# =================================================================
# 6 & 7. Theoretical Grids and Plotting
# =================================================================
frac_core_scan = np.linspace(0.0, 1.0, 9)
ur_scan, ew_scan = [], []
for fc in frac_core_scan:
    fd = 1.0 - fc
    mc = np.array(sps_core_neb.get_mags(tage=TAGE_CORE, bands=["sdss_u","sdss_r"]))
    md = np.array(sps_disk_neb.get_mags(tage=TAGE_DISK, bands=["sdss_u","sdss_r"]))
    fu = fc * 10**(-0.4*mc[0]) + fd * 10**(-0.4*md[0])
    fr = fc * 10**(-0.4*mc[1]) + fd * 10**(-0.4*md[1])
    ur_scan.append(-2.5 * np.log10(max(fu/max(fr, 1e-30), 1e-30)))
    sp_tmp = fc * (spec_cn / norm_cn) + fd * (spec_dn / norm_dn)
    ew_scan.append(compute_ew_simple(wave_rest, sp_tmp, HA_REST))

DARK_BG, PANEL_BG = "#090912", "#0d0d20"
C_HC, C_DISK, C_MIX, C_STD, C_OBS, C_GRID, C_TXT, C_EVANS = "#ff8a65", "#4dd0e1", "#ce93d8", "#80cbc4", "#ffd54f", "#1e1e3a", "#e8eaf6", "#a5d6a7"

plt.rcParams.update({"figure.facecolor": DARK_BG, "text.color": C_TXT, "axes.labelcolor": C_TXT,
                     "xtick.color": C_TXT, "ytick.color": C_TXT, "axes.edgecolor": "#3d3d6b", "font.family": "monospace"})

fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.38, top=0.92, bottom=0.06, left=0.06, right=0.97)

# Panel A: Synthetic Spectrum
axA = fig.add_subplot(gs[0, :], facecolor=PANEL_BG)
mask_obs = (wave_c > 3000) & (wave_c < 22000)
norm_ref = np.percentile(spec_mix[mask_obs], 95) or 1.0
axA.plot(wave_c[mask_obs] / 1e4, spec_mix[mask_obs] / norm_ref, color=C_MIX, lw=2.0, label=f"Mixed (Core {FRAC_CORE:.0%} + Disk {FRAC_DISK:.0%})")
axA.plot(wave_c[mask_obs] / 1e4, spec_mix_neb[mask_obs] / norm_ref, color=C_OBS, lw=0.9, alpha=0.65, label="Mixed + Nebular")
axA.plot(wave_c[mask_obs] / 1e4, (spec_c[mask_obs]/norm_c)*FRAC_CORE/norm_ref, color=C_HC, lw=1.2, ls="--", alpha=0.7, label="Core Component")
axA.plot(wave_d[mask_obs] / 1e4, (spec_d[mask_obs]/norm_d)*FRAC_DISK/norm_ref, color=C_DISK, lw=1.2, ls=":", alpha=0.7, label="Disk Component")
axA.set(xlabel="Observed Wavelength [μm]", ylabel="Normalized Flux", title=f"Phase 5 Red Misfit - Synthetic Spectrum (z={ZRED_OBS})")
axA.legend(loc="upper right", framealpha=0.25)

# Panel B: Rest-frame Spectrum
axB = fig.add_subplot(gs[1, 0], facecolor=PANEL_BG)
mask_r = (wave_rest > 3000) & (wave_rest < 8000)
norm_B = np.percentile(spec_mix[mask_r], 98) or 1.0
axB.plot(wave_rest[mask_r], spec_mix[mask_r] / norm_B, color=C_MIX, lw=1.5, label="Mixed Model")
axB.plot(wave_rest[mask_r], spec_mix_neb[mask_r] / norm_B, color=C_OBS, lw=0.9, alpha=0.7, label="Mixed + Nebular")
axB.set(xlabel="Rest-frame Wavelength [A]", title="Rest-frame Spectrum (3000-8000A)")
axB.legend(framealpha=0.25)

# Panel C: Emission EW Bar Chart
axC = fig.add_subplot(gs[1, 1], facecolor=PANEL_BG)
ew_dict = {name: compute_ew_simple(wave_rest, spec_mix_neb, w0) for name, w0 in 
           {"[OII]":3727.1, "Hb":4862.7, "[OIII]":5007.0, "Ha":6564.6}.items()}
axC.barh(list(ew_dict.keys()), list(ew_dict.values()), color=C_MIX, alpha=0.88)
axC.axvline(1.0, color=C_OBS, ls="--", label="EW=1A Limit")
axC.set(xlabel="Equivalent Width [A]", title="Emission Line EW")
axC.legend(framealpha=0.25)

# Panel D: u-r vs Ha EW Diagnostic
axD = fig.add_subplot(gs[1, 2], facecolor=PANEL_BG)
axD.plot(ur_scan, ew_scan, "o-", color=C_MIX, lw=2.0, label="HC-RSFG Theory Track")
axD.scatter([ur_mix], [ew_ha_mix], s=180, color=C_OBS, marker="*", zorder=8, label=f"This Model")
axD.axvline(UR_RED_SEQ, color=C_OBS, ls="--", alpha=0.6)
axD.set(xlabel="u - r [mag]", ylabel="EW(Ha) [A]", title="u-r vs EW(Ha) Diagnostic Diagram")
axD.legend(framealpha=0.25)

# Panel E: Core Fraction vs Colors
axE = fig.add_subplot(gs[2, 0], facecolor=PANEL_BG)
ax_twin = axE.twinx()
axE.plot(frac_core_scan * 100, ur_scan, "o-", color=C_HC, label="u-r (Left)")
ax_twin.plot(frac_core_scan * 100, ew_scan, "s--", color=C_DISK, label="EW(Ha) (Right)")
axE.set(xlabel="HC-RSFG Core Fraction [%]", ylabel="u - r [mag]", title="Core Fraction vs Color/Lines")

# Panel F: SED Comparison
axF = fig.add_subplot(gs[2, 1], facecolor=PANEL_BG)
band_waves = np.array([0.355, 0.468, 0.617, 0.748, 0.893, 1.235, 1.662, 2.159])
flux_mix_arr = 10**(-0.4*(np.array([mag_mix[b] for b in bands_all]) - np.median([mag_mix[b] for b in bands_all])))
axF.plot(band_waves, flux_mix_arr, "o-", color=C_MIX, lw=2.0, label="Mixed Model")
axF.set(xlabel="Observed Wavelength [μm]", ylabel="Normalized Flux", title="Mock SED Fitting")
axF.legend(framealpha=0.25)

# Panel G: Table (Simplified for text output in Plot)
axG = fig.add_subplot(gs[2, 2], facecolor=PANEL_BG)
axG.axis("off")
axG.text(0.1, 0.6, f"Red Misfit Match:\n\nModel u-r = {ur_mix:.2f} (Obs > 2.22)\nModel Ha EW = {ew_ha_mix:.2f}A (Obs 1-15A)", color=C_TXT, fontsize=12)
axG.set_title("Evans+2018 Red Misfit Comparison", fontsize=10)

fig.suptitle("HC-RSFG Theory - Code 2: Red Misfit Simulation (Phase 5 Transition)", fontsize=14, y=0.97)
plt.savefig("HC_RSFG_Code2_RedMisfit.png", dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"\n=================================================================")
print(f" ✅ Code 2 Completed - Plot saved as: HC_RSFG_Code2_RedMisfit.png")
print("=================================================================")
