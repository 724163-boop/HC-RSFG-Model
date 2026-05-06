"""
HC-RSFG Theory Code 1: Little Red Dot (LRD) State Simulation
Phase 4: HC-RSFG Full Activation Phase (z=6.0, tage=0.9 Gyr)
Reference: Hviding et al. 2025 (RUBIES Survey), Labbe et al. 2023
"""

import numpy as np
import warnings
warnings.filterwarnings("ignore")
import fsps

# =================================================================
# 1. Parameter Definitions
# =================================================================
Z = 6.0
TAGE = 0.9
LOGZSOL = -1.0
DUST2 = 0.08
TAU = 0.10
GAS_LOGU = -3.0

IMF_TYPE = 2
IMF1 = 1.8
IMF2 = 1.8
IMF3 = 3.0
IMF_UMAX = 1.5

LRD_THRESH = 1.0
BANDS = ["jwst_f115w", "jwst_f150w", "jwst_f200w", "jwst_f277w", "jwst_f356w", "jwst_f444w"]

# =================================================================
# 2. FSPS Instance Creation
# =================================================================
def make_sps(neb=False, imf_type=IMF_TYPE, logzsol=LOGZSOL, dust2=DUST2, zred=Z, tage=TAGE, imf_umax=IMF_UMAX):
    kw = dict(compute_vega_mags=False, zcontinuous=1)
    if neb:
        kw["add_neb_emission"] = True
    sps = fsps.StellarPopulation(**kw)
    sps.params["imf_type"] = imf_type
    sps.params["imf_upper_limit"] = imf_umax
    sps.params["imf1"] = IMF1
    sps.params["imf2"] = IMF2
    sps.params["imf3"] = IMF3
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
print(f" IMF: Kroupa Custom (imf1={IMF1}, imf2={IMF2}, imf3={IMF3}, imf_upper_limit={IMF_UMAX})")
print("=================================================================")

sps_hc = make_sps(neb=False)
sps_neb = make_sps(neb=True)

# =================================================================
# 3. Band Color Calculation
# =================================================================
mags_hc = sps_hc.get_mags(tage=TAGE, bands=BANDS)

def color(mags, b1, b2):
    return mags[BANDS.index(b1)] - mags[BANDS.index(b2)]

c200_444_hc = color(mags_hc, "jwst_f200w", "jwst_f444w")
c277_444_hc = color(mags_hc, "jwst_f277w", "jwst_f444w")
c115_200_hc = color(mags_hc, "jwst_f115w", "jwst_f200w")

print("\n[ JWST Band Colors (HC-RSFG) ]")
print(f" {'Color':<15} {'HC-RSFG':>10} {'Status'}")
print("-" * 40)
lrd_flag = " LRD!" if c200_444_hc > LRD_THRESH else ""
print(f" {'F115W-F200W':<15} {c115_200_hc:>+10.3f}")
print(f" {'F200W-F444W':<15} {c200_444_hc:>+10.3f} {lrd_flag}")
print(f" {'F277W-F444W':<15} {c277_444_hc:>+10.3f}")

# =================================================================
# 4. Emission Line Equivalent Width (EW) Calculation
# =================================================================
wave_neb, spec_neb = sps_neb.get_spectrum(tage=TAGE)
wave_rest_neb = wave_neb / (1 + Z)

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
    if m_line.sum() == 0 or m_cont.sum() == 0:
        return 0.0
    peak = np.max(spec[m_line])
    cont = np.median(spec[m_cont])
    if cont <= 0:
        return 0.0
    return (peak - cont) / cont * (2 * dw_line)

print(f"\n[ Emission Line EW (logzsol={LOGZSOL}, gas_logu={GAS_LOGU}) ]")
print(f" {'Line':<15} {'EW [A]':>8}   {'Classification'}")
print("-" * 45)
for name, w0 in EMISSION_LINES.items():
    ew = compute_ew(wave_rest_neb, spec_neb, w0)
    if ew > 2.0:
        tag = "Strong (AGN-like)"
    elif ew > 0.5:
        tag = "Moderate"
    elif ew > 0.05:
        tag = "Weak (Narrow LRD)"
    else:
        tag = "Undetected"
    print(f" {name:<15} {ew:>8.3f}   {tag}")

# =================================================================
# 5. Metallicity Scan
# =================================================================
print(f"\n[ Metallicity Scan (F200W-F444W vs logzsol) ]")
print(f" {'logzsol':>8} {'F200W-F444W':>12} {'LRD?':>6}")
print("-" * 32)
logzsol_grid = np.array([-3.0, -2.5, -2.0, -1.5, -1.3, -1.0, -0.7, -0.5, -0.3, 0.0])
for lz in logzsol_grid:
    sp = make_sps(neb=False, logzsol=lz)
    m = sp.get_mags(tage=TAGE, bands=BANDS)
    c = color(m, "jwst_f200w", "jwst_f444w")
    flag = "YES" if c > LRD_THRESH else "NO"
    marker = " <--" if abs(lz - LOGZSOL) < 0.01 else ""
    print(f" {lz:>+8.1f} {c:>+12.3f} {flag:>6}{marker}")

# =================================================================
# 6. Dust Scan
# =================================================================
print(f"\n[ Dust Scan (F200W-F444W vs Av) ]")
print(f" {'dust2':>6} {'Av [mag]':>10} {'F200W-F444W':>12} {'LRD?':>6}")
print("-" * 40)
dust2_grid = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0])
for d2 in dust2_grid:
    sp = make_sps(neb=False, logzsol=LOGZSOL, dust2=d2)
    m = sp.get_mags(tage=TAGE, bands=BANDS)
    c = color(m, "jwst_f200w", "jwst_f444w")
    av = d2 * 3.1
    flag = "YES" if c > LRD_THRESH else "NO"
    marker = " <--" if abs(d2 - DUST2) < 0.01 else ""
    print(f" {d2:>6.2f} {av:>10.2f} {c:>+12.3f} {flag:>6}{marker}")

print("\n=================================================================")
print(f" ✅ Code 1 Completed")
print("=================================================================")
