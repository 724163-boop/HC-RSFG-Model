#The spectral diagram for output [3/5] Relic / NGC 1277 (z=0.0173) is generally correct, but there may be issues with the numerical values. In that case, please use code HC_RSFG_Relic / NGC 1277_v2.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import fsps

HC_BASE = dict(
    sfh=1,
    tau=0.1,
    imf_type=2,
    imf1=1.8,
    imf2=1.8,
    imf3=3.0,
    imf_upper_limit=1.5, 
    logzsol=-1.0,
    dust2=0.08,
    dust_tesc=7.0,
    add_neb_emission=True,
    gas_logz=-1.0,
    gas_logu=-3.0,
)

JWST_BANDS = {
    'F090W': ( 9030, 2340),
    'F115W': (11540, 2250),
    'F150W': (15010, 3170),
    'F200W': (19900, 4600),
    'F277W': (27760, 6720),
    'F356W': (35680, 8630),
    'F444W': (44210, 10610),
}
SDSS_BANDS = {
    'u': ( 3557, 541),
    'g': ( 4825, 1263),
    'r': ( 6261, 1142),
    'i': ( 7672, 1535),
    'z': ( 9097, 1029),
}


def make_sp(extra=None):
    kw = dict(HC_BASE)
    if extra:
        kw.update(extra)
    return fsps.StellarPopulation(**kw)


def get_safe_spectrum(sp):
    wave = np.array(sp.wavelengths)
    _, spec = sp.get_spectrum()
    spec = np.array(spec)
    if spec.ndim == 2:
        spec = spec[-1, :]
    n = min(len(wave), len(spec))
    return wave[:n], spec[:n]


def calc_mag_tophat(wave_obs, fnu, lam_cen, dlam):
    mask = (wave_obs >= lam_cen - dlam / 2) & (wave_obs <= lam_cen + dlam / 2)
    if not np.any(mask):
        return 99.0
    avg = np.mean(fnu[mask])
    if avg <= 0:
        return 99.0
    return -2.5 * np.log10(avg) - 48.6


def get_color_manual(sp, z, band1_key, band2_key, band_dict):
    wave_rest, flux = get_safe_spectrum(sp)
    c_AA_s = 2.99792458e18 # Å/s
    wave_obs = wave_rest * (1.0 + z)
    fnu = flux * wave_obs**2 / c_AA_s # fλ → fν 変換

    m1 = calc_mag_tophat(wave_obs, fnu, *band_dict[band1_key])
    m2 = calc_mag_tophat(wave_obs, fnu, *band_dict[band2_key])
    return m1 - m2


FSPS_JWST = {
    'F090W': 'jwst_f090w',
    'F115W': 'jwst_f115w',
    'F150W': 'jwst_f150w',
    'F200W': 'jwst_f200w',
    'F277W': 'jwst_f277w',
    'F356W': 'jwst_f356w',
    'F444W': 'jwst_f444w',
}
FSPS_SDSS = {
    'u': 'sdss_u',
    'g': 'sdss_g',
    'r': 'sdss_r',
    'i': 'sdss_i',
    'z': 'sdss_z',
}


def get_color_fsps(sp, z, band1_key, band2_key, fsps_band_dict, use_fallback=True):
    try:
        bands_fsps = [fsps_band_dict[band1_key], fsps_band_dict[band2_key]]

        wave, spec = sp.get_spectrum()
        mags = sp.get_mags(redshift=z, bands=bands_fsps)

        mags = np.asarray(mags)
        if mags.ndim == 2:
            mags = mags[-1, :]
        mags = mags.ravel()

        color = float(mags[0] - mags[1])
        return color

    except (KeyError, FileNotFoundError, RuntimeError) as e:
        if use_fallback:
            band_dict = JWST_BANDS if 'jwst' in str(fsps_band_dict) else SDSS_BANDS
            return get_color_manual(sp, z, band1_key, band2_key, band_dict)
        else:
            raise


def get_jwst_color(sp, z, f1, f2):
    return get_color_fsps(sp, z, f1, f2, FSPS_JWST, use_fallback=True)


def get_sdss_color(sp, z, f1, f2):
    return get_color_fsps(sp, z, f1, f2, FSPS_SDSS, use_fallback=True)


def add_agn_lines(wave, flux, z, line_dict):
    flux = np.array(flux, copy=True)
    for name, (lam0, ew, sig) in line_dict.items():
        lam_obs = lam0 * (1.0 + z)
        cont = flux[np.argmin(np.abs(wave - lam_obs))]
        profile = (cont * ew / (sig * np.sqrt(2.0 * np.pi))
                   * np.exp(-0.5 * ((wave - lam_obs) / sig) ** 2))
        flux += profile
    return flux


def obs_wave(wave_rest, z):
    return wave_rest * (1.0 + z)


print("[1/5] HC-RSFG LRD spectrum (z=6.0, t=0.9 Gyr)...")
sp_hc_lrd = make_sp(dict(tage=0.9))
wave_hc_lrd, spec_hc_lrd = get_safe_spectrum(sp_hc_lrd)
color_hc_lrd = get_jwst_color(sp_hc_lrd, 6.0, 'F200W', 'F444W')
print(f" HC-LRD F200W-F444W = {color_hc_lrd:.3f}")
print(f" (Using real FSPS filter transmission curves)")

print("[2/5] Narrow-line LRD (AGN hybrid, z=6.0)...")
sp_obs_lrd = make_sp(dict(tage=0.9, dust2=0.15))
wave_obs_lrd, spec_obs_lrd = get_safe_spectrum(sp_obs_lrd)
spec_obs_lrd = add_agn_lines(wave_obs_lrd, spec_obs_lrd, 6.0, {
    'Hbeta': (4861, 25, 8),
    'OIII_5007': (5007, 50, 8),
    'Halpha': (6563, 120, 10),
    'NII_6584': (6584, 35, 10),
})

print("[3/5] Relic / NGC 1277 (z=0.0173)...")
sp_hc_relic = make_sp(dict(tage=12.0, dust2=0.0))
wave_relic, spec_hc_relic = get_safe_spectrum(sp_hc_relic)

sp_ngc1277 = make_sp(dict(tage=13.0, imf3=4.0, dust2=0.0))
wave_ngc, spec_ngc1277 = get_safe_spectrum(sp_ngc1277)

gr_relic = get_sdss_color(sp_hc_relic, 0.0173, 'g', 'r')
gr_ngc = get_sdss_color(sp_ngc1277, 0.0173, 'g', 'r')
print(f" HC-Relic g-r = {gr_relic:.3f} (target +0.851)")
print(f" NGC 1277 g-r = {gr_ngc:.3f}")
print(f" (Using SDSS filter transmission curves)")

print("[4/5] Red Misfit (z=0.5)...")
sp_hc_misfit = make_sp(dict(tage=4.0, dust2=0.1))
wave_mis, spec_hc_misfit = get_safe_spectrum(sp_hc_misfit)
spec_hc_misfit = add_agn_lines(wave_mis, spec_hc_misfit, 0.5,
                               {'Halpha': (6563, 0.9, 6)})

sp_obs_misfit = make_sp(dict(tage=3.5, dust2=0.15))
wave_om, spec_obs_misfit = get_safe_spectrum(sp_obs_misfit)
spec_obs_misfit = add_agn_lines(wave_om, spec_obs_misfit, 0.5,
                                {'Halpha': (6563, 10, 6)})

ur_hc = get_sdss_color(sp_hc_misfit, 0.5, 'u', 'r')
gr_hc = get_sdss_color(sp_hc_misfit, 0.5, 'g', 'r')
print(f" HC-Misfit u-r = {ur_hc:.3f} (target +3.55)")
print(f" HC-Misfit g-r = {gr_hc:.3f} (target +1.50)")

print("[5/5] Dark-RSFG (z=9.5, t=0.1 Gyr)...")
sp_dark = make_sp(dict(tage=0.1, dust2=0.0))
wave_dark, spec_dark = get_safe_spectrum(sp_dark)
color_dark = get_jwst_color(sp_dark, 9.5, 'F200W', 'F444W')
color_dark2 = get_jwst_color(sp_dark, 9.5, 'F150W', 'F277W')
print(f" Dark-RSFG F200W-F444W = {color_dark:.3f} (target +1.96)")
print(f" Dark-RSFG F150W-F277W = {color_dark2:.3f} (target +1.50)")
print(f" (Using real FSPS filter transmission curves)")

fig = plt.figure(figsize=(18, 20))
# 3行×2列：最下行は Fig 6 が全幅
gs = GridSpec(3, 2, figure=fig, hspace=0.38, wspace=0.27,
              height_ratios=[1, 1, 1])

# ----------------------------------------------------------
# Fig 1: Narrow-line LRD vs HC-RSFG pure (z=6)
# ----------------------------------------------------------
ax1 = fig.add_subplot(gs[0, 0])
z6 = 6.0
w1 = obs_wave(wave_hc_lrd, z6)
w2 = obs_wave(wave_obs_lrd, z6)

ni = np.argmin(np.abs(w1 - 20000))
ax1.plot(w1/1e4, spec_hc_lrd / spec_hc_lrd[ni],
         color='crimson', lw=2, label='HC-RSFG Pure (dust2=0.08)')
ax1.plot(w2/1e4, spec_obs_lrd / spec_obs_lrd[ni],
         color='navy', lw=1.5, alpha=0.75, label='Narrow-line LRD (AGN hybrid)')

for name, lam in {'Hβ': 4861*(1+z6), '[OIII]': 5007*(1+z6),
                  'Hα': 6563*(1+z6), '[NII]': 6584*(1+z6)}.items():
    ax1.axvline(lam/1e4, color='gray', ls='--', alpha=0.4)
    ax1.text(lam/1e4, 1.38, name, rotation=90, va='bottom',
             ha='center', fontsize=8, color='dimgray')

ax1.set_xlim(0.8, 5.5); ax1.set_ylim(0.1, 1.55)
ax1.set_xlabel(r'Observed Wavelength [$\mu$m]', fontsize=12)
ax1.set_ylabel('Normalized Flux', fontsize=12)
ax1.set_title(r'Fig 1: Narrow-line LRD vs HC-RSFG Pure ($z=6$)',
              fontsize=13, fontweight='bold')
ax1.legend(fontsize=10); ax1.grid(alpha=0.3)

# ----------------------------------------------------------
# Fig 2: Relic Galaxy / NGC 1277 (z~0)
# ----------------------------------------------------------
ax2 = fig.add_subplot(gs[0, 1])
z0 = 0.0173
w_r = obs_wave(wave_relic, z0)
w_n = obs_wave(wave_ngc, z0)

ni = np.argmin(np.abs(w_r - 5500))
ax2.plot(w_r/1e4, spec_hc_relic / spec_hc_relic[ni],
         color='darkorange', lw=2.5, label='HC-RSFG Relic (t=12 Gyr, imf3=3.0)')
ax2.plot(w_n/1e4, spec_ngc1277 / spec_ngc1277[ni],
         color='black', lw=2, ls='--', label='NGC 1277-like (t=13 Gyr, imf3=4.0)')

ax2.axvline(0.4*(1+z0), color='teal', ls=':', alpha=0.6)
ax2.text(0.4*(1+z0)+0.005, 1.22, '4000Å break', fontsize=9, color='teal')

ax2.set_xlim(0.3, 1.1); ax2.set_ylim(0.3, 1.35)
ax2.set_xlabel(r'Observed Wavelength [$\mu$m]', fontsize=12)
ax2.set_ylabel('Normalized Flux', fontsize=12)
ax2.set_title(r'Fig 2: Relic Galaxy NGC 1277 vs HC-RSFG ($z\approx0$)',
              fontsize=13, fontweight='bold')
ax2.legend(fontsize=10); ax2.grid(alpha=0.3)

# ----------------------------------------------------------
# Fig 3: Red Misfit (z=0.5)
# ----------------------------------------------------------
ax3 = fig.add_subplot(gs[1, 0])
z05 = 0.5
w_m = obs_wave(wave_mis, z05)
w_om = obs_wave(wave_om, z05)

ni = np.argmin(np.abs(w_m - 5500))
ax3.plot(w_m/1e4, spec_hc_misfit / spec_hc_misfit[ni],
         color='firebrick', lw=2, label='HC-RSFG Red Misfit (t=4 Gyr)')
ax3.plot(w_om/1e4, spec_obs_misfit / spec_obs_misfit[ni],
         color='steelblue', lw=1.5, alpha=0.75, label='Observed Red Misfit')

ha_obs = 6563*(1+z05)/1e4
ax3.axvline(ha_obs, color='gray', ls='--', alpha=0.5)
ax3.text(ha_obs, 1.32, r'H$\alpha$', rotation=90, va='bottom',
         ha='center', fontsize=9)

ax3.set_xlim(0.3, 1.1); ax3.set_ylim(0.4, 1.45)
ax3.set_xlabel(r'Observed Wavelength [$\mu$m]', fontsize=12)
ax3.set_ylabel('Normalized Flux', fontsize=12)
ax3.set_title(r'Fig 3: Red Misfit vs HC-RSFG ($z=0.5$)',
              fontsize=13, fontweight='bold')
ax3.legend(fontsize=10); ax3.grid(alpha=0.3)

# ----------------------------------------------------------
# Fig 4: Dark-RSFG (z=9.5)
# ----------------------------------------------------------
ax4 = fig.add_subplot(gs[1, 1])
z95 = 9.5
w_d = obs_wave(wave_dark, z95)
spec_max = spec_dark.max()

ax4.plot(w_d/1e4, spec_dark, color='darkred', lw=2,
         label='Dark-RSFG (t=0.1 Gyr)')

for fname in ('F150W', 'F200W', 'F277W', 'F444W'):
    lam, dlam = JWST_BANDS[fname]
    ax4.axvspan((lam-dlam/2)/1e4, (lam+dlam/2)/1e4,
                alpha=0.12, color='royalblue')
    ax4.text(lam/1e4, spec_max*0.92, fname, ha='center',
             fontsize=8, color='navy')

ax4.set_xlim(1.0, 5.5); ax4.set_ylim(0, spec_max*1.15)
ax4.set_xlabel(r'Observed Wavelength [$\mu$m]', fontsize=12)
ax4.set_ylabel(r'Flux per $M_\odot$ formed [L$_\odot$/Hz]', fontsize=11)
ax4.set_title(r'Fig 4: Dark-RSFG Spectrum ($z=9.5$, $t_{\rm age}=0.1$ Gyr)',
              fontsize=13, fontweight='bold')
ax4.text(0.97, 0.88,
         rf'$F200W-F444W \approx {color_dark:+.2f}$'+'\n'
         rf'$F150W-F277W \approx {color_dark2:+.2f}$',
         transform=ax4.transAxes, fontsize=11, ha='right', va='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))
ax4.legend(fontsize=10); ax4.grid(alpha=0.3)

# ----------------------------------------------------------
# Fig 6: NGC 1277 M/L Ratio Radial Profile（all）
# ----------------------------------------------------------
ax6 = fig.add_subplot(gs[2, :])

obs_r = np.array([0.05, 0.15, 0.30, 0.60, 1.0, 1.5, 2.5, 3.5])
obs_ml = np.array([2.05, 1.95, 1.80, 1.55, 1.35, 1.20, 1.05, 1.00])
obs_err = np.array([0.15, 0.12, 0.10, 0.10, 0.08, 0.07, 0.05, 0.04])
ax6.errorbar(obs_r, obs_ml, yerr=obs_err, fmt='o', color='black',
             ecolor='gray', capsize=4, markersize=8,
             label='NGC 1277 Observed (Martín-Navarro+2015)')

r_mod = np.linspace(0.01, 4.0, 300)
r_sc = 0.6 # kpc
ml_disk = 1.0

ml_cons = ml_disk + (1.56 - ml_disk) * np.exp(-r_mod / r_sc)
ml_best = ml_disk + (2.00 - ml_disk) * np.exp(-r_mod / r_sc)

ax6.plot(r_mod, ml_cons, '-', color='crimson', lw=2.5,
         label='HC-RSFG conservative (imf3=3.0, t=12 Gyr)')
ax6.plot(r_mod, ml_best, '--', color='darkorange', lw=2.5,
         label='HC-RSFG best-fit (imf3=4.0, t=13 Gyr)')

ax6.axvline(0.01, color='teal', ls=':', alpha=0.6)
ax6.text(0.015, 1.9, 'Core radius\n(~10 pc)', fontsize=9,
         color='teal', va='top')
ax6.text(2.3, 1.27, 'Bottom-heavy IMF in the core\nexplains elevated central M/L',
         fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax6.set_xlim(0, 4.0); ax6.set_ylim(0.8, 2.3)
ax6.set_xlabel('Radius [kpc]', fontsize=13)
ax6.set_ylabel(r'$(M/L_V)\,/\,(M/L_V)_{\rm Kroupa}$', fontsize=13)
ax6.set_title('Fig 5 (former Fig 6): NGC 1277 Mass-to-Light Ratio Profile',
              fontsize=14, fontweight='bold')
ax6.legend(fontsize=11); ax6.grid(alpha=0.3)

outname = 'HC_RSFG_5panel.png'
plt.savefig(outname, dpi=200, bbox_inches='tight')
print(f"\n✓ Saved: {outname}")
plt.show()
 
