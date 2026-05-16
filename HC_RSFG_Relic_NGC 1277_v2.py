# -*- coding: utf-8 -*-
"""Fig 3 only: Relic Galaxy NGC 1277 vs HC-RSFG (z≈0)"""

import numpy as np
import matplotlib.pyplot as plt
import fsps

HC_BASE = dict(
    sfh=1, tau=0.1, imf_type=2, imf1=1.8, imf2=1.8,
    imf_upper_limit=1.5, dust_tesc=7.0,
    add_neb_emission=True, gas_logz=-1.0, gas_logu=-3.0,
)

# --- HC-RSFG Relic ---
sp_relic = fsps.StellarPopulation(
    **HC_BASE, logzsol=-1.0, dust2=0.0, tage=12.0, imf3=3.0
)
wave_r, spec_r = sp_relic.get_spectrum(tage=12.0)
wave_r = np.array(wave_r); spec_r = np.array(spec_r)
if spec_r.ndim == 2:
    spec_r = spec_r[-1, :]
mags_r = sp_relic.get_mags(redshift=0.0173, bands=['sdss_g', 'sdss_r'], tage=12.0)
mags_r = np.asarray(mags_r).ravel()
gr_relic = float(mags_r[0] - mags_r[1])
print(f"HC-Relic g-r = {gr_relic:.3f} (target +0.851)")

# --- NGC 1277-like ---
sp_ngc = fsps.StellarPopulation(
    **HC_BASE, logzsol=-1.0, dust2=0.0, tage=13.0, imf3=4.0
)
wave_n, spec_n = sp_ngc.get_spectrum(tage=13.0)
wave_n = np.array(wave_n); spec_n = np.array(spec_n)
if spec_n.ndim == 2:
    spec_n = spec_n[-1, :]
mags_n = sp_ngc.get_mags(redshift=0.0173, bands=['sdss_g', 'sdss_r'], tage=13.0)
mags_n = np.asarray(mags_n).ravel()
gr_ngc = float(mags_n[0] - mags_n[1])
print(f"NGC 1277 g-r = {gr_ngc:.3f}")

fig, ax = plt.subplots(figsize=(9, 6))

z0 = 0.0173
w_r = wave_r * (1 + z0)
w_n = wave_n * (1 + z0)

idx_r = np.argmin(np.abs(w_r - 5500))
idx_n = np.argmin(np.abs(w_n - 5500))
norm_r = spec_r[idx_r]
norm_n = spec_n[idx_n]

ax.plot(w_r/1e4, spec_r/norm_r, color='darkorange', lw=2.5,
        label=f'HC-RSFG Relic (t=12 Gyr, imf3=3.0)  g-r={gr_relic:.3f}')
ax.plot(w_n/1e4, spec_n/norm_n, color='black', lw=2, ls='--',
        label=f'NGC 1277-like (t=13 Gyr, imf3=4.0)  g-r={gr_ngc:.3f}')

ax.set_xlim(0.3, 1.1)
ax.set_ylim(0.3, 1.3)
ax.set_xlabel(r'Observed Wavelength [${\mu}$m]', fontsize=12)
ax.set_ylabel(r'Normalized Flux', fontsize=12)
ax.set_title(r'Fig 2: Relic Galaxy NGC 1277 vs HC-RSFG ($z{\approx}0$)', fontsize=13, fontweight='bold')
ax.legend(loc='upper right', fontsize=10)
ax.grid(alpha=0.3)
ax.axvline(0.4*(1+z0), color='teal', ls=':', alpha=0.5)
ax.text(0.4*(1+z0), 1.25, '4000Å break', fontsize=9, color='teal')

plt.tight_layout()
plt.savefig('fig3_relic_ngc1277.png', dpi=200, bbox_inches='tight')
plt.show()
