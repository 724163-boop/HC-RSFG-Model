#!/usr/bin/env python3
"""
HC‑RSFG → LRD & AGN Co‑emergence Timeline (BUGFIX)
====================================================
Tracks gas accumulation, star formation (RSFG), and BH growth
from seed to low‑mass AGN.
"""

import numpy as np

# ============================================================
# 1. Model parameters
# ============================================================
Mdot_cold = 20.0 # Msun/yr (total cold stream inflow)
f_conc = 0.3 # fraction reaching inner core
Mdot_core = Mdot_cold * f_conc # 6 Msun/yr

t_onset = 220e6 # yr (time to reach critical density)
t_end = 500e6 # yr (simulation end)

tau_SF = 0.1e9 # yr (e‑folding time for gas → stars)
M_BH_seed = 100.0 # Msun
epsilon = 0.1
t_Salpeter = 4.5e7 * (epsilon / 0.1) / (1.0 - epsilon)
f_Edd = 3.0 # super‑Eddington (realistic for early BH growth)

M_AGN_target = 1e5 # Msun — to power narrow lines

# ============================================================
# 2. Time integration (fixed)
# ============================================================
dt = 1e5 # 0.1 Myr steps — fine enough
time = np.arange(0, t_end + dt, dt)
n = len(time)

M_gas = np.zeros(n)
M_star = np.zeros(n)
M_BH = np.zeros(n)

# Initial conditions
M_gas[0] = 0.0
M_star[0] = 0.0
M_BH[0] = M_BH_seed

# Flags
reached_AGN = False
t_AGN = np.nan

for i in range(1, n):
    t = time[i]
    Mg = M_gas[i-1]
    Ms = M_star[i-1]
    Mb = M_BH[i-1]
    
    # ---- Continuous cold stream inflow ----------
    dMg_in = Mdot_core * dt
    
    # ---- Star formation (only if t >= t_onset) ----
    if t >= t_onset and Mg > 0:
        dMs = (Mg / tau_SF) * dt
    else:
        dMs = 0.0
    
    # ---- Black hole accretion (only if t >= t_onset) ----
    if t >= t_onset and Mb > 0:
        Mdot_Edd = Mb / t_Salpeter * f_Edd
        dMb = min(Mdot_Edd * dt, (Mg + dMg_in) * 0.1) # ≤ 10% of available gas
    else:
        dMb = 0.0
    
    # Ensure we don't consume more gas than available
    gas_consumed = dMs + dMb
    available = Mg + dMg_in
    if gas_consumed > available:
        factor = available / gas_consumed
        dMs *= factor
        dMb *= factor
    
    # Update
    M_gas[i] = Mg + dMg_in - dMs - dMb
    M_star[i] = Ms + dMs
    M_BH[i] = Mb + dMb
    
    # Check AGN threshold
    if not reached_AGN and M_BH[i] >= M_AGN_target:
        reached_AGN = True
        t_AGN = t

# ============================================================
# 3. Output
# ============================================================
print("=" * 60)
print("HC‑RSFG → LRD & AGN CO‑EMERGENCE (corrected)")
print("=" * 60)
print(f" Gas accumulation time : {t_onset/1e6:.0f} Myr")
print(f" LRD phase begins : ~{t_onset/1e6:.0f} Myr")

if reached_AGN:
    delay = t_AGN - t_onset
    print(f" AGN reaches 10⁵ M⊙ : {t_AGN/1e6:.1f} Myr")
    print(f" Delay after LRD onset : {delay/1e6:.1f} Myr")
    if delay < 50e6:
        print(" ✅ AGN appears almost simultaneously with LRD.")
        print(" → Hybrid model (red continuum + narrow lines) is CO‑TEMPORAL.")
    else:
        print(" ✅ AGN appears, but with significant delay.")
else:
    print(" ❌ AGN did NOT reach 10⁵ M⊙ within sim time.")
    print(" → Try higher f_Edd, larger seed, or additional gas inflow.")

print(f"\n Final BH mass : {M_BH[-1]:.2e} M⊙")
print(f" Final stellar mass : {M_star[-1]:.2e} M⊙")
print(f" Remaining gas mass : {M_gas[-1]:.2e} M⊙")
print(f" Star formation efficiency : {M_star[-1] / (Mdot_core * t_onset):.2f}")
print("\n" + "=" * 60)
print("Note: Analytic model. Full hydrodynamic")
print("simulations are needed for confirmation.")
print("=" * 60)
//little worng
