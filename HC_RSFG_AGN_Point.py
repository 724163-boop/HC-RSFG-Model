#!/usr/bin/env python3
"""
HC-RSFG → AGN Growth & Feedback → Core Quenching Simulation
=============================================================
Takes the HC‑RSFG core parameters and calculates:
  1. Black hole growth (seed → low‑mass AGN)
  2. AGN feedback: can the energy blow away the core gas?
No FSPS required – uses analytic accretion & feedback equations.
"""

import numpy as np

# ============================================================
# 1. HC‑RSFG Core Parameters (from earlier gas accumulation)
# ============================================================
M_core_gas = 1.26e9 # Msun (total gas mass in core)
R_core_pc = 10.0 # pc (inner core radius)
R_core_cm = R_core_pc * 3.086e18 # cm
G = 6.674e-8 # cm³ g⁻¹ s⁻²
Msun = 1.989e33 # g
c = 2.998e10 # cm/s

# ============================================================
# 2. Black Hole Seed & Growth
# ============================================================
M_seed = 100.0 # Msun (Pop III remnant or direct collapse seed)

# Time available for growth: from z~7.5 to z~5 (approx. 0.5 Gyr)
t_growth_yr = 5.0e8 # years

# Eddington parameters
epsilon = 0.1 # radiative efficiency
t_salpeter = 4.5e7 * (epsilon / 0.1) / (1.0 - epsilon) # ~5e7 yr

# Option 1: Eddington-limited growth
M_BH_eddington = M_seed * np.exp(t_growth_yr / t_salpeter)

# Option 2: Super‑Eddington growth (f_Edd > 1)
f_Edd = 3.0 # super‑Eddington factor
t_salpeter_super = t_salpeter / f_Edd
M_BH_super = M_seed * np.exp(t_growth_yr / t_salpeter_super)

print("=" * 60)
print("HC‑RSFG → AGN GROWTH SIMULATION")
print("=" * 60)
print(f" Core gas mass : {M_core_gas:.2e} Msun")
print(f" Core radius : {R_core_pc} pc")
print(f" BH seed mass : {M_seed} Msun")
print(f" Available time : {t_growth_yr/1e6:.0f} Myr")
print(f" Eddington growth → M_BH = {M_BH_eddington:.1e} Msun")
print(f" Super‑Eddington (x{f_Edd}) → M_BH = {M_BH_super:.1e} Msun")
print(f" Narrow-line LRD observed BH mass ≈ 1e5–1e6 Msun")
if M_BH_super >= 1e5:
    print(" ✅ Super‑Eddington growth can reach observed BH masses.")
else:
    print(" ⚠ Requires higher seed or more efficient accretion.")

# Use the super‑Eddington BH mass for feedback calculation
M_BH = M_BH_super

# ============================================================
# 3. AGN Feedback: Can the AGN blow away the core gas?
# ============================================================
# Energy released by accretion over a feedback timescale
# Use the gas mass that will be accreted during the BH growth
# Total accreted mass: M_BH - M_seed ≈ M_BH
M_acc = M_BH - M_seed
if M_acc < 0:
    M_acc = M_BH

E_AGN = epsilon * M_acc * Msun * c**2 # erg

# Binding energy of the remaining gas (assuming uniform sphere)
# The gas that hasn't been accreted yet is M_core_gas - M_acc (if any)
# If M_acc > M_core_gas, then all gas was consumed by accretion
M_gas_remaining = max(0.0, M_core_gas - M_acc)
# Binding energy: U_bind = (3/5) * G * M_total * M_gas / R
# M_total includes the growing BH? Actually, we consider the gas bound to the core (including BH and stars).
# For simplicity, use M_core_gas (the initial) as the gravitating mass that binds the remaining gas.
# More precise: M_total = M_BH + M_core_gas (stellar mass is negligible at this early stage).
M_total = M_BH + M_core_gas
U_bind = (3.0/5.0) * G * (M_total * Msun) * (M_gas_remaining * Msun) / R_core_cm # erg

print("\n" + "=" * 60)
print("AGN FEEDBACK – CAN IT UNBIND THE CORE GAS?")
print("=" * 60)
print(f" BH mass : {M_BH:.1e} Msun")
print(f" Accreted mass : {M_acc:.1e} Msun")
print(f" Remaining gas mass : {M_gas_remaining:.1e} Msun")
print(f" AGN total energy : {E_AGN:.2e} erg")
print(f" Gas binding energy : {U_bind:.2e} erg")
print(f" Energy ratio E_AGN/U_bind : {E_AGN/U_bind:.1f}")
if E_AGN > U_bind:
    print(" ✅ AGN energy exceeds binding energy → gas can be expelled.")
    print(" → Core becomes gas‑poor, star formation shuts down.")
    print(" → System may become extremely faint or invisible (Pathway 4).")
else:
    print(" ❌ AGN energy insufficient to blow away core gas.")
    print(" → Gas remains, star formation continues or BH grows further.")

print("\n" + "=" * 60)
print("DISCLAIMER: This is a simplified analytic model.")
print("Full hydrodynamic simulations are required to confirm.")
print("=" * 60)
