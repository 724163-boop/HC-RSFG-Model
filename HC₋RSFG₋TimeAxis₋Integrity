import numpy as np

# Physical constants
Msun       = 1.989e33        # g, solar mass
pc_to_cm   = 3.086e18        # cm per parsec
m_p        = 1.673e-24       # g, proton mass
mu         = 1.22            # mean molecular weight
sec_per_yr = 365.25 * 24 * 3600

# Parameters (observation- and theory-based)
R_inner_pc = 10.0            # inner core radius [pc]
n_initial  = 1e4             # initial density [cm^-3]  (after Pop III SNe)
n_crit     = 1e7             # critical density [cm^-3]  (RSFG onset condition)
Mdot_total = 20.0            # cold stream inflow rate [Msun/yr] (Dekel+09)
f_conc     = 0.3             # fraction of inflow reaching the inner core

# Volume and mass calculations
V_inner    = (4./3.) * np.pi * (R_inner_pc * pc_to_cm)**3
M_init     = mu * m_p * n_initial * V_inner   # initial gas mass
M_target   = mu * m_p * n_crit    * V_inner   # target gas mass
dM         = M_target - M_init                # additional mass needed

# Net mass inflow rate into the inner core
Mdot_inner = Mdot_total * f_conc * Msun       # [g/yr]
t_acc_yr   = dM / Mdot_inner                  # accumulation time [yr]
t_acc_Gyr  = t_acc_yr / 1e9                   # accumulation time [Gyr]

print(f"Time to reach critical density: {t_acc_Gyr:.2f} Gyr = {t_acc_Gyr*1000:.0f} Myr")
print(f"Starting from z=10 (cosmic age ~0.48 Gyr) → RSFG onset at z≈7.5 (age ~0.69 Gyr)")
