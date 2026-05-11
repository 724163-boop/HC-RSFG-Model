import numpy as np

G   = 6.67430e-8   # cm^3 g^-1 s^-2
k_B = 1.380649e-16 # erg K^-1
m_p = 1.67262192e-24 # g
mu  = 1.22
h   = 6.62607015e-27 # erg s
c   = 2.99792458e10  # cm s^-1


n           = 1e7
T_gas       = 23.0
Z_sun_ratio = 1e-4


rho    = mu * m_p * n
t_ff   = np.sqrt(3.0 * np.pi / (32.0 * G * rho))
t_ff_yr = t_ff / (365.25 * 24 * 3600)


A_CII, gamma_CII, E_CII_kB = 2.3e-6, 7.0e-10, 91.25
n_crit_CII = A_CII / gamma_CII
nu_CII     = c / 157.74e-4
XC         = 2.69e-4 * Z_sun_ratio
boltz_CII  = np.exp(-E_CII_kB / T_gas)
L_CII = A_CII * h * nu_CII * (4./2.) * XC * boltz_CII \
        / (1. + n_crit_CII/n + (4./2.)*boltz_CII)


A_OI, gamma_OI, E_OI_kB = 8.9e-5, 4.0e-10, 227.72
n_crit_OI = A_OI / gamma_OI
nu_OI     = c / 63.18e-4
XO        = 4.90e-4 * Z_sun_ratio
boltz_OI  = np.exp(-E_OI_kB / T_gas)
L_OI = A_OI * h * nu_OI * (3./5.) * XO * boltz_OI \
       / (1. + n_crit_OI/n + (3./5.)*boltz_OI)


Lambda_total = n**2 * (L_CII + L_OI)
thermal      = 1.5 * n * k_B * T_gas
t_cool_yr    = (thermal / Lambda_total) / (365.25 * 24 * 3600)

print(f"t_ff  = {t_ff_yr:.2e} yr")
print(f"t_cool = {t_cool_yr:.2e} yr")
print(f"t_cool / t_ff = {t_cool_yr/t_ff_yr:.2e}")

if t_cool_yr < t_ff_yr:
    print("=> 冷却優勢。RSFGモードは物理的に発動可能。")
