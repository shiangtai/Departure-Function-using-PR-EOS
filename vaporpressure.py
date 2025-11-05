#calculate vapor pressure based on the Peng-Robinson equation of state
#by STLin 2025 (stlin@ntu.edu.tw)

import numpy as np

def fugacity(T, P, Tc, Pc, omega):
    R =8.314  # J/mol K, universal gas constant
    Tr = T / Tc  # reduced temperature
    Pr = P / Pc  # reduced pressure

    # Peng-Robinson parameters

    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - Tr**0.5))**2

    a = 0.45724 * (R**2 * Tc**2) *alpha / Pc # in unit of Pa m3/mol2
    b = 0.07780 * (R * Tc) / Pc  # in unit of m3/mol

    A = a * P / (R**2 * T**2)
    B = b * P / (R * T)

    # Calculate the compressibility factor Z using the cubic equation of state
    coeffs = [1, -(1 - B), A - 3*B**2 - 2*B, -(A*B - B**2 - B**3)]
    roots = np.roots(coeffs)
    Zv = max(roots[roots.imag == 0].real)  # Select the largest real root
    Zl = min(roots[roots.imag == 0].real)  # Select the smallest real root

 
    # Departure function for Gibbs free energy (G)
    Z=Zv  # for vapor phase fugacity
    delta_G1_RT= (Z - 1) - np.log(Z - B) - A / (2 * np.sqrt(2) * B) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)) 
    fugv = np.exp(delta_G1_RT)*P
    Z=Zl  # for liquid phase fugacity
    delta_G1_RT= (Z - 1) - np.log(Z - B) - A / (2 * np.sqrt(2) * B) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)) 
    fugl = np.exp(delta_G1_RT)*P
    return fugv, fugl


T=298.15 # Temperature in K
P=1e4  # Guess for vapor pressure in Pa
Tc=647.3 # Critical temperature in K for water
Pc=2.2e7 # Critical pressure in Pa for water
omega=0.344 # Acentric factor for water


fugv,fugl=fugacity(T, P, Tc, Pc, omega)
while abs(1 - fugl/fugv) > 1e-6:
    P *= fugl / fugv
    fugv, fugl = fugacity(T, P, Tc, Pc, omega)


print("Fugacity Calculation using Peng-Robinson EOS by STLin")
print(f"Temperature (T): {T} K, Pressure (P): {P} Pa")
print(f"Critical Temperature (Tc): {Tc} K, Critical Pressure (Pc): {Pc} Pa, Acentric Factor (omega): {omega}")
print(f"Calculated Fugacity: (Vapor) {fugv} Pa, (Liquid) {fugl} Pa")
