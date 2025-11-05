#calculate fugacity based on the Peng-Robinson equation of state
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
    #print(a,b)
    A = a * P / (R**2 * T**2)
    B = b * P / (R * T)

    # Calculate the compressibility factor Z using the cubic equation of state
    coeffs = [1, -(1 - B), A - 3*B**2 - 2*B, -(A*B - B**2 - B**3)]
    roots = np.roots(coeffs)
    Z = max(roots[roots.imag == 0].real)  # Select the largest real root

 
    # Departure function for Gibbs free energy (G)
    delta_G1_RT= (Z - 1) - np.log(Z - B) - A / (2 * np.sqrt(2) * B) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)) 
    fugacity = np.exp(delta_G1_RT)*P
    
    return fugacity


T=373.15 # Temperature in K
P=1e6  # Pressure in Pa
Tc=425.2 # Critical temperature in K for butane
Pc=3.8e6 # Critical pressure in Pa for butane
omega=0.193 # Acentric factor for butane


fu=fugacity(T, P, Tc, Pc, omega)

print("Fugacity Calculation using Peng-Robinson EOS by STLin")
print(f"Temperature (T): {T} K, Pressure (P): {P} Pa")
print(f"Critical Temperature (Tc): {Tc} K, Critical Pressure (Pc): {Pc} Pa, Acentric Factor (omega): {omega}")
print(f"Calculated Fugacity: {fu} Pa")

