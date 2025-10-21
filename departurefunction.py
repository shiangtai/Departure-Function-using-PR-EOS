#calculate departure function based on the Peng-Robinson equation of state
#by STLin 2025 (stlin@ntu.edu.tw)

import numpy as np

def departure_function(T, P, Tc, Pc, omega):
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
    Z = max(roots[roots.imag == 0].real)  # Select the largest real root

    dadT = -0.45724 * (R**2 * Tc**2) / Pc * kappa * (alpha/T/Tc)**0.5

    # Departure function for enthalpy (H) and entropy (S)
    delta_H = R * T * (Z - 1) + ( (T*dadT -a) / (2*np.sqrt(2) * b)) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))
    delta_S = R * np.log(Z - B) + (dadT / (2 * np.sqrt(2) * b)) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))

    return delta_H, delta_S, Z, A, B, coeffs, a, b, dadT

def idealgas(T, P, Cp):
    # Ideal gas enthalpy and entropy changes from reference state
    R = 8.314  # J/mol K
    Tref = 298.15  # Reference temperature in K
    Pref = 1e5  # Reference pressure in Pa
    
    dH_ideal = Cp[0]*(T - Tref) + Cp[1]/2*(T**2 - Tref**2) + Cp[2]/3*(T**3 - Tref**3) + Cp[3]/4*(T**4 - Tref**4)
    dS_ideal = Cp[0]*np.log(T/Tref) + Cp[1]*(T - Tref) + Cp[2]/2*(T**2 - Tref**2) + Cp[3]/3*(T**3 - Tref**3) - R*np.log(P/Pref)
    return dH_ideal, dS_ideal



T=600 # Temperature in K
P=1e7  # Pressure in Pa
Tc=154.6 # Critical temperature in K for a specific substance
Pc=5.046e6 # Critical pressure in Pa for a specific substance
omega=0.021 # Acentric factor for a specific substance
Cp_coeffs=[25.46, 1.519e-2, -0.715e-5, 1.311e-9] #Cp coefficients  for O2

dH,dS,Z,A,B,coeffs,a,b,dadT=departure_function(T, P, Tc, Pc, omega)
dH_ideal, dS_ideal = idealgas(T, P, Cp_coeffs)

print("Departure Function Calculation using Peng-Robinson EOS by STLin")
print(f"Temperature (T): {T} K, Pressure (P): {P} Pa")
print(f"Critical Temperature (Tc): {Tc} K, Critical Pressure (Pc): {Pc} Pa, Acentric Factor (omega): {omega}")
print(f"Calculated a: {a} Pa m3/mol2, b: {b} m3/mol, dadT: {dadT} Pa m3/(mol2 K)")
print(f"Peng-Robinson EOS Parameters: A={A}, B={B}")
print(f"Peng-Robinson EOS Coefficients: {coeffs}")
print(f"Compressibility Factor (Z): {Z}")
print(f"Departure Enthalpy (dH): {dH} J/mol")
print(f"Departure Entropy (dS): {dS} J/(mol·K)")
print(f"Ideal Gas Enthalpy Change (dH_ideal): {dH_ideal} J/mol")
print(f"Ideal Gas Entropy Change (dS_ideal): {dS_ideal} J/(mol·K)")

H=dH + dH_ideal
S=dS + dS_ideal
print(f"Total Enthalpy (H): {H} J/mol")
print(f"Total Entropy (S): {S} J/(mol·K)")
