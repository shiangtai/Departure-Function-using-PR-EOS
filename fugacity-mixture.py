#calculate fugacity based on the Peng-Robinson equation of state for mixtures
#by STLin 2025 (stlin@ntu.edu.tw)

import numpy as np

def fugacity(T, P, Tc, Pc, omega, kij, molx):
    R =8.314  # J/mol K, universal gas constant
    Tr = T / Tc  # reduced temperature
    Pr = P / Pc  # reduced pressure

    # Peng-Robinson parameters for each component
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - Tr**0.5))**2

    a = 0.45724 * (R**2 * Tc**2) *alpha / Pc # in unit of Pa m3/mol2
    b = 0.07780 * (R * Tc) / Pc  # in unit of m3/mol
    #print(a,b)

    # Mixing rules for a and b
    aij = np.sqrt(np.outer(a, a)) * (1 - kij)
    am = np.sum(molx[:, None] * molx[None, :] * aij)
    bm = np.sum(molx * b)    

    #print(am,bm)
    # Dimensionless parameters
    Aij = aij * P / (R**2 * T**2)
    Bi = b * P / (R * T)
    Am = am * P / (R**2 * T**2)
    Bm = bm * P / (R * T)

    # Calculate the compressibility factor Z using the cubic equation of state
    coeffs = [1, -(1 - Bm), Am - 3*Bm**2 - 2*Bm, -(Am*Bm - Bm**2 - Bm**3)]
    roots = np.roots(coeffs)
    Z = max(roots[roots.imag == 0].real)  # Select the largest real root

    # Departure function for Gibbs free energy (G)
    lnphi = (Bi/Bm) * (Z - 1) - np.log(Z - Bm) - (Am/(2*np.sqrt(2)*Bm)) * (2*(Aij @ molx)/Am - Bi/Bm) * np.log((Z + (1 + np.sqrt(2))*Bm) / (Z + (1 - np.sqrt(2))*Bm))
    fugacity = np.exp(lnphi)*P*molx
    #print("Fugacity:", fugacity)
    return fugacity


T=373.15 # Temperature in K
P=5e6  # Pressure in Pa
Tc=np.array([469.6, 562.1]) # list of critical temperatures in K for each component
Pc=np.array([3.374e6, 4.894e6]) # list of critical pressures in Pa for each component
omega=np.array([0.251, 0.212]) # list of acentric factors for each component
kij=np.array([[0.0, 0.018],[0.018, 0.0]]) # binary interaction parameters
molx=np.array([0.5, 0.5]) # mole fractions of each component

fu=fugacity(T, P, Tc, Pc, omega, kij, molx)

print("Fugacity Calculation using Peng-Robinson EOS by STLin")
print(f"Temperature (T): {T} K, Pressure (P): {P} Pa")
for i in range(len(molx)):
    print(f"Component {i+1}: Tc={Tc[i]} K, Pc={Pc[i]:.2e} Pa, omega={omega[i]}, mole fraction={molx[i]}, fugacity={fu[i]:.2e} Pa")

fua=[]
fub=[]
for x in np.linspace(0.0, 1.0, 11):
    molx=np.array([x, 1-x])
    fu=fugacity(T, P, Tc, Pc, omega, kij, molx)
    fua.append(fu[0])
    fub.append(fu[1])

import matplotlib.pyplot as plt
plt.plot(np.linspace(0.0, 1.0, 11), fua, label='Component 1 Fugacity')
plt.plot(np.linspace(0.0, 1.0, 11), fub, label='Component 2 Fugacity')
plt.xlabel('Mole Fraction of Component 1')
plt.ylabel('Fugacity (Pa)')
plt.title('Fugacity vs Mole Fraction at T=373.15 K, P=5 MPa')
plt.legend()
plt.grid()
plt.show()


