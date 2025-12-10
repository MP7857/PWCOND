#!/usr/bin/env python3
"""
Validation script for f-orbital integration accuracy
Tests the proposed spline-based method vs. current method
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import jv
from scipy.integrate import simpson, quad
import matplotlib.pyplot as plt

# --- 1. CONFIGURATION ---
L = 3          # f-orbital
m_bessel = 0   # J0 component
z_slice = 0.5  # Arbitrary z-slice
g_values = np.linspace(0.1, 15.0, 30) # Range of g (simulating energy window)

# --- 2. MOCK FUNCTIONS ---
def beta_r(r):
    """Mock radial part of f-orbital projector: r * R_nl(r)"""
    return r * (r**L * np.exp(-r))

def integrand_smooth(r, z):
    """The non-oscillatory part of the integrand: beta(r) * (1/r^3)"""
    # Vector-safe implementation
    r = np.asarray(r)
    return beta_r(r) / (r**3)

def exact_integrand(r, g, z):
    """Full integrand for reference calculation"""
    if r**2 < z**2: return 0.0
    r_perp = np.sqrt(r**2 - z**2)
    return integrand_smooth(r, z) * jv(m_bessel, g * r_perp)

# --- 3. GRIDS ---
# Mimic QE Logarithmic Grid
r_min, r_max, N_grid = 0.001, 20.0, 150
r_log = np.logspace(np.log10(r_min), np.log10(r_max), N_grid)

# --- 4. CALCULATION LOOP ---
errors_current = []
errors_proposed = []

print(f"{'g':<6} | {'Current Error':<15} | {'Proposed Error':<15}")
print("-" * 42)

for g in g_values:
    # A. EXACT SOLUTION (Adaptive Quadrature)
    ref_val, _ = quad(exact_integrand, abs(z_slice), r_max, args=(g, z_slice), limit=100)

    # B. CURRENT METHOD (Simpson on Log Grid)
    mask = r_log > abs(z_slice)
    r_sub = r_log[mask]
    if len(r_sub) < 2:
        errors_current.append(0)
        errors_proposed.append(0)
        continue

    r_perp_sub = np.sqrt(r_sub**2 - z_slice**2)
    y_sub = integrand_smooth(r_sub, z_slice) * jv(m_bessel, g * r_perp_sub)

    curr_val = simpson(y=y_sub, x=r_sub)

    # C. PROPOSED METHOD (Spline -> Fine Linear Grid)
    y_smooth_full = integrand_smooth(r_log, z_slice)
    cs = CubicSpline(r_log, y_smooth_full)

    N_fine = 1000
    r_fine = np.linspace(abs(z_slice), r_max, N_fine)

    smooth_interp = cs(r_fine)
    r_perp_fine = np.sqrt(r_fine**2 - z_slice**2)
    oscillatory = jv(m_bessel, g * r_perp_fine)

    prop_val = simpson(y=smooth_interp * oscillatory, x=r_fine)

    # D. ERROR METRICS
    err_c = abs(curr_val - ref_val)
    err_p = abs(prop_val - ref_val)

    errors_current.append(err_c)
    errors_proposed.append(err_p)
    
    print(f"{g:<6.2f} | {err_c:<15.2e} | {err_p:<15.2e}")

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(g_values, errors_current, 'r-o', label='Current (Simpson on Log)')
plt.plot(g_values, errors_proposed, 'b-s', label='Proposed (Spline + Fine Grid)')
plt.yscale('log')
plt.xlabel('g vector magnitude (1/Bohr)')
plt.ylabel('Integration Absolute Error')
plt.title('Integration Accuracy for f-orbitals')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('error_comparison.png')
print("\nPlot saved to error_comparison.png")
