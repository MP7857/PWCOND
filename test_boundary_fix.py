#!/usr/bin/env python3
"""
Test script to verify the boundary condition fix improves accuracy.
This demonstrates that proper endpoint handling is critical for accuracy.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import jv
from scipy.integrate import simpson, quad
import matplotlib.pyplot as plt

# Configuration
L = 3
m_bessel = 0
z_slice = 0.5
g_values = np.linspace(0.1, 15.0, 30)

def beta_r(r):
    """Mock radial part of f-orbital projector"""
    return r * (r**L * np.exp(-r))

def integrand_smooth(r, z):
    """The non-oscillatory part: beta(r) * (1/r^3)"""
    r = np.asarray(r)
    return beta_r(r) / (r**3)

def exact_integrand(r, g, z):
    """Full integrand for reference calculation"""
    if r**2 < z**2 - 1e-12: return 0.0
    if abs(r**2 - z**2) < 1e-12:
        r_perp = 0.0
    else:
        r_perp = np.sqrt(r**2 - z**2)
    return integrand_smooth(r, z) * jv(m_bessel, g * r_perp)

# Grids
r_min, r_max, N_grid = 0.001, 20.0, 150
r_log = np.logspace(np.log10(r_min), np.log10(r_max), N_grid)

def integrate_with_boundary_fix(r_grid, smooth_func, z, g, m):
    """
    Improved method with proper boundary condition.
    Interpolates from (|z|, smooth_func(|z|)) instead of clamping at r(iz).
    """
    z_abs = abs(z)
    
    # Find first grid point > |z|
    mask = r_grid > z_abs
    if not mask.any():
        return 0.0
    
    iz = np.where(mask)[0][0]
    r_sub = r_grid[iz:]
    
    # Compute smooth function values on the grid
    y_smooth = smooth_func(r_sub)
    
    # Compute value at |z| by linear interpolation
    if iz > 0:
        # Interpolate between r[iz-1] and r[iz]
        r0, r1 = r_grid[iz-1], r_grid[iz]
        y0, y1 = smooth_func(r0), smooth_func(r1)
        val_at_z = y0 + (y1 - y0) * (z_abs - r0) / (r1 - r0)
    else:
        # Extrapolate if z < r[0]
        val_at_z = smooth_func(r_grid[0]) * (z_abs / r_grid[0])
    
    # Create fine grid starting exactly at |z|
    n_fine = max(500, int(g * (r_grid[-1] - z_abs) * 5.0))
    n_fine = min(n_fine, 5000)
    if n_fine % 2 == 0:
        n_fine += 1
    
    r_fine = np.linspace(z_abs, r_grid[-1], n_fine)
    
    # Interpolate smooth part onto fine grid
    # Need to handle the gap [z_abs, r[iz]]
    y_fine = np.zeros(n_fine)
    for i, r_val in enumerate(r_fine):
        if r_val < r_grid[iz]:
            # Interpolate from (z_abs, val_at_z) to (r[iz], y_smooth[0])
            y_fine[i] = val_at_z + (y_smooth[0] - val_at_z) * (r_val - z_abs) / (r_grid[iz] - z_abs)
        else:
            # Standard interpolation on the grid
            y_fine[i] = np.interp(r_val, r_sub, y_smooth)
    
    # Compute Bessel on fine grid
    r_perp_fine = np.sqrt(np.maximum(0, r_fine**2 - z_abs**2))
    oscillatory = jv(m, g * r_perp_fine)
    
    # Integrate
    result = simpson(y=y_fine * oscillatory, x=r_fine)
    return result

# Test both methods
errors_old = []
errors_new = []

print(f"{'g':<6} | {'Old Error':<15} | {'New Error':<15} | {'Improvement':<12}")
print("-" * 65)

for g in g_values:
    # Reference solution
    ref_val, _ = quad(exact_integrand, abs(z_slice), r_max, args=(g, z_slice), limit=100)
    
    # Old method (with clamping bug)
    mask = r_log > abs(z_slice)
    r_sub = r_log[mask]
    if len(r_sub) < 2:
        errors_old.append(0)
        errors_new.append(0)
        continue
    
    r_perp_sub = np.sqrt(r_sub**2 - z_slice**2)
    y_sub = integrand_smooth(r_sub, z_slice) * jv(m_bessel, g * r_perp_sub)
    old_val = simpson(y=y_sub, x=r_sub)
    
    # New method (with boundary fix)
    new_val = integrate_with_boundary_fix(r_log, 
                                          lambda r: integrand_smooth(r, z_slice),
                                          z_slice, g, m_bessel)
    
    err_old = abs(old_val - ref_val)
    err_new = abs(new_val - ref_val)
    
    errors_old.append(err_old)
    errors_new.append(err_new)
    
    if err_old > 0:
        improvement = err_old / err_new if err_new > 0 else np.inf
    else:
        improvement = 1.0
    
    print(f"{g:<6.2f} | {err_old:<15.2e} | {err_new:<15.2e} | {improvement:<12.1f}x")

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(g_values, errors_old, 'r-o', label='Old (clamping bug)', markersize=4)
plt.plot(g_values, errors_new, 'b-s', label='New (boundary fix)', markersize=4)
plt.yscale('log')
plt.xlabel('g vector magnitude (1/Bohr)')
plt.ylabel('Integration Absolute Error')
plt.title('Integration Accuracy: Before and After Boundary Fix')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('boundary_fix_comparison.png', dpi=150)
print("\nPlot saved to boundary_fix_comparison.png")
print(f"\nMedian improvement factor: {np.median([e1/e2 if e2 > 0 else 1 for e1, e2 in zip(errors_old, errors_new) if e1 > 0]):.1f}x")
