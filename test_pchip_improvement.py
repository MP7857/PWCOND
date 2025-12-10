#!/usr/bin/env python3
"""
Test script to validate PCHIP/cubic spline improvements for f-orbital integration
"""

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
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

def integrate_with_cubic_spline(r_grid, smooth_func, z, g, m, n_fine_factor=40):
    """
    Method using cubic spline interpolation with high fine-grid density
    """
    z_abs = abs(z)
    
    mask = r_grid > z_abs
    if not mask.any():
        return 0.0
    
    iz = np.where(mask)[0][0]
    r_sub = r_grid[iz:]
    
    # Compute smooth function values
    y_smooth = smooth_func(r_sub)
    
    # Compute value at |z|
    if iz > 0:
        r0, r1 = r_grid[iz-1], r_grid[iz]
        y0, y1 = smooth_func(r0), smooth_func(r1)
        val_at_z = y0 + (y1 - y0) * (z_abs - r0) / (r1 - r0)
    else:
        val_at_z = smooth_func(r_grid[0]) * (z_abs / r_grid[0])
    
    # Create cubic spline with boundary point
    r_nodes = np.concatenate([[z_abs], r_sub])
    y_nodes = np.concatenate([[val_at_z], y_smooth])
    cs = CubicSpline(r_nodes, y_nodes)
    
    # High-density fine grid (40 points per Bessel period)
    n_fine = max(500, int(n_fine_factor * g * (r_grid[-1] - z_abs)))
    n_fine = min(n_fine, 10000)
    if n_fine % 2 == 0:
        n_fine += 1
    
    r_fine = np.linspace(z_abs, r_grid[-1], n_fine)
    
    # Interpolate and compute Bessel
    y_fine = cs(r_fine)
    r_perp_fine = np.sqrt(np.maximum(0, r_fine**2 - z_abs**2))
    oscillatory = jv(m, g * r_perp_fine)
    
    # Integrate
    result = simpson(y=y_fine * oscillatory, x=r_fine)
    return result

def integrate_with_pchip(r_grid, smooth_func, z, g, m, n_fine_factor=40):
    """
    Method using PCHIP (monotonic) interpolation
    """
    z_abs = abs(z)
    
    mask = r_grid > z_abs
    if not mask.any():
        return 0.0
    
    iz = np.where(mask)[0][0]
    r_sub = r_grid[iz:]
    
    # Compute smooth function values
    y_smooth = smooth_func(r_sub)
    
    # Compute value at |z|
    if iz > 0:
        r0, r1 = r_grid[iz-1], r_grid[iz]
        y0, y1 = smooth_func(r0), smooth_func(r1)
        val_at_z = y0 + (y1 - y0) * (z_abs - r0) / (r1 - r0)
    else:
        val_at_z = smooth_func(r_grid[0]) * (z_abs / r_grid[0])
    
    # Create PCHIP interpolator with boundary point
    r_nodes = np.concatenate([[z_abs], r_sub])
    y_nodes = np.concatenate([[val_at_z], y_smooth])
    pchip = PchipInterpolator(r_nodes, y_nodes)
    
    # High-density fine grid
    n_fine = max(500, int(n_fine_factor * g * (r_grid[-1] - z_abs)))
    n_fine = min(n_fine, 10000)
    if n_fine % 2 == 0:
        n_fine += 1
    
    r_fine = np.linspace(z_abs, r_grid[-1], n_fine)
    
    # Interpolate and compute Bessel
    y_fine = pchip(r_fine)
    r_perp_fine = np.sqrt(np.maximum(0, r_fine**2 - z_abs**2))
    oscillatory = jv(m, g * r_perp_fine)
    
    # Integrate
    result = simpson(y=y_fine * oscillatory, x=r_fine)
    return result

# Test both methods
errors_linear = []  # Previous method (linear interp, 5x factor)
errors_cubic = []   # New: cubic spline, 40x factor
errors_pchip = []   # New: PCHIP, 40x factor

print(f"{'g':<6} | {'Linear (old)':<15} | {'Cubic (new)':<15} | {'PCHIP (new)':<15}")
print("-" * 68)

for g in g_values:
    # Reference solution
    ref_val, _ = quad(exact_integrand, abs(z_slice), r_max, args=(g, z_slice), limit=100)
    
    # Linear interpolation (old method - factor 5)
    linear_val = integrate_with_cubic_spline(r_log, 
                                            lambda r: integrand_smooth(r, z_slice),
                                            z_slice, g, m_bessel, n_fine_factor=5)
    
    # Cubic spline (new method - factor 40)
    cubic_val = integrate_with_cubic_spline(r_log, 
                                           lambda r: integrand_smooth(r, z_slice),
                                           z_slice, g, m_bessel, n_fine_factor=40)
    
    # PCHIP (new method - factor 40)
    pchip_val = integrate_with_pchip(r_log, 
                                     lambda r: integrand_smooth(r, z_slice),
                                     z_slice, g, m_bessel, n_fine_factor=40)
    
    err_linear = abs(linear_val - ref_val)
    err_cubic = abs(cubic_val - ref_val)
    err_pchip = abs(pchip_val - ref_val)
    
    errors_linear.append(err_linear)
    errors_cubic.append(err_cubic)
    errors_pchip.append(err_pchip)
    
    print(f"{g:<6.2f} | {err_linear:<15.2e} | {err_cubic:<15.2e} | {err_pchip:<15.2e}")

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(g_values, errors_linear, 'r-o', label='Linear interp (5x)', markersize=4)
plt.plot(g_values, errors_cubic, 'b-s', label='Cubic spline (40x)', markersize=4)
plt.plot(g_values, errors_pchip, 'g-^', label='PCHIP (40x)', markersize=4)
plt.yscale('log')
plt.xlabel('g vector magnitude (1/Bohr)')
plt.ylabel('Integration Absolute Error')
plt.title('Integration Accuracy: PCHIP/Cubic Spline vs Linear')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pchip_spline_comparison.png', dpi=150)
print("\nPlot saved to pchip_spline_comparison.png")

median_improve_cubic = np.median([e1/e2 if e2 > 0 else 1 for e1, e2 in zip(errors_linear, errors_cubic) if e1 > 0])
median_improve_pchip = np.median([e1/e2 if e2 > 0 else 1 for e1, e2 in zip(errors_linear, errors_pchip) if e1 > 0])

print(f"\nMedian improvement (cubic spline): {median_improve_cubic:.1f}x")
print(f"Median improvement (PCHIP): {median_improve_pchip:.1f}x")
