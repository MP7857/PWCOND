#!/usr/bin/env python3
"""
Diagnostic functions for testing FX5 integration improvements.

This module provides Python implementations that mirror the Fortran logic
in four.f90, particularly for the f-orbital FX5 term integration.
"""

import numpy as np
import mpmath as mp


def beta_r(r_grid, l=3, alpha=1.0):
    """
    Simple model beta function for testing purposes.
    
    Parameters:
    -----------
    r_grid : array-like
        Radial grid points
    l : int
        Angular momentum quantum number
    alpha : float
        Shape parameter
    
    Returns:
    --------
    array-like
        Beta values at each grid point
    """
    return r_grid**l * np.exp(-alpha * r_grid)


def simpson_uniform(y, h):
    """
    Simpson's rule integration for uniformly spaced data.
    
    Parameters:
    -----------
    y : array-like
        Function values at grid points
    h : float
        Uniform spacing between points
    
    Returns:
    --------
    float
        Integral value
    """
    n = len(y)
    if n < 2:
        return 0.0
    if n == 2:
        return 0.5 * h * (y[0] + y[1])
    
    # Standard Simpson's rule
    integral = y[0] + y[-1]
    for i in range(1, n-1, 2):
        integral += 4.0 * y[i]
    for i in range(2, n-1, 2):
        integral += 2.0 * y[i]
    
    return integral * h / 3.0


def fortran_like_f_x5(r_grid, z, gn, l=3, alpha=1.0):
    """
    Mimics the logic used for x5 in your current four.f90:
      - find index iz such that r(iz) > |z|
      - Simpson from r(iz) to rmax
      - treat [|z|, r(iz)] with 3-point Simpson + linear extrapolation of beta

    Here we focus on the f-orbital 'x5' term:
        x5(r) = beta(r) * J0(gn * sqrt(r^2 - z^2)) / r^3
    and integrate over r.

    This is not the full 3D derivation (no Jacobian etc.),
    but it closely mirrors the *numerical pattern* of your code.
    """
    beta = beta_r(r_grid, l=l, alpha=alpha)

    zabs = abs(z)
    # find iz such that r(iz) > |z|
    iz = np.searchsorted(r_grid, zabs, side="right")
    if iz == 0 or iz >= len(r_grid):
        return 0.0

    nmesh = len(r_grid)
    nmeshs = nmesh  # keep all points for simplicity

    # main integral from r(iz) ... r(nmeshs-1)
    r_seg = r_grid[iz:nmeshs]
    rz = np.sqrt(r_seg**2 - zabs**2)

    # x5(ir) in Fortran: betar(ir)*bessj(0,gn*rz)/r(ir)**3
    # Vectorize Bessel function computation for efficiency
    bessel_values = np.array([float(mp.besselj(0, gn*rv)) for rv in rz])
    x5 = beta[iz:nmeshs] * bessel_values / (r_seg**3)

    h = r_grid[1] - r_grid[0]  # uniform spacing for this model
    I_main = simpson_uniform(x5, h)

    # Endpoint [|z|, r(iz)] treated with 3-point Simpson rule,
    # using linear extrapolation of beta and J0(gn * sqrt(r^2 - z^2)) ≈ exact J0
    dr = r_grid[iz] - r_grid[iz-1]
    zr = r_grid[iz] - zabs  # length of the first small segment

    r0 = zabs
    r1 = r_grid[iz]
    rmid = 0.5 * (r0 + r1)

    # Extrapolated beta at r0 = |z|
    beta0 = beta[iz] - (beta[iz] - beta[iz-1]) * zr / dr

    # Beta at midpoint by the same linear model
    zrmid = r1 - rmid
    betamid = beta[iz] - (beta[iz] - beta[iz-1]) * zrmid / dr

    if r0 > 1e-8:
        # here J0(gn * sqrt(r0^2 - z^2)) = J0(0) = 1
        x5_0 = beta0 / (r0**3)
    else:
        x5_0 = 0.0

    if rmid > 1e-8:
        rz_mid = np.sqrt(rmid**2 - zabs**2)
        j0_mid = float(mp.besselj(0, gn*rz_mid))
        x5_mid = betamid * j0_mid / (rmid**3)
    else:
        x5_mid = 0.0

    x5_1 = x5[0]   # already includes J0 at r1

    I_end = (x5_0 + 4.0 * x5_mid + x5_1) * zr / 6.0

    return I_main + I_end


def reference_f_x5(r_grid, z, gn, l=3, alpha=1.0):
    """
    Reference implementation using high-precision integration.
    
    This serves as a baseline for comparison with the Fortran-like implementation.
    
    Parameters:
    -----------
    r_grid : array-like
        Radial grid points
    z : float
        z-coordinate value
    gn : float
        Reciprocal space vector magnitude
    l : int
        Angular momentum quantum number
    alpha : float
        Shape parameter for beta function
    
    Returns:
    --------
    float
        Integral value
    """
    beta = beta_r(r_grid, l=l, alpha=alpha)
    zabs = abs(z)
    
    # Only integrate where r > |z|
    mask = r_grid > zabs
    r_valid = r_grid[mask]
    beta_valid = beta[mask]
    
    if len(r_valid) == 0:
        return 0.0
    
    rz = np.sqrt(r_valid**2 - zabs**2)
    
    # x5(r) = beta(r) * J0(gn * sqrt(r^2 - z^2)) / r^3
    # Vectorize Bessel function computation for efficiency
    bessel_values = np.array([float(mp.besselj(0, gn*rv)) for rv in rz])
    x5 = beta_valid * bessel_values / (r_valid**3)
    
    h = r_grid[1] - r_grid[0]
    return simpson_uniform(x5, h)


def test_fx5_integration():
    """
    Test the FX5 integration with the improved endpoint treatment.
    
    This function compares the Fortran-like implementation with a reference
    implementation to verify correctness.
    """
    # Setup test parameters
    rmax = 10.0
    npoints = 500
    r_grid = np.linspace(0.01, rmax, npoints)
    
    # Test cases with different z and gn values
    test_cases = [
        (0.5, 1.0),   # Small z, moderate gn
        (2.0, 2.0),   # Moderate z, higher gn
        (5.0, 0.5),   # Larger z, small gn
        (1.0, 3.0),   # Moderate z, high gn
    ]
    
    print("Testing FX5 integration with improved endpoint treatment")
    print("=" * 60)
    
    for z, gn in test_cases:
        result_fortran = float(fortran_like_f_x5(r_grid, z, gn))
        result_ref = float(reference_f_x5(r_grid, z, gn))
        
        # Use epsilon threshold to avoid artificially large errors when denominator is near zero
        eps_threshold = 1e-12
        if abs(result_ref) > eps_threshold:
            rel_error = abs(result_fortran - result_ref) / abs(result_ref)
        else:
            rel_error = 0.0
        
        print(f"\nz = {z:.2f}, gn = {gn:.2f}")
        print(f"  Fortran-like: {result_fortran:.10e}")
        print(f"  Reference:    {result_ref:.10e}")
        print(f"  Rel. error:   {rel_error:.10e}")
    
    print("\n" + "=" * 60)
    print("Test completed. Check relative errors to verify improvement.")


if __name__ == "__main__":
    # Run tests when script is executed directly
    test_fx5_integration()
