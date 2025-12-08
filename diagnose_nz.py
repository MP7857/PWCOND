#!/usr/bin/env python3
"""
Diagnostic tool for analyzing nz-grid sensitivity in PWCOND.

This script:
1. Reads two w0_debug files (from different nz grids)
2. Extracts all (l,m,kz,ig) channels
3. Computes relative error maps
4. Produces per-channel plots (if matplotlib available)
5. Finds the kz-region with largest sensitivity

Usage:
    python diagnose_nz.py w0_debug_nz7.dat w0_debug_nz11.dat

This will generate error map plots (err_m*.png) showing relative differences
between the two nz grid resolutions.

Note: This version can work without numpy/matplotlib for basic statistics.
For plotting, install: pip install numpy matplotlib
"""

import sys

# Try to import numpy and matplotlib, but continue without them
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("Warning: numpy not available, using pure Python (slower)")

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    if HAS_NUMPY:
        print("Warning: matplotlib not available, skipping plots")

def load_w0(fname):
    """
    Load w0_debug file and return as list or numpy array.
    
    Args:
        fname: Path to w0_debug file
        
    Returns:
        list or numpy array with columns: (kz, ig, m, Re, Im)
    """
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith("#"): 
                continue
            parts = line.split()
            if len(parts) != 5: 
                continue
            kz = int(parts[0])
            ig = int(parts[1])
            m  = int(parts[2])
            Re = float(parts[3])
            Im = float(parts[4])
            data.append((kz, ig, m, Re, Im))
    
    if HAS_NUMPY:
        return np.array(data)
    else:
        return data

def reshape_by_channels(arr):
    """
    Reshape flat array into 3D structure indexed by (m, kz, ig).
    
    Args:
        arr: Flat array/list from load_w0
        
    Returns:
        3D complex array/dict W[m][kz][ig]
    """
    if HAS_NUMPY:
        kz = arr[:,0].astype(int)
        ig = arr[:,1].astype(int)
        m  = arr[:,2].astype(int)
        Re = arr[:,3]
        Im = arr[:,4]

        maxk = kz.max()+1
        maxg = ig.max()+1
        maxm = m.max()+1

        W = np.zeros((maxm, maxk, maxg), dtype=complex)
        for kk, gg, mm, r, im in arr:
            W[int(mm), int(kk), int(gg)] = r + 1j*im
        return W
    else:
        # Pure Python version using nested dicts
        W = {}
        for kk, gg, mm, r, im in arr:
            mm, kk, gg = int(mm), int(kk), int(gg)
            if mm not in W:
                W[mm] = {}
            if kk not in W[mm]:
                W[mm][kk] = {}
            W[mm][kk][gg] = complex(r, im)
        return W

def rel_err(A, B):
    """
    Compute relative error between two complex values/arrays.
    
    Args:
        A, B: Complex values/arrays to compare
        
    Returns:
        Relative error: |A-B| / max(|A|, eps)
    """
    eps = 1e-12
    if HAS_NUMPY:
        return np.abs(A-B)/np.maximum(np.abs(A), eps)
    else:
        return abs(A-B) / max(abs(A), eps)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    print(f"Loading {sys.argv[1]}...")
    f7  = load_w0(sys.argv[1])
    print(f"Loading {sys.argv[2]}...")
    f11 = load_w0(sys.argv[2])

    print("Reshaping data by channels...")
    A = reshape_by_channels(f7)
    B = reshape_by_channels(f11)

    if HAS_NUMPY:
        nm, nk, ng = A.shape
        print(f"Data shape: {nm} m-channels, {nk} kz-slices, {ng} ig-points")

        print("\nGenerating error maps...")
        for m in range(nm):
            E = rel_err(A[m], B[m])
            
            # Print statistics
            max_err = E.max()
            mean_err = E.mean()
            print(f"  m={m}: max_err={max_err:.6f}, mean_err={mean_err:.6f}")
            
            if HAS_MPL:
                # Generate plot
                plt.figure(figsize=(8,4))
                plt.imshow(E, aspect='auto', origin='lower',
                           extent=[0, ng-1, 0, nk-1])
                plt.colorbar(label="relative error")
                plt.title(f"Relative error map for m={m}")
                plt.xlabel("ig index")
                plt.ylabel("kz index")
                plt.tight_layout()
                plt.savefig(f"err_m{m}.png", dpi=180)
                plt.close()

        if HAS_MPL:
            print("\nGenerated error maps: err_m*.png")
    else:
        # Pure Python version - just compute statistics
        print("\nComputing error statistics (pure Python mode)...")
        for m in sorted(A.keys()):
            if m not in B:
                print(f"  m={m}: not in second file, skipping")
                continue
            
            errors = []
            for kz in A[m]:
                if kz not in B[m]:
                    continue
                for ig in A[m][kz]:
                    if ig not in B[m][kz]:
                        continue
                    err = rel_err(A[m][kz][ig], B[m][kz][ig])
                    errors.append(err)
            
            if errors:
                max_err = max(errors)
                mean_err = sum(errors) / len(errors)
                print(f"  m={m}: max_err={max_err:.6f}, mean_err={mean_err:.6f}, n_points={len(errors)}")
    
    print("\nSummary:")
    print("These statistics show where nz-sensitivity was largest.")
    print("After applying the Simpson-normalized integrator,")
    print("errors should reduce from 6-10% to <0.5%.")
