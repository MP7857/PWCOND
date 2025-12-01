#!/usr/bin/env python3
"""
Analysis script for w0_debug.dat output from PWCOND four.f90

This script analyzes the w0 (Fourier transform of beta functions) values
to check sign correctness and understand the role of z position.

Usage:
    python analyze_w0.py [w0_debug.dat]

The script will:
1. Parse the w0_debug.dat file
2. Check sign patterns for each orbital type
3. Analyze z-dependence of w0 values
4. Generate summary statistics and diagnostics
"""

import sys
import re
import numpy as np
from collections import defaultdict

# Orbital type names
ORBITAL_NAMES = {0: 's', 1: 'p', 2: 'd', 3: 'f'}

# Expected m-values for each orbital type (index 1-based as in Fortran)
# m index -> (actual m quantum number, name)
M_VALUES = {
    0: {1: (0, 'm=0')},           # s: 1 component (m=0)
    1: {1: (0, 'pz'), 2: (-1, 'p-x'), 3: (-1, 'p-y')},  # p: 3 components
    2: {1: (0, 'dz2'), 2: (-1, 'd-xz'), 3: (-1, 'd-yz'), 
        4: (2, 'dx2-y2'), 5: (2, 'dxy')},  # d: 5 components
    3: {1: (0, 'fz3'), 2: (1, 'fx(5z2-r2)'), 3: (-1, 'f-y(5z2-r2)'), 
        4: (2, 'fz(x2-y2)'), 5: (-2, 'f-xyz'), 
        6: (3, 'fx(x2-3y2)'), 7: (-3, 'f-y(3x2-y2)')}  # f: 7 components
}

# Number of m-components for each orbital type
M_COUNT = {0: 1, 1: 3, 2: 5, 3: 7}

def parse_w0_debug(filename):
    """Parse w0_debug.dat file and return structured data."""
    data = []
    current_orbital = None
    current_z0 = None
    current_dz = None
    skipped_lines = 0
    
    with open(filename, 'r') as f:
        for line in f:
            original_line = line
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Parse header information
            if line.startswith('# Orbital lb ='):
                match = re.search(r'lb = \s*(\d+)', line)
                if match:
                    current_orbital = int(match.group(1))
            elif line.startswith('# z0 ='):
                match = re.search(r'z0 = \s*([-\d.E+]+)', line, re.IGNORECASE)
                if match:
                    current_z0 = float(match.group(1))
            elif line.startswith('# dz ='):
                match = re.search(r'dz = \s*([-\d.E+]+)', line, re.IGNORECASE)
                if match:
                    current_dz = float(match.group(1))
            elif line.startswith('#'):
                continue
            else:
                # Handle Fortran exponential format (D instead of E)
                line = line.replace('D', 'E').replace('d', 'e')
                
                # Skip corrupted lines (lines that don't start with a number or are too short)
                if not line or not line[0].isdigit() and line[0] not in ' -':
                    skipped_lines += 1
                    continue
                
                # Skip lines that appear to be partial/corrupted (contain multiple records or strange patterns)
                if line.count('E+') + line.count('E-') > 4:  # Too many exponents - likely merged lines
                    skipped_lines += 1
                    continue
                    
                try:
                    # Try whitespace-based parsing first (works for well-formatted data)
                    parts = line.split()
                    if len(parts) >= 7 and current_orbital is not None:
                        # Validate that first 3 parts look like integers
                        kz = int(parts[0])
                        ig = int(parts[1])
                        m_idx = int(parts[2])
                        
                        # Validate reasonable ranges
                        if kz < 1 or kz > 10000 or ig < 1 or ig > 100000:
                            skipped_lines += 1
                            continue
                        
                        zsl = float(parts[3])
                        gn = float(parts[4])
                        re_w0 = float(parts[5])
                        im_w0 = float(parts[6])
                        
                        # Validate m_idx is within expected range for this orbital
                        max_m = M_COUNT.get(current_orbital, 7)
                        if m_idx < 1 or m_idx > max_m:
                            skipped_lines += 1
                            continue  # Skip invalid entries
                        
                        # Skip entries with unreasonably large values (likely parsing errors)
                        if abs(re_w0) > 1e10 or abs(im_w0) > 1e10:
                            skipped_lines += 1
                            continue
                        
                        data.append({
                            'lb': current_orbital,
                            'z0': current_z0,
                            'dz': current_dz,
                            'kz': kz,
                            'ig': ig,
                            'm_idx': m_idx,  # 1-based index
                            'zsl': zsl,
                            'gn': gn,
                            're_w0': re_w0,
                            'im_w0': im_w0,
                            'w0': complex(re_w0, im_w0)
                        })
                except (ValueError, IndexError) as e:
                    skipped_lines += 1
                    continue
    
    if skipped_lines > 0:
        print(f"Note: Skipped {skipped_lines} corrupted/invalid lines during parsing")
    
    return data

def get_m_name(lb, m_idx):
    """Get the name of the m-component for a given orbital and index."""
    if lb in M_VALUES and m_idx in M_VALUES[lb]:
        return M_VALUES[lb][m_idx][1]
    return f'm_idx={m_idx}'

def analyze_signs(data):
    """Analyze sign patterns of w0 values."""
    print("\n" + "="*60)
    print("SIGN ANALYSIS")
    print("="*60)
    
    # Group by orbital type and m_idx
    by_orbital_m = defaultdict(list)
    for d in data:
        key = (d['lb'], d['m_idx'])
        by_orbital_m[key].append(d)
    
    for (lb, m_idx), entries in sorted(by_orbital_m.items()):
        orbital_name = ORBITAL_NAMES.get(lb, f'l={lb}')
        m_name = get_m_name(lb, m_idx)
        
        re_vals = [e['re_w0'] for e in entries]
        im_vals = [e['im_w0'] for e in entries]
        
        # Count positive/negative/zero
        re_pos = sum(1 for v in re_vals if v > 1e-12)
        re_neg = sum(1 for v in re_vals if v < -1e-12)
        re_zero = len(re_vals) - re_pos - re_neg
        
        im_pos = sum(1 for v in im_vals if v > 1e-12)
        im_neg = sum(1 for v in im_vals if v < -1e-12)
        im_zero = len(im_vals) - im_pos - im_neg
        
        print(f"\n{orbital_name}-orbital, {m_name} (lb={lb}, m_idx={m_idx}):")
        print(f"  Total entries: {len(entries)}")
        print(f"  Re(w0): {re_pos} positive, {re_neg} negative, {re_zero} ~zero")
        print(f"  Im(w0): {im_pos} positive, {im_neg} negative, {im_zero} ~zero")
        
        if re_vals:
            print(f"  Re(w0) range: [{min(re_vals):.6e}, {max(re_vals):.6e}]")
        if im_vals:
            print(f"  Im(w0) range: [{min(im_vals):.6e}, {max(im_vals):.6e}]")
        
        # Check for sign consistency based on z
        z_positive = [e for e in entries if e['zsl'] > 1e-10]
        z_negative = [e for e in entries if e['zsl'] < -1e-10]
        
        if z_positive and z_negative:
            re_sign_zpos = np.sign([e['re_w0'] for e in z_positive if abs(e['re_w0']) > 1e-12])
            re_sign_zneg = np.sign([e['re_w0'] for e in z_negative if abs(e['re_w0']) > 1e-12])
            
            if len(re_sign_zpos) > 0 and len(re_sign_zneg) > 0:
                if np.all(re_sign_zpos == re_sign_zpos[0]) and np.all(re_sign_zneg == re_sign_zneg[0]):
                    if re_sign_zpos[0] == re_sign_zneg[0]:
                        print(f"  Sign pattern: SAME sign for z>0 and z<0 (EVEN in z)")
                    else:
                        print(f"  Sign pattern: OPPOSITE signs for z>0 and z<0 (ODD in z)")

def analyze_z_dependence(data):
    """Analyze how w0 depends on z position."""
    print("\n" + "="*60)
    print("Z-DEPENDENCE ANALYSIS")
    print("="*60)
    
    # Group by orbital type, m_idx, and g
    by_orbital_m_g = defaultdict(list)
    for d in data:
        key = (d['lb'], d['m_idx'], d['ig'])
        by_orbital_m_g[key].append(d)
    
    # For each group, analyze z-dependence
    for (lb, m_idx, ig), entries in sorted(by_orbital_m_g.items()):
        if len(entries) < 2:
            continue
            
        orbital_name = ORBITAL_NAMES.get(lb, f'l={lb}')
        m_name = get_m_name(lb, m_idx)
        
        # Sort by z
        entries_sorted = sorted(entries, key=lambda e: e['zsl'])
        z_vals = [e['zsl'] for e in entries_sorted]
        w0_vals = [e['w0'] for e in entries_sorted]
        
        # Only print first g-vector for each orbital/m combination
        if ig == 1:
            print(f"\n{orbital_name}-orbital, {m_name} (lb={lb}, m_idx={m_idx}), ig={ig}:")
            print(f"  gn = {entries[0]['gn']:.6e}")
            print(f"  z range: [{min(z_vals):.6f}, {max(z_vals):.6f}]")
            
            # Check symmetry properties
            # For even functions: w0(z) = w0(-z)
            # For odd functions: w0(z) = -w0(-z)
            
            # Find pairs of z and -z values
            tolerance = 1e-6
            symmetric_pairs = []
            for i, e1 in enumerate(entries_sorted):
                for e2 in entries_sorted[i+1:]:
                    if abs(e1['zsl'] + e2['zsl']) < tolerance:
                        symmetric_pairs.append((e1, e2))
            
            if symmetric_pairs:
                even_count = 0
                odd_count = 0
                for e1, e2 in symmetric_pairs:
                    mag = max(abs(e1['w0']), abs(e2['w0']), 1e-12)
                    # Check if even: w0(z) ≈ w0(-z)
                    if abs(e1['w0'] - e2['w0']) < 1e-6 * mag:
                        even_count += 1
                    # Check if odd: w0(z) ≈ -w0(-z)
                    elif abs(e1['w0'] + e2['w0']) < 1e-6 * mag:
                        odd_count += 1
                
                if even_count > odd_count and even_count > 0:
                    print(f"  Symmetry: EVEN in z (w0(z) = w0(-z))")
                elif odd_count > even_count and odd_count > 0:
                    print(f"  Symmetry: ODD in z (w0(z) = -w0(-z))")
                elif even_count > 0 or odd_count > 0:
                    print(f"  Symmetry: MIXED (even:{even_count}, odd:{odd_count})")
            
            # Show a few sample values
            print(f"  Sample values:")
            for e in entries_sorted[:3]:
                print(f"    z={e['zsl']:+.6f}: w0 = {e['re_w0']:+.6e} + {e['im_w0']:+.6e}i")
            if len(entries_sorted) > 6:
                print(f"    ...")
            for e in entries_sorted[-3:]:
                print(f"    z={e['zsl']:+.6f}: w0 = {e['re_w0']:+.6e} + {e['im_w0']:+.6e}i")

def analyze_g_dependence(data):
    """Analyze how w0 depends on g-vector magnitude."""
    print("\n" + "="*60)
    print("G-VECTOR DEPENDENCE ANALYSIS")
    print("="*60)
    
    # Group by orbital type and m_idx, fixed z (middle z point)
    by_orbital_m = defaultdict(list)
    for d in data:
        key = (d['lb'], d['m_idx'])
        by_orbital_m[key].append(d)
    
    for (lb, m_idx), entries in sorted(by_orbital_m.items()):
        orbital_name = ORBITAL_NAMES.get(lb, f'l={lb}')
        m_name = get_m_name(lb, m_idx)
        
        # Get unique z values
        z_vals = sorted(set(e['zsl'] for e in entries))
        if not z_vals:
            continue
            
        # Pick middle z value
        mid_z = z_vals[len(z_vals)//2]
        
        # Filter entries at this z
        at_mid_z = [e for e in entries if abs(e['zsl'] - mid_z) < 1e-8]
        
        # Sort by g
        at_mid_z_sorted = sorted(at_mid_z, key=lambda e: e['gn'])
        
        if at_mid_z_sorted:
            print(f"\n{orbital_name}-orbital, {m_name} at z={mid_z:.6f}:")
            print(f"  |w0| vs gn:")
            
            prev_gn = -1
            for e in at_mid_z_sorted[:10]:  # First 10 g-values
                if abs(e['gn'] - prev_gn) < 1e-8:
                    continue
                prev_gn = e['gn']
                w0_mag = abs(e['w0'])
                print(f"    gn={e['gn']:.6e}: |w0|={w0_mag:.6e}")

def check_expected_signs(data):
    """Check if signs match theoretical expectations."""
    print("\n" + "="*60)
    print("SIGN CORRECTNESS CHECK (THEORETICAL)")
    print("="*60)
    
    print("""
Expected sign patterns based on spherical harmonics:

For s-orbital (lb=0, m_idx=1):
  - w0 should be REAL and EVEN in z
  
For p-orbital (lb=1):
  - m_idx=1 (pz): REAL, ODD in z (contains factor z)
  - m_idx=2,3 (px,py): IMAGINARY, EVEN in z
  
For d-orbital (lb=2):
  - m_idx=1 (dz²): REAL, EVEN in z (contains z²)
  - m_idx=2,3 (dxz,dyz): IMAGINARY, ODD in z (contains z)
  - m_idx=4,5 (dx²-y²,dxy): REAL, EVEN in z
  
For f-orbital (lb=3):
  - m_idx=1 (fz³): REAL, ODD in z (contains z³, z)
  - m_idx=2,3: IMAGINARY, EVEN in z (contains z²)
  - m_idx=4,5: REAL, ODD in z (contains z)
  - m_idx=6,7: IMAGINARY, EVEN in z
""")
    
    # Group and check
    by_orbital_m = defaultdict(list)
    for d in data:
        key = (d['lb'], d['m_idx'])
        by_orbital_m[key].append(d)
    
    issues = []
    
    for (lb, m_idx), entries in sorted(by_orbital_m.items()):
        orbital_name = ORBITAL_NAMES.get(lb, f'l={lb}')
        m_name = get_m_name(lb, m_idx)
        
        re_vals = [e['re_w0'] for e in entries if abs(e['re_w0']) > 1e-15]
        im_vals = [e['im_w0'] for e in entries if abs(e['im_w0']) > 1e-15]
        
        is_real = len(re_vals) > 0 and len(im_vals) == 0
        is_imag = len(im_vals) > 0 and len(re_vals) == 0
        is_complex = len(re_vals) > 0 and len(im_vals) > 0
        
        # Determine expected behavior based on orbital type
        expected_type = "UNKNOWN"
        if lb == 0:  # s-orbital
            expected_type = "REAL"
            if not is_real and is_complex:
                # Check if imaginary part is negligible
                re_mag = max(abs(v) for v in re_vals) if re_vals else 0
                im_mag = max(abs(v) for v in im_vals) if im_vals else 0
                if im_mag < 1e-6 * re_mag:
                    is_real = True
                    is_complex = False
        elif lb == 1:  # p-orbital
            if m_idx == 1:  # pz
                expected_type = "REAL"
            elif m_idx in [2, 3]:  # px, py
                expected_type = "IMAGINARY"
        elif lb == 2:  # d-orbital
            if m_idx == 1:  # dz²
                expected_type = "REAL"
            elif m_idx in [2, 3]:  # dxz, dyz
                expected_type = "IMAGINARY"
            elif m_idx in [4, 5]:  # dx²-y², dxy
                expected_type = "REAL"
        elif lb == 3:  # f-orbital
            if m_idx == 1:  # fz³
                expected_type = "REAL"
            elif m_idx in [2, 3]:
                expected_type = "IMAGINARY"
            elif m_idx in [4, 5]:
                expected_type = "REAL"
            elif m_idx in [6, 7]:
                expected_type = "IMAGINARY"
        
        actual = "REAL" if is_real else ("IMAGINARY" if is_imag else "COMPLEX")
        status = "✓" if actual == expected_type else "✗"
        
        print(f"{status} {orbital_name}-orbital {m_name} (m_idx={m_idx}): {actual} (expected: {expected_type})")
        
        if actual != expected_type:
            issues.append(f"{orbital_name}-orbital {m_name}: Expected {expected_type}, got {actual}")
    
    if issues:
        print(f"\n⚠️  {len(issues)} POTENTIAL SIGN ISSUES DETECTED")
    else:
        print("\n✓ All sign patterns match expectations")

def find_anomalies(data, threshold=1e-10):
    """Find specific grid points with anomalous values in the 'wrong' component.
    
    For REAL orbitals, finds entries with large imaginary parts.
    For IMAGINARY orbitals, finds entries with large real parts.
    """
    print("\n" + "="*60)
    print("ANOMALY DETECTION")
    print("="*60)
    print(f"Threshold for 'anomalous': |value| > {threshold:.0e}")
    
    # Expected types based on orbital and m_idx
    expected_types = {
        (0, 1): 'REAL',    # s-orbital m=0
        (1, 1): 'REAL',    # p-orbital pz
        (1, 2): 'IMAGINARY',  # p-orbital px
        (1, 3): 'IMAGINARY',  # p-orbital py
        (2, 1): 'REAL',    # d-orbital dz²
        (2, 2): 'IMAGINARY',  # d-orbital dxz
        (2, 3): 'IMAGINARY',  # d-orbital dyz
        (2, 4): 'REAL',    # d-orbital dx²-y²
        (2, 5): 'REAL',    # d-orbital dxy
        (3, 1): 'REAL',    # f-orbital m=1
        (3, 2): 'IMAGINARY',  # f-orbital m=2
        (3, 3): 'IMAGINARY',  # f-orbital m=3
        (3, 4): 'REAL',    # f-orbital m=4
        (3, 5): 'REAL',    # f-orbital m=5
        (3, 6): 'IMAGINARY',  # f-orbital m=6
        (3, 7): 'IMAGINARY',  # f-orbital m=7
    }
    
    anomalies_by_orbital = defaultdict(list)
    
    for d in data:
        key = (d['lb'], d['m_idx'])
        expected = expected_types.get(key, 'UNKNOWN')
        
        if expected == 'REAL' and abs(d['im_w0']) > threshold:
            anomalies_by_orbital[key].append({
                'type': 'unexpected_imaginary',
                'kz': d['kz'],
                'ig': d['ig'],
                'z': d['zsl'],
                'gn': d['gn'],
                'value': d['im_w0'],
                're_w0': d['re_w0']
            })
        elif expected == 'IMAGINARY' and abs(d['re_w0']) > threshold:
            anomalies_by_orbital[key].append({
                'type': 'unexpected_real',
                'kz': d['kz'],
                'ig': d['ig'],
                'z': d['zsl'],
                'gn': d['gn'],
                'value': d['re_w0'],
                'im_w0': d['im_w0']
            })
    
    total_anomalies = sum(len(v) for v in anomalies_by_orbital.values())
    print(f"\nTotal anomalies found: {total_anomalies}")
    
    if total_anomalies == 0:
        print("✓ No anomalies detected")
        return
    
    for (lb, m_idx), anomalies in sorted(anomalies_by_orbital.items()):
        if not anomalies:
            continue
            
        orbital_name = ORBITAL_NAMES.get(lb, f'l={lb}')
        m_name = get_m_name(lb, m_idx)
        expected = expected_types.get((lb, m_idx), 'UNKNOWN')
        
        print(f"\n{orbital_name}-orbital {m_name} (expected {expected}):")
        print(f"  {len(anomalies)} anomalous entries")
        
        # Find max anomaly
        max_anomaly = max(anomalies, key=lambda x: abs(x['value']))
        print(f"  Max anomalous value: {max_anomaly['value']:.6e}")
        
        # Show details for top 5 anomalies
        sorted_anomalies = sorted(anomalies, key=lambda x: abs(x['value']), reverse=True)[:5]
        print(f"  Top anomalous points:")
        for a in sorted_anomalies:
            if a['type'] == 'unexpected_imaginary':
                print(f"    kz={a['kz']:4d}, ig={a['ig']:5d}, z={a['z']:+.6f}, gn={a['gn']:.6e}: "
                      f"Im={a['value']:+.6e} (Re={a['re_w0']:+.6e})")
            else:
                print(f"    kz={a['kz']:4d}, ig={a['ig']:5d}, z={a['z']:+.6f}, gn={a['gn']:.6e}: "
                      f"Re={a['value']:+.6e} (Im={a['im_w0']:+.6e})")
        
        # Check if anomalies cluster around specific conditions
        gn_zero = [a for a in anomalies if abs(a['gn']) < 1e-8]
        z_zero = [a for a in anomalies if abs(a['z']) < 1e-8]
        
        if gn_zero:
            print(f"  ⚠️  {len(gn_zero)} anomalies at gn≈0 (Γ-point)")
        if z_zero:
            print(f"  ⚠️  {len(z_zero)} anomalies at z≈0")

def main():
    # Default filename
    filename = 'w0_debug.dat'
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    
    print(f"Analyzing: {filename}")
    print("="*60)
    
    try:
        data = parse_w0_debug(filename)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        print("Please run a PWCOND calculation first to generate the debug output.")
        sys.exit(1)
    
    if not data:
        print("No data found in file. Check if the file format is correct.")
        sys.exit(1)
    
    print(f"Loaded {len(data)} data points")
    
    # Get summary of orbital types found
    orbitals_found = set(d['lb'] for d in data)
    print(f"Orbital types found: {[ORBITAL_NAMES.get(lb, f'l={lb}') for lb in sorted(orbitals_found)]}")
    
    # Run analyses
    analyze_signs(data)
    analyze_z_dependence(data)
    analyze_g_dependence(data)
    check_expected_signs(data)
    find_anomalies(data)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)

if __name__ == '__main__':
    main()
