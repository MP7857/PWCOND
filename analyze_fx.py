#!/usr/bin/env python3
"""
Utility script to analyze fx_full.dat diagnostic output from PWCOND.

This script can:
1. Verify the format of fx_full.dat files
2. Compare two fx_full.dat files (e.g., from nz=7 and nz=11)
3. Generate basic statistics

Usage:
    python analyze_fx.py verify fx_full.dat
    python analyze_fx.py compare fx_full_nz7.dat fx_full_nz11.dat
    python analyze_fx.py stats fx_full.dat

Note: This script uses only Python standard library (no numpy required)
"""

import sys
from collections import defaultdict


def read_fx_full(filename):
    """
    Read fx_full.dat file and return structured data.
    
    Returns:
        list: List of tuples (lb, kz, ig, ign, fx1, fx2, fx3, fx4, fx5, fx6, fx7)
    """
    try:
        records = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 11:
                    continue
                try:
                    record = (
                        int(parts[0]),    # lb
                        int(parts[1]),    # kz
                        int(parts[2]),    # ig
                        int(parts[3]),    # ign
                        float(parts[4]),  # fx1
                        float(parts[5]),  # fx2
                        float(parts[6]),  # fx3
                        float(parts[7]),  # fx4
                        float(parts[8]),  # fx5
                        float(parts[9]),  # fx6
                        float(parts[10])  # fx7
                    )
                    records.append(record)
                except ValueError as e:
                    print(f"Warning: Skipping malformed line: {line}")
                    continue
        
        if not records:
            print(f"Warning: {filename} contains no valid data")
            return None
        
        return records
    except FileNotFoundError:
        print(f"Error: File {filename} not found")
        return None
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None


def verify_file(filename):
    """Verify the format and content of fx_full.dat file."""
    print(f"\nVerifying {filename}...")
    
    records = read_fx_full(filename)
    if records is None:
        return False
    
    n_records = len(records)
    print(f"  Total records: {n_records}")
    
    # Collect statistics
    lb_values = set()
    kz_by_lb = defaultdict(set)
    ig_values = set()
    ign_values = set()
    
    for rec in records:
        lb, kz, ig, ign = rec[0], rec[1], rec[2], rec[3]
        lb_values.add(lb)
        kz_by_lb[lb].add(kz)
        ig_values.add(ig)
        ign_values.add(ign)
    
    unique_lb = sorted(lb_values)
    print(f"  Orbital types (lb): {unique_lb}")
    
    # Check z-slices per lb
    for lb in unique_lb:
        kz_set = kz_by_lb[lb]
        print(f"    lb={lb}: {len(kz_set)} z-slices (kz: {min(kz_set)}-{max(kz_set)})")
    
    # Check ig and ign ranges
    print(f"  ig range: {min(ig_values)}-{max(ig_values)}")
    print(f"  ign range: {min(ign_values)}-{max(ign_values)}")
    
    # Check for expected patterns
    for lb in unique_lb:
        lb_records = [r for r in records if r[0] == lb]
        if lb == 0:
            # s-orbital: only fx1 should be non-zero
            nz_fx1 = sum(1 for r in lb_records if abs(r[4]) > 1e-15)
            nz_fx2 = sum(1 for r in lb_records if abs(r[5]) > 1e-15)
            print(f"  lb=0: fx1 non-zero: {nz_fx1}, fx2 non-zero: {nz_fx2} (expect 0)")
        elif lb == 1:
            # p-orbital: fx1, fx2 should be non-zero
            nz_fx1 = sum(1 for r in lb_records if abs(r[4]) > 1e-15)
            nz_fx2 = sum(1 for r in lb_records if abs(r[5]) > 1e-15)
            nz_fx3 = sum(1 for r in lb_records if abs(r[6]) > 1e-15)
            print(f"  lb=1: fx1,fx2 non-zero: {nz_fx1},{nz_fx2}, fx3 non-zero: {nz_fx3} (expect 0)")
        elif lb == 2:
            # d-orbital: fx1-fx4 should be non-zero
            nz_fx4 = sum(1 for r in lb_records if abs(r[7]) > 1e-15)
            nz_fx5 = sum(1 for r in lb_records if abs(r[8]) > 1e-15)
            print(f"  lb=2: fx4 non-zero: {nz_fx4}, fx5 non-zero: {nz_fx5} (expect 0)")
        elif lb == 3:
            # f-orbital: fx1-fx6 should be non-zero
            nz_fx6 = sum(1 for r in lb_records if abs(r[9]) > 1e-15)
            nz_fx7 = sum(1 for r in lb_records if abs(r[10]) > 1e-15)
            print(f"  lb=3: fx6 non-zero: {nz_fx6}, fx7 non-zero: {nz_fx7} (expect 0)")
    
    print(f"\n{filename} verification complete.")
    return True


def compare_files(file1, file2):
    """Compare two fx_full.dat files."""
    print(f"\nComparing {file1} and {file2}...")
    
    records1 = read_fx_full(file1)
    records2 = read_fx_full(file2)
    
    if records1 is None or records2 is None:
        return False
    
    print(f"  {file1}: {len(records1)} records")
    print(f"  {file2}: {len(records2)} records")
    
    # Find common lb values
    lb1 = set(r[0] for r in records1)
    lb2 = set(r[0] for r in records2)
    common_lb = sorted(lb1 & lb2)
    
    print(f"  Common orbital types: {common_lb}")
    
    for lb in common_lb:
        recs1 = [r for r in records1 if r[0] == lb]
        recs2 = [r for r in records2 if r[0] == lb]
        
        kz1 = set(r[1] for r in recs1)
        kz2 = set(r[1] for r in recs2)
        
        print(f"\n  lb={lb}:")
        print(f"    {file1}: {len(kz1)} z-slices")
        print(f"    {file2}: {len(kz2)} z-slices")
        
        # For meaningful comparison, you would need to interpolate or
        # find corresponding slices. This is left for more advanced analysis.
    
    print(f"\nComparison complete. Note: Detailed difference analysis requires")
    print(f"interpolation between different nz grids.")
    return True


def print_stats(filename):
    """Print statistics for fx_full.dat file."""
    print(f"\nStatistics for {filename}...")
    
    records = read_fx_full(filename)
    if records is None:
        return False
    
    lb_values = sorted(set(r[0] for r in records))
    
    for lb in lb_values:
        lb_records = [r for r in records if r[0] == lb]
        print(f"\nlb={lb}:")
        
        for i in range(4, 11):  # fx1 (index 4) through fx7 (index 10)
            fx_num = i - 3
            fx_values = [r[i] for r in lb_records]
            non_zero = [v for v in fx_values if abs(v) > 1e-15]
            
            if non_zero:
                mean = sum(non_zero) / len(non_zero)
                variance = sum((v - mean)**2 for v in non_zero) / len(non_zero)
                std = variance ** 0.5
                
                print(f"  fx{fx_num}: min={min(non_zero):.6e}, max={max(non_zero):.6e}, "
                      f"mean={mean:.6e}, std={std:.6e} ({len(non_zero)} non-zero)")
    
    return True


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    
    command = sys.argv[1]
    
    if command == 'verify' and len(sys.argv) == 3:
        verify_file(sys.argv[2])
    elif command == 'compare' and len(sys.argv) == 4:
        compare_files(sys.argv[2], sys.argv[3])
    elif command == 'stats' and len(sys.argv) == 3:
        print_stats(sys.argv[2])
    else:
        print(__doc__)
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
