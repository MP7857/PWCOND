#!/usr/bin/env python3
"""
Analyze w0 norms across different nz1 values to validate nz-stability.

Usage: python analyze_w0_norms_vs_nz.py w0_norms_nz7.dat w0_norms_nz11.dat [...]
"""

import sys
import numpy as np
from collections import defaultdict

def load_norms(fname, nz_label):
    """Load norm data from a w0_norms.dat file."""
    data = []
    with open(fname, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or not line.strip():
                continue
            parts = line.split()
            nz1 = int(parts[0])
            lb  = int(parts[1])
            m   = int(parts[2])
            norm = float(parts[3])
            data.append((nz_label, nz1, lb, m, norm))
    return data

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analyze_w0_norms_vs_nz.py w0_norms_nz7.dat w0_norms_nz11.dat [w0_norms_nzXX.dat ...]")
        print()
        print("This script compares w0 norms across different nz1 values to validate")
        print("the nz-stability of the Simpson integration implementation.")
        sys.exit(1)

    all_data = []
    for fname in sys.argv[1:]:
        # Try to extract nz-label from file name (optional)
        # e.g. w0_norms_nz7.dat -> nz_label=7
        label = None
        for tok in fname.replace('.', '_').split('_'):
            if tok.lower().startswith('nz'):
                try:
                    label = int(tok[2:])
                except ValueError:
                    pass
        if label is None:
            label = 0
        all_data.extend(load_norms(fname, label))

    # group by (lb,m)
    by_lm = defaultdict(list)
    for nz_label, nz1, lb, m, norm in all_data:
        by_lm[(lb, m)].append((nz_label, nz1, norm))

    print("Per-(l,m) norm comparison across nz runs")
    print("  lb  m   nz_label   nz1     norm")
    print("  ----------------------------------------")
    for (lb, m), vals in sorted(by_lm.items()):
        for nz_label, nz1, norm in sorted(vals):
            print(f"  {lb:2d} {m:2d}   {nz_label:8d}  {nz1:4d}  {norm:12.5e}")
        
        # Compute relative change if we have multiple nz values
        if len(vals) > 1:
            sorted_vals = sorted(vals)
            ref_norm = sorted_vals[0][2]
            for i, (nz_label, nz1, norm) in enumerate(sorted_vals[1:], 1):
                if ref_norm > 1e-15:
                    rel_change = abs(norm - ref_norm) / ref_norm * 100
                    print(f"       -> relative change from first: {rel_change:6.2f}%")
                elif abs(norm) > 1e-15:
                    print(f"       -> relative change from first: undefined (ref_norm~0, norm={norm:.3e})")
                else:
                    print(f"       -> relative change from first: both norms ~0")
        print()

    print("\nSummary:")
    print("For nz-stable integration, relative changes should be < 5-10%")
    print("Large changes (>10%) indicate nz-dependent behavior in that channel")
