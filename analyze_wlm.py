#!/usr/bin/env python3
"""
WLM_SUMMARY Output Parser and Plotter

This script parses the WLM_SUMMARY output from PWCOND's four.f90 module
and generates analysis plots for tunneling conductance calculations.

It produces three types of plots:
1. κ(E) - Decay constant vs energy (identifies tunneling resonances)
2. ⟨g²⟩(E) - Average transverse momentum vs energy
3. Orbital contributions - Bar chart showing which orbitals dominate tunneling

Usage:
    python3 analyze_wlm.py [output_file]

If no file is specified, it defaults to "out.txt"

Output:
    - kappa_vs_E.png: Tunneling decay constant plot
    - g2_vs_E.png: Transverse momentum plot  
    - orbital_contrib_best_E.png: Orbital decomposition at best tunneling energy
"""

import re
import sys
import matplotlib.pyplot as plt
from collections import defaultdict

def get_float(label, line, default=None):
    """Extract a floating point value from a labeled field in a line."""
    m = re.search(rf'{label}\s*=\s*([+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?)', line)
    return float(m.group(1)) if m else default

def get_int(label, line, default=None):
    """Extract an integer value from a labeled field in a line."""
    m = re.search(rf'{label}\s*=\s*([+-]?\d+)', line)
    return int(m.group(1)) if m else default

def parse_wlm_output(filename):
    """
    Parse WLM_SUMMARY output file.
    
    Handles both single-line and two-line (WLM_SUMMARY + WLM_SUMMARY_CONT) formats.
    
    Returns:
        state_data: List of dicts with MODE=B:STATE information
        state_lm_data: List of dicts with MODE=B:STATE_LM information
    """
    state_data = []      # list of dicts for MODE=B:STATE
    state_lm_data = []   # list of dicts for MODE=B:STATE_LM
    
    with open(filename, "r") as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        if "WLM_SUMMARY" not in line:
            i += 1
            continue
        
        # Parse MODE=B:STATE (two-line or single-line format)
        if "MODE=B:STATE" in line and "STATE_LM" not in line:
            # Line 1: energy, k-points, state index, kappa values
            E = float(line.split()[2])
            k1k2 = re.search(r'k1,k2=\s*([+-]?\d+\.\d*)\s+([+-]?\d+\.\d*)', line)
            k1 = float(k1k2.group(1)) if k1k2 else 0.0
            k2 = float(k1k2.group(2)) if k1k2 else 0.0
            
            n = get_int("n", line)
            kappa_bohr = get_float("kappa_bohr", line)
            kappa_ang = get_float("kappa_ang", line)
            
            # Check if g2 data is on the same line (legacy single-line format)
            g2_bohr = get_float("g2_bohr", line)
            g2_ang = get_float("g2_ang", line)
            norm = get_float("norm", line)
            
            # Or on the next line (two-line continuation format)
            if g2_bohr is None and i + 1 < len(lines) and "WLM_SUMMARY_CONT" in lines[i + 1]:
                cont_line = lines[i + 1]
                g2_bohr = get_float("g2_bohr", cont_line)
                g2_ang = get_float("g2_ang", cont_line)
                norm = get_float("norm", cont_line)
                i += 1  # Skip continuation line in next iteration
            
            state_data.append({
                "E": E, "k1": k1, "k2": k2, "n": n,
                "kappa_bohr": kappa_bohr, "kappa_ang": kappa_ang,
                "g2_bohr": g2_bohr, "g2_ang": g2_ang,
                "norm": norm
            })
        
        # Parse MODE=B:STATE_LM (two-line or single-line format)
        elif "MODE=B:STATE_LM" in line:
            # Line 1: energy, k-points, state index, l, m
            E = float(line.split()[2])
            k1k2 = re.search(r'k1,k2=\s*([+-]?\d+\.\d*)\s+([+-]?\d+\.\d*)', line)
            k1 = float(k1k2.group(1)) if k1k2 else 0.0
            k2 = float(k1k2.group(2)) if k1k2 else 0.0
            
            n = get_int("n", line)
            l = get_int("l", line)
            m = get_int("m", line)
            
            # Check if g2 data is on the same line (legacy single-line format)
            g2_bohr = get_float("g2_bohr", line)
            g2_ang = get_float("g2_ang", line)
            norm_lm = get_float("norm_lm", line)
            
            # Or on the next line (two-line continuation format)
            if g2_bohr is None and i + 1 < len(lines) and "WLM_SUMMARY_CONT" in lines[i + 1]:
                cont_line = lines[i + 1]
                g2_bohr = get_float("g2_bohr", cont_line)
                g2_ang = get_float("g2_ang", cont_line)
                norm_lm = get_float("norm_lm", cont_line)
                i += 1  # Skip continuation line in next iteration
            
            state_lm_data.append({
                "E": E, "k1": k1, "k2": k2, "n": n,
                "l": l, "m": m,
                "g2_bohr": g2_bohr, "g2_ang": g2_ang,
                "norm_lm": norm_lm
            })
        
        i += 1
    
    return state_data, state_lm_data

def plot_kappa_vs_energy(state_data, output_file="kappa_vs_E.png"):
    """Plot decay constant κ vs energy."""
    if not state_data:
        print("Warning: No state data available for κ(E) plot")
        return
    
    state_data = sorted(state_data, key=lambda d: d["E"])
    
    E_list = [d["E"] for d in state_data if d["kappa_ang"] is not None]
    kappa_list = [d["kappa_ang"] for d in state_data if d["kappa_ang"] is not None]
    
    if not E_list:
        print("Warning: No valid kappa data to plot")
        return
    
    plt.figure(figsize=(8, 6))
    plt.plot(E_list, kappa_list, marker="o", linewidth=2, markersize=6)
    plt.xlabel("Energy (eV)", fontsize=12)
    plt.ylabel(r"$\kappa$ (Å$^{-1}$)", fontsize=12)
    plt.title("Slowest-decaying CBS state: κ(E)", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()

def plot_g2_vs_energy(state_data, output_file="g2_vs_E.png"):
    """Plot average transverse momentum squared ⟨g²⟩ vs energy."""
    if not state_data:
        print("Warning: No state data available for ⟨g²⟩(E) plot")
        return
    
    state_data = sorted(state_data, key=lambda d: d["E"])
    
    E_list = [d["E"] for d in state_data if d["g2_ang"] is not None]
    g2_list = [d["g2_ang"] for d in state_data if d["g2_ang"] is not None]
    
    if not E_list:
        print("Warning: No valid g2 data to plot")
        return
    
    plt.figure(figsize=(8, 6))
    plt.plot(E_list, g2_list, marker="o", linewidth=2, markersize=6, color='green')
    plt.xlabel("Energy (eV)", fontsize=12)
    plt.ylabel(r"$\langle g^2 \rangle$ (Å$^{-2}$)", fontsize=12)
    plt.title("Slowest-decaying CBS state: ⟨g²⟩(E)", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()

def plot_orbital_contributions(state_data, state_lm_data, output_file="orbital_contrib_best_E.png"):
    """
    Plot orbital contributions at the best tunneling energy.
    
    Best energy = where κ is minimum (slowest decay = strongest tunneling).
    """
    if not state_data:
        print("Warning: No state data available for orbital contribution plot")
        return
    
    if not state_lm_data:
        print("Warning: No STATE_LM data available for orbital contribution plot")
        return
    
    # Find energy where kappa is minimum
    valid_states = [d for d in state_data if d["kappa_ang"] is not None]
    if not valid_states:
        print("Warning: No valid kappa data to find best tunneling energy")
        return
    
    best_state = min(valid_states, key=lambda d: d["kappa_ang"])
    E_best = best_state["E"]
    n_best = best_state["n"]
    kappa_best = best_state["kappa_ang"]
    
    print(f"\nBest tunneling state:")
    print(f"  Energy: {E_best:.3f} eV")
    print(f"  State n: {n_best}")
    print(f"  κ: {kappa_best:.4f} Å⁻¹")
    
    # Filter LM data for that energy and state (tolerance for floating-point)
    tol_E = 1e-6
    lm_best = [d for d in state_lm_data 
               if abs(d["E"] - E_best) < tol_E and d["n"] == n_best and d["norm_lm"] is not None]
    
    if not lm_best:
        print(f"Warning: No STATE_LM data found for E={E_best:.3f} eV, n={n_best}")
        return
    
    # Group by (l,m)
    by_lm = defaultdict(float)
    for d in lm_best:
        by_lm[(d["l"], d["m"])] += d["norm_lm"]
    
    # Normalize
    total = sum(by_lm.values())
    lm_labels = []
    lm_weights = []
    for (l, m), wt in sorted(by_lm.items()):
        lm_labels.append(f"(l={l},m={m:+d})")
        lm_weights.append(wt / total if total > 0 else 0.0)
    
    if not lm_labels:
        print("Warning: No orbital data to plot")
        return
    
    # Create bar chart
    plt.figure(figsize=(10, 6))
    colors = plt.cm.viridis([i/len(lm_labels) for i in range(len(lm_labels))])
    plt.bar(range(len(lm_labels)), lm_weights, color=colors)
    plt.xticks(range(len(lm_labels)), lm_labels, rotation=45, ha="right")
    plt.ylabel("Normalized weight", fontsize=12)
    plt.title(f"Orbital contributions at E={E_best:.3f} eV (state n={n_best}, κ={kappa_best:.4f} Å⁻¹)", 
              fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()
    
    # Print summary
    print("\nOrbital contributions (normalized):")
    for i, (label, weight) in enumerate(zip(lm_labels, lm_weights)):
        print(f"  {label}: {weight:.4f} ({weight*100:.1f}%)")

def print_help():
    """Print help message."""
    help_text = """
WLM_SUMMARY Output Parser and Plotter

Usage:
    python3 analyze_wlm.py [options] [output_file]

Arguments:
    output_file    Path to PWCOND output file (default: out.txt)

Options:
    -h, --help     Show this help message and exit
    -v, --version  Show version information

Examples:
    python3 analyze_wlm.py                    # Use default out.txt
    python3 analyze_wlm.py my_output.txt      # Parse my_output.txt
    python3 analyze_wlm.py --help             # Show this help

Output:
    - kappa_vs_E.png: Tunneling decay constant plot
    - g2_vs_E.png: Transverse momentum plot
    - orbital_contrib_best_E.png: Orbital decomposition

For more information, see README_WLM_ANALYSIS.md
"""
    print(help_text)

def main():
    """Main analysis routine."""
    # Parse command line arguments
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg in ['-h', '--help']:
            print_help()
            return 0
        elif arg in ['-v', '--version']:
            print("analyze_wlm.py version 1.0")
            print("Part of PWCOND quantum transport package")
            return 0
        else:
            filename = arg
    else:
        filename = "out.txt"
    
    print(f"Parsing WLM_SUMMARY output from: {filename}")
    print("=" * 70)
    
    try:
        state_data, state_lm_data = parse_wlm_output(filename)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        print(f"Usage: python3 {sys.argv[0]} [output_file]")
        return 1
    except Exception as e:
        print(f"Error parsing file: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print(f"Found {len(state_data)} MODE=B:STATE entries")
    print(f"Found {len(state_lm_data)} MODE=B:STATE_LM entries")
    print("=" * 70)
    
    if not state_data and not state_lm_data:
        print("\nWarning: No WLM_SUMMARY data found in file")
        print("Make sure the file contains WLM_SUMMARY output from PWCOND")
        return 1
    
    # Generate plots
    print("\nGenerating plots...")
    print("-" * 70)
    
    if state_data:
        plot_kappa_vs_energy(state_data)
        plot_g2_vs_energy(state_data)
    
    if state_data and state_lm_data:
        plot_orbital_contributions(state_data, state_lm_data)
    
    print("-" * 70)
    print("\nAnalysis complete!")
    print("\nGenerated files:")
    print("  - kappa_vs_E.png: Tunneling decay constant vs energy")
    print("  - g2_vs_E.png: Average transverse momentum vs energy")
    print("  - orbital_contrib_best_E.png: Orbital decomposition at best tunneling energy")
    print("\nInterpretation guide:")
    print("  • Deep minima in κ(E) → dominant tunneling channels")
    print("  • Low ⟨g²⟩(E) → more normal-incidence-like, better coupling")
    print("  • Orbital bars → which atomic orbitals carry the tunneling current")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
