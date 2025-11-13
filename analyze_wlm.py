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
import csv

def get_float(label, line, default=None):
    """Extract a floating point value from a labeled field in a line."""
    m = re.search(rf'{label}\s*=\s*([+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?)', line)
    return float(m.group(1)) if m else default

def get_int(label, line, default=None):
    """Extract an integer value from a labeled field in a line."""
    m = re.search(rf'{label}\s*=\s*([+-]?\d+)', line)
    return int(m.group(1)) if m else default

def orbital_name(l, m):
    """
    Map (l, m) quantum numbers to physical orbital names.
    
    Based on the ordering in four.f90:
    - l=0: s
    - l=1: p_z (m=0), p_x (m=-1), p_y (m=+1)
    - l=2: d_z²−1 (m=0), d_xz (m=-1), d_yz (m=+1), d_x²−y² (m=+2), d_xy (m=-2)
    - l=3: f_z(5z²−3r²) (m=0), f_x(5z²−r²) (m=-1), f_y(5z²−r²) (m=+1),
            f_z(x²−y²) (m=-2), f_xyz (m=+2), f_x(x²−3y²) (m=-3), f_y(3x²−y²) (m=+3)
    
    Args:
        l: Angular momentum quantum number
        m: Magnetic quantum number
    
    Returns:
        Human-readable orbital name string
    """
    if l == 0:
        return "s"
    elif l == 1:
        if m == 0:
            return "p_z"
        elif m == -1:
            return "p_x"
        elif m == 1:
            return "p_y"
    elif l == 2:
        if m == 0:
            return "d_z²−1"
        elif m == -1:
            return "d_xz"
        elif m == 1:
            return "d_yz"
        elif m == 2:
            return "d_x²−y²"
        elif m == -2:
            return "d_xy"
    elif l == 3:
        if m == 0:
            return "f_z(5z²−3r²)"
        elif m == -1:
            return "f_x(5z²−r²)"
        elif m == 1:
            return "f_y(5z²−r²)"
        elif m == -2:
            return "f_z(x²−y²)"
        elif m == 2:
            return "f_xyz"
        elif m == -3:
            return "f_x(x²−3y²)"
        elif m == 3:
            return "f_y(3x²−y²)"
    
    # Fallback for undefined cases
    return f"(l={l},m={m:+d})"

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
        lm_labels.append(orbital_name(l, m))
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

def generate_tables(state_data, state_lm_data, output_prefix="wlm_tables"):
    """
    Generate comprehensive tables for tunneling state analysis.
    
    Creates:
    - Table 1: Top tunneling states per energy (CSV format)
    - Table 2: Orbital character decomposition per state (CSV format)
    - Text summary for inclusion in papers
    
    Args:
        state_data: List of dicts with MODE=B:STATE information
        state_lm_data: List of dicts with MODE=B:STATE_LM information
        output_prefix: Prefix for output filenames
    """
    if not state_data or not state_lm_data:
        print("Warning: Insufficient data for table generation")
        return
    
    # Group data by energy
    states_by_energy = defaultdict(list)
    lm_by_energy_state = defaultdict(lambda: defaultdict(list))
    
    for state in state_data:
        if state["kappa_ang"] is not None and state["n"] is not None:
            E = state["E"]
            states_by_energy[E].append(state)
    
    for lm_state in state_lm_data:
        if lm_state["norm_lm"] is not None:
            E = lm_state["E"]
            n = lm_state["n"]
            lm_by_energy_state[E][n].append(lm_state)
    
    # Table 1: Top tunneling states per energy
    table1_filename = f"{output_prefix}_top_states.csv"
    print(f"\nGenerating Table 1: {table1_filename}")
    
    with open(table1_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Energy_eV', 'State_n', 'kappa_Bohr', 'kappa_Ang', 
                        'g2_Bohr2', 'g2_Ang2', 'norm', 'dominant_orbitals'])
        
        for E in sorted(states_by_energy.keys()):
            # Sort states by kappa (ascending) and take top 3
            states = sorted(states_by_energy[E], key=lambda s: s["kappa_ang"])[:3]
            
            for state in states:
                n = state["n"]
                
                # Get dominant orbitals for this state
                lm_states = lm_by_energy_state[E][n]
                if lm_states:
                    # Compute total weight and fractional contributions
                    total_weight = sum(lm["norm_lm"] for lm in lm_states)
                    
                    if total_weight > 0:
                        # Sort by weight descending, take top 3
                        lm_sorted = sorted(lm_states, 
                                         key=lambda x: x["norm_lm"], 
                                         reverse=True)[:3]
                        
                        orbital_strs = []
                        for lm in lm_sorted:
                            frac = lm["norm_lm"] / total_weight
                            if frac > 0.05:  # Only show if > 5%
                                name = orbital_name(lm["l"], lm["m"])
                                orbital_strs.append(f"{name}({frac*100:.0f}%)")
                        
                        dominant = ", ".join(orbital_strs) if orbital_strs else "mixed"
                    else:
                        dominant = "N/A"
                else:
                    dominant = "N/A"
                
                writer.writerow([
                    f"{E:.3f}",
                    n,
                    f"{state['kappa_bohr']:.4f}",
                    f"{state['kappa_ang']:.4f}",
                    f"{state['g2_bohr']:.2e}",
                    f"{state['g2_ang']:.2e}",
                    f"{state['norm']:.2e}",
                    dominant
                ])
    
    print(f"Saved: {table1_filename}")
    
    # Table 2: Detailed orbital character per state
    table2_filename = f"{output_prefix}_orbital_character.csv"
    print(f"\nGenerating Table 2: {table2_filename}")
    
    with open(table2_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Energy_eV', 'State_n', 'kappa_Ang', 'l', 'm', 
                        'orbital_name', 'norm_lm', 'fractional_weight', 'g2_Ang2'])
        
        for E in sorted(states_by_energy.keys()):
            states = sorted(states_by_energy[E], key=lambda s: s["kappa_ang"])[:3]
            
            for state in states:
                n = state["n"]
                kappa_ang = state["kappa_ang"]
                
                lm_states = lm_by_energy_state[E][n]
                if not lm_states:
                    continue
                
                # Compute total weight
                total_weight = sum(lm["norm_lm"] for lm in lm_states)
                
                # Sort by contribution
                lm_sorted = sorted(lm_states, 
                                 key=lambda x: x["norm_lm"], 
                                 reverse=True)
                
                for lm in lm_sorted:
                    frac = lm["norm_lm"] / total_weight if total_weight > 0 else 0
                    if frac > 0.01:  # Only include if > 1%
                        name = orbital_name(lm["l"], lm["m"])
                        writer.writerow([
                            f"{E:.3f}",
                            n,
                            f"{kappa_ang:.4f}",
                            lm["l"],
                            lm["m"],
                            name,
                            f"{lm['norm_lm']:.3e}",
                            f"{frac:.4f}",
                            f"{lm['g2_ang']:.2e}"
                        ])
    
    print(f"Saved: {table2_filename}")
    
    # Generate text summary
    summary_filename = f"{output_prefix}_summary.txt"
    print(f"\nGenerating text summary: {summary_filename}")
    
    with open(summary_filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("TUNNELING STATE ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("This summary identifies the dominant tunneling channels and their\n")
        f.write("orbital character at each energy. States are ranked by decay constant\n")
        f.write("κ (smaller κ = slower decay = better tunneling).\n\n")
        
        for E in sorted(states_by_energy.keys()):
            f.write("-" * 80 + "\n")
            f.write(f"Energy: {E:.3f} eV\n")
            f.write("-" * 80 + "\n")
            
            states = sorted(states_by_energy[E], key=lambda s: s["kappa_ang"])[:3]
            
            for rank, state in enumerate(states, 1):
                n = state["n"]
                kappa = state["kappa_ang"]
                g2 = state["g2_ang"]
                
                f.write(f"\nRank {rank}: State n={n}\n")
                f.write(f"  κ = {kappa:.4f} Å⁻¹\n")
                f.write(f"  ⟨g²⟩ = {g2:.2f} Å⁻²\n")
                
                lm_states = lm_by_energy_state[E][n]
                if lm_states:
                    total_weight = sum(lm["norm_lm"] for lm in lm_states)
                    
                    if total_weight > 0:
                        # Group by l and compute total per l
                        by_l = defaultdict(float)
                        for lm in lm_states:
                            by_l[lm["l"]] += lm["norm_lm"] / total_weight
                        
                        f.write("  Orbital character:\n")
                        l_names = {0: "s", 1: "p", 2: "d", 3: "f"}
                        for l in sorted(by_l.keys()):
                            frac = by_l[l]
                            f.write(f"    {l_names.get(l, f'l={l}')}: {frac*100:.1f}%\n")
                        
                        # Show top 3 individual (l,m) contributions
                        lm_sorted = sorted(lm_states, 
                                         key=lambda x: x["norm_lm"], 
                                         reverse=True)[:3]
                        
                        f.write("  Dominant (l,m) channels:\n")
                        for lm in lm_sorted:
                            frac = lm["norm_lm"] / total_weight
                            if frac > 0.05:
                                name = orbital_name(lm["l"], lm["m"])
                                f.write(f"    {name}: {frac*100:.1f}%\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("END OF SUMMARY\n")
        f.write("=" * 80 + "\n")
    
    print(f"Saved: {summary_filename}")
    print("\nTable generation complete!")
    print(f"\nFiles created:")
    print(f"  • {table1_filename} - Top tunneling states per energy")
    print(f"  • {table2_filename} - Detailed orbital character breakdown")
    print(f"  • {summary_filename} - Human-readable summary for papers")

def plot_orbital_evolution(state_data, state_lm_data, output_file="orbital_evolution_vs_E.png"):
    """
    Plot how the dominant orbital character evolves with energy.
    
    For each energy, takes the slowest-decaying state and plots the
    fractional contribution of s, p, d, f orbitals vs energy.
    """
    if not state_data or not state_lm_data:
        print("Warning: Insufficient data for orbital evolution plot")
        return
    
    # Group data by energy
    states_by_energy = defaultdict(list)
    lm_by_energy_state = defaultdict(lambda: defaultdict(list))
    
    for state in state_data:
        if state["kappa_ang"] is not None:
            states_by_energy[state["E"]].append(state)
    
    for lm_state in state_lm_data:
        if lm_state["norm_lm"] is not None:
            lm_by_energy_state[lm_state["E"]][lm_state["n"]].append(lm_state)
    
    # For each energy, get the slowest-decaying state and its orbital character
    energies = []
    s_frac = []
    p_frac = []
    d_frac = []
    f_frac = []
    
    for E in sorted(states_by_energy.keys()):
        states = sorted(states_by_energy[E], key=lambda s: s["kappa_ang"])
        if not states:
            continue
        
        best_state = states[0]
        n = best_state["n"]
        
        lm_states = lm_by_energy_state[E][n]
        if not lm_states:
            continue
        
        # Compute total weight and fractions per l
        total_weight = sum(lm["norm_lm"] for lm in lm_states)
        if total_weight <= 0:
            continue
        
        by_l = defaultdict(float)
        for lm in lm_states:
            by_l[lm["l"]] += lm["norm_lm"] / total_weight
        
        energies.append(E)
        s_frac.append(by_l.get(0, 0.0))
        p_frac.append(by_l.get(1, 0.0))
        d_frac.append(by_l.get(2, 0.0))
        f_frac.append(by_l.get(3, 0.0))
    
    if not energies:
        print("Warning: No valid data for orbital evolution plot")
        return
    
    # Create stacked area plot
    plt.figure(figsize=(10, 6))
    plt.stackplot(energies, s_frac, p_frac, d_frac, f_frac,
                  labels=['s', 'p', 'd', 'f'],
                  alpha=0.8,
                  colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
    
    plt.xlabel("Energy (eV)", fontsize=12)
    plt.ylabel("Fractional orbital character", fontsize=12)
    plt.title("Orbital character evolution of slowest-decaying state vs Energy", fontsize=13)
    plt.legend(loc='upper right', fontsize=11)
    plt.grid(True, alpha=0.3, axis='y')
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.close()

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
        plot_orbital_evolution(state_data, state_lm_data)
    
    # Generate comprehensive tables
    print("\n" + "=" * 70)
    print("Generating comprehensive analysis tables...")
    print("=" * 70)
    
    if state_data and state_lm_data:
        generate_tables(state_data, state_lm_data)
    
    print("-" * 70)
    print("\nAnalysis complete!")
    print("\nGenerated files:")
    print("  Plots:")
    print("    - kappa_vs_E.png: Tunneling decay constant vs energy")
    print("    - g2_vs_E.png: Average transverse momentum vs energy")
    print("    - orbital_contrib_best_E.png: Orbital decomposition at best tunneling energy")
    print("    - orbital_evolution_vs_E.png: Orbital character evolution with energy")
    print("  Tables:")
    print("    - wlm_tables_top_states.csv: Top tunneling states per energy")
    print("    - wlm_tables_orbital_character.csv: Detailed orbital character breakdown")
    print("    - wlm_tables_summary.txt: Human-readable summary for papers")
    print("\nInterpretation guide:")
    print("  • Deep minima in κ(E) → dominant tunneling channels")
    print("  • Low ⟨g²⟩(E) → more normal-incidence-like, better coupling")
    print("  • Orbital bars/evolution → which atomic orbitals carry the tunneling current")
    print("  • Tables provide quantitative orbital character for paper writing")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
