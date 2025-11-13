# WLM_SUMMARY Analysis Tool

This directory contains tools for analyzing the `WLM_SUMMARY` output from the PWCOND quantum transport code.

## Overview

The `analyze_wlm.py` script parses WLM_SUMMARY output and generates visualization plots to help understand:

1. **Where tunneling happens** - by identifying energies with minimum decay constant κ(E)
2. **Tunneling channel characteristics** - via transverse momentum ⟨g²⟩(E) analysis
3. **Orbital contributions** - which atomic orbitals (s, p, d, f) carry the tunneling current

## Requirements

- Python 3.x
- matplotlib library

Install matplotlib if needed:
```bash
pip3 install matplotlib
```

## Usage

### Basic Usage

Run the script on your PWCOND output file:

```bash
python3 analyze_wlm.py out.txt
```

If no filename is provided, it defaults to `out.txt`:

```bash
python3 analyze_wlm.py
```

### Output Files

The script generates three PNG files:

1. **kappa_vs_E.png** - Decay constant κ vs energy
   - Deep minima indicate strong tunneling resonances
   - The minimum κ corresponds to the dominant tunneling channel

2. **g2_vs_E.png** - Average transverse momentum squared ⟨g²⟩ vs energy
   - Low ⟨g²⟩ indicates more normal-incidence-like states
   - States with both low κ and low ⟨g²⟩ couple most efficiently

3. **orbital_contrib_best_E.png** - Orbital decomposition bar chart
   - Shows which atomic orbitals dominate at the best tunneling energy
   - Helps identify the physical character of the tunneling state

### Example Output

```
Parsing WLM_SUMMARY output from: out.txt
======================================================================
Found 10 MODE=B:STATE entries
Found 40 MODE=B:STATE_LM entries
======================================================================

Generating plots...
----------------------------------------------------------------------
Saved: kappa_vs_E.png
Saved: g2_vs_E.png

Best tunneling state:
  Energy: 5.400 eV
  State n: 7
  κ: 0.1864 Å⁻¹
Saved: orbital_contrib_best_E.png

Orbital contributions (normalized):
  (l=0,m=+0): 0.7658 (76.6%)    ← s-orbital dominates
  (l=1,m=-1): 0.1595 (16.0%)    ← p_x orbital
  (l=1,m=+0): 0.0638 (6.4%)     ← p_z orbital
  (l=1,m=+1): 0.0108 (1.1%)     ← p_y orbital
----------------------------------------------------------------------

Analysis complete!
```

## Interpreting the Results

### κ(E) Curve (Decay Constant)

The decay constant κ measures how quickly the wavefunction decays into the tunneling barrier:

- **Low κ** (sharp minima) → slow decay → **strong tunneling**
- **High κ** (broad plateaus) → fast decay → weak tunneling
- Local minima identify the dominant tunneling resonances

### ⟨g²⟩(E) Curve (Transverse Momentum)

The average transverse momentum squared characterizes the angular distribution:

- **Low ⟨g²⟩** → state is more normal-incidence-like (small transverse component)
- **High ⟨g²⟩** → state has large oblique momentum components
- States with low ⟨g²⟩ couple more efficiently to electrodes at k∥≈0

### Orbital Contributions

The bar chart shows the relative contribution of each (l,m) channel:

| Angular Momentum | Orbital Character | Physical Interpretation |
|------------------|-------------------|------------------------|
| l=0, m=0 | s-like | Spherically symmetric |
| l=1, m=0 | p_z | Aligned with transport direction |
| l=1, m=±1 | p_x, p_y | Perpendicular to transport |
| l=2, m=... | d-orbitals | Various d-orbital characters |
| l=3, m=... | f-orbitals | Various f-orbital characters |

High contributions from p_z or d_{z²-r²} orbitals indicate tunneling aligned with the transport direction.

## Technical Details

### WLM_SUMMARY Output Format

The script supports both legacy single-line and modern two-line output formats.

#### Modern Two-Line Format (Current)

The script parses three types of output from `four.f90`:

#### MODE=B:STATE (State-resolved)
```
WLM_SUMMARY MODE=B:STATE    5.400 k1,k2=  0.000000  0.000000 n=   7 kappa_bohr=   0.0987 kappa_ang=   0.1864
WLM_SUMMARY_CONT g2_bohr=  1.12345E+00 g2_ang=  4.01195E+00 norm= 1.567E-01
```

#### MODE=B:STATE_LM (State and orbital-resolved)
```
WLM_SUMMARY MODE=B:STATE_LM    5.400 k1,k2=  0.000000  0.000000 n=   7 l=  0 m=  0
WLM_SUMMARY_CONT g2_bohr=  1.10000E+00 g2_ang=  3.92821E+00 norm_lm= 1.200E-01
```

#### Legacy Single-Line Format

The script also supports older single-line format where all data is on one line:

```
WLM_SUMMARY MODE=B:STATE    5.400 k1,k2=  0.000000  0.000000 n=   7 kappa_bohr=   0.0987 kappa_ang=   0.1864 g2_bohr=  1.12345E+00 g2_ang=  4.01195E+00 norm= 1.567E-01
```

### Data Fields

**MODE=B:STATE** provides:
- Energy (eV)
- Transverse k-point (k1, k2)
- State index (n)
- Decay constant κ in Bohr⁻¹ and Å⁻¹
- Average transverse momentum ⟨g²⟩ in Bohr⁻² and Å⁻²
- State weight (norm)

**MODE=B:STATE_LM** provides:
- All of the above, plus:
- Angular momentum quantum numbers (l, m)
- Orbital-resolved ⟨g²⟩
- Orbital-resolved weight (norm_lm)

### Algorithm

1. **Parse** WLM_SUMMARY lines handling the two-line format
2. **Find best tunneling energy** by locating minimum κ(E)
3. **Extract orbital data** for that energy
4. **Normalize weights** and generate bar chart
5. **Plot** κ(E) and ⟨g²⟩(E) curves

## Troubleshooting

**No WLM_SUMMARY data found:**
- Ensure your PWCOND calculation includes Complex Band Structure (CBS) analysis
- Check that `four.f90` is being called with CBS eigenvector data available

**Missing plots:**
- Verify matplotlib is installed: `pip3 install matplotlib`
- Check file permissions in the output directory

**Incorrect parsing:**
- The script expects the two-line format (WLM_SUMMARY + WLM_SUMMARY_CONT)
- Verify your output matches the format shown above

## Citation

If you use this analysis tool in your research, please cite:

- The PWCOND quantum transport package
- Relevant methodology papers for complex band structure analysis

## References

For more information on the physical interpretation of these curves, see:

- M. Pourfath, "The Non-Equilibrium Green's Function Method for Nanoscale Device Simulation" (2014)
- A. Smogunov et al., "Ballistic conductance of magnetic Co and Ni nanowires" Phys. Rev. B 70, 045417 (2004)

## Contact

For issues or questions about this tool, please open an issue in the repository.
