# WLM_SUMMARY Analysis Tool

This directory contains tools for analyzing the `WLM_SUMMARY` output from the PWCOND quantum transport code.

## Overview

The `analyze_wlm.py` script parses WLM_SUMMARY output and generates visualization plots and comprehensive tables to help understand:

1. **Where tunneling happens** - by identifying energies with minimum decay constant κ(E)
2. **Tunneling channel characteristics** - via transverse momentum ⟨g²⟩(E) analysis
3. **Orbital contributions** - which atomic orbitals (s, p, d, f) carry the tunneling current
4. **Quantitative tables** - for direct use in research papers with clear physics statements

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

The script generates both visualization plots and analysis tables:

#### Plots (PNG files)

1. **kappa_vs_E.png** - Decay constant κ vs energy
   - Deep minima indicate strong tunneling resonances
   - The minimum κ corresponds to the dominant tunneling channel

2. **g2_vs_E.png** - Average transverse momentum squared ⟨g²⟩ vs energy
   - Low ⟨g²⟩ indicates more normal-incidence-like states
   - States with both low κ and low ⟨g²⟩ couple most efficiently

3. **orbital_contrib_best_E.png** - Orbital decomposition bar chart
   - Shows which atomic orbitals dominate at the best tunneling energy
   - Helps identify the physical character of the tunneling state

4. **orbital_evolution_vs_E.png** - Orbital character evolution with energy (best state)
   - Stacked area plot showing how s, p, d, f contributions change with energy
   - Shows the orbital character of the single best (slowest-decaying) state at each energy

5. **orbital_vs_energy_weighted.png** - Weighted orbital contributions vs energy (NEW)
   - Weighted-average orbital fractions across ALL tunneling-relevant states per energy
   - Uses filtered states (small κ, reasonable norm) weighted by importance
   - Provides more comprehensive view than single-state analysis
   - Particularly useful when multiple states contribute to tunneling at a given energy

#### Tables (CSV and TXT files)

1. **wlm_tables_top_states.csv** - Top tunneling states per energy
   - Energy, state index, decay constant κ, transverse momentum ⟨g²⟩
   - Dominant orbital characters with percentages
   - Directly usable for Table 1 in research papers

2. **wlm_tables_orbital_character.csv** - Detailed orbital character breakdown
   - Complete (l,m) decomposition for each tunneling state
   - Fractional weights showing which orbitals contribute
   - Provides data for Table 2 in research papers

3. **wlm_tables_summary.txt** - Human-readable summary
   - Formatted text suitable for copying into paper drafts
   - Ranks states by decay constant with orbital character analysis
   - Provides clear physics statements about tunneling

4. **wlm_orbital_vs_energy.csv** - Weighted orbital fractions per energy (NEW)
   - Energy-resolved s, p, d, f contributions computed via weighted averaging
   - Columns: E_eV, frac_s, frac_p, frac_d, frac_f
   - Each row represents the collective orbital character at that energy
   - Can be plotted directly or imported for further analysis

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
  s: 0.7658 (76.6%)    ← s-orbital dominates
  p_x: 0.1595 (16.0%)  ← p_x orbital
  p_z: 0.0638 (6.4%)   ← p_z orbital
  p_y: 0.0108 (1.1%)   ← p_y orbital
Saved: orbital_evolution_vs_E.png

======================================================================
Generating comprehensive analysis tables...
======================================================================

Generating Table 1: wlm_tables_top_states.csv
Saved: wlm_tables_top_states.csv

Generating Table 2: wlm_tables_orbital_character.csv
Saved: wlm_tables_orbital_character.csv

Generating text summary: wlm_tables_summary.txt
Saved: wlm_tables_summary.txt
----------------------------------------------------------------------

Analysis complete!
```

#### Sample Table 1: Top Tunneling States (wlm_tables_top_states.csv)

```csv
Energy_eV,State_n,kappa_Bohr,kappa_Ang,g2_Bohr2,g2_Ang2,norm,dominant_orbitals
5.400,7,0.0987,0.1864,1.12e+00,4.01e+00,1.57e-01,"s(77%), p_x(16%), p_z(6%)"
6.900,5,0.0450,0.0851,2.05e+01,7.32e+01,1.23e+00,"d_z²−1(54%), f_z(5z²−3r²)(31%), p_z(15%)"
```

#### Sample Table 2: Orbital Character (wlm_tables_orbital_character.csv)

```csv
Energy_eV,State_n,kappa_Ang,l,m,orbital_name,norm_lm,fractional_weight,g2_Ang2
5.400,7,0.1864,0,0,s,1.200e-01,0.7658,3.93e+00
5.400,7,0.1864,1,-1,p_x,2.500e-02,0.1595,4.11e+00
6.900,5,0.0851,2,0,d_z²−1,5.500e-01,0.5392,6.43e+01
6.900,5,0.0851,3,0,f_z(5z²−3r²),3.200e-01,0.3137,6.61e+01
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

The analysis identifies which atomic orbitals carry the tunneling current:

| Orbital | Symbol | Physical Interpretation |
|---------|--------|------------------------|
| l=0, m=0 | **s** | Spherically symmetric |
| l=1, m=0 | **p_z** | Aligned with transport direction |
| l=1, m=±1 | **p_x, p_y** | Perpendicular to transport |
| l=2, m=0 | **d_z²−1** | Axial d-orbital |
| l=2, m=±1 | **d_xz, d_yz** | Off-axis d-orbitals |
| l=2, m=±2 | **d_xy, d_x²−y²** | Planar d-orbitals |
| l=3, m=0 | **f_z(5z²−3r²)** | Axial f-orbital |
| l=3, m=±1 | **f_x(5z²−r²), f_y(5z²−r²)** | Off-axis f-orbitals |
| l=3, m=±2,±3 | **f_xyz, etc.** | Various f-orbital characters |

High contributions from **p_z** or **d_z²−1** orbitals indicate tunneling aligned with the transport direction, which typically leads to higher transmission coefficients.

### Weighted Orbital Analysis (NEW)

The weighted orbital analysis provides a more comprehensive picture of tunneling character at each energy by considering **all** tunneling-relevant states, not just the single best state.

#### How it works:

1. **Filter states** by thresholds:
   - κ (decay constant) ≤ 0.5 Bohr⁻¹ (configurable)
   - norm ≥ 10⁻³ (configurable)
   
2. **Weight each state** by its importance:
   - Default: uses state norm as weight
   - Alternatively: could use exp(-2κL) for barrier thickness L
   
3. **Compute weighted average** of orbital fractions:
   ```
   F_l(E) = Σ_j [weight_j × frac_l^(j)] / Σ_j weight_j
   ```
   where j sums over all filtered states at energy E

#### When to use weighted vs single-state analysis:

| Scenario | Best Analysis |
|----------|---------------|
| One dominant state per energy | Either approach works similarly |
| Multiple states contribute | **Weighted** gives complete picture |
| Near resonances/crossings | **Weighted** captures all contributions |
| Far from resonances | Single-state may be sufficient |

#### Output:

- **CSV**: `wlm_orbital_vs_energy.csv` with columns E_eV, frac_s, frac_p, frac_d, frac_f
- **Plot**: `orbital_vs_energy_weighted.png` showing stacked area curves

This provides the data structure suggested in the problem statement for combining orbital contributions across states at each energy.

### Using Tables in Research Papers

#### Table 1: Top Tunneling States

This table provides a quick overview of the dominant tunneling channels:

- **Use in papers**: "At E = 5.4 eV, the slowest-decaying state (n=7, κ=0.186 Å⁻¹) is dominated by s-orbital character (77%), with smaller p contributions."
- Shows which states dominate at each energy
- Provides quantitative orbital percentages for clear physics statements

#### Table 2: Detailed Orbital Character

This table gives complete (l,m) decomposition:

- **Use in papers**: "The dominant tunneling channel at E = 6.9 eV exhibits mixed d-f character, with d_z²−1 (54%) and f_z(5z²−3r²) (31%) being the primary contributors."
- Allows precise identification of specific orbital contributions
- Supports detailed analysis of tunneling mechanisms

#### Text Summary

The generated `wlm_tables_summary.txt` provides ready-to-use text:

- Can be directly adapted for Methods or Results sections
- Includes all key physics statements
- Properly formatted with correct notation

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

## Making Physics Statements from the Analysis

### Generic Template for Papers

Based on the analysis output, you can write statements like:

> "At each energy, we rank the complex band structure (CBS) evanescent states by their decay constant κ = |Im k_z|. The states with the smallest κ dominate the tunneling current. For each such state we compute the transverse momentum moment ⟨g_⊥²⟩ and orbital character using the CBS plane-wave amplitudes.
>
> We find that near E = 6.9 eV the slowest-decaying channel (κ = 0.085 Å⁻¹) is predominantly d_z²−1 (54%) and f_z(5z²−3r²) (31%) like, with smaller p_z contributions (15%). At lower energies (E = 5.4 eV), the dominant channel becomes more s-like (77%), with smaller p-orbital character."

### Workflow for Paper Writing

1. **Run the analysis**: `python3 analyze_wlm.py out.txt`

2. **Examine Table 1** (`wlm_tables_top_states.csv`):
   - Identify which energies have the best tunneling (lowest κ)
   - Note the dominant orbital characters

3. **Examine Table 2** (`wlm_tables_orbital_character.csv`):
   - Get precise percentages for each orbital contribution
   - Identify trends in orbital character vs energy

4. **Examine weighted orbital CSV** (`wlm_orbital_vs_energy.csv`) (NEW):
   - Get energy-dependent orbital fractions averaged across all relevant states
   - Use for comprehensive statements about tunneling character vs energy
   - Plot or analyze trends in collective orbital contributions

5. **Use the text summary** (`wlm_tables_summary.txt`):
   - Copy relevant sections into your paper
   - Adapt the language to match your writing style

6. **Include plots in your paper**:
   - κ(E) plot shows where tunneling is strongest
   - Orbital evolution plot shows energy-dependent character
   - Weighted orbital plot shows comprehensive multi-state analysis
   - Orbital bar chart illustrates dominant channels

### Example Physics Statements

**From Table 1 analysis**:
> "The slowest-decaying state at E = 6.9 eV (κ = 0.085 Å⁻¹) exhibits mixed d-f character, dominated by d_z²−1 (54%) and f_z(5z²−3r²) (31%) orbitals."

**From orbital evolution plot**:
> "As energy increases from 5.4 eV to 6.9 eV, the dominant tunneling channel transitions from predominantly s-character to mixed d-f character, indicating a change in the electronic structure of the evanescent states."

**From weighted orbital analysis (NEW)**:
> "Using weighted-average orbital decomposition across all tunneling-relevant states (κ < 0.5 Bohr⁻¹), we find that at E = 5.4 eV the collective tunneling character is 65% s-orbital and 30% p-orbital, while at E = 6.9 eV it shifts to 48% d-orbital and 28% f-orbital character. This weighted analysis accounts for contributions from multiple evanescent channels and provides a comprehensive view of the orbital composition of the tunneling current."

**From combined κ(E) and orbital analysis**:
> "The minimum decay constant occurs at E = 6.9 eV where the tunneling channel is dominated by d_z²−1 and f_z(5z²−3r²) orbitals, both aligned with the transport direction, which explains the enhanced transmission at this energy."

## References

For more information on the physical interpretation of these curves, see:

- M. Pourfath, "The Non-Equilibrium Green's Function Method for Nanoscale Device Simulation" (2014)
- A. Smogunov et al., "Ballistic conductance of magnetic Co and Ni nanowires" Phys. Rev. B 70, 045417 (2004)

## Contact

For issues or questions about this tool, please open an issue in the repository.
