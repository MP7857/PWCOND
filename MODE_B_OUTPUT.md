# Mode B Output Filtering Guide

## Overview

The Mode B output provides state-resolved and energy-dependent transverse momentum analysis (`⟨g²⟩`) for Complex Band Structure (CBS) eigenvectors. This is essential for understanding tunneling transport in ballistic conductance calculations.

## What is Mode B?

Mode B computes the energy- and state-dependent weighted average of transverse momentum squared:

```
⟨g²⟩^(n) = Σ_ig g²(ig) |C^(n)(ig)|² / Σ_ig |C^(n)(ig)|²
```

where:
- `C^(n)(ig)` are the CBS eigenvector plane-wave components for state `n`
- `g` is the transverse momentum vector (perpendicular to transport direction)
- The sum is over all transverse plane-wave components `ig`

## Mode A vs Mode B

- **Mode A (LM)**: Energy-independent baseline showing orbital-projector `⟨g²⟩` per `(l,m)` channel
  - Printed only once at first energy/k-point
  - Shows intrinsic orbital symmetry effects (e.g., m≠0 suppression)
  
- **Mode B (STATE)**: Energy-dependent `⟨g²⟩` for actual CBS evanescent modes
  - Varies with energy as states evolve
  - Filtered to show only physically relevant (slowest-decaying) states

## Output Filtering Parameters

All filtering parameters are defined in `four.f90` in the `compute_mode_b_g2` subroutine:

### `NKEEP` (default: 8)
Number of slowest-decaying states to keep at each energy. States with smallest Im(k_z) = κ are most relevant for tunneling.

```fortran
INTEGER, PARAMETER :: NKEEP = 8
```

**Recommendation**: 3-8 states are typically sufficient. Increase if studying many competing channels.

### `KAPPA_MAX` (default: 2.0 Bohr⁻¹)
Maximum decay constant threshold. Only states with κ ≤ KAPPA_MAX are printed. Set ≤0 to disable.

```fortran
REAL(DP), PARAMETER :: KAPPA_MAX = 2.0_DP
```

**Recommendation**: 
- For vacuum tunneling: 0.5-1.0 Bohr⁻¹
- For barrier tunneling: 1.0-2.0 Bohr⁻¹
- Disable (≤0) to see all NKEEP states regardless of decay

### `MIN_WT` (default: 1.0E-6)
Minimum state weight (norm Σ|C^(n)(ig)|²) threshold. States below this are skipped as numerically irrelevant.

```fortran
REAL(DP), PARAMETER :: MIN_WT = 1.0E-6_DP
```

**Recommendation**: Keep at 1e-6 to 1e-8. Lower values risk numerical noise.

### `ENERGY_STRIDE` (default: 1)
Print Mode B every N-th energy point. Use >1 to reduce output for dense energy scans.

```fortran
INTEGER, PARAMETER :: ENERGY_STRIDE = 1
```

**Recommendation**: 
- For dense scans (>100 energies): try 2-5
- For coarse scans (<50 energies): keep at 1

### `ENABLE_CSV` (default: .FALSE.)
Enable CSV output to separate file `wlm_mode_b.csv` for analysis-friendly data.

```fortran
LOGICAL, PARAMETER :: ENABLE_CSV = .FALSE.
```

**CSV Format**:
```
MODE,TYPE,METRIC,E_eV,k1,k2,n,Re_kz,Im_kz,kappa_Bohrm1,g2_Bohrm2
MODE,B,STATE,6.700,0.000000,0.000000,    12,  0.12345678,  0.23456789,  0.4567,  1.23456E+00
```

## Unit Conversions

### Output Units
- **κ (kappa)**: Bohr⁻¹
- **Re(k_z), Im(k_z)**: Bohr⁻¹
- **⟨g²⟩**: Bohr⁻²

### Conversion Factors
To convert for plotting:
- **Bohr⁻² to Å⁻²**: multiply by (1 Bohr / 0.529177 Å)² ≈ **3.57106**
- **Bohr⁻¹ to Å⁻¹**: multiply by (1 Bohr / 0.529177 Å) ≈ **1.88973**

### Internal Representation
The code uses these conversions:
- `kvall` stores k-values in units of 2π/a
- `gper` stores G-vectors in reciprocal lattice units
- `tpiba = 2π/a` [Bohr⁻¹] converts to physical units
- `κ = |Im(kvall)| × tpiba` gives decay constant in Bohr⁻¹

## Output Interpretation

### Standard Output Format
```
WLM_SUMMARY MODE=B:STATE 6.700 k1,k2= 0.000000  0.000000 n=   12 kappa=  0.4567 g2= 1.23456E+00 units:kappa=Bohr^-1 g2=Bohr^-2
```

Fields:
- `E`: Energy in eV
- `k1,k2`: Transverse k-point coordinates
- `n`: CBS state index (matches CBS_KZ output)
- `kappa`: Decay constant κ = |Im(k_z)| in Bohr⁻¹
- `g2`: ⟨g²⟩ in Bohr⁻²

### Physical Interpretation
- **Smaller κ**: Slower decay → more relevant for long-range tunneling
- **Larger ⟨g²⟩**: More transverse momentum → typically faster decay
- **Expected correlation**: κ should generally increase with ⟨g²⟩

## Sanity Checks

Before using Mode B data for analysis, verify:

1. **Energy variation**: Compare two different energies. For the same state index (or state with smallest Im(k_z)), ⟨g²⟩ should differ across energies.

2. **Decay correlation**: At fixed energy, larger κ should tend to correlate with larger ⟨g²⟩.

3. **Normalization**: The state weight (internally computed) should be O(1) if eigenvectors are properly normalized. States with weight << 1 may have significant local/projector character not captured in ⟨g²⟩.

4. **State indexing**: Verify `n` matches the `CBS_KZ` output block for the same energy.

## Common Use Cases

### Reduce Output Volume
For a dense energy scan producing too much output:
```fortran
INTEGER, PARAMETER :: NKEEP = 5           ! Fewer states
INTEGER, PARAMETER :: ENERGY_STRIDE = 3   ! Every 3rd energy
REAL(DP), PARAMETER :: KAPPA_MAX = 1.0_DP ! Tighter threshold
```

### Analysis-Ready Data
For post-processing and plotting:
```fortran
LOGICAL, PARAMETER :: ENABLE_CSV = .TRUE. ! Export to CSV
INTEGER, PARAMETER :: NKEEP = 8           ! Keep more states
INTEGER, PARAMETER :: ENERGY_STRIDE = 1   ! All energies
```

### Focus on Dominant Channel
To track only the most relevant mode:
```fortran
INTEGER, PARAMETER :: NKEEP = 1           ! Single state
REAL(DP), PARAMETER :: KAPPA_MAX = 0.5_DP ! Very slow decay only
```

## Troubleshooting

### No Mode B Output
- Check that CBS eigenvector data is available (`cbs_vec_l_ready = .TRUE.`)
- Mode B only prints after CBS calculation completes
- Verify `ENERGY_STRIDE` isn't skipping all energies

### Too Much Output
- Increase `ENERGY_STRIDE` (e.g., 2-5)
- Decrease `NKEEP` (e.g., 3-5)
- Reduce `KAPPA_MAX` (e.g., 0.5-1.0 Bohr⁻¹)
- Enable `ENABLE_CSV` and redirect CSV to separate file

### Unexpected ⟨g²⟩ Values
- Check unit conversions (output is in Bohr⁻²)
- Verify states have sufficient plane-wave weight (check internally computed norm)
- Compare with Mode A baseline for orbital character

### CSV File Not Created
- Check write permissions in working directory
- Look for warning message: "Warning: Failed to open CSV file, CSV output disabled"
- Verify `ENABLE_CSV = .TRUE.`

## References

This implementation follows recommendations from the problem statement addressing:
1. Output volume control via selective filtering
2. Energy stride for dense scans
3. CSV format for analysis
4. Unit consistency verification
5. Documentation of filtering strategy

For more details on the CBS method and tunneling transport, see the main PWCOND documentation.
