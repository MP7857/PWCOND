# Per-Channel Integrand Diagnostics for PWCOND

## Summary

This repository now includes diagnostic capabilities for analyzing per-channel integrand differences in the `four.f90` module of PWCOND (Quantum Espresso's ballistic conductance module).

## What Was Implemented

### 1. New Diagnostic Subroutine

A new subroutine `write_fx_full` has been added to `four.f90` that outputs per-slice integrands (fx1-fx7) for all orbital types:

- **s-orbitals (lb=0)**: Outputs fx1
- **p-orbitals (lb=1)**: Outputs fx1, fx2
- **d-orbitals (lb=2)**: Outputs fx1, fx2, fx3, fx4
- **f-orbitals (lb=3)**: Outputs fx1, fx2, fx3, fx4, fx5, fx6

The diagnostic is called for every (lb, kz, ig, ign) combination during the Fourier transform computation.

### 2. Output File

The diagnostic writes to `fx_full.dat` with the following format:

```
# lb kz ig ign  fx1 fx2 fx3 fx4 fx5 fx6 fx7
    0    1    1    1   1.234567890123E+00   0.000000000000E+00   ...
```

Where:
- `lb`: Orbital angular momentum quantum number (0-3)
- `kz`: z-slice index (1 to nz1)
- `ig`: Global g-perp point index
- `ign`: g-perp shell index
- `fx1-fx7`: Integrand components in scientific notation (es20.12)

### 3. Analysis Utility

A Python utility script `analyze_fx.py` is provided for:
- **Verifying** fx_full.dat file format and content
- **Comparing** two fx_full.dat files (e.g., from different nz grids)
- **Computing** basic statistics on the integrand values

No external dependencies required - uses only Python standard library.

## Usage

### Step 1: Run PWCOND with nz1=7

```bash
# Configure your PWCOND input with nz1=7
pwcond.x < input_nz7.in > output_nz7.out

# Backup the diagnostic output
mv fx_full.dat fx_full_nz7.dat
```

### Step 2: Run PWCOND with nz1=11

```bash
# Configure your PWCOND input with nz1=11
pwcond.x < input_nz11.in > output_nz11.out

# Backup the diagnostic output
mv fx_full.dat fx_full_nz11.dat
```

### Step 3: Analyze Results

```bash
# Verify file format
python3 analyze_fx.py verify fx_full_nz7.dat
python3 analyze_fx.py verify fx_full_nz11.dat

# Compare the two files
python3 analyze_fx.py compare fx_full_nz7.dat fx_full_nz11.dat

# Get statistics
python3 analyze_fx.py stats fx_full_nz7.dat
```

## Analysis Capabilities

With the raw fx integrand data from both nz=7 and nz=11 calculations, you can now:

1. **Plot per-channel difference maps** for all (l,m) between nz=7 and nz=11
2. **Compute relative differences** vs z-slice index
3. **Identify radial regions** contributing to integration errors
4. **Generate integrand shape distortion plots**
5. **Perform radial contribution attribution** analysis
6. **Propose minimal-impact fixes** based on identified error sources

## Files Added/Modified

### Modified Files
- `four.f90`: Added `write_fx_full` subroutine and diagnostic calls for all orbital types

### New Files
- `DIAGNOSTIC_OUTPUT.md`: Detailed technical documentation
- `analyze_fx.py`: Python utility for analyzing diagnostic output
- `README_FX_DIAGNOSTICS.md`: This file - overview and usage guide

## Technical Details

### Parallelization
The diagnostic is MPI-safe - only the I/O node writes to prevent race conditions:
```fortran
if (.not. ionode) return
```

### File I/O
- First call creates new file with header
- Subsequent calls append data
- File is closed after each write for crash safety
- Uses unit 333 (ensure no conflicts in your setup)

### Performance Impact
The diagnostic adds minimal overhead:
- Simple file I/O operations
- Only active on I/O node in parallel runs
- No computational overhead, just data recording

## Next Steps

After collecting `fx_full_nz7.dat` and `fx_full_nz11.dat`:

1. **Upload both files** to your analysis environment
2. **Run full diagnostic pipeline** to:
   - Generate per-channel heatmaps
   - Identify error sources (near-core, far-tail, stencil, etc.)
   - Quantify relative nz sensitivity
3. **Propose and implement fix** based on analysis:
   - Regularized evaluation near axis
   - Hybrid integration scheme
   - Analytic continuation for small rz
   - Targeted fixes for unstable channels

## Support

For detailed information about the implementation, see:
- `DIAGNOSTIC_OUTPUT.md` - Complete technical documentation
- `four.f90` - Source code with inline comments
- `analyze_fx.py` - Analysis utility with built-in help

## References

This diagnostic capability was added to support analysis of nz-grid sensitivity in PWCOND's Fourier transform calculations, specifically for debugging integrand differences between different z-slice discretizations.
