# Diagnostic Output for PWCOND fx Analysis

## Overview

This document describes the diagnostic output feature added to `four.f90` for analyzing per-channel integrand differences between different nz grid resolutions.

## What Was Added

### New Subroutine: `write_fx_full`

A new subroutine has been added at the end of `four.f90` that writes per-slice integrands (fx1-fx7) to a diagnostic file.

```fortran
subroutine write_fx_full(lb, kz, ig, ign, f1, f2, f3, f4, f5, f6, f7)
```

**Parameters:**
- `lb`: Orbital angular momentum quantum number (0=s, 1=p, 2=d, 3=f)
- `kz`: z-slice index (1 to nz1)
- `ig`: Global g-perp point index
- `ign`: g-perp shell index
- `f1-f7`: The seven fx integrand components

### Integration Points

The subroutine is called after computing w0 arrays for each orbital type:

1. **lb=0 (s-orbital)**: Outputs fx1 only (others are 0)
2. **lb=1 (p-orbital)**: Outputs fx1, fx2 (others are 0)
3. **lb=2 (d-orbital)**: Outputs fx1, fx2, fx3, fx4 (others are 0)
4. **lb=3 (f-orbital)**: Outputs fx1, fx2, fx3, fx4, fx5, fx6 (fx7 is 0)

## Output File Format

The diagnostic data is written to `fx_full.dat` in the working directory.

### File Structure

```
# lb kz ig ign  fx1 fx2 fx3 fx4 fx5 fx6 fx7
    0    1    1    1   1.234567890123E+00   0.000000000000E+00   ...
    0    2    1    1   2.345678901234E+00   0.000000000000E+00   ...
    ...
```

### Column Descriptions

| Column | Type | Description |
|--------|------|-------------|
| lb | integer | Orbital angular momentum (0-3) |
| kz | integer | z-slice index (1 to nz1) |
| ig | integer | Global g-perp point index |
| ign | integer | g-perp shell index |
| fx1-fx7 | real(dp) | Integrand components (es20.12 format) |

## Usage Instructions

### Step 1: Run PWCOND with nz1=7

```bash
# Set nz1=7 in your PWCOND input file
pwcond.x < input_nz7.in > output_nz7.out

# Save the diagnostic output
mv fx_full.dat fx_full_nz7.dat
```

### Step 2: Run PWCOND with nz1=11

```bash
# Set nz1=11 in your PWCOND input file
pwcond.x < input_nz11.in > output_nz11.out

# Save the diagnostic output
mv fx_full.dat fx_full_nz11.dat
```

### Step 3: Analyze the Results

You now have the raw fx integrand data for both grid resolutions:
- `fx_full_nz7.dat` - Data from nz1=7 calculation
- `fx_full_nz11.dat` - Data from nz1=11 calculation

These files can be used to:
1. Plot per-channel difference maps for all (l,m) between nz=7 and nz=11
2. Compute relative differences vs z-slice index
3. Identify radial regions contributing to errors
4. Generate integrand shape distortion plots
5. Perform radial contribution attribution analysis

## Implementation Details

### Parallelization

The subroutine uses the `ionode` check to ensure only the I/O node writes to the file, preventing race conditions in parallel runs:

```fortran
if (.not. ionode) return
```

### File Handling

- On first call: Creates new file with header
- Subsequent calls: Appends to existing file
- Uses unit 333 (ensure this doesn't conflict with other I/O in your setup)

### Thread Safety

The `first` variable is declared as `save` to maintain state across calls:

```fortran
logical, save :: first = .true.
```

## Notes

1. The output file can become quite large for fine grids or many orbitals
2. All values are written in scientific notation with 12 decimal places
3. Unused fx components are written as 0.0 for consistency
4. The file is closed after each write to ensure data persistence in case of crashes

## Future Enhancements

Potential improvements for future versions:
1. Option to output only specific orbital types
2. Binary output format for reduced file size
3. Option to limit output to specific ign values
4. Statistical summaries computed on-the-fly

## References

See the original `four.f90` header comments for information about the bidimensional Fourier transform implementation.
