# w0 Norm Diagnostics for nz-Stability Validation

This directory includes diagnostic tools to validate the nz-stability improvements from Simpson integration.

## Diagnostic Files Generated

When running pwcond with the diagnostic code enabled, the following files are produced:

### `w0_norms.dat`
Compact file with global norms per (l,m) channel:
- **Format**: `nz1  lb  m  norm`
- **Content**: Sum of |w0(kz,ig,m)|² over all kz,ig for each projector call
- **Purpose**: Track how total norm changes with nz1

### `w0_full.dat` (optional)
Full w0 arrays for detailed analysis:
- **Format**: `kz  ig  lb  m  Re  Im`
- **Content**: Complete w0(kz,ig,m) values for all slices and G-vectors
- **Purpose**: Inspect shape changes in w0 as nz increases
- **Note**: Only enable for diagnostic runs (uncomment call in four.f90)

## Usage

### 1. Run with different nz1 values

```bash
# Run with nz1=7
pwcond.x < input_nz7.in > output_nz7.out
mv w0_norms.dat w0_norms_nz7.dat

# Run with nz1=11
pwcond.x < input_nz11.in > output_nz11.out
mv w0_norms.dat w0_norms_nz11.dat
```

### 2. Compare norms across nz values

```bash
python analyze_w0_norms_vs_nz.py w0_norms_nz7.dat w0_norms_nz11.dat
```

This will show:
- Norm values for each (l,m) channel at different nz1
- Relative changes between nz1 values
- Summary indicating if nz-stability is achieved

### 3. Expected Results

For successful nz-stable integration:
- **Relative changes < 5-10%** between nz1=7 and nz1=11
- All orbital types (s, p, d, f) show similar stability
- No systematic drift in any particular m-channel

Large changes (>10%) in specific channels indicate nz-dependent behavior that may need investigation.

## Implementation Details

### Code in four.f90

**`accumulate_w0_norms`**: 
- Called automatically at end of each `four()` invocation
- Computes full norm over all kz,ig for each m-channel
- Appends results to `w0_norms.dat`
- Enabled by default

**`write_w0_full`**: 
- Optional detailed dump of all w0 values
- Disabled by default (uncomment call in `four()` to enable)
- Useful for debugging specific channel behaviors

### Python Analysis Script

`analyze_w0_norms_vs_nz.py`:
- Loads norm data from multiple nz runs
- Groups by (l,m) channel
- Computes relative changes
- Provides summary assessment

## Testing Workflow

1. **Baseline**: Run with original (pre-Simpson) code, compare nz1=7 vs nz1=11
   - Should show significant differences (multi-meV scale)

2. **With Simpson**: Run with Simpson integration, compare nz1=7 vs nz1=11
   - Should show minimal differences (<5% in norms)

3. **Multiple nz values**: Run nz1 = 7, 9, 11, 13 to confirm monotonic convergence

4. **Channel analysis**: If any (l,m) shows large variation, investigate with `w0_full.dat`

## Notes

- Diagnostic output is only written by the I/O node (ionode check in subroutines)
- File unit 320 is used for `w0_norms.dat`
- File unit 321 is used for `w0_full.dat`
- The `save :: first` variable ensures headers are written only once
- m values are stored as m-1 to match Python 0-based indexing conventions
