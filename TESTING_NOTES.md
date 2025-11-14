# Testing Notes for Multi-Orbital STATE_LM Fix

## What Changed

The code has been modified to support **full orbital-resolved tunneling analysis**. Previously, only one angular momentum channel (typically `l=0`, s-orbital) appeared in `MODE=B:STATE_LM` output because only the first projector's data was used. Now, **all** angular momentum channels (s, p, d, f) from all projectors contribute to the analysis.

## Files Modified

1. **four.f90**: 
   - Extended `mode_b_guard` module with accumulation infrastructure
   - Split `compute_mode_b_g2` into `compute_mode_b_state` (MODE=B:STATE only)
   - Added accumulation and flushing logic for STATE_LM
   
2. **scatter_forw.f90**: 
   - Added call to flush accumulated STATE_LM data after all projectors processed

## How to Test

### 1. Build the Code

```bash
cd /path/to/PWCOND
make clean
make
```

If the build fails, check:
- Module dependencies are correct
- No syntax errors in the modified files
- The parent Quantum ESPRESSO installation is properly configured

### 2. Run a Test Calculation

Run a PWCOND calculation with CBS (Complex Band Structure) analysis enabled. The calculation should include atoms with different orbital characters (s, p, d, f).

### 3. Check the Output

#### Before the Fix
You would see output like:
```
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  0   0  ...
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      2  0   0  ...
WLM_CSV,STATE_LM   -0.350  0.000000  0.000000      1  0   0  ...
```
Notice: **ALL entries have `l=0` (s-orbital only)**

#### After the Fix
You should see output like:
```
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  0   0  ...
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  1  -1  ...
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  1   0  ...
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  1   1  ...
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  2  -2  ...
WLM_CSV,STATE_LM   -0.400  0.000000  0.000000      1  2  -1  ...
...
```
Notice: **Entries now include `l=1` (p), `l=2` (d), `l=3` (f) channels**

### 4. Verify Analysis

Run the analysis script:
```bash
python3 analyze_wlm.py output_file.txt
```

The analysis should now show:
- Multiple orbital types (s, p, d, f) contributing to tunneling at each energy
- Orbital contributions that vary realistically with energy
- No artificial dominance of s-orbitals at all energies

## Expected Behavior

### MODE=B:STATE Output
- **Unchanged**: Still prints once per (ik, ien) pair
- Shows decay constant κ and average transverse momentum ⟨g²⟩ for each state

### MODE=B:STATE_LM Output
- **Changed**: Now includes **all** angular momentum channels (l=0,1,2,3) present in the system
- Each (n, l, m) combination shows its contribution to that tunneling state
- Multiple l values should appear at each energy (not just l=0)

## Debugging

If the output still only shows l=0:

1. **Check that projectors with l>0 exist**: 
   - Look for atoms with p, d, or f character in your system
   - Verify the pseudopotentials include these orbitals

2. **Check accumulation is working**:
   - Add debug print statements in `accumulate_wlm2` to verify it's being called for different l values
   - Check that `accum_lm_count` has non-zero values for l>0

3. **Check flushing is happening**:
   - Add debug print statement in `flush_state_lm_analysis` to verify it's being called
   - Verify `accum_has_data` is TRUE when flushing

## Compatibility

The output format is **fully compatible** with the existing `analyze_wlm.py` script:
- Two-line format preserved (WLM_SUMMARY + WLM_SUMMARY_CONT)
- CSV format unchanged (WLM_CSV,STATE_LM)
- Column order and units remain the same

## Performance Impact

- **Minimal**: Accumulation arrays are small (typically <1 MB)
- **One-time cost**: Flushing happens once per (ik, ien) pair
- **No repeated computation**: Each l value is analyzed only once using accumulated data

## Rollback

If issues arise, you can revert to the previous version:
```bash
git checkout be604e9  # or the appropriate commit before the change
```

## Contact

For issues or questions, please open an issue in the repository with:
- Full build log
- Sample output showing the issue
- System information (pseudopotentials, atoms, etc.)
