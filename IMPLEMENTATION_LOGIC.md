# Implementation Logic Flow

## Before Fix (Original Behavior)

```
For each (ik, ien):
  For each projector (iorb):
    Call four(ik, ien, lb=tblm(3,iorb))
      |
      |--> if should_print_mode_b(ik, ien):  # TRUE only for first projector
      |      compute_mode_b_g2(lb, ...)
      |        |--> Print MODE=B:STATE (all states)
      |        |--> Print MODE=B:STATE_LM (only for this lb)
      |
      |--> should_print_mode_b(ik, ien):  # FALSE for remaining projectors
             (Nothing printed)
```

**Result**: Only one `lb` value per (ik, ien) appears in STATE_LM output.

---

## After Fix (New Behavior)

```
For each (ik, ien):
  For each projector (iorb):
    Call four(ik, ien, lb=tblm(3,iorb))
      |
      |--> if should_print_mode_b(ik, ien):  # TRUE only for first projector
      |      init_wlm2_accumulator(...)      # Initialize accumulators
      |      compute_mode_b_state(...)        # Print MODE=B:STATE
      |
      |--> accumulate_wlm2(w0, lb)          # ALL projectors contribute
      |      accum_Wlm2(:,:,lb) += |w0|^2    # Accumulate for this lb
      |      accum_lm_count(:,lb) += 1
      |
  After all projectors done:
    Call flush_state_lm_analysis()
      |--> For each l in [0,1,2,3]:
      |      if accum_lm_count(:,l) > 0:
      |        normalize: Wlm2 / count (average over atoms)
      |        compute_state_lm_for_l(Wlm2, l)
      |          |--> Print MODE=B:STATE_LM for all states at this l
```

**Result**: All `l` values (0,1,2,3) with projectors appear in STATE_LM output.

---

## Key Data Structures

### Accumulation Arrays (in mode_b_guard module)
```fortran
! Dimensions
accum_ngper          ! Number of G-vectors
accum_lmax = 3       ! Max l (0=s, 1=p, 2=d, 3=f)
max_mdim = 7         ! Max 2*l+1 (for f-orbitals)

! Arrays
accum_Wlm2(ig, m_idx, l)     ! Accumulated z-integrated weights
  - ig: 1..accum_ngper        ! G-vector index
  - m_idx: 1..2*l+1           ! Magnetic quantum number index
  - l: 0..3                   ! Angular momentum

accum_lm_count(m_idx, l)     ! Count of contributions per (m,l)
accum_gper(2, ngper)         ! Copy of G-vector coordinates
accum_energy, accum_xyk      ! Energy and k-point for this accumulation
```

### Flow Through Subroutines

1. **init_wlm2_accumulator**: Called once when first projector is encountered
   - Allocates arrays if needed
   - Resets accumulators to zero
   - Stores energy, k-point, G-vectors

2. **accumulate_wlm2**: Called for each projector
   - Computes Wlm2 = sum_z |w0(z,ig,m)|^2 for this projector
   - Adds to accum_Wlm2(ig, m_idx, lb)
   - Increments accum_lm_count(m_idx, lb)

3. **flush_state_lm_analysis**: Called after all projectors
   - For each l with data:
     - Normalize: Wlm2 /= count (average over atoms)
     - Call compute_state_lm_for_l(Wlm2, l)
       - Compute state sorting (by κ)
       - For top states:
         - For each m:
           - Compute separable g²_lm
           - Print WLM_SUMMARY and WLM_CSV lines

---

## Example Output Sequence

```
# First energy, first k-point
Four() called for iorb=1, lb=0 (s-orbital)
  |--> should_print_mode_b = TRUE
  |--> init_wlm2_accumulator()
  |--> compute_mode_b_state()
         [prints MODE=B:STATE for states 1,2,3,...]
  |--> accumulate_wlm2(lb=0)
         accum_Wlm2(:,:,0) += |w0|^2

Four() called for iorb=2, lb=1 (p-orbital)
  |--> should_print_mode_b = FALSE
  |--> accumulate_wlm2(lb=1)
         accum_Wlm2(:,:,1) += |w0|^2

Four() called for iorb=3, lb=2 (d-orbital)
  |--> should_print_mode_b = FALSE
  |--> accumulate_wlm2(lb=2)
         accum_Wlm2(:,:,2) += |w0|^2

After all slabs/projectors:
  flush_state_lm_analysis()
    |--> Process l=0: compute_state_lm_for_l(Wlm2_s, 0)
           [prints MODE=B:STATE_LM for states 1,2,3,... with l=0, m=0]
    |--> Process l=1: compute_state_lm_for_l(Wlm2_p, 1)
           [prints MODE=B:STATE_LM for states 1,2,3,... with l=1, m=-1,0,1]
    |--> Process l=2: compute_state_lm_for_l(Wlm2_d, 2)
           [prints MODE=B:STATE_LM for states 1,2,3,... with l=2, m=-2,-1,0,1,2]
```

---

## Comparison of Output

### Before Fix (Only l=0)
```
WLM_SUMMARY MODE=B:STATE    -0.400 k1,k2= 0.0 0.0 n=1 kappa_bohr=0.15 kappa_ang=0.28
WLM_SUMMARY_CONT g2_bohr=2.34E+00 g2_ang=8.35E+00 norm=1.2E-01

WLM_SUMMARY MODE=B:STATE_LM -0.400 k1,k2= 0.0 0.0 n=1 l=0 m=0
WLM_SUMMARY_CONT g2_bohr=2.10E+00 g2_ang=7.50E+00 norm_lm=8.5E-02
WLM_CSV,STATE_LM -0.400 0.0 0.0 1 0 0 2.10E+00 7.50E+00 8.5E-02
```

### After Fix (All l values)
```
WLM_SUMMARY MODE=B:STATE    -0.400 k1,k2= 0.0 0.0 n=1 kappa_bohr=0.15 kappa_ang=0.28
WLM_SUMMARY_CONT g2_bohr=2.34E+00 g2_ang=8.35E+00 norm=1.2E-01

WLM_SUMMARY MODE=B:STATE_LM -0.400 k1,k2= 0.0 0.0 n=1 l=0 m=0
WLM_SUMMARY_CONT g2_bohr=2.10E+00 g2_ang=7.50E+00 norm_lm=8.5E-02
WLM_CSV,STATE_LM -0.400 0.0 0.0 1 0 0 2.10E+00 7.50E+00 8.5E-02

WLM_SUMMARY MODE=B:STATE_LM -0.400 k1,k2= 0.0 0.0 n=1 l=1 m=-1
WLM_SUMMARY_CONT g2_bohr=3.45E+00 g2_ang=1.23E+01 norm_lm=2.1E-02
WLM_CSV,STATE_LM -0.400 0.0 0.0 1 1 -1 3.45E+00 1.23E+01 2.1E-02

WLM_SUMMARY MODE=B:STATE_LM -0.400 k1,k2= 0.0 0.0 n=1 l=1 m=0
WLM_SUMMARY_CONT g2_bohr=1.89E+00 g2_ang=6.75E+00 norm_lm=1.3E-02
WLM_CSV,STATE_LM -0.400 0.0 0.0 1 1 0 1.89E+00 6.75E+00 1.3E-02

... (more l=1,2,3 entries)
```

**Key difference**: Multiple l values, realistic orbital character distribution.
