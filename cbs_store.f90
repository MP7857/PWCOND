!
! Copyright (C) 2025
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module to store CBS eigenvector data for Mode B calculations
!
MODULE cbs_store
  USE kinds, ONLY: DP
  IMPLICIT NONE
  SAVE
  !
  ! Dimensions
  INTEGER :: nz1_m = 0      ! number of subslabs (kz grid points)
  INTEGER :: ngper_m = 0    ! number of perpendicular G vectors
  INTEGER :: nstl_m = 0     ! number of left-moving states
  INTEGER :: nchanl_m = 0   ! number of propagating channels
  !
  ! CBS eigenvectors over the (kz, ig) plane-wave grid
  ! cbs_vec_l(kz, ig, n) where n runs from 1 to (nstl_m + nchanl_m)
  ! This includes both evanescent and propagating states
  COMPLEX(DP), ALLOCATABLE :: cbs_vec_l(:,:,:)
  !
  ! Flag indicating if eigenvector data is ready
  LOGICAL :: cbs_vec_l_ready = .FALSE.
  !
CONTAINS
  !
  ! Helper subroutine to allocate storage
  SUBROUTINE allocate_cbs_vec(nz1, ngper, nstl, nchanl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nz1, ngper, nstl, nchanl
    !
    IF (ALLOCATED(cbs_vec_l)) DEALLOCATE(cbs_vec_l)
    !
    nz1_m = nz1
    ngper_m = ngper
    nstl_m = nstl
    nchanl_m = nchanl
    !
    ALLOCATE(cbs_vec_l(nz1, ngper, nstl+nchanl))
    cbs_vec_l = (0.0_DP, 0.0_DP)
    cbs_vec_l_ready = .FALSE.
    !
  END SUBROUTINE allocate_cbs_vec
  !
  ! Helper subroutine to deallocate storage
  SUBROUTINE deallocate_cbs_vec()
    IMPLICIT NONE
    !
    IF (ALLOCATED(cbs_vec_l)) DEALLOCATE(cbs_vec_l)
    cbs_vec_l_ready = .FALSE.
    nz1_m = 0
    ngper_m = 0
    nstl_m = 0
    nchanl_m = 0
    !
  END SUBROUTINE deallocate_cbs_vec
  !
END MODULE cbs_store
