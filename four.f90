!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Modified by M. Pourfath (2025)
! Extended to include f-orbital (l = 3) calculations.
!
subroutine four(w0, z0, dz, tblm, taunew, r, rab, betar, ik, ien)
!
! This routine computes the bidimensional fourier transform of the
! beta function. It has been implemented for s, p, d, f-orbitals.
!
!   w0(z,g,m)=1/S * \int w(r) \exp{-ig r_\perp} dr_\perp
!   where w(r) - beta function of the alpha's orbital.
!
!   (see Gradshtein "Tables of integrals")
! For a fixed l it computes w0 for all m.
!
! The order of spherical harmonics used:
!             s ;
!             p_z, p_{-x}, p_{-y} ;
!             d_{z^2-1}, d_{-xz}, d_{-yz}, d_{x^2-y^2}, d_{xy}
!             f_{z(5z^2-3r^2)}, f_{x(5z^2-r^2)}, f_{y(5z^2-r^2)},
!             f_{z(x^2-y^2)}, f_{xyz}, f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}
!
! input:  tblm   -  array characterizing the orbital.
!         taunew -  coordinates and radius of the orbital.
!         z0     -  the initial z
!         dz     -  the slab width
!
! output: w0(z, g, m), where
!                      z0< z <z0+dz
!                      g - 2D g-vector
!
  USE kinds, ONLY: DP
  USE constants, ONLY : tpi, fpi
  USE radial_grids, only : ndmx
  USE cell_base, ONLY : alat, tpiba
  USE cond, ONLY : sarea, nz1, ngper, gper, ninsh, gnsh, ngpsh, earr, xyk

implicit none

  integer :: kz, ig, ign, igphi, &
             indexr, iz, lb, ir, nmesh, nmeshs, tblm(4), m_idx, mdim, &
             ik, ien, m_val, abs_m
  real(DP), parameter :: eps=1.d-8
  complex(DP), parameter :: cim=(0.d0, 1.d0)
  real(DP) :: gn, s1, s2, s3, s4, cs, sn, cs2, sn2, cs3, sn3, rz, dz1, zr, &
                   dr, z0, dz,  bessj, taunew(4), r(ndmx),         &
                   rab(ndmx), betar(ndmx), &
                   sum_w2, sum_g2w2, g2_avg, gmag2, w2
  real(DP), allocatable :: x1(:), x2(:), x3(:), x4(:), x5(:), x6(:)
  real(DP), allocatable :: fx1(:), fx2(:), fx3(:), fx4(:), fx5(:), fx6(:), zsl(:)
  complex(DP) :: w0(nz1, ngper, 7)
  complex(DP), allocatable :: wadd(:,:), wadd2(:,:), wadd3(:,:)
  complex(DP) :: t1, t2, t3, t4, t5, t6, t7, wa1, wa2, wa3


  allocate( x1(0:ndmx) )
  allocate( x2(0:ndmx) )
  allocate( x3(0:ndmx) )
  allocate( x4(0:ndmx) )
  allocate( x5(0:ndmx) )
  allocate( x6(0:ndmx) )
  allocate( fx1( nz1 ) )
  allocate( fx2( nz1 ) )
  allocate( fx3( nz1 ) )
  allocate( fx4( nz1 ) )
  allocate( fx5( nz1 ) )
  allocate( fx6( nz1 ) )
  allocate( zsl( nz1) )
  allocate( wadd( nz1, ngper ) )
  allocate( wadd2( nz1, ngper ) )
  allocate( wadd3( nz1, ngper ) )

  lb = tblm(3)
  nmesh=indexr(taunew(4)*alat,ndmx,r)
  dz1=dz/nz1
  zsl(1)=(z0+dz1*0.5d0-taunew(3))*alat
  do kz = 2, nz1
    zsl(kz) = zsl(kz-1)+dz1*alat
  enddo


  ig=0
  do ign=1, ngpsh

     gn=gnsh(ign)
     do kz=1, nz1
       if (abs(zsl(kz))+eps.le.taunew(4)*alat) then
         iz=indexr(zsl(kz),nmesh,r)
         if ((nmesh-iz)/2*2.eq.nmesh-iz) then
            nmeshs=nmesh
         else
            nmeshs=nmesh+1
         endif
         do ir=iz, nmeshs
            rz=sqrt(r(ir)**2-zsl(kz)**2)
            if (lb.eq.0) then
               x1(ir)=betar(ir)*bessj(0,gn*rz)
            elseif (lb.eq.1) then
               x1(ir)=betar(ir)*bessj(1,gn*rz)/r(ir)*rz
               x2(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)
            elseif (lb.eq.2) then
               x1(ir)=betar(ir)*bessj(2,gn*rz)*rz**2/r(ir)**2
               x2(ir)=betar(ir)*bessj(1,gn*rz)*rz/r(ir)**2
               x3(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)**2
               x4(ir)=betar(ir)*bessj(0,gn*rz)
            elseif (lb.eq.3) then
               x1(ir)=betar(ir)*bessj(3,gn*rz)*rz**3/r(ir)**3
               x2(ir)=betar(ir)*bessj(2,gn*rz)*rz**2/r(ir)**3
               x3(ir)=betar(ir)*bessj(1,gn*rz)*rz/r(ir)**3
               x4(ir)=betar(ir)*bessj(1,gn*rz)*rz**3/r(ir)**3
               x5(ir)=betar(ir)*bessj(0,gn*rz)/r(ir)**3
               x6(ir)=betar(ir)*bessj(0,gn*rz)*rz**2/r(ir)**3
            else
               call errore ('four','ls not programmed ',1)
            endif
         enddo
         call simpson(nmeshs-iz+1,x1(iz),rab(iz),fx1(kz))
         if (iz.eq.1) then
            dr=r(iz)
         else
            dr=r(iz)-r(iz-1)
         endif
         zr=r(iz)-abs(zsl(kz))
         if (lb.eq.0) then
            if (iz.eq.1) then
               x1(iz-1)=betar(iz)-betar(iz)/dr*zr
            else
               x1(iz-1)=betar(iz)-(betar(iz)-betar(iz-1))/dr*zr
            endif
            fx1(kz)=fx1(kz)+(x1(iz-1)+x1(iz))*0.5d0*zr
         else
            fx1(kz)=fx1(kz)+x1(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x2(iz),rab(iz),fx2(kz))
         endif
         if (lb.eq.1) then
            if(iz.eq.1) then
              x2(iz-1)=0.d0
            else
              x2(iz-1)=(betar(iz)-(betar(iz)-   &
                        betar(iz-1))/dr*zr)/abs(zsl(kz))
            endif
            fx2(kz)=fx2(kz)+(x2(iz-1)+x2(iz))*0.5d0*zr
         endif
         if (lb.eq.2) then
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            if(iz.eq.1) then
               x3(iz-1)=0.d0
               x4(iz-1)=0.d0
            else
               x3(iz-1)=(betar(iz)-(betar(iz)-   &
                         betar(iz-1))/dr*zr)/abs(zsl(kz))**2
               x4(iz-1)=betar(iz)-(betar(iz)-       &
                         betar(iz-1))/dr*zr
            endif
            fx3(kz)=fx3(kz)+(x3(iz-1)+x3(iz))*0.5d0*zr
            fx4(kz)=fx4(kz)+(x4(iz-1)+x4(iz))*0.5d0*zr
         elseif (lb.eq.3) then
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            fx3(kz)=fx3(kz)+x3(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            fx4(kz)=fx4(kz)+x4(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x5(iz),rab(iz),fx5(kz))
            call simpson(nmeshs-iz+1,x6(iz),rab(iz),fx6(kz))
            if(iz.eq.1) then
               x5(iz-1)=0.d0
            else
               x5(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr)/(abs(zsl(kz))**3)
            endif
            x6(iz-1)=0.d0
            fx5(kz)=fx5(kz)+(x5(iz-1)+x5(iz))*0.5d0*zr
            fx6(kz)=fx6(kz)+(x6(iz-1)+x6(iz))*0.5d0*zr
         endif
       else
          fx1(kz)=0.d0
          fx2(kz)=0.d0
          fx3(kz)=0.d0
          fx4(kz)=0.d0
          fx5(kz)=0.d0
          fx6(kz)=0.d0
       endif
     enddo
     do igphi=1, ninsh(ign)
        ig=ig+1
        if (gn.gt.eps) then
          cs=gper(1,ig)*tpiba/gn
          sn=gper(2,ig)*tpiba/gn
        else
          cs=0.d0
          sn=0.d0
        endif
        cs2=cs**2-sn**2
        sn2=2*cs*sn
        cs3 = cs * (4.d0*cs**2 - 3.d0)
        sn3 = sn * (3.d0 - 4.d0*sn**2)

        do kz=1, nz1
            if (lb.eq.0) then
               w0(kz,ig,1)=fx1(kz)
            elseif (lb.eq.1) then
               w0(kz,ig,2)=cs*fx1(kz)
               w0(kz,ig,1)=fx2(kz)
               w0(kz,ig,3)=sn*fx1(kz)
            elseif (lb.eq.2) then
               w0(kz,ig,5)=sn2*fx1(kz)
               w0(kz,ig,2)=cs*fx2(kz)
               w0(kz,ig,1)=fx3(kz)
               w0(kz,ig,3)=sn*fx2(kz)
               w0(kz,ig,4)=cs2*fx1(kz)
               wadd(kz,ig)=fx4(kz)
            elseif (lb.eq.3) then
               w0(kz,ig,1)=fx5(kz)
               wadd(kz,ig)=fx6(kz)
               w0(kz,ig,2)=cs*fx3(kz)
               wadd2(kz,ig)=cs*fx4(kz)
               w0(kz,ig,3)=sn*fx3(kz)
               wadd3(kz,ig)=sn*fx4(kz)
               w0(kz,ig,4)=cs2*fx2(kz)
               w0(kz,ig,5)=sn2*fx2(kz)
               w0(kz,ig,6)=cs3*fx1(kz)
               w0(kz,ig,7)=sn3*fx1(kz)
            endif
        enddo
     enddo

  enddo

  if (lb.eq.0) then
     s1=tpi/sarea/sqrt(fpi)
  elseif (lb.eq.1) then
     s1=tpi/sarea*sqrt(3.d0/fpi)
  elseif (lb.eq.2) then
     s1=-tpi/2.d0/sarea*sqrt(15.d0/fpi)
     s2=tpi/sarea*sqrt(5.d0/tpi/8.d0)
  elseif (lb.eq.3) then
     s1 = tpi/sarea * sqrt(35.d0/(fpi*8.d0))
     s2 = tpi/sarea * sqrt(105.d0/(fpi*4.d0))
     s3 = tpi/sarea * sqrt(21.d0/(fpi*8.d0))
     s4 = tpi/sarea * sqrt(7.d0/(fpi*4.d0))
  endif
  do ig=1, ngper
    do kz=1, nz1
      if (lb.eq.0) then
        w0(kz,ig,1)=s1*w0(kz,ig,1)
      elseif (lb.eq.1) then
        w0(kz,ig,2)=cim*s1*w0(kz,ig,2)
        w0(kz,ig,1)=s1*zsl(kz)*w0(kz,ig,1)
        w0(kz,ig,3)=cim*s1*w0(kz,ig,3)
      elseif (lb.eq.2) then
        w0(kz,ig,5)=s1*w0(kz,ig,5)
        w0(kz,ig,2)=-2.d0*cim*s1*zsl(kz)*w0(kz,ig,2)
        w0(kz,ig,1)=3.d0*zsl(kz)**2*s2*w0(kz,ig,1)-s2*wadd(kz,ig)
        w0(kz,ig,3)=-2.d0*cim*s1*zsl(kz)*w0(kz,ig,3)
        w0(kz,ig,4)=s1*w0(kz,ig,4)
      elseif (lb.eq.3) then
        t1=w0(kz,ig,1);wa1=wadd(kz,ig)
        t2=w0(kz,ig,2);wa2=wadd2(kz,ig)
        t3=w0(kz,ig,3);wa3=wadd3(kz,ig)
        t4=w0(kz,ig,4);t5=w0(kz,ig,5)
        t6=w0(kz,ig,6);t7=w0(kz,ig,7)
        w0(kz,ig,1)=s4*(2.d0*zsl(kz)**3*t1-3.d0*zsl(kz)*wa1)
        w0(kz,ig,2)=-cim*s3*(4.d0*zsl(kz)**2*t2-wa2)
        w0(kz,ig,3)=-cim*s3*(4.d0*zsl(kz)**2*t3-wa3)
        w0(kz,ig,4)=s2*zsl(kz)*t4
        w0(kz,ig,5)=s2*zsl(kz)*t5
        w0(kz,ig,6)=-cim*s1*t6
        w0(kz,ig,7)=-cim*s1*t7
      endif
    enddo
  enddo

!
! Compute and print WLM_SUMMARY: ⟨g²⟩ per (l,m) channel (Mode A)
! Mode A is energy-independent (projector-only), print once per k-point
!
  mdim = 2*lb + 1
  ! Only print Mode A once per k-point (at first energy)
  if (ien.eq.1) then
    do m_idx = 1, mdim
      sum_w2   = 0.0_dp
      sum_g2w2 = 0.0_dp
      do kz = 1, nz1
        do ig = 1, ngper
          gmag2 = (gper(1,ig)*tpiba)**2 + (gper(2,ig)*tpiba)**2
          w2    = real(w0(kz,ig,m_idx))**2 + aimag(w0(kz,ig,m_idx))**2
          sum_w2   = sum_w2   + w2
          sum_g2w2 = sum_g2w2 + gmag2 * w2
        enddo
      enddo
      if (sum_w2 > 0.0_dp) then
        g2_avg = sum_g2w2 / sum_w2
      else
        g2_avg = -1.0_dp
      endif
      ! Compute m quantum number from channel index
      m_val = m_idx - 1 - lb
      abs_m = abs(m_val)
      ! Print clearly labeled as energy-independent baseline
      ! Convert to Angstrom units: 1 Bohr = 0.529177 Å, so 1 Bohr^-2 = 3.5711 Å^-2
      write(*,'(A,1x,A,1x,A,2F10.6,1x,A,I2,1x,A,I2,1x,A,ES12.5,1x,A,ES12.5,1x,A)') &
           'WLM_SUMMARY', 'MODE=A:LM(baseline)', 'k1,k2=', xyk(1,ik), xyk(2,ik), &
           'l=', lb, 'm=', m_val, 'g2_bohr=', g2_avg, 'g2_ang=', g2_avg*3.5711_dp, &
           'units:Bohr^-2,Ang^-2 (energy-independent)'
    enddo
  endif

!
! Mode B: State-resolved and state-channel-resolved ⟨g²⟩
!
  call compute_mode_b_g2(w0, nz1, ngper, lb, gper, tpiba, earr(ien), xyk(:,ik))


  deallocate(x1)
  deallocate(x2)
  deallocate(x3)
  deallocate(x4)
  deallocate(x5)
  deallocate(x6)
  deallocate(fx1)
  deallocate(fx2)
  deallocate(fx3)
  deallocate(fx4)
  deallocate(fx5)
  deallocate(fx6)
  deallocate(zsl)
  deallocate(wadd)
  deallocate(wadd2)
  deallocate(wadd3)

  return

CONTAINS

!-----------------------------------------------------------------------
subroutine compute_mode_b_g2(w0, nz1, ngper, lb, gper, tpiba, energy, xyk)
!-----------------------------------------------------------------------
!
! Computes and prints Mode B: state-resolved and state-channel-resolved ⟨g²⟩
! Only executes if CBS eigenvector data is available.
!
! Mode B provides energy- and state-dependent transverse momentum analysis:
!
! MODE=B:STATE - State-resolved ⟨g²⟩:
!   ⟨g²⟩^(n) = Σ_ig g²(ig) |C^(n)(ig)|² / Σ_ig |C^(n)(ig)|²
!   where C^(n)(ig) are the CBS eigenvector components for state n
!
! MODE=B:STATE_LM - State and (l,m)-channel-resolved ⟨g²⟩ (separable approx):
!   ⟨g²⟩^(n)_lm = Σ_ig g²(ig) |C^(n)(ig)|² (Σ_z |w0_lm(z,ig)|²) / norm
!   where w0_lm is the Fourier transform of the beta function for orbital (l,m)
!
! This separable approximation weights each CBS state by its overlap with
! different orbital characters, helping identify which states carry which
! orbital character at different energies.
!
! Output is filtered to keep only physically relevant states (slowest decay).
!
  USE kinds, ONLY: DP
  USE cond, ONLY: nz1_m, ngper_m, nstl_m, nchanl_m, ntot_m, cbs_vec_l, cbs_vec_l_ready, kvall
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nz1, ngper, lb
  REAL(DP), INTENT(IN) :: gper(2,ngper), tpiba, energy, xyk(2)
  COMPLEX(DP), INTENT(IN) :: w0(nz1, ngper, 7)
  !
  ! Output policy parameters
  INTEGER, PARAMETER :: NKEEP = 16          ! Keep slowest-decaying N states
  REAL(DP), PARAMETER :: KAPPA_MAX = 5.0_DP ! Max kappa in Bohr^-1; <=0 to disable
  REAL(DP), PARAMETER :: MIN_WT = 1.0E-8_DP ! Min state weight threshold
  !
  ! State sorting structure
  TYPE state_row
    REAL(DP) :: kappa       ! Decay constant in Bohr^-1
    INTEGER :: idx          ! Original state index
    REAL(DP) :: g2          ! State-resolved g2
    REAL(DP) :: wt          ! State weight (norm)
  END TYPE state_row
  !
  INTEGER :: n, kz, ig, m_idx, mdim, m_val, j, kept, i, ncomp
  REAL(DP) :: gmag2, sum_w2, sum_g2w2, g2_avg_state, g2_avg_lm_state
  REAL(DP) :: c2, w2lm, w2, kappa_j, wt_j, kappa_ang, g2_ang
  REAL(DP), PARAMETER :: bohr_to_ang = 0.529177_DP  ! Bohr to Angstrom conversion
  REAL(DP), PARAMETER :: bohr2_to_ang2 = 3.5711_DP  ! Bohr^-2 to Ang^-2
  COMPLEX(DP) :: c
  LOGICAL :: is_finite
  TYPE(state_row), ALLOCATABLE :: rows(:)
  TYPE(state_row) :: key
  !
  ! Check if CBS eigenvector data is available
  IF (.NOT. cbs_vec_l_ready) THEN
    ! No CBS data available yet - this is expected on first calls
    RETURN
  ENDIF
  !
  ! Mode B.A: State-resolved ⟨g²⟩ with filtering and sorting
  ! First pass: compute g2 and weights for all states
  !
  ALLOCATE(rows(ntot_m))
  !
  ! Determine number of components (currently always 1, but future-proofing for spinors)
  ncomp = SIZE(cbs_vec_l, 1)
  !
  DO n = 1, ntot_m
    sum_w2   = 0.0_DP
    sum_g2w2 = 0.0_DP
    !
    ! Loop over kz (which is 1 for boundary values) and ig
    DO kz = 1, nz1_m
      DO ig = 1, MIN(ngper_m, ngper)
        gmag2 = (gper(1,ig)*tpiba)**2 + (gper(2,ig)*tpiba)**2
        !
        ! Sum over components for spinor/multi-component safety
        c2 = 0.0_DP
        DO i = 1, ncomp
          c = cbs_vec_l(i, ig, n)
          c2 = c2 + REAL(c*CONJG(c), DP)
        ENDDO
        !
        ! Check for finite values
        is_finite = (gmag2 < 1.0E30_DP) .AND. (c2 < 1.0E30_DP) .AND. &
                    (gmag2 > -1.0E30_DP) .AND. (c2 > -1.0E30_DP)
        !
        IF (is_finite) THEN
          sum_w2   = sum_w2   + c2
          sum_g2w2 = sum_g2w2 + gmag2 * c2
        ENDIF
      ENDDO
    ENDDO
    !
    IF (sum_w2 > 0.0_DP) THEN
      g2_avg_state = sum_g2w2 / sum_w2
    ELSE
      g2_avg_state = -1.0_DP
    ENDIF
    !
    ! Compute kappa (decay constant) in Bohr^-1
    IF (n <= SIZE(kvall)) THEN
      kappa_j = ABS(AIMAG(kvall(n))) * tpiba
    ELSE
      kappa_j = 1.0E30_DP  ! Invalid/missing state
    ENDIF
    !
    ! Store state information
    rows(n)%kappa = kappa_j
    rows(n)%idx   = n
    rows(n)%g2    = g2_avg_state
    rows(n)%wt    = sum_w2
  ENDDO
  !
  ! Sort states by kappa (ascending = slowest decay first)
  ! Simple insertion sort
  DO i = 2, ntot_m
    key = rows(i)
    j = i - 1
    DO WHILE (j >= 1)
      IF (rows(j)%kappa <= key%kappa) EXIT
      rows(j+1) = rows(j)
      j = j - 1
    ENDDO
    rows(j+1) = key
  ENDDO
  !
  ! Print only the NKEEP slowest-decaying states that pass filters
  kept = 0
  DO j = 1, ntot_m
    ! Apply filters
    IF (KAPPA_MAX > 0.0_DP) THEN
      IF (rows(j)%kappa > KAPPA_MAX) CYCLE
    ENDIF
    IF (rows(j)%wt < MIN_WT) CYCLE
    IF (rows(j)%g2 < 0.0_DP) CYCLE  ! Skip invalid states
    !
    kept = kept + 1
    IF (kept > NKEEP) EXIT
    !
    ! Convert to Angstrom units
    kappa_ang = rows(j)%kappa / bohr_to_ang  ! Bohr^-1 to Ang^-1
    g2_ang = rows(j)%g2 * bohr2_to_ang2      ! Bohr^-2 to Ang^-2
    !
    ! Print in compact format with kappa, g2, and normalization information
    ! Format: MODE | Energy(eV) | k1,k2 | state_n | kappa(Bohr^-1) | kappa(Ang^-1) | g2(Bohr^-2) | g2(Ang^-2) | norm
    WRITE(*,'(A,1x,A,1x,F8.3,1x,A,2F10.6,1x,A,I4,1x,A,F8.4,1x,A,F8.4,1x,A,ES12.5,1x,A,ES12.5,1x,A,ES10.3)') &
      'WLM_SUMMARY', 'MODE=B:STATE', energy, 'k1,k2=', xyk(1), xyk(2), &
      'n=', rows(j)%idx, 'kappa_bohr=', rows(j)%kappa, 'kappa_ang=', kappa_ang, &
      'g2_bohr=', rows(j)%g2, 'g2_ang=', g2_ang, 'norm=', rows(j)%wt
  ENDDO
  !
  DEALLOCATE(rows)
  !
  ! Mode B.B: State and (l,m)-resolved ⟨g²⟩ (commented out to reduce output)
  ! Uncomment if detailed orbital character analysis is needed
  !
  ! mdim = 2*lb + 1
  ! DO n = 1, MIN(NKEEP, ntot_m)
  !   DO m_idx = 1, mdim
  !     ... (full code available but commented out)
  !   ENDDO
  ! ENDDO
  !
  RETURN

END SUBROUTINE compute_mode_b_g2

end subroutine four
!
!-----------------------------------------------------------------------
function indexr(zz, ndim, r)
  USE kinds, only : DP
  implicit none

  integer :: iz, ndim, indexr
  real(DP) :: zz, r(ndim)
!
!     abs(zz)<r(indexr)
!
  iz = 1
  do while(r(iz).le.abs(zz)+1.d-10)
    iz=iz+1
  enddo
  indexr=iz
  return
end function indexr
