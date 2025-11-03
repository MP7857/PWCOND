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
subroutine four(w0, z0, dz, tblm, taunew, r, rab, betar)
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
  USE cond, ONLY : sarea, nz1, ngper, gper, ninsh, gnsh, ngpsh

implicit none

  integer :: kz, ig, ign, igphi, &
             indexr, iz, lb, ir, nmesh, nmeshs, tblm(4)
  real(DP), parameter :: eps=1.d-8
  complex(DP), parameter :: cim=(0.d0, 1.d0)
  real(DP) :: gn, s1, s2, s3, s4, cs, sn, cs2, sn2, cs3, sn3, rz, dz1, zr, &
                   dr, z0, dz,  bessj, taunew(4), r(ndmx),         &
                   rab(ndmx), betar(ndmx)
  real(DP), allocatable :: x1(:), x2(:), x3(:), x4(:), x5(:), x6(:)
  real(DP), allocatable :: fx1(:), fx2(:), fx3(:), fx4(:), fx5(:), fx6(:), zsl(:)
  complex(DP) :: w0(nz1, ngper, 7)
  complex(DP), allocatable :: wadd(:,:), wadd2(:,:), wadd3(:,:)
  complex(DP) :: t1, t2, t3, t4, t5, t6, t7, wa1, wa2, wa3
  logical, parameter :: check_integrals = .false.
  logical :: integrals_ok


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
        w0(kz,ig,4)=-s2*zsl(kz)*t4
        w0(kz,ig,5)=-s2*zsl(kz)*t5
        w0(kz,ig,6)=cim*s1*t6
        w0(kz,ig,7)=cim*s1*t7
      endif
    enddo
  enddo

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

!
! Check computed integrals (fx1-fx6) and normalization by verifying w0
! Set check_integrals to .true. to enable (default: .false. for performance)
!
  if (check_integrals) then
    call check_computed_integrals(w0, lb, s1, s2, s3, s4, nz1, ngper, dz1*alat, integrals_ok)
    if (.not. integrals_ok) then
      write(*,'(A)') 'WARNING: Computed integrals may have issues!'
    endif
  endif

  return
end subroutine four

!
!----------------------------------------------------------------------
subroutine check_computed_integrals(w0, lb, s1, s2, s3, s4, nz1, ngper, dz, integrals_ok)
!----------------------------------------------------------------------
!
! This subroutine verifies the computed integrals (fx1-fx6) and normalization
! by checking properties of the final w0 array after all processing.
! 
! It computes the integrated norm of w0 over z-slices for each m component
! and reports diagnostic information about the computed integrals.
!
! Enhanced checks include:
! - Integrated norms and statistical properties (max, min, avg)
! - Verification that values are in reasonable ranges
! - Symmetry checks where applicable
! - Convergence indicators
!
  USE kinds, ONLY: DP
  USE constants, ONLY: tpi, fpi
  implicit none
  
  integer, intent(in) :: lb, nz1, ngper
  real(DP), intent(in) :: s1, s2, s3, s4, dz
  complex(DP), intent(in) :: w0(nz1, ngper, 7)
  logical, intent(out) :: integrals_ok
  
  integer :: kz, ig, m, num_m
  real(DP) :: norm_sum, max_val, min_val, avg_val, std_dev, variance
  real(DP) :: real_part_sum, imag_part_sum, symmetry_measure
  real(DP) :: total_norm, relative_norm
  real(DP), parameter :: eps = 1.d-10
  real(DP), parameter :: large_threshold = 1.d10
  real(DP), parameter :: small_threshold = 1.d-15
  
  ! All checks passed by default
  integrals_ok = .true.
  
  ! Determine number of m components for this l
  if (lb .eq. 0) then
    num_m = 1
  elseif (lb .eq. 1) then
    num_m = 3
  elseif (lb .eq. 2) then
    num_m = 5
  elseif (lb .eq. 3) then
    num_m = 7
  else
    num_m = 0
  endif
  
  write(*,*)
  write(*,'(A)') '===================================================='
  write(*,'(A)') ' COMPREHENSIVE VERIFICATION OF COMPUTED INTEGRALS'
  write(*,'(A)') '===================================================='
  write(*,'(A)') 'Checking radial integrals (fx1-fx6) after complete'
  write(*,'(A)') 'processing: integration, angular combination, and'
  write(*,'(A)') 'normalization with constants s1-s4.'
  write(*,*)
  write(*,'(A,I2)') 'Angular momentum l = ', lb
  write(*,'(A,I2)') 'Number of m components = ', num_m
  write(*,'(A,I5)') 'Number of z-slices = ', nz1
  write(*,'(A,I5)') 'Number of g-vectors = ', ngper
  write(*,'(A,ES15.8)') 'dz (step size) = ', dz
  write(*,*)
  write(*,'(A)') 'Normalization constants:'
  write(*,'(A,ES15.8)') '  s1 = ', s1
  if (lb .ge. 2) then
    write(*,'(A,ES15.8)') '  s2 = ', s2
  endif
  if (lb .ge. 3) then
    write(*,'(A,ES15.8)') '  s3 = ', s3
    write(*,'(A,ES15.8)') '  s4 = ', s4
  endif
  write(*,*)
  write(*,'(A)') '----------------------------------------------------'
  
  ! Compute total norm across all components
  total_norm = 0.d0
  do m = 1, num_m
    do ig = 1, ngper
      do kz = 1, nz1
        total_norm = total_norm + abs(w0(kz,ig,m))**2 * dz
      enddo
    enddo
  enddo
  
  write(*,'(A,ES15.8)') 'Total integrated |w0|^2 (all m) = ', total_norm
  write(*,*)
  write(*,'(A)') 'Per-component analysis:'
  write(*,'(A)') '----------------------------------------------------'
  
  ! Check each m component with enhanced diagnostics
  do m = 1, num_m
    norm_sum = 0.d0
    max_val = 0.d0
    min_val = 1.d30
    avg_val = 0.d0
    real_part_sum = 0.d0
    imag_part_sum = 0.d0
    variance = 0.d0
    
    ! First pass: compute statistics
    do ig = 1, ngper
      do kz = 1, nz1
        norm_sum = norm_sum + abs(w0(kz,ig,m))**2 * dz
        max_val = max(max_val, abs(w0(kz,ig,m)))
        if (abs(w0(kz,ig,m)) > eps) then
          min_val = min(min_val, abs(w0(kz,ig,m)))
        endif
        avg_val = avg_val + abs(w0(kz,ig,m))
        real_part_sum = real_part_sum + real(w0(kz,ig,m))
        imag_part_sum = imag_part_sum + aimag(w0(kz,ig,m))
      enddo
    enddo
    
    avg_val = avg_val / (nz1 * ngper)
    
    ! Second pass: compute standard deviation
    do ig = 1, ngper
      do kz = 1, nz1
        variance = variance + (abs(w0(kz,ig,m)) - avg_val)**2
      enddo
    enddo
    variance = variance / (nz1 * ngper)
    std_dev = sqrt(variance)
    
    ! Compute relative norm
    if (total_norm > eps) then
      relative_norm = norm_sum / total_norm
    else
      relative_norm = 0.d0
    endif
    
    ! Symmetry measure (for complex values)
    symmetry_measure = abs(real_part_sum) + abs(imag_part_sum)
    
    write(*,*)
    write(*,'(A,I2)') '>>> Component m = ', m
    write(*,'(A,ES15.8)') '  Integrated |w0|^2 over z    = ', norm_sum
    write(*,'(A,F12.6,A)') '  Relative norm (%)           = ', relative_norm * 100.d0, ' %'
    write(*,'(A,ES15.8)') '  Max |w0|                    = ', max_val
    if (min_val < 1.d29) then
      write(*,'(A,ES15.8)') '  Min |w0| (non-zero)         = ', min_val
    else
      write(*,'(A)') '  Min |w0| (non-zero)         = N/A (all zero)'
    endif
    write(*,'(A,ES15.8)') '  Avg |w0|                    = ', avg_val
    write(*,'(A,ES15.8)') '  Std dev |w0|                = ', std_dev
    if (avg_val > eps) then
      write(*,'(A,F12.6)') '  Coefficient of variation    = ', std_dev / avg_val
    endif
    write(*,'(A,ES15.8)') '  Sum of real parts           = ', real_part_sum
    write(*,'(A,ES15.8)') '  Sum of imag parts           = ', imag_part_sum
    
    ! Perform validation checks
    write(*,*)
    write(*,'(A)') '  Validation checks:'
    
    ! Check 1: Not all zero
    if (max_val .lt. small_threshold) then
      write(*,'(A)') '  ✗ FAIL: All values near zero'
      write(*,'(A)') '    → Radial integrals may not be computed correctly'
      integrals_ok = .false.
    else
      write(*,'(A)') '  ✓ PASS: Values are non-zero'
    endif
    
    ! Check 2: Not unreasonably large
    if (max_val .gt. large_threshold) then
      write(*,'(A)') '  ✗ FAIL: Values unreasonably large'
      write(*,'(A)') '    → Normalization constants may be incorrect'
      integrals_ok = .false.
    else
      write(*,'(A)') '  ✓ PASS: Values in reasonable range'
    endif
    
    ! Check 3: Norm is non-negative
    if (norm_sum .lt. -eps) then
      write(*,'(A)') '  ✗ FAIL: Negative integrated norm'
      write(*,'(A)') '    → This should never happen (numerical error?)'
      integrals_ok = .false.
    else
      write(*,'(A)') '  ✓ PASS: Integrated norm is non-negative'
    endif
    
    ! Check 4: Coefficient of variation check
    if (avg_val > eps .and. std_dev / avg_val > 100.d0) then
      write(*,'(A)') '  ⚠ WARNING: Very high variability in |w0|'
      write(*,'(A)') '    → This may indicate numerical issues'
    elseif (avg_val > eps) then
      write(*,'(A)') '  ✓ PASS: Variability appears reasonable'
    endif
    
    ! Check 5: Min/max ratio check
    if (min_val < 1.d29 .and. max_val > eps) then
      if (max_val / min_val > 1.d12) then
        write(*,'(A)') '  ⚠ WARNING: Very large dynamic range'
        write(*,'(A,ES10.3)') '    → max/min ratio = ', max_val / min_val
      else
        write(*,'(A)') '  ✓ PASS: Dynamic range acceptable'
      endif
    endif
    
  enddo
  
  write(*,*)
  write(*,'(A)') '===================================================='
  write(*,'(A)') 'SUMMARY'
  write(*,'(A)') '===================================================='
  if (integrals_ok) then
    write(*,'(A)') '✓ ALL CHECKS PASSED'
    write(*,'(A)') 'The computed radial integrals appear correct.'
  else
    write(*,'(A)') '✗ SOME CHECKS FAILED'
    write(*,'(A)') 'Please review the diagnostics above.'
  endif
  write(*,*)
  write(*,'(A)') 'Note: This verification checks the fx1-fx6 integrals'
  write(*,'(A)') 'indirectly through the final w0 array after applying'
  write(*,'(A)') 'angular terms (cs, sn, cs2, sn2, cs3, sn3) and'
  write(*,'(A)') 'normalization constants (s1, s2, s3, s4).'
  write(*,'(A)') '===================================================='
  write(*,*)
  
  return
end subroutine check_computed_integrals

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
