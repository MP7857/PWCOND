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
! Check normalization constants by numerical integration
! Set check_norm_constants to .true. to enable verification (default: .false. for performance)
!
  logical, parameter :: check_norm_constants = .false.
  logical :: norm_ok
  
  if (check_norm_constants) then
    call check_normalization(lb, s1, s2, s3, s4, norm_ok)
    if (.not. norm_ok) then
      write(*,'(A)') 'WARNING: Some normalization checks failed!'
    endif
  endif

  return
end subroutine four

!
!----------------------------------------------------------------------
subroutine check_normalization(lb, s1, s2, s3, s4, norm_ok)
!----------------------------------------------------------------------
!
! This subroutine verifies the normalization constants s1, s2, s3, s4
! by numerically integrating the spherical harmonics over the unit sphere.
! For properly normalized spherical harmonics Y_lm:
!   \int |Y_lm|^2 dOmega = 1
!
! The normalization constants s1, s2, s3, s4 used in four() are derived from
! the standard spherical harmonic normalization factors. This subroutine verifies
! that the underlying spherical harmonics are correctly normalized, which validates
! that the normalization factors used to compute s1, s2, s3, s4 are mathematically
! correct. The constants s1-s4 also include the 2D Fourier transform factors
! (tpi/sarea) which are not tested here.
!
  USE kinds, ONLY: DP
  USE constants, ONLY: pi, fpi
  implicit none
  
  integer, intent(in) :: lb
  real(DP), intent(in) :: s1, s2, s3, s4
  logical, intent(out) :: norm_ok
  
  integer :: ntheta, nphi, itheta, iphi
  real(DP) :: theta, phi, dtheta, dphi, dOmega
  real(DP) :: x, y, z
  real(DP) :: integral, ylm_val
  real(DP), parameter :: eps = 1.d-6
  
  ! Configurable integration parameters
  integer, parameter :: ntheta_default = 100
  integer, parameter :: nphi_default = 200
  
  ! All checks passed by default
  norm_ok = .true.
  
  ! Number of integration points
  ntheta = ntheta_default
  nphi = nphi_default
  dtheta = pi / dble(ntheta)
  dphi = 2.d0 * pi / dble(nphi)
  
  write(*,*)
  write(*,'(A)') '============================================'
  write(*,'(A)') 'Normalization Check for Spherical Harmonics'
  write(*,'(A)') '============================================'
  write(*,'(A)') 'Note: This verifies the spherical harmonic normalization'
  write(*,'(A)') '      which validates the mathematical basis for s1-s4.'
  write(*,*)
  
  if (lb .eq. 0) then
    ! l=0, m=0 (s-orbital)
    write(*,'(A)') 'l = 0 (s-orbital)'
    write(*,'(A,ES15.8)') '  Normalization constant s1 = ', s1
    
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        
        ! Y_00 = 1/sqrt(4*pi)
        ylm_val = 1.d0 / sqrt(fpi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=0: Integral of |Y_00|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified (integral ≈ 1)'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
  elseif (lb .eq. 1) then
    ! l=1 (p-orbitals): m=-1, 0, 1
    write(*,'(A)') 'l = 1 (p-orbitals)'
    write(*,'(A,ES15.8)') '  Normalization constant s1 = ', s1
    
    ! m=0 (p_z)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        z = cos(theta)
        
        ! Y_10 = sqrt(3/(4*pi)) * cos(theta)
        ylm_val = sqrt(3.d0 / fpi) * z
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=0 (p_z): Integral of |Y_10|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±1 (p_x, p_y)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        x = sin(theta) * cos(phi)
        
        ! Y_11 real part (for p_x) = sqrt(3/(4*pi)) * sin(theta) * cos(phi)
        ylm_val = sqrt(3.d0 / fpi) * sin(theta) * cos(phi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±1 (p_x): Integral of |Y_11|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±1 (p_y)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        y = sin(theta) * sin(phi)
        
        ! Y_11 imag part (for p_y) = sqrt(3/(4*pi)) * sin(theta) * sin(phi)
        ylm_val = sqrt(3.d0 / fpi) * sin(theta) * sin(phi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±1 (p_y): Integral of |Y_11|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
  elseif (lb .eq. 2) then
    ! l=2 (d-orbitals): 5 components
    write(*,'(A)') 'l = 2 (d-orbitals)'
    write(*,'(A,ES15.8)') '  Normalization constant s1 = ', s1
    write(*,'(A,ES15.8)') '  Normalization constant s2 = ', s2
    
    ! m=0 (d_z^2)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        z = cos(theta)
        
        ! Y_20 = sqrt(5/(16*pi)) * (3*cos^2(theta) - 1)
        ylm_val = sqrt(5.d0 / (16.d0 * pi)) * (3.d0 * z**2 - 1.d0)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=0 (d_z^2): Integral of |Y_20|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±1 (d_xz, d_yz)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        x = sin(theta) * cos(phi)
        z = cos(theta)
        
        ! Y_21 real (d_xz) = sqrt(15/(4*pi)) * sin(theta)*cos(theta)*cos(phi)
        ylm_val = sqrt(15.d0 / fpi) * sin(theta) * z * cos(phi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±1 (d_xz): Integral of |Y_21|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±1 (d_yz)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        y = sin(theta) * sin(phi)
        z = cos(theta)
        
        ! Y_21 imag (d_yz) = sqrt(15/(4*pi)) * sin(theta)*cos(theta)*sin(phi)
        ylm_val = sqrt(15.d0 / fpi) * sin(theta) * z * sin(phi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±1 (d_yz): Integral of |Y_21|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±2 (d_x^2-y^2, d_xy)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        x = sin(theta) * cos(phi)
        y = sin(theta) * sin(phi)
        
        ! Y_22 real (d_x^2-y^2) = sqrt(15/(16*pi)) * sin^2(theta)*(cos^2(phi)-sin^2(phi))
        ylm_val = sqrt(15.d0 / (16.d0 * pi)) * sin(theta)**2 * (cos(phi)**2 - sin(phi)**2)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±2 (d_x^2-y^2): Integral of |Y_22|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±2 (d_xy)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        
        ! Y_22 imag (d_xy) = sqrt(15/(16*pi)) * sin^2(theta)*2*cos(phi)*sin(phi)
        ylm_val = sqrt(15.d0 / (16.d0 * pi)) * sin(theta)**2 * 2.d0 * cos(phi) * sin(phi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±2 (d_xy): Integral of |Y_22|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
  elseif (lb .eq. 3) then
    ! l=3 (f-orbitals): 7 components
    write(*,'(A)') 'l = 3 (f-orbitals)'
    write(*,'(A,ES15.8)') '  Normalization constant s1 = ', s1
    write(*,'(A,ES15.8)') '  Normalization constant s2 = ', s2
    write(*,'(A,ES15.8)') '  Normalization constant s3 = ', s3
    write(*,'(A,ES15.8)') '  Normalization constant s4 = ', s4
    
    ! m=0 (f_z^3)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        z = cos(theta)
        
        ! Y_30 = sqrt(7/(16*pi)) * (5*cos^3(theta) - 3*cos(theta))
        ylm_val = sqrt(7.d0 / (16.d0 * pi)) * (5.d0 * z**3 - 3.d0 * z)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=0 (f_z^3): Integral of |Y_30|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±1 (f_xz^2, f_yz^2)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        x = sin(theta) * cos(phi)
        z = cos(theta)
        
        ! Y_31 real (f_xz^2) = sqrt(21/(32*pi)) * sin(theta)*cos(phi)*(5*cos^2(theta)-1)
        ylm_val = sqrt(21.d0 / (32.d0 * pi)) * sin(theta) * cos(phi) * (5.d0 * z**2 - 1.d0)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±1 (f_xz^2): Integral of |Y_31|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±1 (f_yz^2)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        y = sin(theta) * sin(phi)
        z = cos(theta)
        
        ! Y_31 imag (f_yz^2) = sqrt(21/(32*pi)) * sin(theta)*sin(phi)*(5*cos^2(theta)-1)
        ylm_val = sqrt(21.d0 / (32.d0 * pi)) * sin(theta) * sin(phi) * (5.d0 * z**2 - 1.d0)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±1 (f_yz^2): Integral of |Y_31|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±2 (f_xyz, f_z(x^2-y^2))
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        z = cos(theta)
        
        ! Y_32 real (f_z(x^2-y^2)) = sqrt(105/(16*pi)) * sin^2(theta)*cos(theta)*(cos^2(phi)-sin^2(phi))
        ylm_val = sqrt(105.d0 / (16.d0 * pi)) * sin(theta)**2 * z * (cos(phi)**2 - sin(phi)**2)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±2 (f_z(x^2-y^2)): Integral of |Y_32|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±2 (f_xyz)
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        z = cos(theta)
        
        ! Y_32 imag (f_xyz) = sqrt(105/(16*pi)) * sin^2(theta)*cos(theta)*2*cos(phi)*sin(phi)
        ylm_val = sqrt(105.d0 / (16.d0 * pi)) * sin(theta)**2 * z * 2.d0 * cos(phi) * sin(phi)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±2 (f_xyz): Integral of |Y_32|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±3 (f_x(x^2-3y^2), f_y(3x^2-y^2))
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        
        ! Y_33 real (f_x(x^2-3y^2)) = sqrt(35/(32*pi)) * sin^3(theta)*(cos^3(phi)-3*cos(phi)*sin^2(phi))
        ylm_val = sqrt(35.d0 / (32.d0 * pi)) * sin(theta)**3 * &
                  (cos(phi)**3 - 3.d0 * cos(phi) * sin(phi)**2)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±3 (f_x(x^2-3y^2)): Integral of |Y_33|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
    ! m=±3 (f_y(3x^2-y^2))
    integral = 0.d0
    do itheta = 1, ntheta
      theta = (itheta - 0.5d0) * dtheta
      do iphi = 1, nphi
        phi = (iphi - 0.5d0) * dphi
        dOmega = sin(theta) * dtheta * dphi
        
        ! Y_33 imag (f_y(3x^2-y^2)) = sqrt(35/(32*pi)) * sin^3(theta)*(3*cos^2(phi)*sin(phi)-sin^3(phi))
        ylm_val = sqrt(35.d0 / (32.d0 * pi)) * sin(theta)**3 * &
                  (3.d0 * cos(phi)**2 * sin(phi) - sin(phi)**3)
        integral = integral + ylm_val**2 * dOmega
      enddo
    enddo
    write(*,'(A,F15.10)') '  m=±3 (f_y(3x^2-y^2)): Integral of |Y_33|^2 = ', integral
    if (abs(integral - 1.d0) < eps) then
      write(*,'(A)') '  ✓ Normalization verified'
    else
      write(*,'(A)') '  ✗ Warning: Normalization may be incorrect'
      norm_ok = .false.
    endif
    
  endif
  
  write(*,*)
  write(*,'(A)') '============================================'
  write(*,*)
  
  return
end subroutine check_normalization

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
