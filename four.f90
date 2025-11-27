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
! Refined segment-based quadrature for l=3 added (2025):
!   - Introduces 1/r^3 regularization via eps_r_f for numerical robustness.
!   - Segment-based integration with adaptive Gauss-Legendre for high ewind.
!   - Controlled by refine_f_quadrature parameter (default .false.).
!   - See four.tex for analytic derivation.
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
!             f_{z(5z^2-3r^2)}, f_{x(5z^2-r^2)}, f_{-y(5z^2-r^2)},
!             f_{z(x^2-y^2)}, f_{-xyz}, f_{x(x^2-3y^2)}, f_{-y(3x^2-y^2)}
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
  !
  ! Parameters for refined f-channel (l=3) quadrature:
  !   refine_f_quadrature: .true. enables segment-based integration for l=3
  !   eps_r_f: regularization epsilon for 1/r^3 singularity (Bohr units)
  !            replaces 1/r^3 by (r^2 + eps_r_f^2)^(-3/2)
  !   dphi_smooth_f: Bessel phase threshold; segments with dphi <= this value
  !                  use simple trapezoid, else use 4-point Gauss-Legendre
  !   (See four.tex for derivation and rationale.)
  !
  logical, parameter :: refine_f_quadrature = .false.
  real(DP), parameter :: eps_r_f = 1.d-10
  real(DP), parameter :: dphi_smooth_f = 0.5d0
  !
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
            !
            ! f-channel (lb=3) integration:
            ! When refine_f_quadrature = .true., use segment-based integration
            ! with adaptive Gauss-Legendre for improved numerical robustness at
            ! high ewind. Otherwise, use the original Simpson integration.
            !
            if (refine_f_quadrature) then
               !
               ! Segment-based integration for f-channel (l=3).
               ! Starts integration at r >= |z| + eps_r_f for regularization.
               ! Dispatches to trapezoid or Gauss-Legendre based on Bessel phase.
               !
               call integrate_f_segments(gn, zsl(kz), r, nmeshs, iz, betar, &
                                         fx1(kz), fx2(kz), fx3(kz), fx4(kz), &
                                         fx5(kz), fx6(kz))
               !
               ! Add endpoint corrections to match original Simpson behavior.
               ! These handle the partial integral near r ≈ |z| that the segment-based
               ! integration skips. Critical for numerical consistency with original code.
               !
               fx1(kz)=fx1(kz)+x1(iz)*0.5d0*zr
               fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
               fx3(kz)=fx3(kz)+x3(iz)*0.5d0*zr
               fx4(kz)=fx4(kz)+x4(iz)*0.5d0*zr
               if(iz.eq.1) then
                  x5(iz-1)=0.d0
               else
                  x5(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr) / &
                           (abs(zsl(kz))**3 + eps_r_f)
               endif
               x6(iz-1)=0.d0
               fx5(kz)=fx5(kz)+(x5(iz-1)+x5(iz))*0.5d0*zr
               fx6(kz)=fx6(kz)+(x6(iz-1)+x6(iz))*0.5d0*zr
            else
               !
               ! Original Simpson-based integration for f-channel (lb=3)
               ! with regularized 1/r^3 singularity at small zsl.
               !
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
                  ! Use eps_r_f for regularization of 1/r^3 singularity at small zsl
                  x5(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr) / &
                           (abs(zsl(kz))**3 + eps_r_f)
               endif
               x6(iz-1)=0.d0
               fx5(kz)=fx5(kz)+(x5(iz-1)+x5(iz))*0.5d0*zr
               fx6(kz)=fx6(kz)+(x6(iz-1)+x6(iz))*0.5d0*zr
            endif
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
        w0(kz,ig,2)=cim*s3*(4.d0*zsl(kz)**2*t2-wa2)
        w0(kz,ig,3)=cim*s3*(4.d0*zsl(kz)**2*t3-wa3)
        w0(kz,ig,4)=-s2*zsl(kz)*t4
        w0(kz,ig,5)=-s2*zsl(kz)*t5
        w0(kz,ig,6)=-cim*s1*t6
        w0(kz,ig,7)=-cim*s1*t7
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

  return

contains

  !--------------------------------------------------------------------
  ! integrate_f_segments: Drives segment-based integration for f-channel (l=3)
  !
  ! For each radial segment [r(i), r(i+1)] starting at r >= |z| + eps_r_f,
  ! computes the Bessel phase increment and dispatches to either simple
  ! trapezoid (slow phase variation) or 4-point Gauss-Legendre (fast phase).
  !
  ! Uses regularized 1/r^3 denominators: (r^2 + eps_r_f^2)^(-3/2)
  ! See four.tex for analytic derivation.
  !--------------------------------------------------------------------
  subroutine integrate_f_segments(gn_val, z_val, r_arr, nmeshs_val, iz_val, &
                                  betar_arr, fx1_out, fx2_out, fx3_out, &
                                  fx4_out, fx5_out, fx6_out)
    implicit none
    real(DP), intent(in) :: gn_val, z_val
    real(DP), intent(in) :: r_arr(ndmx), betar_arr(ndmx)
    integer, intent(in) :: nmeshs_val, iz_val
    real(DP), intent(out) :: fx1_out, fx2_out, fx3_out, fx4_out, fx5_out, fx6_out

    integer :: i_seg, i_start
    real(DP) :: z_abs, r_lower
    real(DP) :: seg_fx1, seg_fx2, seg_fx3, seg_fx4, seg_fx5, seg_fx6

    ! Initialize output accumulators
    fx1_out = 0.d0
    fx2_out = 0.d0
    fx3_out = 0.d0
    fx4_out = 0.d0
    fx5_out = 0.d0
    fx6_out = 0.d0

    z_abs = abs(z_val)
    r_lower = z_abs + eps_r_f

    ! Find first grid point such that r(i) >= r_lower
    i_start = iz_val
    do while (i_start <= nmeshs_val .and. r_arr(i_start) < r_lower)
       i_start = i_start + 1
    enddo

    ! If no valid starting point, return zero integrals
    if (i_start > nmeshs_val) return

    ! Loop over segments [r(i), r(i+1)]
    do i_seg = i_start, nmeshs_val - 1
       call integrate_segment_f(i_seg, i_seg+1, z_val, gn_val, r_arr, betar_arr, &
                                seg_fx1, seg_fx2, seg_fx3, seg_fx4, seg_fx5, seg_fx6)
       fx1_out = fx1_out + seg_fx1
       fx2_out = fx2_out + seg_fx2
       fx3_out = fx3_out + seg_fx3
       fx4_out = fx4_out + seg_fx4
       fx5_out = fx5_out + seg_fx5
       fx6_out = fx6_out + seg_fx6
    enddo

  end subroutine integrate_f_segments

  !--------------------------------------------------------------------
  ! integrate_segment_f: Computes integral over one segment [r(i1), r(i2)]
  !
  ! Computes Bessel phase increment dphi = gn * |rz2 - rz1|.
  ! If dphi <= dphi_smooth_f, uses simple trapezoid rule.
  ! Otherwise, uses 4-point Gauss-Legendre quadrature.
  !--------------------------------------------------------------------
  subroutine integrate_segment_f(i1, i2, z_val, gn_val, r_arr, betar_arr, &
                                 fx1_seg, fx2_seg, fx3_seg, fx4_seg, &
                                 fx5_seg, fx6_seg)
    implicit none
    integer, intent(in) :: i1, i2
    real(DP), intent(in) :: z_val, gn_val
    real(DP), intent(in) :: r_arr(ndmx), betar_arr(ndmx)
    real(DP), intent(out) :: fx1_seg, fx2_seg, fx3_seg, fx4_seg, fx5_seg, fx6_seg

    real(DP) :: r1, r2, z_sq, rz1, rz2, dphi

    r1 = r_arr(i1)
    r2 = r_arr(i2)
    z_sq = z_val**2

    ! Compute rz at segment endpoints with regularization
    rz1 = sqrt(max(r1**2 - z_sq, 0.d0))
    rz2 = sqrt(max(r2**2 - z_sq, 0.d0))

    ! Compute Bessel phase increment
    dphi = gn_val * abs(rz2 - rz1)

    if (dphi <= dphi_smooth_f) then
       ! Use simple trapezoid rule for smooth (slow phase) segments
       call simple_segment_f(r1, r2, z_val, gn_val, betar_arr(i1), betar_arr(i2), &
                             fx1_seg, fx2_seg, fx3_seg, fx4_seg, fx5_seg, fx6_seg)
    else
       ! Use 4-point Gauss-Legendre for rapidly oscillating segments
       call gauss4_segment_f(r1, r2, z_val, gn_val, betar_arr(i1), betar_arr(i2), &
                             fx1_seg, fx2_seg, fx3_seg, fx4_seg, fx5_seg, fx6_seg)
    endif

  end subroutine integrate_segment_f

  !--------------------------------------------------------------------
  ! simple_segment_f: Simple trapezoid rule for segment [r1, r2]
  !
  ! Evaluates integrands at endpoints and uses trapezoidal approximation:
  !   integral ≈ 0.5 * (f(r1) + f(r2)) * (r2 - r1)
  !
  ! Uses regularized denominator: (r^2 + eps_r_f^2)^(-3/2) instead of 1/r^3
  !--------------------------------------------------------------------
  subroutine simple_segment_f(r1, r2, z_val, gn_val, betar1, betar2, &
                              fx1_seg, fx2_seg, fx3_seg, fx4_seg, &
                              fx5_seg, fx6_seg)
    implicit none
    real(DP), intent(in) :: r1, r2, z_val, gn_val, betar1, betar2
    real(DP), intent(out) :: fx1_seg, fx2_seg, fx3_seg, fx4_seg, fx5_seg, fx6_seg

    real(DP) :: z_sq, dr
    real(DP) :: rz1, rz2, denom1, denom2
    real(DP) :: f1_x1, f1_x2, f1_x3, f1_x4, f1_x5, f1_x6
    real(DP) :: f2_x1, f2_x2, f2_x3, f2_x4, f2_x5, f2_x6

    z_sq = z_val**2
    dr = r2 - r1

    ! Endpoint 1: r1
    rz1 = sqrt(max(r1**2 - z_sq, 0.d0))
    ! Regularized denominator: (r1^2 + eps_r_f^2)^(-3/2)
    denom1 = (r1**2 + eps_r_f**2)**1.5d0

    ! Integrands at r1 matching four.tex definitions for lb=3:
    !   x1 ~ betar * J3 * rz^3 / r^3
    !   x2 ~ betar * J2 * rz^2 / r^3
    !   x3 ~ betar * J1 * rz   / r^3
    !   x4 ~ betar * J1 * rz^3 / r^3
    !   x5 ~ betar * J0        / r^3
    !   x6 ~ betar * J0 * rz^2 / r^3
    f1_x1 = betar1 * bessj(3, gn_val*rz1) * rz1**3 / denom1
    f1_x2 = betar1 * bessj(2, gn_val*rz1) * rz1**2 / denom1
    f1_x3 = betar1 * bessj(1, gn_val*rz1) * rz1    / denom1
    f1_x4 = betar1 * bessj(1, gn_val*rz1) * rz1**3 / denom1
    f1_x5 = betar1 * bessj(0, gn_val*rz1)          / denom1
    f1_x6 = betar1 * bessj(0, gn_val*rz1) * rz1**2 / denom1

    ! Endpoint 2: r2
    rz2 = sqrt(max(r2**2 - z_sq, 0.d0))
    denom2 = (r2**2 + eps_r_f**2)**1.5d0

    f2_x1 = betar2 * bessj(3, gn_val*rz2) * rz2**3 / denom2
    f2_x2 = betar2 * bessj(2, gn_val*rz2) * rz2**2 / denom2
    f2_x3 = betar2 * bessj(1, gn_val*rz2) * rz2    / denom2
    f2_x4 = betar2 * bessj(1, gn_val*rz2) * rz2**3 / denom2
    f2_x5 = betar2 * bessj(0, gn_val*rz2)          / denom2
    f2_x6 = betar2 * bessj(0, gn_val*rz2) * rz2**2 / denom2

    ! Trapezoid rule
    fx1_seg = 0.5d0 * (f1_x1 + f2_x1) * dr
    fx2_seg = 0.5d0 * (f1_x2 + f2_x2) * dr
    fx3_seg = 0.5d0 * (f1_x3 + f2_x3) * dr
    fx4_seg = 0.5d0 * (f1_x4 + f2_x4) * dr
    fx5_seg = 0.5d0 * (f1_x5 + f2_x5) * dr
    fx6_seg = 0.5d0 * (f1_x6 + f2_x6) * dr

  end subroutine simple_segment_f

  !--------------------------------------------------------------------
  ! gauss4_segment_f: 4-point Gauss-Legendre quadrature for segment [r1, r2]
  !
  ! Uses standard 4-point Gauss-Legendre nodes and weights on [-1,1],
  ! mapped to [r1, r2]. Linearly interpolates betar between betar1 and betar2.
  !
  ! Uses regularized denominator: (r^2 + eps_r_f^2)^(-3/2) instead of 1/r^3
  !--------------------------------------------------------------------
  subroutine gauss4_segment_f(r1, r2, z_val, gn_val, betar1, betar2, &
                              fx1_seg, fx2_seg, fx3_seg, fx4_seg, &
                              fx5_seg, fx6_seg)
    implicit none
    real(DP), intent(in) :: r1, r2, z_val, gn_val, betar1, betar2
    real(DP), intent(out) :: fx1_seg, fx2_seg, fx3_seg, fx4_seg, fx5_seg, fx6_seg

    ! 4-point Gauss-Legendre nodes and weights on [-1, 1]
    ! Reference: Abramowitz & Stegun, Table 25.4 (n=4)
    real(DP), parameter :: xi(4) = (/ -0.861136311594053d0, &
                                      -0.339981043584856d0, &
                                       0.339981043584856d0, &
                                       0.861136311594053d0 /)
    real(DP), parameter :: wt(4) = (/  0.347854845137454d0, &
                                       0.652145154862546d0, &
                                       0.652145154862546d0, &
                                       0.347854845137454d0 /)

    integer :: k
    real(DP) :: z_sq, dr, r_mid, half_dr
    real(DP) :: r_sub, t_interp, betar_sub, rz_sub, denom_sub
    real(DP) :: f_x1, f_x2, f_x3, f_x4, f_x5, f_x6

    z_sq = z_val**2
    dr = r2 - r1
    r_mid = 0.5d0 * (r1 + r2)
    half_dr = 0.5d0 * dr

    ! Initialize accumulators
    fx1_seg = 0.d0
    fx2_seg = 0.d0
    fx3_seg = 0.d0
    fx4_seg = 0.d0
    fx5_seg = 0.d0
    fx6_seg = 0.d0

    ! Guard against zero-width segments
    if (abs(dr) < 1.d-15) return

    ! Loop over 4 Gauss-Legendre points
    do k = 1, 4
       ! Map xi(k) from [-1,1] to r in [r1, r2]
       r_sub = r_mid + half_dr * xi(k)

       ! Linear interpolation of betar: t_interp in [0,1]
       t_interp = (r_sub - r1) / dr
       betar_sub = betar1 + t_interp * (betar2 - betar1)

       ! Compute rz at this sub-point
       rz_sub = sqrt(max(r_sub**2 - z_sq, 0.d0))

       ! Regularized denominator
       denom_sub = (r_sub**2 + eps_r_f**2)**1.5d0

       ! Evaluate integrands (same as four.tex definitions for lb=3)
       f_x1 = betar_sub * bessj(3, gn_val*rz_sub) * rz_sub**3 / denom_sub
       f_x2 = betar_sub * bessj(2, gn_val*rz_sub) * rz_sub**2 / denom_sub
       f_x3 = betar_sub * bessj(1, gn_val*rz_sub) * rz_sub    / denom_sub
       f_x4 = betar_sub * bessj(1, gn_val*rz_sub) * rz_sub**3 / denom_sub
       f_x5 = betar_sub * bessj(0, gn_val*rz_sub)             / denom_sub
       f_x6 = betar_sub * bessj(0, gn_val*rz_sub) * rz_sub**2 / denom_sub

       ! Accumulate with Gauss-Legendre weights
       fx1_seg = fx1_seg + wt(k) * f_x1
       fx2_seg = fx2_seg + wt(k) * f_x2
       fx3_seg = fx3_seg + wt(k) * f_x3
       fx4_seg = fx4_seg + wt(k) * f_x4
       fx5_seg = fx5_seg + wt(k) * f_x5
       fx6_seg = fx6_seg + wt(k) * f_x6
    enddo

    ! Scale by (r2-r1)/2 for the mapping from [-1,1] to [r1,r2]
    fx1_seg = fx1_seg * half_dr
    fx2_seg = fx2_seg * half_dr
    fx3_seg = fx3_seg * half_dr
    fx4_seg = fx4_seg * half_dr
    fx5_seg = fx5_seg * half_dr
    fx6_seg = fx6_seg * half_dr

  end subroutine gauss4_segment_f

end subroutine four

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
