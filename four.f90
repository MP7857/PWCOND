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
! This version uses only cubic-spline interpolation on the fine grid
! for all f-orbital integrals. The monotone PCHIP interpolation is
! retained in the source for reference but is not invoked.
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
!             f_{z(5z^2-3r^2)}, f_{-x(5z^2-r^2)}, f_{-y(5z^2-r^2)},
!             f_{z(x^2-y^2)}, f_{xyz}, f_{-x(x^2-3y^2)}, f_{-y(3x^2-y^2)}
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
                   rab(ndmx), betar(ndmx), val_beta_z, val_z_x5
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
               ! Store SMOOTH parts only (without Bessel functions)
               ! Bessel functions will be evaluated on fine grid later
               x1(ir)=betar(ir)*rz**3/r(ir)**3
               x2(ir)=betar(ir)*rz**2/r(ir)**3
               x3(ir)=betar(ir)*rz/r(ir)**3
               x4(ir)=betar(ir)*rz**3/r(ir)**3
               x5(ir)=betar(ir)/r(ir)**3
               x6(ir)=betar(ir)*rz**2/r(ir)**3
            else
               call errore ('four','ls not programmed ',1)
            endif
         enddo
         if (lb.ne.3) then
            ! Standard integration for lb=0,1,2
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
            ! Use robust fine-grid integration for f-orbitals
            ! This addresses accuracy issues at high g-values by:
            ! 1. Interpolating smooth part onto fine grid
            ! 2. Evaluating Bessel functions on fine grid
            ! 3. Integrating the product
            
            ! --- Calculate Smooth Integrand Value at r = |z| ---
            ! Most terms (x1,x2,x3,x4,x6) depend on r_perp = sqrt(r^2-z^2).
            ! At r=|z|, r_perp=0, so these terms are exactly 0.0.
            ! Only x5 (which is beta(r)/r^3) is non-zero at |z|.
            
            ! Calculate x5 at |z| by interpolating beta(r)
            zr = abs(zsl(kz))
            
            ! Linear interpolation of beta to |z|
            if (iz > 1) then
               ! Interpolate between iz-1 and iz using standard form
               dr = r(iz) - r(iz-1)
               val_beta_z = betar(iz-1) + (betar(iz)-betar(iz-1)) * (zr-r(iz-1)) / dr
            else
               ! Extrapolate if zr < r(1) - assumes linear behavior near origin
               val_beta_z = betar(1) * (zr / r(1))
            endif
            
            ! Calculate x5 start value (beta(z) / z^3)
            if (zr > 1.0d-9) then
               val_z_x5 = val_beta_z / (zr**3)
            else
               val_z_x5 = 0.d0
            endif
            
            ! --- Call Robust Integration with Start Value ---
            ! Use cubic spline interpolation for all f-orbital integrals
            ! (PCHIP is retained in source for reference but not invoked)
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x1, zsl(kz), gn, 3, 0.d0, fx1(kz), .false.)
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x2, zsl(kz), gn, 2, 0.d0, fx2(kz), .false.)
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x3, zsl(kz), gn, 1, 0.d0, fx3(kz), .false.)
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x4, zsl(kz), gn, 1, 0.d0, fx4(kz), .false.)
            
            ! Pass the calculated non-zero start value for x5
            ! Using cubic spline (not PCHIP) for simplicity
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x5, zsl(kz), gn, 0, val_z_x5, fx5(kz), .false.)
            
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x6, zsl(kz), gn, 0, 0.d0, fx6(kz), .false.)
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

!-----------------------------------------------------------------------
subroutine integrate_fine_from_arrays(npts, iz, r, x_smooth, z, g, m, val_at_z, result, use_pchip)
  !-----------------------------------------------------------------------
  ! Interpolates smooth integrand from x_smooth array onto a fine linear 
  ! grid, evaluates Bessel function J_m on that grid, and integrates.
  ! 
  ! This addresses integration accuracy issues for highly oscillatory
  ! integrands (f-orbitals at high g-values) by sampling oscillations
  ! on a fine grid rather than the coarse logarithmic grid.
  !
  ! Input:
  !   npts: number of points from iz to end
  !   iz: starting index in arrays
  !   r: radial grid  
  !   x_smooth: smooth part of integrand (from iz onwards)
  !   z: z-coordinate of slice
  !   g: magnitude of g-vector
  !   m: Bessel function order
  !   val_at_z: Value of smooth function at r=|z|
  !   use_pchip: kept for compatibility but always uses cubic spline
  ! Output:
  !   result: integrated value
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE radial_grids, only : ndmx
  implicit none
  integer, intent(in) :: npts, iz, m
  real(DP), intent(in) :: r(ndmx), x_smooth(0:ndmx), z, g, val_at_z
  logical, intent(in) :: use_pchip
  real(DP), intent(out) :: result
  
  ! DECLARE BESSJ AS EXTERNAL TO AVOID ARRAY CONFUSION
  real(DP), external :: bessj

  integer :: i, n_fine, ir, nmesh_end
  real(DP) :: zabs, r_start, r_end, dr_fine, r_curr, r_perp
  real(DP) :: val_smooth, val_bess, b_arg
  real(DP), allocatable :: y_fine(:)

  zabs = abs(z)
  nmesh_end = iz + npts - 1
  
  ! Start EXACTLY at |z|
  r_start = zabs
  r_end = r(nmesh_end)
  
  ! Safety check
  if (r_start >= r_end) then
     result = 0.d0
     return
  endif

  ! Define Fine Grid Density
  ! Bessel oscillation period ~ 2*pi/g, so we want at least 20-40 points per period
  ! for high accuracy with cubic/PCHIP interpolation
  ! Target: 40 points/period gives spacing = (2*pi/g)/40 = pi/(20*g)
  ! Therefore: n_fine = range / spacing = range * 20*g/pi ~ 6.37 * g * range
  ! We use factor of 6.4 for ~40 points per Bessel period
  n_fine = max(500, int(g * (r_end - r_start) * 6.4d0))
  n_fine = min(n_fine, 10000) ! Cap to avoid excessive cost
  
  ! Ensure odd number for Simpson's rule
  if (mod(n_fine, 2) .eq. 0) n_fine = n_fine + 1
  
  allocate(y_fine(n_fine))

  dr_fine = (r_end - r_start) / dble(n_fine - 1)

  ! Build Integrand on Fine Grid
  do i = 1, n_fine
     r_curr = r_start + dble(i-1)*dr_fine
     
     ! Use cubic spline interpolation for all f-orbital integrals
     ! (use_pchip flag is ignored; PCHIP code retained for reference)
     call cubic_spline_interp(npts, iz, r, x_smooth, r_curr, zabs, val_at_z, val_smooth)
     
     ! Compute Bessel function
     r_perp = 0.d0
     if (r_curr > zabs) r_perp = sqrt(r_curr**2 - zabs**2)
     
     b_arg = g * r_perp
     val_bess = bessj(m, b_arg)
     
     ! Combine
     y_fine(i) = val_smooth * val_bess
  enddo

  ! Simpson's Rule Integration on uniform grid
  result = 0.d0
  if (n_fine >= 3) then
     result = y_fine(1) + y_fine(n_fine)
     do i = 2, n_fine - 1, 2
        result = result + 4.d0 * y_fine(i)
     enddo
     do i = 3, n_fine - 2, 2
        result = result + 2.d0 * y_fine(i)
     enddo
     result = result * dr_fine / 3.d0
  endif

  deallocate(y_fine)

end subroutine integrate_fine_from_arrays

!-----------------------------------------------------------------------
subroutine lin_interp_arrays(npts, iz, x, y, x_eval, z_abs, val_at_z, y_eval)
  !-----------------------------------------------------------------------
  ! Linear interpolation for arrays starting at index iz
  ! Input added: z_abs (start of physical range), val_at_z (value at z_abs)
  ! Interpolates y values defined on x grid from iz to iz+npts-1
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE radial_grids, only : ndmx
  implicit none
  integer, intent(in) :: npts, iz
  real(DP), intent(in) :: x(ndmx), y(0:ndmx), x_eval, z_abs, val_at_z
  real(DP), intent(out) :: y_eval
  integer :: i, iend

  iend = iz + npts - 1
  
  ! CASE 1: x_eval is in the "gap" between |z| and the first grid point r(iz)
  if (x_eval < x(iz)) then
      ! Interpolate between (z_abs, val_at_z) and (x(iz), y(iz))
      ! Add tolerance check to prevent division by near-zero values
      if (x(iz) - z_abs > 1.0d-12) then
         y_eval = val_at_z + (y(iz) - val_at_z) * (x_eval - z_abs) / (x(iz) - z_abs)
      else
         y_eval = val_at_z
      endif
      return
  endif

  ! CASE 2: x_eval is beyond the grid (clamp)
  if (x_eval >= x(iend)) then
     y_eval = y(iend)
     return
  endif

  ! CASE 3: Standard Linear search within the grid
  do i = iz, iend-1
     if (x_eval >= x(i) .and. x_eval <= x(i+1)) then
        y_eval = y(i) + (y(i+1)-y(i)) * (x_eval - x(i)) / (x(i+1)-x(i))
        return
     endif
  enddo
  
  y_eval = y(iend)
  
end subroutine lin_interp_arrays

!-----------------------------------------------------------------------
subroutine pchip_interp(npts, iz, x, y, x_eval, z_abs, val_at_z, y_eval)
  !-----------------------------------------------------------------------
  ! PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) interpolation
  ! This is a monotone-preserving cubic interpolation, ideal for functions
  ! with singularities like 1/r^3. Based on Fritsch-Carlson algorithm.
  !
  ! Reference: Fritsch & Carlson (1980), "Monotone Piecewise Cubic Interpolation"
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE radial_grids, only : ndmx
  implicit none
  integer, intent(in) :: npts, iz
  real(DP), intent(in) :: x(ndmx), y(0:ndmx), x_eval, z_abs, val_at_z
  real(DP), intent(out) :: y_eval
  
  integer :: i, iend, n_nodes
  real(DP), allocatable :: x_nodes(:), y_nodes(:), d(:), delta_k(:)
  real(DP) :: h, t, h00, h10, h01, h11, h1, h2, w1, w2
  real(DP) :: alpha, beta, tau
  
  iend = iz + npts - 1
  
  ! Handle boundary case: in gap [z_abs, x(iz)]
  if (x_eval < x(iz)) then
     if (x(iz) - z_abs > 1.0d-12) then
        y_eval = val_at_z + (y(iz) - val_at_z) * (x_eval - z_abs) / (x(iz) - z_abs)
     else
        y_eval = val_at_z
     endif
     return
  endif
  
  ! Beyond grid - clamp
  if (x_eval >= x(iend)) then
     y_eval = y(iend)
     return
  endif
  
  ! Build node list including boundary point
  n_nodes = npts + 1
  allocate(x_nodes(n_nodes), y_nodes(n_nodes), d(n_nodes))
  
  x_nodes(1) = z_abs
  y_nodes(1) = val_at_z
  do i = 2, n_nodes
     x_nodes(i) = x(iz + i - 2)
     y_nodes(i) = y(iz + i - 2)
  enddo
  
  ! Compute PCHIP tangents using Fritsch-Carlson method
  ! Step 1: Compute secant slopes (deltas)
  allocate(delta_k(n_nodes-1))
  do i = 1, n_nodes - 1
     delta_k(i) = (y_nodes(i+1) - y_nodes(i)) / (x_nodes(i+1) - x_nodes(i))
  enddo
  
  ! Step 2: Estimate derivatives at interior points
  d(1) = delta_k(1)  ! One-sided at start
  d(n_nodes) = delta_k(n_nodes-1)  ! One-sided at end
  
  do i = 2, n_nodes - 1
     ! Check for monotonicity - if slopes change sign or either is zero, set d=0
     if (delta_k(i-1) * delta_k(i) <= 0.d0 .or. abs(delta_k(i-1)) < 1.0d-12 .or. abs(delta_k(i)) < 1.0d-12) then
        d(i) = 0.d0
     else
        ! Weighted harmonic mean (safe since both deltas are non-zero and same sign)
        h1 = x_nodes(i) - x_nodes(i-1)
        h2 = x_nodes(i+1) - x_nodes(i)
        w1 = 2.d0*h2 + h1
        w2 = h2 + 2.d0*h1
        d(i) = (w1 + w2) / (w1/delta_k(i-1) + w2/delta_k(i))
     endif
  enddo
  
  ! Step 3: Adjust derivatives to ensure monotonicity in each interval
  ! This uses the Fritsch-Carlson constraint: alpha^2 + beta^2 <= 9
  ! where alpha = d(i)/delta(i) and beta = d(i+1)/delta(i)
  do i = 1, n_nodes - 1
     if (abs(delta_k(i)) < 1.0d-12) then
        d(i) = 0.d0
        d(i+1) = 0.d0
     else
        alpha = d(i) / delta_k(i)
        beta = d(i+1) / delta_k(i)
        ! Check if (alpha,beta) is outside the monotonicity region
        ! The factor of 9 comes from the Fritsch-Carlson PCHIP algorithm
        if (alpha**2 + beta**2 > 9.d0) then
           tau = 3.d0 / sqrt(alpha**2 + beta**2)
           d(i) = tau * alpha * delta_k(i)
           d(i+1) = tau * beta * delta_k(i)
        endif
     endif
  enddo
  
  deallocate(delta_k)
  
  ! Step 3: Find interval and evaluate cubic Hermite polynomial
  do i = 1, n_nodes - 1
     if (x_eval >= x_nodes(i) .and. x_eval <= x_nodes(i+1)) then
        h = x_nodes(i+1) - x_nodes(i)
        t = (x_eval - x_nodes(i)) / h
        
        ! Hermite basis functions
        h00 = (1.d0 + 2.d0*t) * (1.d0 - t)**2
        h10 = t * (1.d0 - t)**2
        h01 = t**2 * (3.d0 - 2.d0*t)
        h11 = t**2 * (t - 1.d0)
        
        y_eval = h00*y_nodes(i) + h10*h*d(i) + h01*y_nodes(i+1) + h11*h*d(i+1)
        
        deallocate(x_nodes, y_nodes, d)
        return
     endif
  enddo
  
  ! Fallback
  y_eval = y(iend)
  deallocate(x_nodes, y_nodes, d)
  
end subroutine pchip_interp

!-----------------------------------------------------------------------
subroutine cubic_spline_interp(npts, iz, x, y, x_eval, z_abs, val_at_z, y_eval)
  !-----------------------------------------------------------------------
  ! Cubic spline interpolation (natural spline with zero second derivative
  ! at boundaries). Better than linear for smooth functions.
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE radial_grids, only : ndmx
  implicit none
  integer, intent(in) :: npts, iz
  real(DP), intent(in) :: x(ndmx), y(0:ndmx), x_eval, z_abs, val_at_z
  real(DP), intent(out) :: y_eval
  
  integer :: i, j, iend, n_nodes
  real(DP), allocatable :: x_nodes(:), y_nodes(:), a(:), b(:), c(:), d_coef(:), h(:), alpha(:), l(:), mu(:), z_vec(:)
  real(DP) :: t, dx
  
  iend = iz + npts - 1
  
  ! Handle boundary case
  if (x_eval < x(iz)) then
     if (x(iz) - z_abs > 1.0d-12) then
        y_eval = val_at_z + (y(iz) - val_at_z) * (x_eval - z_abs) / (x(iz) - z_abs)
     else
        y_eval = val_at_z
     endif
     return
  endif
  
  if (x_eval >= x(iend)) then
     y_eval = y(iend)
     return
  endif
  
  ! Build node list
  n_nodes = npts + 1
  allocate(x_nodes(n_nodes), y_nodes(n_nodes))
  allocate(a(n_nodes), b(n_nodes-1), c(n_nodes-1), d_coef(n_nodes-1), h(n_nodes-1))
  allocate(alpha(n_nodes-1), l(n_nodes), mu(n_nodes), z_vec(n_nodes))
  
  x_nodes(1) = z_abs
  y_nodes(1) = val_at_z
  do i = 2, n_nodes
     x_nodes(i) = x(iz + i - 2)
     y_nodes(i) = y(iz + i - 2)
  enddo
  
  ! Natural cubic spline algorithm
  do i = 1, n_nodes
     a(i) = y_nodes(i)
  enddo
  
  do i = 1, n_nodes - 1
     h(i) = x_nodes(i+1) - x_nodes(i)
  enddo
  
  do i = 2, n_nodes - 1
     alpha(i) = 3.d0*(a(i+1)-a(i))/h(i) - 3.d0*(a(i)-a(i-1))/h(i-1)
  enddo
  
  ! Solve tridiagonal system
  l(1) = 1.d0
  mu(1) = 0.d0
  z_vec(1) = 0.d0
  
  do i = 2, n_nodes - 1
     l(i) = 2.d0*(x_nodes(i+1)-x_nodes(i-1)) - h(i-1)*mu(i-1)
     mu(i) = h(i)/l(i)
     z_vec(i) = (alpha(i) - h(i-1)*z_vec(i-1))/l(i)
  enddo
  
  l(n_nodes) = 1.d0
  z_vec(n_nodes) = 0.d0
  c(n_nodes-1) = 0.d0
  
  do j = n_nodes - 1, 1, -1
     c(j) = z_vec(j) - mu(j)*c(j+1)
     b(j) = (a(j+1)-a(j))/h(j) - h(j)*(c(j+1)+2.d0*c(j))/3.d0
     d_coef(j) = (c(j+1)-c(j))/(3.d0*h(j))
  enddo
  
  ! Evaluate at x_eval
  do i = 1, n_nodes - 1
     if (x_eval >= x_nodes(i) .and. x_eval <= x_nodes(i+1)) then
        dx = x_eval - x_nodes(i)
        y_eval = a(i) + b(i)*dx + c(i)*dx**2 + d_coef(i)*dx**3
        
        deallocate(x_nodes, y_nodes, a, b, c, d_coef, h, alpha, l, mu, z_vec)
        return
     endif
  enddo
  
  ! Fallback
  y_eval = y(iend)
  deallocate(x_nodes, y_nodes, a, b, c, d_coef, h, alpha, l, mu, z_vec)
  
end subroutine cubic_spline_interp
