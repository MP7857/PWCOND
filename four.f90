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
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x1, zsl(kz), gn, 3, fx1(kz))
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x2, zsl(kz), gn, 2, fx2(kz))
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x3, zsl(kz), gn, 1, fx3(kz))
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x4, zsl(kz), gn, 1, fx4(kz))
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x5, zsl(kz), gn, 0, fx5(kz))
            call integrate_fine_from_arrays(nmeshs-iz+1, iz, r, x6, zsl(kz), gn, 0, fx6(kz))
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
subroutine integrate_fine_from_arrays(npts, iz, r, x_smooth, z, g, m, result)
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
  ! Output:
  !   result: integrated value
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE radial_grids, only : ndmx
  implicit none
  integer, intent(in) :: npts, iz, m
  real(DP), intent(in) :: r(ndmx), x_smooth(0:ndmx), z, g
  real(DP), intent(out) :: result

  integer :: i, n_fine, ir, nmesh_end
  real(DP) :: zabs, r_start, r_end, dr_fine, r_curr, r_perp
  real(DP) :: val_smooth, val_bess, b_arg, bessj
  real(DP), allocatable :: y_fine(:)

  zabs = abs(z)
  nmesh_end = iz + npts - 1
  ! Start from abs(z) or first grid point, whichever is larger
  ! This ensures we don't integrate in the unphysical region r < |z|
  r_start = max(zabs, r(1))
  r_end = r(nmesh_end)
  
  ! Safety check
  if (r_start >= r_end) then
     result = 0.d0
     return
  endif

  ! Define Fine Grid Density
  ! Bessel oscillation period ~ 2*pi/g, so we want at least ~10 points per period
  ! This gives grid spacing ~ pi/(5*g), hence n_fine ~ g * range / (pi/(5*g)) ~ 5*g*range/pi
  ! Using factor of 5 (instead of ~5/pi) provides extra safety margin
  n_fine = max(500, int(g * (r_end - r_start) * 5.0d0))
  n_fine = min(n_fine, 5000) ! Cap to avoid excessive cost
  
  ! Ensure odd number for Simpson's rule
  if (mod(n_fine, 2) .eq. 0) n_fine = n_fine + 1
  
  allocate(y_fine(n_fine))

  dr_fine = (r_end - r_start) / dble(n_fine - 1)

  ! Build Integrand on Fine Grid
  do i = 1, n_fine
     r_curr = r_start + dble(i-1)*dr_fine
     
     ! Interpolate smooth part from x_smooth array
     call lin_interp_arrays(npts, iz, r, x_smooth, r_curr, val_smooth)
     
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
subroutine lin_interp_arrays(npts, iz, x, y, x_eval, y_eval)
  !-----------------------------------------------------------------------
  ! Linear interpolation for arrays starting at index iz
  ! Interpolates y values defined on x grid from iz to iz+npts-1
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  USE radial_grids, only : ndmx
  implicit none
  integer, intent(in) :: npts, iz
  real(DP), intent(in) :: x(ndmx), y(0:ndmx), x_eval
  real(DP), intent(out) :: y_eval
  integer :: i, iend

  iend = iz + npts - 1
  
  if (x_eval <= x(iz)) then
     y_eval = y(iz)
     return
  endif
  if (x_eval >= x(iend)) then
     y_eval = y(iend)
     return
  endif

  ! Linear search
  do i = iz, iend-1
     if (x_eval >= x(i) .and. x_eval <= x(i+1)) then
        y_eval = y(i) + (y(i+1)-y(i)) * (x_eval - x(i)) / (x(i+1)-x(i))
        return
     endif
  enddo
  
  ! Fallback
  y_eval = y(iend)
  
end subroutine lin_interp_arrays
