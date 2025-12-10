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
                   rab(ndmx), betar(ndmx), val_z_x5, betar_z
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
               ! Store smooth part without Bessel for fine-grid interpolation
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
         if (iz.eq.1) then
            dr=r(iz)
         else
            dr=r(iz)-r(iz-1)
         endif
         zr=r(iz)-abs(zsl(kz))
         if (lb.ne.3) then
            call simpson(nmeshs-iz+1,x1(iz),rab(iz),fx1(kz))
         endif
         if (lb.eq.0) then
            if (iz.eq.1) then
               x1(iz-1)=betar(iz)-betar(iz)/dr*zr
            else
               x1(iz-1)=betar(iz)-(betar(iz)-betar(iz-1))/dr*zr
            endif
            fx1(kz)=fx1(kz)+(x1(iz-1)+x1(iz))*0.5d0*zr
         elseif (lb.ne.3) then
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
            ! Use fine-grid integration for f-orbitals with improved accuracy
            ! Compute boundary value for singular channel (x5)
            if(iz.eq.1) then
               betar_z = betar(iz)
            else
               betar_z = betar(iz) - (betar(iz)-betar(iz-1))/dr*zr
            endif
            val_z_x5 = betar_z / (abs(zsl(kz))**3)
            
            ! Call integrate_fine for each channel with appropriate m-values
            ! x1: m=3, x2: m=2, x3: m=1, x4: m=1, x5: m=0, x6: m=0
            call integrate_fine(nmeshs-iz+1, iz, r, x1, zsl(kz), gn, 3, 0.d0, fx1(kz))
            call integrate_fine(nmeshs-iz+1, iz, r, x2, zsl(kz), gn, 2, 0.d0, fx2(kz))
            call integrate_fine(nmeshs-iz+1, iz, r, x3, zsl(kz), gn, 1, 0.d0, fx3(kz))
            call integrate_fine(nmeshs-iz+1, iz, r, x4, zsl(kz), gn, 1, 0.d0, fx4(kz))
            call integrate_fine(nmeshs-iz+1, iz, r, x5, zsl(kz), gn, 0, val_z_x5, fx5(kz))
            call integrate_fine(nmeshs-iz+1, iz, r, x6, zsl(kz), gn, 0, 0.d0, fx6(kz))
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
subroutine integrate_fine(n_coarse, iz, r_coarse, y_smooth, z_val, &
                          gn, m_bessel, val_z, result)
!-----------------------------------------------------------------------
!
! This subroutine performs fine-grid integration for f-orbital (l=3)
! calculations with improved accuracy.
!
! Input:
!   n_coarse   : number of coarse grid points
!   iz         : starting index in coarse grid
!   r_coarse   : radial grid points (coarse)
!   y_smooth   : smooth integrand values on coarse grid (without Bessel)
!   z_val      : z-coordinate value
!   gn         : g-vector magnitude
!   m_bessel   : Bessel function order (0, 1, 2, or 3)
!   val_z      : boundary value at |z|
!
! Output:
!   result     : integral value
!
  USE kinds, ONLY: DP
  USE radial_grids, only : ndmx
  implicit none
  
  integer, intent(in) :: n_coarse, iz, m_bessel
  real(DP), intent(in) :: r_coarse(ndmx), y_smooth(0:ndmx), z_val, gn, val_z
  real(DP), intent(out) :: result
  
  integer :: i, n_fine, idx
  real(DP) :: r_start, r_end, dr_fine, r_perp, bessj
  real(DP), allocatable :: r_fine(:), y_fine(:), integrand(:)
  real(DP), allocatable :: a(:), b(:), c(:), d(:)
  real(DP) :: h, dx, t, y0, y1, y2
  integer :: j
  real(DP), parameter :: eps = 1.d-10
  
  ! Define starting and ending points
  r_start = abs(z_val)
  r_end = r_coarse(iz + n_coarse - 1)
  
  ! Build dense linear grid (ensure odd number for Simpson's rule)
  n_fine = max(1001, 2*n_coarse + 1)
  if (mod(n_fine, 2) == 0) n_fine = n_fine + 1
  
  allocate(r_fine(n_fine))
  allocate(y_fine(n_fine))
  allocate(integrand(n_fine))
  
  ! Create fine linear grid
  dr_fine = (r_end - r_start) / real(n_fine - 1, DP)
  do i = 1, n_fine
    r_fine(i) = r_start + real(i - 1, DP) * dr_fine
  enddo
  
  ! Natural cubic spline interpolation
  ! Allocate spline coefficients
  allocate(a(n_coarse))
  allocate(b(n_coarse))
  allocate(c(n_coarse))
  allocate(d(n_coarse))
  
  ! Build spline from coarse grid
  call build_natural_spline(n_coarse, r_coarse(iz:iz+n_coarse-1), &
                            y_smooth(iz:iz+n_coarse-1), a, b, c, d)
  
  ! Interpolate onto fine grid
  do i = 1, n_fine
    call eval_spline(n_coarse, r_coarse(iz:iz+n_coarse-1), &
                     a, b, c, d, r_fine(i), y_fine(i))
  enddo
  
  ! Multiply by Bessel function and compute integrand
  do i = 1, n_fine
    ! Guard against roundoff: r_perp = sqrt(r^2 - z^2)
    if (r_fine(i)**2 > z_val**2 + eps) then
      r_perp = sqrt(r_fine(i)**2 - z_val**2)
    else
      r_perp = 0.d0
    endif
    
    ! Multiply smooth part by appropriate Bessel function
    integrand(i) = y_fine(i) * bessj(m_bessel, gn * r_perp)
  enddo
  
  ! Handle boundary value at r = |z|
  integrand(1) = val_z
  
  ! Apply composite Simpson's rule
  result = 0.d0
  do i = 1, n_fine - 2, 2
    result = result + (integrand(i) + 4.d0*integrand(i+1) + integrand(i+2))
  enddo
  result = result * dr_fine / 3.d0
  
  deallocate(r_fine, y_fine, integrand)
  deallocate(a, b, c, d)
  
end subroutine integrate_fine

!-----------------------------------------------------------------------
subroutine build_natural_spline(n, x, y, a, b, c, d)
!-----------------------------------------------------------------------
!
! Build natural cubic spline coefficients
!
  USE kinds, ONLY: DP
  implicit none
  
  integer, intent(in) :: n
  real(DP), intent(in) :: x(n), y(n)
  real(DP), intent(out) :: a(n), b(n), c(n), d(n)
  
  integer :: i
  real(DP), allocatable :: h(:), alpha(:), l(:), mu(:), z(:)
  
  allocate(h(n-1), alpha(n-1), l(n), mu(n-1), z(n))
  
  ! Copy y values to a
  do i = 1, n
    a(i) = y(i)
  enddo
  
  ! Compute h values
  do i = 1, n-1
    h(i) = x(i+1) - x(i)
  enddo
  
  ! Compute alpha values
  do i = 2, n-1
    alpha(i) = (3.d0/h(i)) * (a(i+1) - a(i)) - &
               (3.d0/h(i-1)) * (a(i) - a(i-1))
  enddo
  
  ! Solve tridiagonal system (natural spline: c(1)=c(n)=0)
  l(1) = 1.d0
  mu(1) = 0.d0
  z(1) = 0.d0
  
  do i = 2, n-1
    l(i) = 2.d0 * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1)
    mu(i) = h(i) / l(i)
    z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
  enddo
  
  l(n) = 1.d0
  z(n) = 0.d0
  c(n) = 0.d0
  
  ! Back substitution
  do i = n-1, 1, -1
    c(i) = z(i) - mu(i) * c(i+1)
    b(i) = (a(i+1) - a(i)) / h(i) - h(i) * (c(i+1) + 2.d0*c(i)) / 3.d0
    d(i) = (c(i+1) - c(i)) / (3.d0 * h(i))
  enddo
  
  deallocate(h, alpha, l, mu, z)
  
end subroutine build_natural_spline

!-----------------------------------------------------------------------
subroutine eval_spline(n, x, a, b, c, d, x_eval, y_eval)
!-----------------------------------------------------------------------
!
! Evaluate cubic spline at a given point
!
  USE kinds, ONLY: DP
  implicit none
  
  integer, intent(in) :: n
  real(DP), intent(in) :: x(n), a(n), b(n), c(n), d(n), x_eval
  real(DP), intent(out) :: y_eval
  
  integer :: i
  real(DP) :: dx
  
  ! Find interval containing x_eval
  i = 1
  do while (i < n .and. x_eval > x(i+1))
    i = i + 1
  enddo
  
  ! Evaluate spline
  if (i < n) then
    dx = x_eval - x(i)
    y_eval = a(i) + b(i)*dx + c(i)*dx**2 + d(i)*dx**3
  else
    ! Extrapolate using last interval
    dx = x_eval - x(n-1)
    y_eval = a(n-1) + b(n-1)*dx + c(n-1)*dx**2 + d(n-1)*dx**3
  endif
  
end subroutine eval_spline
