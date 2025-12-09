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
                   rab(ndmx), betar(ndmx), fx1w, fx2w, fx3w, fx4w, fx5w, fx6w
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
            ! Integrate r >= r(iz) with Simpson (unchanged)
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            fx3(kz)=fx3(kz)+x3(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            fx4(kz)=fx4(kz)+x4(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x5(iz),rab(iz),fx5(kz))
            call simpson(nmeshs-iz+1,x6(iz),rab(iz),fx6(kz))

            ! NEW: more accurate wedge [|z|, r(iz)] for all fx1..fx6
            call f_wedge_integral_gauss( abs(zsl(kz)), r(iz), gn, nmesh, r, betar, &
                                         fx1w, fx2w, fx3w, fx4w, fx5w, fx6w )

            fx1(kz) = fx1(kz) + fx1w
            fx2(kz) = fx2(kz) + fx2w
            fx3(kz) = fx3(kz) + fx3w
            fx4(kz) = fx4(kz) + fx4w
            fx5(kz) = fx5(kz) + fx5w
            fx6(kz) = fx6(kz) + fx6w
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
! Quadratic interpolation of betar(rq) on the radial grid.
! Uses three nearby points when possible, otherwise falls back to linear.
subroutine interp_betar_quad(rq, nmesh, r, betar, bq)
  USE kinds, only : DP
  implicit none
  integer,  intent(in)  :: nmesh
  real(DP), intent(in)  :: rq, r(nmesh), betar(nmesh)
  real(DP), intent(out) :: bq

  integer :: i, i1, i2, i3
  real(DP) :: x1,x2,x3, y1,y2,y3
  real(DP) :: L1,L2,L3, denom12, denom13, denom23

  ! If outside grid, clamp
  if (rq <= r(1)) then
     bq = betar(1)
     return
  elseif (rq >= r(nmesh)) then
     bq = betar(nmesh)
     return
  endif

  ! Find interval [i, i+1] such that r(i) <= rq <= r(i+1)
  i = 1
  do while (r(i+1) < rq .and. i < nmesh-1)
     i = i + 1
  enddo

  ! Choose three points for quadratic interpolation
  if (i == 1) then
     i1 = 1;       i2 = 2;       i3 = 3
  elseif (i == nmesh-1) then
     i1 = nmesh-2; i2 = nmesh-1; i3 = nmesh
  else
     i1 = i-1;     i2 = i;       i3 = i+1
  endif

  x1 = r(i1); y1 = betar(i1)
  x2 = r(i2); y2 = betar(i2)
  x3 = r(i3); y3 = betar(i3)

  denom12 = (x1 - x2)
  denom13 = (x1 - x3)
  denom23 = (x2 - x3)

  ! Guard against degeneracy: fall back to linear if needed
  if (abs(denom12*denom13*denom23) < 1.d-14) then
     ! simple linear between i and i+1
     bq = y2 + (y3 - y2) * (rq - x2) / (x3 - x2)
     return
  endif

  L1 = (rq - x2)*(rq - x3) / (denom12*denom13)
  L2 = (rq - x1)*(rq - x3) / ((x2 - x1)*denom23)
  L3 = (rq - x1)*(rq - x2) / ((x3 - x1)*(x3 - x2))

  bq = y1*L1 + y2*L2 + y3*L3

  return
end subroutine interp_betar_quad

!-----------------------------------------------------------------------
! More accurate treatment of the wedge r in [|z|, r_iz] for l=3.
! Returns contributions to fx1..fx6 from this wedge.
subroutine f_wedge_integral_gauss( zabs, r_hi, gn, nmesh, r, betar,    &
                                   fx1w, fx2w, fx3w, fx4w, fx5w, fx6w )
  USE kinds, only : DP
  implicit none
  real(DP), intent(in)  :: zabs, r_hi, gn
  integer,  intent(in)  :: nmesh
  real(DP), intent(in)  :: r(nmesh), betar(nmesh)
  real(DP), intent(out) :: fx1w, fx2w, fx3w, fx4w, fx5w, fx6w

  real(DP), parameter :: eps = 1.d-12
  real(DP) :: dr, t1, t2, w1, w2, rq, rzq, bj, bessj
  real(DP) :: x1,x2,x3,x4,x5,x6

  fx1w = 0.d0
  fx2w = 0.d0
  fx3w = 0.d0
  fx4w = 0.d0
  fx5w = 0.d0
  fx6w = 0.d0

  dr = r_hi - zabs
  if (dr <= eps) return

  ! 2-point Gauss-Legendre on [0,1]: nodes and weights
  t1 = 0.5d0 - sqrt(3.d0)/6.d0
  t2 = 0.5d0 + sqrt(3.d0)/6.d0
  w1 = 0.5d0
  w2 = 0.5d0

  ! First node
  rq = zabs + t1*dr
  call interp_betar_quad( rq, nmesh, r, betar, bj )
  if (rq > zabs) then
     rzq = sqrt( max( rq*rq - zabs*zabs, 0.d0 ) )
  else
     rzq = 0.d0
  endif
  if (rzq > eps) then
     x1 = bj * bessj(3, gn*rzq) * rzq**3 / rq**3
     x2 = bj * bessj(2, gn*rzq) * rzq**2 / rq**3
     x3 = bj * bessj(1, gn*rzq) * rzq     / rq**3
     x4 = bj * bessj(1, gn*rzq) * rzq**3 / rq**3
     x5 = bj * bessj(0, gn*rzq)          / rq**3
     x6 = bj * bessj(0, gn*rzq) * rzq**2 / rq**3

     fx1w = fx1w + w1*x1
     fx2w = fx2w + w1*x2
     fx3w = fx3w + w1*x3
     fx4w = fx4w + w1*x4
     fx5w = fx5w + w1*x5
     fx6w = fx6w + w1*x6
  endif

  ! Second node
  rq = zabs + t2*dr
  call interp_betar_quad( rq, nmesh, r, betar, bj )
  if (rq > zabs) then
     rzq = sqrt( max( rq*rq - zabs*zabs, 0.d0 ) )
  else
     rzq = 0.d0
  endif
  if (rzq > eps) then
     x1 = bj * bessj(3, gn*rzq) * rzq**3 / rq**3
     x2 = bj * bessj(2, gn*rzq) * rzq**2 / rq**3
     x3 = bj * bessj(1, gn*rzq) * rzq     / rq**3
     x4 = bj * bessj(1, gn*rzq) * rzq**3 / rq**3
     x5 = bj * bessj(0, gn*rzq)          / rq**3
     x6 = bj * bessj(0, gn*rzq) * rzq**2 / rq**3

     fx1w = fx1w + w2*x1
     fx2w = fx2w + w2*x2
     fx3w = fx3w + w2*x3
     fx4w = fx4w + w2*x4
     fx5w = fx5w + w2*x5
     fx6w = fx6w + w2*x6
  endif

  ! Scale by interval length
  fx1w = fx1w * dr
  fx2w = fx2w * dr
  fx3w = fx3w * dr
  fx4w = fx4w * dr
  fx5w = fx5w * dr
  fx6w = fx6w * dr

  return
end subroutine f_wedge_integral_gauss
