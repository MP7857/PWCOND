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
            ! For p,d keep old wedge; for f we will override below
            if (lb.ne.3) then
               fx1(kz) = fx1(kz) + x1(iz)*0.5d0*zr
            endif
            call simpson( nmeshs-iz+1, x2(iz), rab(iz), fx2(kz) )
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
            ! --- f: high-order wedge treatment + standard Simpson bulk ---

            ! Bulk part (r >= r(iz)) as before
            fx2(kz) = fx2(kz) + x2(iz)*0.5d0*zr
            call simpson( nmeshs-iz+1, x3(iz), rab(iz), fx3(kz) )
            fx3(kz) = fx3(kz) + x3(iz)*0.5d0*zr
            call simpson( nmeshs-iz+1, x4(iz), rab(iz), fx4(kz) )
            fx4(kz) = fx4(kz) + x4(iz)*0.5d0*zr
            call simpson( nmeshs-iz+1, x5(iz), rab(iz), fx5(kz) )
            call simpson( nmeshs-iz+1, x6(iz), rab(iz), fx6(kz) )

            ! High-order Gauss–Legendre integration over the wedge r ∈ [|z|, r(iz)]
            call wedge_gl_f( gn, abs(zsl(kz)), betar, r, ndmx, iz,  &
                             x1(0), x2(0), x3(0), x4(0), x5(0), x6(0) )

            fx1(kz) = fx1(kz) + x1(0)
            fx2(kz) = fx2(kz) + x2(0)
            fx3(kz) = fx3(kz) + x3(0)
            fx4(kz) = fx4(kz) + x4(0)
            fx5(kz) = fx5(kz) + x5(0)
            fx6(kz) = fx6(kz) + x6(0)

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

subroutine wedge_gl_f( gn, zabs, betar, r, ndmx, iz, fx1w, fx2w, fx3w, fx4w, fx5w, fx6w )
  !
  ! High-order Gauss–Legendre integration of the f-projector wedge:
  !   r ∈ [ |z|, r(iz) ], with r(iz-1) <= |z| < r(iz)
  !
  ! For lb = 3 we need the radial pieces corresponding to:
  !   x1 ~ betar * J3(gn*ρ) * ρ^3 / r^3
  !   x2 ~ betar * J2(gn*ρ) * ρ^2 / r^3
  !   x3 ~ betar * J1(gn*ρ) * ρ   / r^3
  !   x4 ~ betar * J1(gn*ρ) * ρ^3 / r^3
  !   x5 ~ betar * J0(gn*ρ)       / r^3
  !   x6 ~ betar * J0(gn*ρ) * ρ^2 / r^3
  !
  USE kinds,     ONLY : DP
  implicit none

  integer,  intent(in)  :: ndmx, iz
  real(DP), intent(in)  :: gn, zabs
  real(DP), intent(in)  :: betar(ndmx), r(ndmx)
  real(DP), intent(out) :: fx1w, fx2w, fx3w, fx4w, fx5w, fx6w

  ! 4-point Gauss–Legendre nodes/weights on [0,1]
  real(DP), parameter :: t1 = 0.8611363115940525752_DP
  real(DP), parameter :: t2 = 0.3399810435848562648_DP
  real(DP), parameter :: w1 = 0.3478548451374538574_DP
  real(DP), parameter :: w2 = 0.6521451548625461426_DP

  real(DP) :: s(4), w(4)
  real(DP) :: rlow, rhigh, dr, rj, rho, bj0, bj1, bj2, bj3
  real(DP) :: br, r_im1, br_im1
  real(DP) :: alpha
  real(DP) :: bessj
  integer  :: j

  external bessj

  ! Initialize outputs
  fx1w = 0.d0
  fx2w = 0.d0
  fx3w = 0.d0
  fx4w = 0.d0
  fx5w = 0.d0
  fx6w = 0.d0

  ! If gn = 0, the f-channel contributions vanish anyway
  if (gn .eq. 0.d0) return

  rlow  = zabs
  rhigh = r(iz)
  dr    = rhigh - rlow
  if (dr <= 0.d0) return

  ! Map 4-point GL from [-1,1] to [0,1]
  s(1) = 0.5d0*(1.d0 - t1)
  s(2) = 0.5d0*(1.d0 - t2)
  s(3) = 0.5d0*(1.d0 + t2)
  s(4) = 0.5d0*(1.d0 + t1)
  w(1) = 0.5d0*w1
  w(2) = 0.5d0*w2
  w(3) = 0.5d0*w2
  w(4) = 0.5d0*w1

  do j = 1, 4
     rj = rlow + s(j)*dr
     if (rj <= zabs) cycle

     ! Interpolate betar(rj) between r(iz-1) and r(iz)
     if (iz > 1) then
        r_im1  = r(iz-1)
        br_im1 = betar(iz-1)
        if (r(iz) > r_im1) then
           alpha = (rj - r_im1) / (r(iz) - r_im1)
           if (alpha < 0.d0) alpha = 0.d0
           if (alpha > 1.d0) alpha = 1.d0
           br = br_im1 + alpha*(betar(iz) - br_im1)
        else
           br = betar(iz)
        endif
     else
        ! Near the origin: assume regular f-like behavior ∝ r^3
        if (r(iz) > 0.d0) then
           br = betar(iz) * (rj / r(iz))**3
        else
           br = 0.d0
        endif
     endif

     rho = sqrt( max( rj*rj - zabs*zabs, 0.d0 ) )
     if (rho <= 0.d0) cycle

     bj0 = bessj(0, gn*rho)
     bj1 = bessj(1, gn*rho)
     bj2 = bessj(2, gn*rho)
     bj3 = bessj(3, gn*rho)

     ! f-integrands consistent with x1..x6 definitions in four.f90
     fx1w = fx1w + w(j)*dr * br * bj3 * (rho**3) / (rj**3)
     fx2w = fx2w + w(j)*dr * br * bj2 * (rho**2) / (rj**3)
     fx3w = fx3w + w(j)*dr * br * bj1 * (rho    ) / (rj**3)
     fx4w = fx4w + w(j)*dr * br * bj1 * (rho**3) / (rj**3)
     fx5w = fx5w + w(j)*dr * br * bj0               / (rj**3)
     fx6w = fx6w + w(j)*dr * br * bj0 * (rho**2) / (rj**3)

  enddo

end subroutine wedge_gl_f

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
