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
!             f_{xyz}, f_{z(x^2-y^2)}, f_{x(x^2-3y^2)}, f_{y(3x^2-y^2)}
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
!
! --- f-only numerical safety & logging (minimal, surgical) ---
! thresholds
  real(DP), parameter :: gtol   = 1.0d-12   ! |g_perp|~0: Γ guard for f
  real(DP), parameter :: angtol = 1.0d-10   ! log angle-norm drift when |den-1|>angtol
  real(DP), parameter :: symtol = 1.0d-12   ! snap radius to hex/trig symmetry angles
  real(DP), parameter :: xtol   = 1.0d-10   ! small-x threshold for bessj stabilization
! logging
  logical, save :: flog_opened = .false.
  integer, save :: flog_unit
! angle helpers (f-only section)
  real(DP) :: den, cs0, sn0, bestd, dc, ds, tcs, tsn, d2
  integer  :: isnap, i
  real(DP), parameter :: rt3 = 1.73205080756887729353d0
  real(DP), parameter :: cs_snap(6) = (/ 1.d0,  0.5d0, -0.5d0, -1.d0,  0.5d0, -0.5d0 /)
  real(DP), parameter :: sn_snap(6) = (/ 0.d0,  0.5d0*rt3, 0.5d0*rt3,  0.d0, -0.5d0*rt3, -0.5d0*rt3 /)
! f-radial helpers
  real(DP) :: rz2

  integer :: kz, ig, ign, igphi, &
             indexr, iz, lb, ir, nmesh, nmeshs, tblm(4)
  real(DP), parameter :: eps=1.d-8
  complex(DP), parameter :: cim=(0.d0, 1.d0)
  real(DP) :: gn, s1, s2, s3, s4, cs, sn, cs2, sn2, cs3, sn3, rz, dz1, zr, &
                   dr, z0, dz,  bessj, taunew(4), r(ndmx),         &
                   rab(ndmx), betar(ndmx)
  real(DP), allocatable :: x1(:), x2(:), x3(:), x4(:), x5(:), x6(:)
  real(DP), allocatable :: fx1(:), fx2(:), fx3(:), fx4(:), fx5(:), fx6(:), zsl(:)
  real(DP), allocatable :: gn_ig(:)  ! store gn per ig for Γ-point guard
  complex(DP) :: w0(nz1, ngper, 7)
  complex(DP), allocatable :: wadd(:,:), wadd2(:,:), wadd3(:,:)
  complex(DP) :: t1, t2, t3, t4, t5, t6, t7, wa1, wa2, wa3

  contains
!
! Small-x stabilized cylindrical Bessel for k=0..3 (used in f-only integrands)
  pure real(DP) function jstab(k,x) result(j)
     implicit none
     integer,  intent(in) :: k
     real(DP), intent(in) :: x
     if (x < xtol) then
        select case (k)
        case (0); j = 1.d0                              ! J0 ≈ 1
        case (1); j = 0.5d0*x                          ! J1 ≈ x/2
        case (2); j = 0.125d0*x*x                      ! J2 ≈ x^2/8
        case (3); j = (x*x*x)/48.d0                    ! J3 ≈ x^3/48
        case default; j = 0.d0
        end select
     else
        j = bessj(k,x)
     end if
  end function jstab

  ! ====== start of body ======

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
  allocate( gn_ig( ngper ) )

  ! open unified f-safety log once
  if (.not. flog_opened) then
     open(newunit=flog_unit, file='four_f_safety.log', status='unknown', &
 &        position='append', action='write')
     write(flog_unit,'(a)') '---- four f-safety log ----'
     write(flog_unit,'(a)') 'TAGS: GZERO_GUARD, ANGLE_NORM, ANGLE_SNAP, RZ_DOMAIN, BESSEL_SMALLX'
     write(flog_unit,'(a)') 'FORMAT per tag: problem | correction | location/values'
     flog_opened = .true.
  end if

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
            if (lb.eq.3) then
               ! domain guard: ensure r^2 - z^2 >= 0 (skip violated points)
               rz2 = r(ir)**2 - zsl(kz)**2
               if (rz2 <= 0.d0) then
                  x1(ir)=0.d0; x2(ir)=0.d0; x3(ir)=0.d0
                  x4(ir)=0.d0; x5(ir)=0.d0; x6(ir)=0.d0
                  write(flog_unit,'(a,1x,i8,1x,i8,1x,i8,1x,es12.4,1x,es12.4)') &
      &              'RZ_DOMAIN | set x1..x6=0 | ig,kz,ir,r,z=', ig, kz, ir, r(ir), zsl(kz)
                  cycle
               end if
               rz = sqrt( rz2 )
               ! replace bessj→jstab for stability at tiny x
               if (gn*rz < xtol) then
                  write(flog_unit,'(a,1x,i8,1x,i8,1x,i8,1x,es12.4)') &
      &              'BESSEL_SMALLX | use series Jk | ig,kz,ir,x=', ig, kz, ir, gn*rz
               end if
               x1(ir)=betar(ir)*jstab(3,gn*rz)*rz**3/r(ir)**3
               x2(ir)=betar(ir)*jstab(2,gn*rz)*rz**2/r(ir)**3
               x3(ir)=betar(ir)*jstab(1,gn*rz)*rz/r(ir)**3
               x4(ir)=betar(ir)*jstab(1,gn*rz)*rz**3/r(ir)**3
               x5(ir)=betar(ir)*jstab(0,gn*rz)/r(ir)**3
               x6(ir)=betar(ir)*jstab(0,gn*rz)*rz**2/r(ir)**3
            else
               rz=sqrt(r(ir)**2-zsl(kz)**2)
            endif
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
               ! Already handled above
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
         endif
         if (lb.eq.3) then
            fx2(kz)=fx2(kz)+x2(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x3(iz),rab(iz),fx3(kz))
            call simpson(nmeshs-iz+1,x4(iz),rab(iz),fx4(kz))
            fx3(kz)=fx3(kz)+x3(iz)*0.5d0*zr
            fx4(kz)=fx4(kz)+x4(iz)*0.5d0*zr
            call simpson(nmeshs-iz+1,x5(iz),rab(iz),fx5(kz))
            call simpson(nmeshs-iz+1,x6(iz),rab(iz),fx6(kz))
            if(iz.eq.1) then
               x5(iz-1)=0.d0
               x6(iz-1)=0.d0
            else
               x5(iz-1)=(betar(iz)-(betar(iz)-betar(iz-1))/dr*zr)/(abs(zsl(kz))**3)
               x6(iz-1)=0.d0
            endif
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
        gn_ig(ig) = gn  ! store gn for this ig
        if (gn.gt.eps) then
          cs=gper(1,ig)*tpiba/gn
          sn=gper(2,ig)*tpiba/gn
        else
          cs=0.d0
          sn=0.d0
        endif
        cs2=cs**2-sn**2
        sn2=2*cs*sn
        if (lb.eq.3) then
           ! ---- f-only: normalize (cs,sn) and rebuild 2φ,3φ ----
           den = sqrt(cs*cs + sn*sn)
           if (den > 0.d0) then
              cs0 = cs; sn0 = sn
              cs  = cs/den; sn = sn/den
              if (abs(den-1.d0) > angtol) then
                 write(flog_unit,'(a,1x,i8,1x,es12.4,1x,2(1x,es12.4),1x,es12.4)') &
     &               'ANGLE_NORM | normalized cs,sn | ig,gn,cs0,sn0,den-1=', ig, gn, cs0, sn0, den-1.d0
              end if
           else
              cs = 1.d0; sn = 0.d0
              write(flog_unit,'(a,1x,i8,1x,es12.4)') &
     &            'ANGLE_NORM | undefined angle -> set cs=1,sn=0 | ig,gn=', ig, gn
           end if
           ! optional: snap to nearest hex/trig symmetry if extremely close (LaF3-friendly)
           bestd = 1.d99; isnap = 0
           do i=1,6
              dc = cs - cs_snap(i); ds = sn - sn_snap(i)
              d2 = dc*dc + ds*ds
              if (d2 < bestd) then
                 bestd = d2; isnap = i
              end if
           end do
           if (bestd < symtol) then
              tcs = cs; tsn = sn
              cs  = cs_snap(isnap); sn = sn_snap(isnap)
              write(flog_unit,'(a,1x,i8,1x,es12.4,1x,4(1x,es12.4),1x,i1)') &
     &            'ANGLE_SNAP | snapped to hex angle | ig,gn,bestd,cs0,sn0,cs,sn,isnap=', &
     &             ig, gn, bestd, tcs, tsn, cs, sn, isnap
           end if
           cs2 = cs*cs - sn*sn
           sn2 = 2.d0*cs*sn
           cs3 = cs * (4.d0*cs**2 - 3.d0)
           sn3 = sn * (3.d0 - 4.d0*sn**2)
        endif

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
     s1 = tpi/sarea * sqrt(35.d0/(fpi*4.d0))
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
        ! ---- f-only G≈0 guard: only m=0 survives; zero |m|>0 and short-circuit ----
        if (gn_ig(ig) < gtol) then
           write(flog_unit,'(a,1x,i8,1x,i8,1x,es12.4,1x,es12.4)') &
     &         'GZERO_GUARD | zero |m|>0; assemble m=0 | ig,kz,gn,z=', ig, kz, gn_ig(ig), zsl(kz)
           w0(kz,ig,2)=0.d0
           w0(kz,ig,3)=0.d0
           w0(kz,ig,4)=0.d0
           w0(kz,ig,5)=0.d0
           w0(kz,ig,6)=0.d0
           w0(kz,ig,7)=0.d0
           t1=w0(kz,ig,1);wa1=wadd(kz,ig)
           w0(kz,ig,1)=s4*(2.d0*zsl(kz)**3*t1-3.d0*zsl(kz)*wa1)
           cycle
        end if

        t1=w0(kz,ig,1);wa1=wadd(kz,ig)
        t2=w0(kz,ig,2);wa2=wadd2(kz,ig)
        t3=w0(kz,ig,3);wa3=wadd3(kz,ig)
        t4=w0(kz,ig,4);t5=w0(kz,ig,5)
        t6=w0(kz,ig,6);t7=w0(kz,ig,7)
        w0(kz,ig,1)=s4*(2.d0*zsl(kz)**3*t1-3.d0*zsl(kz)*wa1)
        w0(kz,ig,2)=cim*s3*(4.d0*zsl(kz)**2*t2-wa2)
        w0(kz,ig,3)=cim*s3*(4.d0*zsl(kz)**2*t3-wa3)
        w0(kz,ig,4)=-s2*zsl(kz)*t5
        w0(kz,ig,5)=-s2*zsl(kz)*t4
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
  deallocate(gn_ig)

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
