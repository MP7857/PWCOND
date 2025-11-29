!-----------------------------------------------------------------------
! four.f90
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
! This routine computes the bidimensional Fourier transform of the
! beta function. It has been implemented for s, p, d, f-orbitals.
!
!   w0(z,g,m) = 1/S * ∫ w(r) exp{-i g · r_⊥} d r_⊥
!   where w(r) is the beta function of the alpha's orbital.
!
! For a fixed l it computes w0 for all m.
!
! Order of spherical harmonics used:
!   s ;
!   p_z, p_{-x}, p_{-y} ;
!   d_{z^2-1}, d_{-xz}, d_{-yz}, d_{x^2-y^2}, d_{xy}
!   f_{z(5z^2-3r^2)}, f_{x(5z^2-r^2)}, f_{-y(5z^2-r^2)},
!   f_{z(x^2-y^2)}, f_{-xyz}, f_{x(x^2-3y^2)}, f_{-y(3x^2-y^2)}
!
! Input:
!   tblm   - array characterizing the orbital (QE convention)
!   taunew - (x,y,z,r_c) of the projector center (Bohr units / alat-scaled)
!   z0     - initial z
!   dz     - slab width
!
! Output:
!   w0(z, g, m), where z0 < z < z0 + dz, g is 2D G-vector
!
!-----------------------------------------------------------------------
subroutine four( w0, z0, dz, tblm, taunew, r, rab, betar )

  use kinds,       only : DP
  use constants,   only : tpi, fpi
  use radial_grids,only : ndmx
  use cell_base,   only : alat, tpiba
  use cond,        only : sarea, nz1, ngper, gper, ninsh, gnsh, ngpsh

  implicit none

  ! Arguments
  complex(DP) :: w0(nz1, ngper, 7)
  real(DP)    :: z0, dz
  integer     :: tblm(4)
  real(DP)    :: taunew(4)
  real(DP)    :: r(ndmx), rab(ndmx), betar(ndmx)

  ! Locals
  integer :: kz, ig, ign, igphi
  integer :: indexr, iz, lb, ir, nmesh, nmeshs
  real(DP), parameter :: eps = 1.d-8
  complex(DP), parameter :: cim = (0.d0, 1.d0)

  real(DP) :: gn, s1, s2, s3, s4
  real(DP) :: cs, sn, cs2, sn2, cs3, sn3
  real(DP) :: rz, dz1, zr, dr, bessj
  real(DP) :: zsl_val

  real(DP), allocatable :: x1(:), x2(:), x3(:), x4(:), x5(:), x6(:)
  real(DP), allocatable :: fx1(:), fx2(:), fx3(:), fx4(:), fx5(:), fx6(:)
  real(DP), allocatable :: zsl(:)

  complex(DP), allocatable :: wadd(:,:), wadd2(:,:), wadd3(:,:)
  complex(DP) :: t1, t2, t3, t4, t5, t6, t7, wa1, wa2, wa3


  ! Allocate work arrays
  allocate( x1(0:ndmx) )
  allocate( x2(0:ndmx) )
  allocate( x3(0:ndmx) )
  allocate( x4(0:ndmx) )
  allocate( x5(0:ndmx) )
  allocate( x6(0:ndmx) )

  allocate( fx1(nz1) )
  allocate( fx2(nz1) )
  allocate( fx3(nz1) )
  allocate( fx4(nz1) )
  allocate( fx5(nz1) )
  allocate( fx6(nz1) )

  allocate( zsl(nz1) )

  allocate( wadd(nz1, ngper) )
  allocate( wadd2(nz1, ngper) )
  allocate( wadd3(nz1, ngper) )

  ! Angular momentum
  lb = tblm(3)

  ! Radial mesh up to projector cutoff
  nmesh = indexr(taunew(4)*alat, ndmx, r)

  ! Z-slab discretization
  dz1      = dz / nz1
  zsl(1)   = ( z0 + dz1*0.5d0 - taunew(3) ) * alat

  do kz = 2, nz1
     zsl(kz) = zsl(kz-1) + dz1*alat
  enddo

  !--------------------------------------------------------
  ! Loop over G shells (gnsh): gn = |G_⊥|
  !--------------------------------------------------------
  ig = 0
  do ign = 1, ngpsh

     gn = gnsh(ign)

     !-----------------------------------------------
     ! For each z-slice kz, compute radial integrals
     !-----------------------------------------------
     do kz = 1, nz1

        zsl_val = zsl(kz)

        if ( abs(zsl_val) + eps .le. taunew(4)*alat ) then

           ! Find first radial index with r > |z|
           iz = indexr( zsl_val, nmesh, r )

           ! Ensure (nmeshs - iz + 1) is odd (Simpson)
           if ( ( (nmesh-iz)/2*2 ) .eq. (nmesh-iz) ) then
              nmeshs = nmesh
           else
              nmeshs = nmesh + 1
           endif

           ! Build integrands on [iz..nmeshs]
           do ir = iz, nmeshs
              rz = sqrt( r(ir)**2 - zsl_val**2 )
              select case ( lb )
              case (0)
                 x1(ir) = betar(ir) * bessj(0, gn*rz)
              case (1)
                 x1(ir) = betar(ir) * bessj(1, gn*rz) * rz / r(ir)
                 x2(ir) = betar(ir) * bessj(0, gn*rz) / r(ir)
              case (2)
                 x1(ir) = betar(ir) * bessj(2, gn*rz) * rz**2 / r(ir)**2
                 x2(ir) = betar(ir) * bessj(1, gn*rz) * rz    / r(ir)**2
                 x3(ir) = betar(ir) * bessj(0, gn*rz)         / r(ir)**2
                 x4(ir) = betar(ir) * bessj(0, gn*rz)
              case (3)
                 x1(ir) = betar(ir) * bessj(3, gn*rz) * rz**3 / r(ir)**3
                 x2(ir) = betar(ir) * bessj(2, gn*rz) * rz**2 / r(ir)**3
                 x3(ir) = betar(ir) * bessj(1, gn*rz) * rz    / r(ir)**3
                 x4(ir) = betar(ir) * bessj(1, gn*rz) * rz**3 / r(ir)**3
                 x5(ir) = betar(ir) * bessj(0, gn*rz)         / r(ir)**3
                 x6(ir) = betar(ir) * bessj(0, gn*rz) * rz**2 / r(ir)**3
              case default
                 call errore( 'four', 'l not programmed', 1 )
              end select
           enddo

           ! Main Simpson integral on [r(iz)..r(nmeshs)]
           call simpson( nmeshs-iz+1, x1(iz), rab(iz), fx1(kz) )

           ! Radial spacing at the lower endpoint
           if ( iz .eq. 1 ) then
              dr = r(iz)
           else
              dr = r(iz) - r(iz-1)
           endif

           ! Distance from |z| to r(iz)
           zr = r(iz) - abs(zsl_val)

           if ( lb .eq. 0 ) then

              ! s-channel: extrapolate beta to r = |z| and integrate [|z|, r(iz)]
              if ( iz .eq. 1 ) then
                 x1(iz-1) = betar(iz) - betar(iz)/dr*zr
              else
                 x1(iz-1) = betar(iz) - (betar(iz)-betar(iz-1))/dr*zr
              endif
              fx1(kz) = fx1(kz) + ( x1(iz-1) + x1(iz) ) * 0.5d0 * zr

           else

              ! For l>0: add [|z|,r(iz)] contribution using trapezoid and
              ! keep separate integral of the second integrand (x2, x3, x4 …)
              fx1(kz) = fx1(kz) + x1(iz) * 0.5d0 * zr
              call simpson( nmeshs-iz+1, x2(iz), rab(iz), fx2(kz) )

           endif

           ! p-channel additional endpoint correction
           if ( lb .eq. 1 ) then
              if ( iz .eq. 1 ) then
                 x2(iz-1) = 0.d0
              else
                 x2(iz-1) = ( betar(iz) - (betar(iz)-betar(iz-1))/dr*zr ) / abs(zsl_val)
              endif
              fx2(kz) = fx2(kz) + ( x2(iz-1) + x2(iz) ) * 0.5d0 * zr
           endif

           ! d-channel additional integrals and endpoint correction
           if ( lb .eq. 2 ) then

              fx2(kz) = fx2(kz) + x2(iz) * 0.5d0 * zr
              call simpson( nmeshs-iz+1, x3(iz), rab(iz), fx3(kz) )
              call simpson( nmeshs-iz+1, x4(iz), rab(iz), fx4(kz) )

              if ( iz .eq. 1 ) then
                 x3(iz-1) = 0.d0
                 x4(iz-1) = 0.d0
              else
                 x3(iz-1) = ( betar(iz) - (betar(iz)-betar(iz-1))/dr*zr ) / abs(zsl_val)**2
                 x4(iz-1) =   betar(iz) - (betar(iz)-betar(iz-1))/dr*zr
              endif

              fx3(kz) = fx3(kz) + ( x3(iz-1) + x3(iz) ) * 0.5d0 * zr
              fx4(kz) = fx4(kz) + ( x4(iz-1) + x4(iz) ) * 0.5d0 * zr

           ! f-channel
           elseif ( lb .eq. 3 ) then

              fx2(kz) = fx2(kz) + x2(iz) * 0.5d0 * zr

              call simpson( nmeshs-iz+1, x3(iz), rab(iz), fx3(kz) )
              fx3(kz) = fx3(kz) + x3(iz) * 0.5d0 * zr

              call simpson( nmeshs-iz+1, x4(iz), rab(iz), fx4(kz) )
              fx4(kz) = fx4(kz) + x4(iz) * 0.5d0 * zr

              call simpson( nmeshs-iz+1, x5(iz), rab(iz), fx5(kz) )
              call simpson( nmeshs-iz+1, x6(iz), rab(iz), fx6(kz) )

              if ( iz .eq. 1 ) then
                 x5(iz-1) = 0.d0
              else
                 x5(iz-1) = ( betar(iz) - (betar(iz)-betar(iz-1))/dr*zr ) / (abs(zsl_val)**3)
              endif
              x6(iz-1) = 0.d0

              fx5(kz) = fx5(kz) + ( x5(iz-1) + x5(iz) ) * 0.5d0 * zr
              fx6(kz) = fx6(kz) + ( x6(iz-1) + x6(iz) ) * 0.5d0 * zr

           endif

        else
           ! Outside projector radius: all integrals vanish
           fx1(kz) = 0.d0
           fx2(kz) = 0.d0
           fx3(kz) = 0.d0
           fx4(kz) = 0.d0
           fx5(kz) = 0.d0
           fx6(kz) = 0.d0
        endif

     enddo  ! kz loop

     !---------------------------
     ! φ-dependence: distribute
     !---------------------------
     do igphi = 1, ninsh(ign)

        ig = ig + 1

        if ( gn .gt. eps ) then
           cs = gper(1,ig)*tpiba / gn
           sn = gper(2,ig)*tpiba / gn
        else
           cs = 0.d0
           sn = 0.d0
        endif

        cs2 = cs**2 - sn**2
        sn2 = 2.d0*cs*sn
        cs3 = cs * (4.d0*cs**2 - 3.d0)
        sn3 = sn * (3.d0 - 4.d0*sn**2)

        do kz = 1, nz1
           select case ( lb )
           case (0)
              w0(kz,ig,1) = fx1(kz)

           case (1)
              w0(kz,ig,2) = cs * fx1(kz)
              w0(kz,ig,1) = fx2(kz)
              w0(kz,ig,3) = sn * fx1(kz)

           case (2)
              w0(kz,ig,5) = sn2 * fx1(kz)
              w0(kz,ig,2) = cs  * fx2(kz)
              w0(kz,ig,1) = fx3(kz)
              w0(kz,ig,3) = sn  * fx2(kz)
              w0(kz,ig,4) = cs2 * fx1(kz)
              wadd(kz,ig) = fx4(kz)

           case (3)
              w0(kz,ig,1) = fx5(kz)
              wadd(kz,ig)  = fx6(kz)
              w0(kz,ig,2) = cs * fx3(kz)
              wadd2(kz,ig) = cs * fx4(kz)
              w0(kz,ig,3) = sn * fx3(kz)
              wadd3(kz,ig) = sn * fx4(kz)
              w0(kz,ig,4) = cs2 * fx2(kz)
              w0(kz,ig,5) = sn2 * fx2(kz)
              w0(kz,ig,6) = cs3 * fx1(kz)
              w0(kz,ig,7) = sn3 * fx1(kz)
           end select
        enddo

     enddo  ! igphi

  enddo  ! ign

  !--------------------------------------
  ! Normalization / spherical harmonics
  !--------------------------------------
  if ( lb .eq. 0 ) then
     s1 = tpi / sarea / sqrt(fpi)
  elseif ( lb .eq. 1 ) then
     s1 = tpi / sarea * sqrt(3.d0/fpi)
  elseif ( lb .eq. 2 ) then
     s1 = -tpi / (2.d0*sarea) * sqrt(15.d0/fpi)
     s2 =  tpi / sarea        * sqrt(5.d0/tpi/8.d0)
  elseif ( lb .eq. 3 ) then
     s1 = tpi / sarea * sqrt(35.d0/(fpi*8.d0))
     s2 = tpi / sarea * sqrt(105.d0/(fpi*4.d0))
     s3 = tpi / sarea * sqrt(21.d0/(fpi*8.d0))
     s4 = tpi / sarea * sqrt(7.d0/(fpi*4.d0))
  endif

  do ig = 1, ngper
     do kz = 1, nz1
        select case ( lb )
        case (0)
           w0(kz,ig,1) = s1 * w0(kz,ig,1)

        case (1)
           w0(kz,ig,2) =  cim * s1 * w0(kz,ig,2)
           w0(kz,ig,1) =        s1 * zsl(kz) * w0(kz,ig,1)
           w0(kz,ig,3) =  cim * s1 * w0(kz,ig,3)

        case (2)
           w0(kz,ig,5) = s1 * w0(kz,ig,5)
           w0(kz,ig,2) = -2.d0*cim*s1*zsl(kz) * w0(kz,ig,2)
           w0(kz,ig,1) = 3.d0*zsl(kz)**2*s2*w0(kz,ig,1) - s2*wadd(kz,ig)
           w0(kz,ig,3) = -2.d0*cim*s1*zsl(kz) * w0(kz,ig,3)
           w0(kz,ig,4) = s1 * w0(kz,ig,4)

        case (3)
           t1  = w0(kz,ig,1);  wa1 = wadd(kz,ig)
           t2  = w0(kz,ig,2);  wa2 = wadd2(kz,ig)
           t3  = w0(kz,ig,3);  wa3 = wadd3(kz,ig)
           t4  = w0(kz,ig,4);  t5  = w0(kz,ig,5)
           t6  = w0(kz,ig,6);  t7  = w0(kz,ig,7)

           w0(kz,ig,1) =  s4 * ( 2.d0*zsl(kz)**3*t1 - 3.d0*zsl(kz)*wa1 )
           w0(kz,ig,2) =  cim*s3 * ( 4.d0*zsl(kz)**2*t2 - wa2 )
           w0(kz,ig,3) =  cim*s3 * ( 4.d0*zsl(kz)**2*t3 - wa3 )
           w0(kz,ig,4) = -s2*zsl(kz)*t4
           w0(kz,ig,5) = -s2*zsl(kz)*t5
           w0(kz,ig,6) = -cim*s1*t6
           w0(kz,ig,7) = -cim*s1*t7
        end select
     enddo
  enddo

  ! Deallocate work arrays
  deallocate( x1 )
  deallocate( x2 )
  deallocate( x3 )
  deallocate( x4 )
  deallocate( x5 )
  deallocate( x6 )

  deallocate( fx1 )
  deallocate( fx2 )
  deallocate( fx3 )
  deallocate( fx4 )
  deallocate( fx5 )
  deallocate( fx6 )

  deallocate( zsl )

  deallocate( wadd )
  deallocate( wadd2 )
  deallocate( wadd3 )

  return
end subroutine four
!-----------------------------------------------------------------------

function indexr( zz, ndim, r ) result( idx )
  use kinds, only : DP
  implicit none

  integer :: ndim, idx
  real(DP) :: zz, r(ndim)
  integer :: iz

  ! Find first index such that abs(zz) < r(idx)
  iz = 1
  do while ( r(iz) .le. abs(zz) + 1.d-10 )
     iz = iz + 1
  enddo
  idx = iz

end function indexr
