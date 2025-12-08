! =====================================================================
!  Quantum ESPRESSO – PWCOND
!  nz-stable integration (Simpson-normalized)
!  M. Pourfath, 2025
! =====================================================================

function int1d(fun, zk, dz, dz1, nz1, tpiba, sign)
  USE kinds, only : DP
  implicit none
  integer :: ik, nz1, sign
  real(DP), parameter :: eps = 1.d-8
  real(DP) :: tpi, dz, dz1, tpiba
  real(DP), allocatable :: wz(:)
  real(DP) :: Wsum, scale
  complex(DP), parameter :: cim = (0.d0,1.d0)
  complex(DP) :: zk, fun(nz1), fact, fact0, arg, int1d

  tpi = 8.d0*atan(1.d0)

  ! Simpson/trapezoid weights
  allocate( wz(nz1) )
  call z_simpson_weights(nz1, wz)

  ! Normalization: QE formula assumes flat rule — fix with scale
  Wsum  = sum(wz)
  scale = Wsum / real(nz1,DP)

  int1d = (0.d0,0.d0)
  arg   = sign*tpi*cim*zk*dz1
  fact0 = exp(arg)
  fact  = (1.d0,0.d0)

  do ik = 1, nz1
     int1d = int1d + conjg(fun(ik)) * fact * wz(ik)
     fact   = fact * fact0
  enddo

  ! Apply normalization scale (no global /3)
  int1d = int1d * scale

  if (abs(real(zk))+abs(aimag(zk)) > eps) then
     int1d = -sign*cim * int1d * (1.d0 - exp(-arg)) / (zk*tpiba)
     if (sign.lt.0) int1d = int1d * exp(tpi*cim*zk*dz)
  else
     int1d = int1d * dz1/tpiba * tpi
  endif

  deallocate(wz)
  return
end function int1d
! =====================================================================

function int2d(fun1, fun2, int1, int2, fact1, fact2, zk, dz1, tpiba, nz1)
  USE kinds, only : DP
  USE constants, ONLY : tpi
  implicit none
  integer :: ik, nz1
  real(DP), parameter :: eps = 1.d-8
  real(DP) :: dz1, tpiba
  real(DP), allocatable :: wz(:)
  real(DP) :: Wsum, scale
  complex(DP), parameter :: cim = (0.d0,1.d0), one = (1.d0,0.d0)
  complex(DP) :: fun1(nz1), fun2(nz1)
  complex(DP) :: int1(nz1), int2(nz1)
  complex(DP) :: fact1(nz1), fact2(nz1)
  complex(DP) :: s1, s2, s3, ff, fact, fact0, f1, f2, zk, int2d

  allocate( wz(nz1) )
  call z_simpson_weights(nz1, wz)

  Wsum  = sum(wz)
  scale = Wsum / real(nz1,DP)

  s1 = (0.d0,0.d0)
  s2 = (0.d0,0.d0)
  s3 = (0.d0,0.d0)

  fact  = fact1(1)
  fact0 = fact2(1)

  do ik = 1, nz1
     ff = conjg(fun1(ik))
     s1 = s1 + int1(ik)*ff*fact1(ik)*wz(ik)
     s2 = s2 + int2(ik)*ff*fact2(ik)*wz(ik)
     s3 = s3 + fun2(ik)*ff*wz(ik)
  enddo

  ! Apply normalization
  s1 = s1 * scale
  s2 = s2 * scale
  s3 = s3 * scale

  f1 = cim*zk*dz1*tpi
  f2 = one/(zk*tpiba)**2

  if (abs(f1) > eps) then
     int2d = ((1.d0-fact+f1)*s3*2.d0 + (2.d0-fact-fact0)*(s1+s2)) * f2
  else
     int2d = (s1+s2+s3)*(dz1*tpi/tpiba)**2
  endif

  deallocate(wz)
  return
end function int2d

! =====================================================================

subroutine setint(fun,int1,int2,fact1,fact2,nz1)
  USE kinds, only : DP
  implicit none
  integer :: ik, nz1
  complex(DP) :: fun(nz1), int1(nz1), int2(nz1), fact1(nz1), fact2(nz1)

  int1(1)   = (0.d0,0.d0)
  int2(nz1) = (0.d0,0.d0)

  do ik = 2, nz1
     int1(ik) = int1(ik-1) + fun(ik-1)*fact2(ik-1)
  enddo

  do ik = nz1-1, 1, -1
     int2(ik) = int2(ik+1) + fun(ik+1)*fact1(ik+1)
  enddo
end subroutine setint

! =====================================================================

subroutine z_simpson_weights(nz1, w)
  USE kinds, only : DP
  implicit none
  integer, intent(in) :: nz1
  real(DP), intent(out) :: w(nz1)
  integer :: k

  if (nz1 < 3 .or. mod(nz1,2)==0) then
     if (nz1 == 1) then
        w(1) = 1.d0
     else
        w(1) = 0.5d0
        do k = 2, nz1-1
           w(k) = 1.d0
        enddo
        w(nz1) = 0.5d0
     endif
     return
  endif

  w(1) = 1.d0
  do k = 2, nz1-1
     if (mod(k,2)==0) then
        w(k) = 4.d0
     else
        w(k) = 2.d0
     endif
  enddo
  w(nz1) = 1.d0
end subroutine z_simpson_weights
