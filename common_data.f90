#include "config.h"
#include "configDropStats.h"


module mod_param
  implicit none
  ! Pi
  real, parameter :: pi = ACOS(-1.)
  ! Total number of grid points in X, Y and Z
  integer, parameter :: nxt = 500, nyt = 500, nzt = 500
#ifdef _MULTIPHASE_SURFACTANT_PF_
  ! refinement factor in 3 directions
  integer, parameter :: gridfx=2, gridfy=2, gridfz=2
  integer, parameter :: nxtf=nxt*gridfx, nytf=nyt*gridfy, nztf=nzt*gridfz
#endif
  ! Size of the domain in X, Y and Z
  real, parameter :: lx = 2.0*pi, ly = 2.0*pi, lz = 2.0*pi
  ! Fluid viscosity and density
  real, parameter :: vis = 0.005
  real, parameter :: rho = 1.0
  ! Grid size in X, Y and Z
  real, parameter :: dx  = lx/(1.0*nxt), dy  = ly/(1.0*nyt), dz  = lz/(1.0*nzt)
  ! Inverse of grid size in X, Y and Z
  real, parameter :: dxi = 1.0/dx,       dyi = 1.0/dy,       dzi = 1.0/dz
end module mod_param

module mod_common_data
  use mod_param, only : nxt,nyt,nzt
  use, intrinsic :: ISO_FORTRAN_ENV
  implicit none
  integer :: rank, ntask
  double precision, parameter :: pi=3.14159265358979
  type ragged_array
    !2D array containing vectors of different lengths
    double precision,allocatable::v(:)
    character(len=12) :: iStepChar
  end type ragged_array
#ifdef MPI_WINDOW 
  integer(kind=int8),pointer,dimension(:,:,:) :: top, car
  double precision, pointer, dimension(:,:,:,:) :: vel, nor
  double precision, pointer, dimension(:,:,:) :: kur
#else
  integer(kind=int8),allocatable,dimension(:,:,:) :: top, car
  double precision, allocatable, dimension(:,:,:,:) :: vel, nor
  double precision, allocatable, dimension(:,:,:) :: kur
#endif
#ifdef CURVATURE_PDF
  double precision, pointer, dimension(:,:,:) :: jpdf
  integer, parameter :: nkur=500 , ndim=500
#endif
#ifdef _MULTIPHASE_SURFACTANT_PF_
#ifdef MPI_WINDOW 
  double precision, pointer, dimension(:,:,:) :: psi
#else
  double precision, allocatable, dimension(:,:,:) :: psi
#endif
#endif

contains

#ifdef RAND_SEARCH
REAL FUNCTION RANX (seed)
  !
  ! Used in function forcingIsotropic
  ! **  RANX(0) returns a random floating point number in [0,1[  **
  ! **  initialize with RANX(seed) with a nonzero seed           **
  ! **                                                           **
  ! **  (C) ANZ 1996  (Alain NOULLEZ, anz@obs-nice.fr)           **
  !
  IMPLICIT  NONE
  SAVE      rvs,rjp,rkp
  INTEGER*4 rvs(0:54),lseed
  INTEGER   rjp,rkp,n,seed
  REAL      RANMAX
  PARAMETER (RANMAX = 4294967296.0)
  IF (seed .NE. 0) THEN
     lseed = seed
     rkp = 0
     DO n = 1,55
        lseed = lseed*69069+1
        rvs(rkp) = lseed
        rkp = rkp+21
        IF (rkp .GE. 55) rkp = rkp-55
     ENDDO
     rjp = 24
     DO n = 1,220
        rjp = rjp-1
        IF (rjp .LT. 0) rjp = 54
        rkp = rkp-1
        IF (rkp .LT. 0) rkp = 54
        rvs(rkp) = rvs(rkp)+rvs(rjp)
     ENDDO
  ELSE
     rjp = rjp-1
     IF (rjp .LT. 0) rjp = 54
     rkp = rkp-1
     IF (rkp .LT. 0) rkp = 54
     rvs(rkp) = rvs(rkp)+rvs(rjp)
  ENDIF
  RANX = FLOAT(rvs(rkp))/RANMAX+0.5
  RETURN
END FUNCTION RANX
#endif
#ifdef _MULTIPHASE_SURFACTANT_PF_
subroutine interpolateCC_f2c(varf,var)
  ! Giovanni Soligo 2022
  ! cell-center linear interpolation
  ! need to fill boundaries and halo cells
  use mod_param, only : gridfx,gridfy,gridfz,nxt,nyt,nzt
  implicit none
  double precision, dimension(nxt*gridfx,nyt*gridfy,nzt*gridfz), intent(in) :: varf
  double precision, dimension(nxt,nyt,nzt), intent(out) :: var
  real :: xa,ya,za
  integer :: i,j,k
  integer :: xbound,ybound,zbound
  xa=real(mod(gridfx+1,2))
  ya=real(mod(gridfy+1,2))
  za=real(mod(gridfz+1,2))
  xbound=floor(real(gridfx-1)*0.5+1.0e-5)
  ybound=floor(real(gridfy-1)*0.5+1.0e-5)
  zbound=floor(real(gridfz-1)*0.5+1.0e-5)
  do k=1,nzt
    do j=1,nyt
      do i=1,nxt
        var(i,j,k)=(xa*ya*za*varf(gridfx*i-xbound-1,gridfy*j-ybound-1,gridfz*k-zbound-1)+ &
                       ya*za*varf(gridfx*i-xbound  ,gridfy*j-ybound-1,gridfz*k-zbound-1)+ &
                    xa*   za*varf(gridfx*i-xbound-1,gridfy*j-ybound  ,gridfz*k-zbound-1)+ &
                          za*varf(gridfx*i-xbound  ,gridfy*j-ybound  ,gridfz*k-zbound-1)+ &
                    xa*ya*   varf(gridfx*i-xbound-1,gridfy*j-ybound-1,gridfz*k-zbound  )+ &
                       ya*   varf(gridfx*i-xbound  ,gridfy*j-ybound-1,gridfz*k-zbound  )+ &
                    xa*      varf(gridfx*i-xbound-1,gridfy*j-ybound  ,gridfz*k-zbound  )+ &
                             varf(gridfx*i-xbound  ,gridfy*j-ybound  ,gridfz*k-zbound  ))/ &
                    (xa*ya*za + ya*za + xa*za + za + xa*ya + ya + xa + 1.0)
      enddo
    enddo
  enddo
  return
end subroutine interpolateCC_f2c
#endif
end module mod_common_data
