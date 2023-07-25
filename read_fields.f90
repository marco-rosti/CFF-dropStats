#include "config.h"
#include "configDropStats.h"
module mod_read_fields
implicit none
contains
subroutine read_fields(nstep,error)
use, intrinsic :: ISO_FORTRAN_ENV
use mod_common_data,   only : vel, rank, top, kur, nor, car
use mod_param,   only : nxt, nyt, nzt, dxi,dyi,dzi
use mod_flood_fill,only : flood_fill
#ifdef _MULTIPHASE_SURFACTANT_PF_
use mod_common_data,   only : psi, interpolateCC_f2c
use mod_param,   only : gridfx, gridfy, gridfz
#endif
implicit none
logical :: check
integer(kind=int64), intent(in) :: nstep
integer, intent(out) :: error
integer(kind=int64) :: filesize, expectedSize, dropVol
integer :: i,j,k, ip,jp,kp,im,jm,km
character(len=200) :: namedir, filename
character(len=12) :: iStepChar
character(len=5) :: rankChar
double precision, allocatable, dimension(:,:,:) :: phi
double precision :: modnor
#ifdef _MULTIPHASE_SURFACTANT_PF_
double precision, allocatable, dimension(:,:,:) :: psiFine
#endif
error=0
namedir='../../data/'
write(iStepChar,'(i12.12)') nstep
write(rankChar,'(i5.5)') rank
! read in physical space
#ifdef _MULTIPHASE_NNF_EVP_
filename=trim(namedir)//'dataEVPSolid'//iStepChar
#else
filename=trim(namedir)//'dataVOF'//iStepChar
#endif
allocate(phi(nxt,nyt,nzt))
inquire(file=trim(filename),exist=check,size=filesize)
expectedSize = 1.0*nxt*nyt*nzt*1*8 !IC cast to real then long integer to prevent overflow
if(rank.eq.0) then
  write(6,*) 'rank ',rank,'reading ',trim(filename)
endif
#ifdef _READ_FROM_BUCKET_    
if((check.eqv..true.)) then 
#else
if((check.eqv..true.).and.(filesize.eq.expectedSize)) then ! input file contains vof double precision
#endif
  open(66,file=filename,access='stream',status='old',form='unformatted',convert='little_endian')
  read(66) phi
  close(66,status='keep')
else
  write(6,*) 'rank ',rank,'wrong file ',trim(filename),'exist ',check,'filesize',filesize,'expect ',expectedSize
  error=1
  return
endif

dropVol = 0 !IC
do k=1,nzt
  do j=1,nyt
    do i=1,nxt
#ifdef _MULTIPHASE_NNF_EVP_
      if(phi(i,j,k).lt.0.01d0)then
#else
      if(phi(i,j,k).ge.0.5d0)then
#endif
        top(i,j,k)=1
        dropVol = dropVol+1 !IC
      else
        top(i,j,k)=0
      endif
!!G Soligo 2023 Mar: not parallel, very non-optimal, but read is single core... Kur will be in shared window
      !compute cell-centred normal
      ip=i+1
      jp=j+1
      kp=k+1
      if(ip.gt.nxt) ip=ip-nxt
      if(jp.gt.nyt) jp=jp-nyt
      if(kp.gt.nzt) kp=kp-nzt
      im=i-1
      jm=j-1
      km=k-1
      if(im.lt.1) im=im+nxt
      if(jm.lt.1) jm=jm+nyt
      if(km.lt.1) km=km+nzt
      nor(i,j,k,1)=(phi(ip,j,k)-phi(im,j,k))*(0.5d0*dxi)
      nor(i,j,k,2)=(phi(i,jp,k)-phi(i,jm,k))*(0.5d0*dyi)
      nor(i,j,k,3)=(phi(i,j,kp)-phi(i,j,km))*(0.5d0*dzi)
      modnor=dsqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2)
      !outward pointing normal
      nor(i,j,k,:)=-nor(i,j,k,:)/modnor
    enddo
  enddo
enddo
if(rank.eq.0)then
  open(44,file='./output/dropVol.dat',access='append',form='formatted')
#ifdef _MULTIPHASE_SURFACTANT_PF_
    write(44,'(2i16,2E15.7)') nstep,dropVol,dropVol/(nxt*nyt*nzt*1.0),sum(psi)/(nxt*nyt*nzt*1.0)
#else
    write(44,'(2i16,E15.7)') nstep,dropVol,dropVol/(nxt*nyt*nzt*1.0)
#endif
  close(44,status='keep')
endif
dropVol=0
do k=1,nzt
  do j=1,nyt
    do i=1,nxt
      !find carrier fluid so we can skip it when looking for voids in get_interface. Assume carrier>25%
      if (dropVol.lt.nxt*nyt*nzt/4.and.top(i,j,k).eq.0) then
        car=0
        dropVol=0
        call flood_fill(i,j,k,0,dropVol,car)
        if (rank.eq.0) write(*,*) 'Carrier fraction',1.0*dropVol/nxt/nyt/nzt
      endif
      im=i-1
      jm=j-1
      km=k-1
      if(im.lt.1) im=im+nxt
      if(jm.lt.1) jm=jm+nyt
      if(km.lt.1) km=km+nzt
      ip=i+1
      jp=j+1
      kp=k+1
      if(ip.gt.nxt) ip=ip-nxt
      if(jp.gt.nyt) jp=jp-nyt
      if(kp.gt.nzt) kp=kp-nzt
      !compute curvature in backward direction
      kur(i,j,k)=(nor(ip,j,k,1)-nor(im,j,k,1))*(0.5d0*dxi)+ &
                 (nor(i,jp,k,2)-nor(i,jm,k,2))*(0.5d0*dyi)+ &
                 (nor(i,j,kp,3)-nor(i,j,km,3))*(0.5d0*dzi)
    enddo
  enddo
enddo
deallocate(phi)
#ifdef _MULTIPHASE_SURFACTANT_PF_
filename=trim(namedir)//'dataSURF_PF'//iStepChar
inquire(file=trim(filename),exist=check,size=filesize)
expectedSize = 1.0*nxt*gridfx*nyt*gridfy*nzt*gridfz*1*8 !IC cast to real then long integer to prevent overflow
if(rank.eq.0) then
  write(6,*) 'rank ',rank,'reading ',trim(filename)
endif
#ifdef _READ_FROM_BUCKET_    
if((check.eqv..true.)) then 
#else
if((check.eqv..true.).and.(filesize.eq.expectedSize)) then ! input file contains vof double precision
#endif
allocate(psiFine(nxt*gridfx,nyt*gridfy,nzt*gridfz))
open(67,file=filename,access='stream',status='old',form='unformatted',convert='little_endian')
  read(67) psiFine
close(67,status='keep')
call interpolateCC_f2c(psiFine,psi)
deallocate(psiFine)
else
  write(6,*) 'rank ',rank,'wrong file ',trim(filename),'exist ',check,'filesize',filesize,'expect ',expectedSize
  error=1
  return
endif
#endif
filename=trim(namedir)//'dataUVW'//iStepChar
inquire(file=trim(filename),exist=check,size=filesize)
expectedSize = 1.0*nxt*nyt*nzt*4*8 !IC contains u,v,w,p. Cast to real then long integer to prevent overflow
if(rank.eq.0) then
  write(6,*) 'rank ',rank,'reading ',trim(filename)
endif
#ifdef _READ_FROM_BUCKET_    
if((check.eqv..true.)) then 
#else
if((check.eqv..true.).and.(filesize.eq.expectedSize)) then ! input file contains vof double precision
#endif
open(67,file=filename,access='stream',status='old',form='unformatted',convert='little_endian')
  read(67) vel
close(67,status='keep')
else
  write(6,*) 'rank ',rank,'wrong file ',trim(filename),'exist ',check,'filesize',filesize,'expect ',expectedSize
  error=1
  return
endif
return
end subroutine read_fields
end module mod_read_fields
