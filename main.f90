#include "config.h"
#include "configDropStats.h"
program mass_center
use mpi,             only : mpi_init, mpi_comm_size, mpi_comm_world, mpi_comm_rank, mpi_barrier, mpi_finalize, MPI_SUM
use mpi,             only : MPI_MAX_PROCESSOR_NAME, MPI_Get_processor_name, MPI_IN_PLACE, MPI_INT, mpi_allreduce
use mod_param,       only : nxt, nyt, nzt,Lx,Ly,Lz,dxi,dyi,dzi
use mod_common_data, only : rank, ntask, top, car, vel, ragged_array, kur, nor
use paraview,        only : write_paraview
#ifdef CURVATURE_PDF
use mod_common_data, only : jpdf, nkur, ndim
#endif
use mod_read_fields, only : read_fields
use mod_get_interface, only : get_interface
!use mpiIO
use, intrinsic :: ISO_FORTRAN_ENV
#ifdef MPI_WINDOW
use mpi,             only : MPI_INTEGER1, mpi_address_kind, mpi_info_null, mpi_double_precision
use mpi,             only : mpi_win_fence, mpi_win_free, mpi_comm_type_shared, OMPI_COMM_TYPE_NUMA
use mpi,             only : OMPI_COMM_TYPE_L1CACHE
use, intrinsic :: iso_c_binding
#endif
#ifdef _MULTIPHASE_SURFACTANT_PF_
use mod_common_data, only : psi
use mod_param,       only : gridfx, gridfy, gridfz
#endif
implicit none
integer(8) :: iout0d, iout1d, iout2d, iout3d, i, ii
integer(kind=int64) :: begin,finish,step,delay
integer :: delayRep
integer:: ierr,namelen,fileReadErr
character(len=200) :: filename
character(len=12) :: iStepChar
character(len=MPI_MAX_PROCESSOR_NAME) :: procname
#ifdef CURVATURE_PDF
integer :: ip,jp
character(len=100) :: fmtstr
#endif
#ifdef MPI_WINDOW
integer :: disp_unit, shcomm, nodeSize, nodeRank, winTop, winCar, winVel, winkur, winNor
type(c_ptr) :: baseptr
integer(kind=mpi_address_kind) :: extent, lb, varsize_tot
#ifdef _MULTIPHASE_SURFACTANT_PF_
integer :: winPsi
#endif
#endif
filename='input.in'
open(3,file=filename)
! Simulation initial and final step
read(3,*)
read(3,*) begin, finish
! Simulation initial time
read(3,*)
read(3,*) 
! Simulation initial time-step
read(3,*)
read(3,*) 
! Output frequency
read(3,*)
read(3,*) iout0d, iout1d, iout2d, iout3d, step, delay
close(3)
call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_barrier(mpi_comm_world,ierr)
call MPI_Get_processor_name(procname, namelen,ierr)
if(rank.eq.0) then   
  write(*,'(1x,a)') '----------------------------------------------------------------------'
  write(*,'(1x,a)') '                 starting number of drops calculation                       '
  write(*,'(1x,a,3(i5,a))') 'nxt * nyt * nzt = ',nxt,' * ',nyt,' * ',nzt,' '
  write(*,*) 'Reading fields in physical space'
  write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',begin,' to ',finish,' with step ',step
  write(*,'(1x,a)') '----------------------------------------------------------------------'
  write(*,*)
endif
#ifdef MPI_WINDOW
call mpi_comm_split_type(mpi_comm_world,mpi_comm_type_shared,rank,mpi_info_null,shcomm,ierr)
call MPI_Comm_rank(shcomm, nodeRank, ierr)
call MPI_Comm_size(shcomm, nodeSize, ierr)

!Allocate window for top array
if(nodeRank.eq.0)then
 ! only one rank (rank 0) allocates memory, others get zero memory
 call mpi_type_get_extent(MPI_INTEGER1,lb,extent,ierr)
 disp_unit=int(extent)
else
 extent=0
 disp_unit=1
endif
varsize_tot=nxt*nyt*nzt*extent
call mpi_win_allocate_shared(varsize_tot,disp_unit,mpi_info_null,shcomm,baseptr,winTop,ierr)
! get location of memory segment
if(nodeRank.ne.0) call mpi_win_shared_query(winTop,0,varsize_tot,disp_unit,baseptr,ierr)
call c_f_pointer(baseptr,top,[nxt,nyt,nzt])

call mpi_win_allocate_shared(varsize_tot,disp_unit,mpi_info_null,shcomm,baseptr,winCar,ierr)
! get location of memory segment
if(nodeRank.ne.0) call mpi_win_shared_query(winCar,0,varsize_tot,disp_unit,baseptr,ierr)
call c_f_pointer(baseptr,car,[nxt,nyt,nzt])

!Allocate window for vel array
if(nodeRank.eq.0)then
 ! only one rank (rank 0) allocates memory, others get zero memory
 call mpi_type_get_extent(mpi_double_precision,lb,extent,ierr)
 disp_unit=int(extent)
endif
varsize_tot=nxt*nyt*nzt*3*extent
call mpi_win_allocate_shared(varsize_tot,disp_unit,mpi_info_null,shcomm,baseptr,winVel,ierr)
! get location of memory segment
if(nodeRank.ne.0) call mpi_win_shared_query(winVel,0,varsize_tot,disp_unit,baseptr,ierr)
call c_f_pointer(baseptr,vel,[nxt,nyt,nzt,3])

!Allocate window for nor array
call mpi_win_allocate_shared(varsize_tot,disp_unit,mpi_info_null,shcomm,baseptr,winNor,ierr)
! get location of memory segment
if(nodeRank.ne.0) call mpi_win_shared_query(winNor,0,varsize_tot,disp_unit,baseptr,ierr)
call c_f_pointer(baseptr,nor,[nxt,nyt,nzt,3])

!Allocate window for kur array
varsize_tot=nxt*nyt*nzt*extent
call mpi_win_allocate_shared(varsize_tot,disp_unit,mpi_info_null,shcomm,baseptr,winkur,ierr)
! get location of memory segment
if(nodeRank.ne.0) call mpi_win_shared_query(winkur,0,varsize_tot,disp_unit,baseptr,ierr)
call c_f_pointer(baseptr,kur,[nxt,nyt,nzt])

#ifdef _MULTIPHASE_SURFACTANT_PF_
!Allocate window for psi array
varsize_tot=nxt*nyt*nzt*extent
call mpi_win_allocate_shared(varsize_tot,disp_unit,mpi_info_null,shcomm,baseptr,winPsi,ierr)
! get location of memory segment
if(nodeRank.ne.0)then
  call mpi_win_shared_query(winPsi,0,varsize_tot,disp_unit,baseptr,ierr)
endif
call c_f_pointer(baseptr,psi,[nxt,nyt,nzt])
#endif
#else
allocate(top(nxt,nyt,nzt))
allocate(kur(nxt,nyt,nzt))
allocate(vel(nxt,nyt,nzt,3))
allocate(nor(nxt,nyt,nzt,3))
#ifdef _MULTIPHASE_SURFACTANT_PF_
allocate(psi(nxt,nyt,nzt))
#endif
#endif
#ifdef CURVATURE_PDF
!! initialize jpdf and axes
allocate(jpdf(0:nkur,0:ndim,2))
jpdf=0.0d0
do jp=1,ndim
  ! maximum drop size is one domain side
  jpdf(0,jp,:)=dble(jp-1)*max(Lx,Ly,Lz)/dble(ndim-1)
enddo
do ip=1,nkur
  ! minimum radius of curvature (1/max curvature) is half a grid cell. Curvature can go from -maxkur to +maxcur
  jpdf(ip,0,:)=dble(ip-1)*2.0d0*(2.0d0*max(dxi,dyi,dzi)/dble(nkur-1))-2.0d0*max(dxi,dyi,dzi)
enddo
#endif
! create output files
open(42,file='./output/dropCount.dat',form='formatted',position='append')
close(42,status='keep')

do ii = begin, finish, step
  do delayRep = 0,1
    !read drop positions after a short delay, to find coalesences and breakups
    if(delay.eq.0.and.delayRep.eq.1) cycle
    i = ii + delayRep*delay
    fileReadErr=0
#ifdef MPI_WINDOW
    call mpi_win_fence(0,winTop,ierr)
    call mpi_win_fence(0,winCar,ierr)
    call mpi_win_fence(0,winVel,ierr)
    call mpi_win_fence(0,winkur,ierr)
    call mpi_win_fence(0,winNor,ierr)
#ifdef _MULTIPHASE_SURFACTANT_PF_
    call mpi_win_fence(0,winPsi,ierr)
#endif
    if(nodeRank.eq.0) then
      call read_fields(i,fileReadErr)
    endif
    call mpi_win_fence(0,winTop,ierr)
    call mpi_win_fence(0,winCar,ierr)
    call mpi_win_fence(0,winVel,ierr)
    call mpi_win_fence(0,winkur,ierr)
    call mpi_win_fence(0,winNor,ierr)
#ifdef _MULTIPHASE_SURFACTANT_PF_
    call mpi_win_fence(0,winPsi,ierr)
#endif
#else
    call read_fields(i,fileReadErr)
    write(*,*) 'rank ',rank,'on',trim(procname),'count',i
#endif
    call mpi_allreduce(MPI_IN_PLACE,fileReadErr,1,MPI_INT,MPI_SUM,mpi_comm_world,ierr)
    if(fileReadErr.ne.0) then
      if(rank.eq.0) write(6,*) 'fileReadErr on step',i,'skip to next step'
      cycle
    endif
    call get_interface(i)
    if(rank.eq.0.and.i.eq.finish) then
      write(iStepChar,'(i12.12)') i
      filename='paraview/top'//iStepChar//'.vtk'
      write(6,*) 'save drops from final frame ',trim(filename)
      call write_paraview(top,filename)
    endif
  enddo
enddo
#ifdef CURVATURE_PDF
!! reduce (sum) cumulative JPDF and write to file
  call MPI_Allreduce(MPI_IN_PLACE,jpdf(1:ndim,1:nkur,2),ndim*nkur,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(rank.eq.0)then
   ! write cumulative jpdf
   ! generate format string for writing 
   write(fmtstr,'(a,i0,a)') '(',ndim+1,'(e16.8))'
   open(45,file='./output/jpdf.dat',form='formatted',status='replace')
     do i=0,nkur
       write(45,fmtstr) jpdf(i,:,2)
     enddo
   close(45,status='keep')
  endif
#endif
! free up memory
#ifdef MPI_WINDOW
call mpi_win_fence(0,winTop,ierr)
call mpi_win_free(winTop,ierr)
call mpi_win_fence(0,winCar,ierr)
call mpi_win_free(winCar,ierr)
call mpi_win_fence(0,winvel,ierr)
call mpi_win_free(winvel,ierr)
call mpi_win_fence(0,winkur,ierr)
call mpi_win_free(winkur,ierr)
call mpi_win_fence(0,winnor,ierr)
call mpi_win_free(winnor,ierr)
#ifdef _MULTIPHASE_SURFACTANT_PF_
call mpi_win_fence(0,winPsi,ierr)
call mpi_win_free(winPsi,ierr)
#endif
#else
deallocate(top)
deallocate(car)
deallocate(vel)
deallocate(nor)
deallocate(kur)
#ifdef _MULTIPHASE_SURFACTANT_PF_
deallocate(psi)
#endif
#endif
#ifdef CURVATURE_PDF
deallocate(jpdf)
#endif
!write(*,*) 'rank ',rank,'on',trim(procname),'finished'
call mpi_finalize(ierr)
if (rank.eq.0) write(6,*) 'This is the end'
end program mass_center
