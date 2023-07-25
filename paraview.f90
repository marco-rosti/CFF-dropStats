module paraview
contains
subroutine write_paraview(fld,namefile)
  use, intrinsic :: ISO_FORTRAN_ENV
  use mod_param, only : nxt, nyt, nzt, lx, ly, lz, dx, dy, dz
  implicit none
  integer(kind=int8), dimension(nxt,nyt,nzt), intent(in)::fld
  integer, parameter :: nfields=1
  integer, parameter :: ir1=1, ir2=nxt, jr1=1, jr2=nyt, kr1=1, kr2=nzt
  integer, parameter :: nxr=ir2-ir1+1, nyr=jr2-jr1+1, nzr=kr2-kr1+1
  integer :: i,j,k
  double precision :: x(nxr), y(nyr), z(nzr)
  character(len=100) :: namefile
  character(len=80) :: buffer
  character(len=10) :: lf
  character(len=8) :: str1,str2,str3
  character(len=16) :: str4

  ! end of line character
  lf=achar(10)

  do i=1,nxr
    x(i)=(dble(ir1-1+i-1)+0.5d0)*dx
  enddo
  do j=1,nyr
    y(j)=(dble(jr1-1+j-1)+0.5d0)*dy
  enddo
  do k=1,nzr
    z(k)=(dble(kr1-1+k-1)+0.5d0)*dz
  enddo

    write(*,*) 'make vtk ', namefile
    open(66,file=trim(namefile),status='replace',form='unformatted',access='stream',convert='big_endian')

    ! start writing vtk file
    ! write header
    buffer='# vtk DataFile Version 3.0'//lf
    write(66) trim(buffer)
    buffer='Reduced data file'//lf
    write(66) trim(buffer)
    buffer='BINARY'//lf
    write(66) trim(buffer)
    buffer='DATASET RECTILINEAR_GRID'//lf
    write(66) trim(buffer)

    ! write grid
    write(str1(1:8),'(i8)') nxr
    write(str2(1:8),'(i8)') nyr
    write(str3(1:8),'(i8)') nzr
    buffer='DIMENSIONS '//str1//str2//str3//lf
    write(66) trim(buffer)

    buffer='X_COORDINATES '//str1//'  float'//lf ;
    write(66) trim(buffer)
    do i=1,nxr
      write(66) real(x(i))
    enddo
    buffer='Y_COORDINATES '//str2//'  float'//lf ;
    write(66) trim(buffer)
    do j=1,nyr
      write(66) real(y(j))
    enddo
    buffer='Z_COORDINATES '//str3//'  float'//lf ;
    write(66) trim(buffer)
    do k=1,nzr
      write(66) real(z(k))
    enddo

    ! write content (data format)
    write(str4(1:16),'(i16)') nxr*nyr*nzr
    buffer='POINT_DATA '//str4//lf
    write(66) trim(buffer)

    write(str1(1:8),'(i8)') nfields
    buffer='FIELD FieldData '//str1//lf
    write(66) trim(buffer)
      ! write vof
      write(str4(1:16),'(i16)') nxr*nyr*nzr

      buffer = 'VOF 1 '//str4//' float'//lf
      write(66) trim(buffer)
      do k=1,nzr
        do j=1,nyr
          do i=1,nxr
            write(66) real(fld(i,j,k))
          enddo
        enddo
      enddo
    buffer=lf
    write(66) trim(buffer)
    close(66,status='keep')
  return
end subroutine write_paraview
end module paraview
