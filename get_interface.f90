#include "config.h"
#include "configDropStats.h"
module mod_get_interface
  implicit none
contains
subroutine get_interface(nstep)
  !use mpi
  use mod_common_data, only : pi, rank, top, ntask, vel, ragged_array, kur, nor, car
  use mod_param, only : nxt, nyt, nzt, lx, ly, lz, dx, dy, dz
  use mod_flood_fill,only : flood_fill
  use paraview      ,only : write_paraview
  use, intrinsic :: ISO_FORTRAN_ENV
#ifdef RAND_SEARCH
  use mod_common_data, only : ranx
#endif
#ifdef _MULTIPHASE_SURFACTANT_PF_
  use mod_common_data, only : psi
#endif
#ifdef CURVATURE_PDF
  use mod_common_data, only : jpdf, nkur, ndim
#endif
  implicit none
  integer,parameter :: maxMom=2,statU=41,posiU=42,veloU=43,MoInU=44,counU=45,topoU=46,voidU=47
  logical :: vtkWritten
  !use int8 for non shared arrays to save memory
  integer(kind=int8), dimension(nxt,nyt,nzt) :: topLocal, s_drop, void
  integer(kind=int8), dimension(0:1,0:1,0:1) :: paint, neigh
  integer :: i,j,k,id,jd,kd,ip,jp,kp,iq,jq,kq,direc,iShifted,mom,ierr,last0(3),last1(3),pos(3),ii,jj,error
  integer :: cols,paintIt,faceOnCorner,genus,nVoids,n1Voids,onInt,handles
  double precision, dimension(maxMom,3) :: dropPos, dropVel
  double precision, dimension(3), parameter :: l = (/lx,ly,lz/)
  double precision :: MoI(3,3), work(8), eiVals(3)
  double precision :: invSize, sizeMean, invsizeMean, surfConc, surfConcStdDev, diag, deformation, dropArea, dA
  double precision :: maxNor, kurMean, kurStdDev, voidSzMn, voidSzStDv, kurInv, kurInvSq
  type(ragged_array) :: hist(3) !histogram of drop mass in x, y and z directions
  integer(kind=int64) :: nstep, dropCount, dropSize, remain, faces, edges, vertices, voidSize
  character(len=200) :: fileEnd 
  character(len=100) :: fmtstr
  character(len=12) :: iStepChar
  character(len=5) :: rankChar
#ifdef CURVATURE_PDF
  double precision :: diameter
  integer, dimension(1) :: indkur,inddim
#endif
#ifdef RAND_SEARCH
  integer :: nStabs !number of locations to look for drops
  integer :: iStab
  real :: ir
#else
  integer, dimension(3) :: decStart, decEnd, nl, decCoord
  integer, dimension(3), parameter :: decDims = (/10,10,10/), nt = (/nxt,nyt,nzt/)
#endif
  if(decDims(1)*decDims(2)*decDims(3).ne.ntask) then
    write(*,*) 'number of processes=', ntask,'is not the product of decDims=', decDims
    call mpi_finalize(ierr) 
  endif
  write(iStepChar,'(i12.12)') nstep
  write(rankChar,'(i5.5)') rank
  fileEnd=iStepChar//'rank'//rankChar//'.dat'
  open(statU,file='./output/statDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
  close(statU,status='keep')
  open(posiU,file='./output/posiDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
  close(posiU,status='keep')
  open(veloU,file='./output/veloDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
  close(veloU,status='keep')
  open(MoInU,file='./output/MoInDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
  close(MoInU,status='keep')
  open(topoU,file='./output/topoDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
  close(topoU,status='keep')
  open(voidU,file='./output/voidDropsT'//trim(fileEnd),access='append',form='formatted',status="REPLACE")
  close(voidU,status='keep')
#ifdef _MULTIPHASE_NNF_EVP_
  write(*,*) iStepChar, ' finding unyielded regions using the von Mises criterion'
#endif
  dropCount=0
  sizeMean=0.0
  invsizeMean=0.0
  surfConc=0.0
  surfConcStdDev=0.0
  dropPos=0.0
  dropVel=0.0
  vtkWritten=.false.
  topLocal=top
#ifdef CURVATURE_PDF
  ! zero JPDF for step nstep (do not zero cumulative JPDF)
  jpdf(1:nkur,1:ndim,1)=0.0d0
#endif
  do direc=1,3
    allocate(hist(direc)%v(nt(direc)))
  enddo
  hist(1)%iStepChar='1'
  hist(2)%iStepChar='2'
  hist(3)%iStepChar='3'
#ifdef RAND_SEARCH
  nStabs = 10**4 !nxt*nyt*nzt/10000 !IC 10*3 takes 4 mins
  write(*,*) 'randomly searching ',nStabs,'points'
  ! Random seed
  ir = ranx(313131+rank+nstep)
  do iStab=1,nStabs
        if ( mod(iStab, nStabs/10) .eq.0) then 
          write(*,*) 'step=',iStepChar,' rank=',rank,' iStab=',iStab, 'out of',nStabs
        endif
        id = int(ranx(0)*nxt+1.0)
        if (id.gt.nxt) cycle    !IC skip to next step if ranx returns 1
        jd = int(ranx(0)*nyt+1.0)
        if (jd.gt.nyt) cycle
        kd = int(ranx(0)*nzt+1.0)
        if (kd.gt.nzt) cycle
#else
  !assign a cuboid subset of the domain to the mpi process
  decCoord(1) = rank/decDims(2)/decDims(3)
  decCoord(2) = rank/decDims(3) - decCoord(1)*decDims(2)
  decCoord(3) = rank - decCoord(2)*decDims(3) - decCoord(1)*decDims(3)*decDims(2)
  do direc=1,3
    if(mod(nt(direc),decDims(direc)).eq.0) then 
      nl(direc) = nt(direc)/decDims(direc)
      decStart(direc) = decCoord(direc)*nl(direc)+1
      decEnd(direc) = decStart(direc) + nl(direc)-1
    else 
      write(*,*) 'array size not divisible by decDims in direc', direc
      call mpi_finalize(ierr) 
    endif
  enddo
  do kd=decStart(3), decEnd(3)
    remain=sum(int(topLocal(decStart(1):decEnd(1),&
      decStart(2):decEnd(2),decStart(3):decEnd(3))))
    if(rank.eq.0.and.mod(kd,nl(3)/10).eq.0) write(*,*) 'rank',rank,'start',decStart(3),'kd',kd,'end',decEnd(3),'remain',remain
    do jd=decStart(2), decEnd(2)
      do id=decStart(1), decEnd(1)
#endif     
        if(topLocal(id,jd,kd).gt.0)then
          !find a single drop with the painter algorithm
          s_drop=0
          s_drop(id,jd,kd)=1
          dropSize=0
          call flood_fill(id,jd,kd,1,dropSize,s_drop)
#ifndef RAND_SEARCH      
          !remove the drop from topLocal so it is not found again
          topLocal=topLocal-s_drop
          remain=sum(int(topLocal(decStart(1):decEnd(1),&
            decStart(2):decEnd(2),decStart(3):decEnd(3))))
#endif
          !calculate drop statistics if minimum i j k in the drop lies in the processor's domain
          extremityLoop: do k=1,nzt
            do j=1,nyt
              do i=1,nxt
                if(s_drop(i,j,k).ne.0) exit extremityLoop
              enddo
            enddo
          enddo extremityLoop
          if(i.ge.decStart(1).and.i.le.decEnd(1).and. &
             j.ge.decStart(2).and.j.le.decEnd(2).and. &
             k.ge.decStart(3).and.k.le.decEnd(3)) then
            dropCount = dropCount + 1
            sizeMean = sizeMean + (dropSize - sizeMean) / (1.0*dropCount)
            invSize = 1.0/dropSize
            invsizeMean = invsizeMean + ( invSize - invsizeMean) / (1.0*dropCount)
            dropArea=0.0
            surfConc=0.0
            surfConcStdDev=0.0
            kurMean=0.0d0
            kurStdDev=0.0d0
            kurInv=0.0d0
            kurInvSq=0.0d0
            nVoids=0
            n1Voids=0
            voidSzMn=0.0
            voidSzStDv=0.0
            faces=0
            edges=0
            vertices=0
#ifdef CURVATURE_PDF
            !get diameter, volume is in cell units
            diameter=(6.0d0*dble(dropSize)*dx*dy*dz/pi)**(1.0/3.0)
            ! find closest index for dimension
            inddim=minloc(dabs(jpdf(0,1:ndim,1)-diameter))
#endif
            do direc=1,3
              hist(direc)%v=0.0
            enddo
            dropVel=0.0
            !make histograms of drop mass in each direction so we can avoid overlap with periodic 
            !boundaries in moment of inertia calculations
            do k=1,nzt
              do j=1,nyt
                do i=1,nxt
                  if(s_drop(i,j,k).ne.0) then
                    hist(1)%v(i)=hist(1)%v(i)+invSize
                    hist(2)%v(j)=hist(2)%v(j)+invSize
                    hist(3)%v(k)=hist(3)%v(k)+invSize
                    do direc=1,3
                      do mom = 1,maxMom
                        dropVel(mom,direc) = dropVel(mom,direc) + (vel(i,j,k,direc)**mom)*invSize
                      enddo
                    enddo
                  endif
                enddo
              enddo
            enddo 
            dropPos=0.0
            last0=0
            last1=0
            do direc=1,3
              !find beginning and end of drop in x,y,z directions
              do i=1,nt(direc)
                ip=i+1
                if(ip.gt.nt(direc)) ip=ip-nt(direc)
                if(hist(direc)%v(i).gt.0.5*invSize.and.hist(direc)%v(ip).lt.0.5*invSize) last1(direc)=i
                if(hist(direc)%v(i).lt.0.5*invSize.and.hist(direc)%v(ip).gt.0.5*invSize) last0(direc)=i
              enddo
              !if part of the drop touches i=1, shift that part to other side of domain so drop is contiguous
              do i=1,nt(direc)
                iShifted = i
                if(last1(direc).lt.last0(direc).and.i.le.last1(direc)) iShifted = i+nt(direc)
                do mom = 1,maxMom
                  dropPos(mom,direc) = dropPos(mom,direc) + hist(direc)%v(i)*(iShifted**mom)
                enddo
              enddo 
            enddo
            MoI=0.0
            do k=1,nzt
              pos(3)=k
              if(last1(3).lt.last0(3).and.k.le.last1(3)) pos(3) = k+nzt
              do j=1,nyt
                pos(2)=j
                if(last1(2).lt.last0(2).and.j.le.last1(2)) pos(2) = j+nyt
                do i=1,nxt
                  pos(1)=i
                  if(last1(1).lt.last0(1).and.i.le.last1(1)) pos(1) = i+nxt
                  !use https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
                  !Bunner & Tryggvason JFM 2003 misses out on diag part of MoI definition
                  !don't bother moving dropPos inside domain until after the moment of inertia calculations
                  if(s_drop(i,j,k).ne.0) then
                    diag=0.0
                    do ii=1,3
                      diag = diag + ( pos(ii)-dropPos(1,ii) )**2 *dx**5
                    enddo
                    do ii=1,3
                      do jj=1,3
                        MoI(ii,jj) = MoI(ii,jj) - ( pos(ii)-dropPos(1,ii) )*( pos(jj)-dropPos(1,jj) )*dx**5
                        if(ii.eq.jj) MoI(ii,jj) = MoI(ii,jj) + diag
                      enddo
                    enddo
                  endif
                  !the interface is here if s_drop changes in any of the 6 directions
                  onInt=0
                  do direc=1,3
                    ip = i
                    jp = j
                    kp = k
                    iq = i
                    jq = j
                    kq = k
                    if (direc.eq.1) then
                      ip=i+1
                      if(ip.gt.nxt)ip=ip-nxt
                      iq=i-1
                      if(iq.lt.1)  iq=iq+nxt
                    endif
                    if (direc.eq.2) then
                      jp=j+1
                      if(jp.gt.nyt)jp=jp-nyt
                      jq=j-1
                      if(jq.lt.1)  jq=jq+nyt
                    endif
                    if (direc.eq.3) then
                      kp=k+1
                      if(kp.gt.nzt)kp=kp-nzt
                      kq=k-1
                      if(kq.lt.1)  kq=kq+nzt
                    endif
                    if (s_drop(ip,jp,kp).ne.s_drop(i,j,k)) faces = faces + 1
                    if (s_drop(ip,jp,kp).ne.s_drop(i,j,k).or.s_drop(iq,jq,kq).ne.s_drop(i,j,k)) onInt=1
                  enddo
                  if(onInt.eq.1)then
                    !find direction which is most aligned with interface normal
                    maxNor=0.0
                    dA=0.0
                    do direc=1,3
                      if(abs(nor(i,j,k,direc)).gt.maxNor)then
                        maxNor = abs(nor(i,j,k,direc))
                        if(direc.eq.1)dA=dy*dz
                        if(direc.eq.2)dA=dz*dx
                        if(direc.eq.3)dA=dx*dy
                      endif
                    enddo
                    !area of interface is 
                    dA=dA/(2.0*maxNor)
                    dropArea = dropArea + dA
                    kurMean = kurMean + dA*kur(i,j,k)
                    kurStdDev = kurStdDev + dA*kur(i,j,k)**2
                    kurInv = kurInv + dA/kur(i,j,k)
                    kurInvSq = kurInvSq + dA/kur(i,j,k)**2
#ifdef CURVATURE_PDF 
                    indkur = minloc(dabs(jpdf(1:nkur,0,1)-kur(i,j,k)))
                    jpdf(indkur(1),inddim(1),:) = jpdf(indkur(1),inddim(1),:) + 1.0d0
#endif
#ifdef _MULTIPHASE_SURFACTANT_PF_
                    surfConc = surfConc + dA*psi(i,j,k)
                    surfConcStdDev = surfConcStdDev + dA*psi(i,j,k)**2
#endif
                  endif
                  !find size of voids
                  ip = modulo(i,nxt)+1
                  if (s_drop(i,j,k).eq.1.and.s_drop(ip,j,k).eq.0.and.car(ip,j,k).eq.0) then
                    void=0
                    voidSize=0
                    call flood_fill(ip,j,k,0,voidSize,void)
                    !if (rank.eq.0) 
                    !write(*,*) 'rank',rank,'dropSize',dropSize,'voidSize',voidSize,'sumCar',sum(int(car))
                    do kp=1,nzt
                      do jp=1,nyt
                        do ip=1,nxt
                          if(void(ip,jp,kp).ne.0) car(ip,jp,kp)=1
                        enddo
                      enddo
                    enddo 
                    if (voidSize.gt.1) then
                      nVoids=nVoids+1
                      voidSzMn = voidSzMn + (voidSize - voidSzMn)/nVoids
                      voidSzStDv = voidSzStDv + (voidSize**2 - voidSzStDv)/nVoids
                    else
                      !count voids of size 1 separately, because they likely result from vof coalescence
                      n1Voids=n1Voids+1
                    endif
                    open(voidU,file='./output/voidDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
                      write(voidU,'(5i16)') rank,dropCount,dropSize,n1Voids+nVoids,voidSize
                    close(voidU,status='keep')
                  endif
                  !I Cannon 2023 May
                  !calculate the genus of the drop using Mendoza et al. 2006 Acta Mater 10.1016/j.actamat.2005.10.010.
                  !genus = tunnels - voids 
                  !for a shape made of simple polygons:
                  !genus = 1 - (faces - edges + corners)/2
                  !use a painter algorithm to count the number of contiguous drop corners at this vertex
                  !paint the vof with up to 8 colours surrounding a vertex between cells
                  cols=0
                  neigh=0
                  paint=0
                  do kq=0,1
                    kp = k + kq
                    if(kp.gt.nzt) kp=kp-nzt
                    do jq=0,1
                      jp = j + jq
                      if(jp.gt.nyt) jp=jp-nyt
                      do iq=0,1
                        ip = i + iq
                        if(ip.gt.nxt) ip=ip-nxt
                        neigh(iq,jq,kq)=s_drop(ip,jp,kp)
                        if(s_drop(ip,jp,kp).ne.0) then
                          cols=cols+1
                          paint(iq,jq,kq)=cols
                        endif
                      enddo
                    enddo
                  enddo
                  faceOnCorner=0
                  !no corners if vertex entirely inside drop
                  if(cols.eq.8) cols=0
                  !trivial case with all drop except for one corner
                  if(cols.eq.7) cols=1
                  !if drop only at the corner, there are three faces touching it
                  if(cols.eq.1) faceOnCorner=3
                  !if just two voids on opposite sides of the vertex, we connect the voids using a tunnel with
                  !Euler characteristic = v-e+f = 3-9+6 = 0 = 0-6+6. I choose v=0 to skip the painting loop below
                  if(cols.eq.6) then
                    cols=0
                    faceOnCorner=6
                    do jq=0,1
                      do iq=0,1
                        if(neigh(iq,jq,1).ne.neigh(1-iq,1-jq,0)) then
                          cols=6
                          faceOnCorner=0
                        endif
                      enddo
                    enddo
                  endif
                  if(cols.gt.1) then
                    !in the worst case, we need 4 iterations to paint the 2x2x2 cube by neighbours
                    !loop 4 times instead of recursion
                    do paintIt=1,4
                      do kq=0,1
                        do jq=0,1
                          do iq=0,1
                            if(neigh(iq,jq,kq).gt.0) then
                              do direc=1,3
                                ip=iq
                                jp=jq
                                kp=kq
                                if(direc.eq.1) ip=1-iq
                                if(direc.eq.2) jp=1-jq
                                if(direc.eq.3) kp=1-kq
                                !paint any touching cells the same colour
                                if(neigh(iq,jq,kq).eq.neigh(ip,jp,kp)) then
                                  if(paint(ip,jp,kp).ne.paint(iq,jq,kq)) then
                                    paint(iq,jq,kq) = min(paint(iq,jq,kq), paint(ip,jp,kp))
                                    paint(ip,jp,kp) = min(paint(iq,jq,kq), paint(ip,jp,kp))
                                  endif
                                !if s_drop has a sign change, there is a face and an edge
                                else if(paintIt.eq.1) then
                                  faceOnCorner=faceOnCorner+1
                                endif
                              enddo
                            endif
                          enddo
                        enddo
                      enddo
                    enddo
                    !count how many colours remain, this is equal to the number of drop corners at the vertex
                    cols=0
                    do paintIt=1,8
                      colourLoop: do kq=0,1
                        do jq=0,1
                          do iq=0,1
                            if(paintIt.eq.paint(iq,jq,kq)) then
                              cols=cols+1
                              exit colourLoop
                            endif
                          enddo
                        enddo
                      enddo colourLoop
                    enddo
                  endif
                  vertices = vertices + cols
                  !the number of edges touching the vertex is equal to the number of faces
                  edges=edges+faceOnCorner
                enddo
              enddo
            enddo
            !every edge is between two corners so we double counted the edges
            edges = edges/2
            if(mod(vertices-edges+faces,2).ne.0)then
              write(*,*)'rank',rank,'count',dropCount,'warning: Euler characteristic is odd, v,e,f=',vertices,edges,faces
            endif
            genus = 1-(vertices-edges+faces)/2
            handles=genus+n1Voids+nVoids
            if(.false..and.dropSize.gt.1000000.and.&
              last1(1).gt.last0(1).and.&
              last1(2).gt.last0(2).and.&
              last1(3).gt.last0(3)) then !.not.vtkWritten.and.
              write(*,*) 'myid',rank,'dropCount',dropCount,'dropSize',dropSize,'handles',handles,'voids',n1Voids+nVoids,&
                'voidSzMn',voidSzMn
              write(fmtstr,'(a,a,a,i5.5,a,i5.5,a,i8.8,a,i5.5,a,i5.5,a)') 'paraview/',iStepChar,'rank',rank,'dropCount',&
                dropCount,'size',dropSize,'handles',handles,'voids',n1Voids+nVoids,'.vtk'
              call write_paraview(s_drop,fmtstr)
              vtkWritten=.true.
            endif
            kurMean = kurMean / dropArea
            kurStdDev = dsqrt( kurStdDev/dropArea - kurMean**2 )
            kurInv = kurInv / dropArea
            kurInvSq = kurInvSq / dropArea
            surfConc = surfConc/dropArea
            surfConcStdDev = dsqrt(surfConcStdDev/dropArea - surfConc**2)
            voidSzStDv = dsqrt(voidSzStDv - voidSzMn**2)
            do direc=1,3
              !move drop back inside domain
              if(dropPos(1,direc).gt.nt(direc)) then
                do mom=1,maxMom
                  dropPos(mom,direc)=dropPos(mom,direc)-nt(direc)**mom
                enddo
              endif
              !put in units of simulation domain size
              do mom=1,maxMom
                dropPos(mom,direc)=dropPos(mom,direc) * ( l(direc)/nt(direc) )**mom
              enddo
            enddo
            !Calculate eigenvalues of moment of inertia
            if(last0(1).eq.0.or.last0(2).eq.0.or.last0(3).eq.0) then
              write(*,*) 'rank',rank,'drop',dropCount,' spans all of domain, deformation is undefined'
              deformation=-1.0
            else
              call DSYEV("N","U",3,MoI,3,eiVals,work,8,error)
              if (error.ne.0) write(*,*) 'dropSize', dropSize, 'eiVals', eiVals
              if (eiVals(1).gt.1.0e-8) then
                deformation=sqrt(eiVals(3)/eiVals(1))
              else
                deformation=1.0
              endif
            endif
            open(statU,file='./output/statDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
              write(statU,'(3i16,6E16.8)') rank,dropCount,dropSize,dropArea,surfConc,surfConcStdDev,deformation,kurMean,kurStdDev
            close(statU,status='keep')
            ! generate format string for writing 
            write(fmtstr,'(a,i0,a)') '(3i16,',maxMom*3,'(e16.8))'
            open(posiU,file='./output/posiDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
              write(posiU,fmtstr) rank,dropCount,dropSize,(dropPos(mom,:),mom=1,maxMom)
            close(posiU,status='keep')
            open(veloU,file='./output/veloDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
              write(veloU,fmtstr) rank,dropCount,dropSize,(dropVel(mom,:),mom=1,maxMom)
            close(veloU,status='keep')
            open(MoInU,file='./output/MoInDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
              write(MoInU,'(3i16,9E16.8,6i16,2E16.8)') rank,dropCount,dropSize,(eiVals(i),i=1,3),&
                MoI(1,1),MoI(1,2),MoI(1,3),MoI(2,2),MoI(2,3),MoI(3,3),last0(:),last1(:),kurInv,kurInvSq
            close(MoInU,status='keep')
            open(topoU,file='./output/topoDropsT'//trim(fileEnd),access='append',form='formatted',status='old')
              write(topoU,'(9i16,2E16.8)')rank,dropCount,dropSize,vertices,edges,faces,genus,n1Voids,nVoids,voidSzMn,voidSzStDv
            close(topoU,status='keep')
          endif
        endif
#ifndef RAND_SEARCH    
      enddo
    enddo
#endif
  enddo
  !2022 Oct 17 Ianto Cannon: inverse of invsizeMean gives the mean drop size with random search
  if(dropCount.gt.0) then
    open(counU,file='./output/dropCount.dat',access='append',form='formatted',status='old')
      write(counU,'(4i16,2E15.7)') nstep,rank,dropCount,remain,sizeMean,invsizeMean
    close(counU,status='keep')
  endif
#ifdef CURVATURE_PDF
!! reduce (sum) nstep JPDF and write to file
  jpdf(1:ndim,1:nkur,1)=dble(rank)
  call MPI_Allreduce(MPI_IN_PLACE,jpdf(1:ndim,1:nkur,1),ndim*nkur,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
  if(rank.eq.0)then
   ! write jpdf at nstep
   ! generate format string for writing 
   write(fmtstr,'(a,i0,a)') '(',ndim+1,'(e16.8))'
   open(45,file='./output/jpdf_'//trim(iStepChar)//'.dat',form='formatted',status='replace')
   do i=0,nkur
     write(45,fmtstr) jpdf(i,:,1)
   enddo
   close(45,status='keep')
  endif
#endif
  return
end subroutine get_interface
end module mod_get_interface
