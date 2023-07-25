#include "config.h"
#include "configDropStats.h"
module mod_flood_fill
contains

recursive subroutine flood_fill(i,j,k,a,sz,s_drop)
!I Cannon 2023 May: Painter algorithm. Finds contiguous ones in top. 
!If top(i,j,k) has neighbors also filled, it calls itself on the neighbour
!Routine stops calling itself when all neighbours have been found, and s_drop contains the entire drop
use mod_common_data,  only : top
use mod_param,  only : nxt, nyt, nzt
use, intrinsic :: ISO_FORTRAN_ENV
implicit none
integer, intent(in) :: i,j,k,a
integer(kind=int64),intent(inout) :: sz
integer(kind=int8), intent(out) ::s_drop(nxt,nyt,nzt)
integer :: di,dj,dk,dir,ip,jp,kp
s_drop(i,j,k)=1
!stop recursion if drop is very large. Recursing 500**3 times uses over 3GB
!if(sz.lt.50**3) then
sz=sz+1
do dk = 0,1
  do dj = 0,1
    do di = 0,1
      !for drops search only in x, y and z direction (6 directions). For voids search also diagonals
      if(a.eq.1.and.(di+dj+dk).gt.1) cycle
      do dir=1,2
        kp=k+dk*(-1)**dir
#ifdef _FLOW_TRIPERIODIC_
        if(kp.gt.nzt)kp=kp-nzt
        if(kp.lt.1)  kp=kp+nzt
#else
        if(kp.gt.nzt)cycle
        if(kp.lt.1)  cycle
#endif
        jp=j+dj*(-1)**dir
        if(jp.gt.nyt)jp=jp-nyt
        if(jp.lt.1)  jp=jp+nxt
        ip=i+di*(-1)**dir
        if(ip.gt.nxt)ip=ip-nyt
        if(ip.lt.1)  ip=ip+nxt
        if((top(ip,jp,kp).eq.a).and.(s_drop(ip,jp,kp).eq.0)) call flood_fill(ip,jp,kp,a,sz,s_drop)
      enddo
    enddo
  enddo
enddo
!endif
return
end subroutine flood_fill

recursive subroutine carrier_fill(s_drop,i,j,k)
!I Cannon 2023 May: Painter algorithm. Finds contiguous zeros in top. This is used to find voids inside drops.
!Same as above routine but looks for zeros not ones, and assumes voids are not connected diagonally
use mod_common_data,  only : top
use mod_param,  only : nxt, nyt, nzt
use, intrinsic :: ISO_FORTRAN_ENV
implicit none
integer(kind=int8), intent(out) ::s_drop(nxt,nyt,nzt)
integer :: i,j,k,di,dj,dk,inew,jnew,knew
s_drop(i,j,k)=1
do dk = -1,1
  knew=k+dk
#ifdef _FLOW_TRIPERIODIC_
  if(knew.gt.nzt) knew=knew-nxt
  if(knew.lt.1)   knew=knew+nxt
#else
  if(knew.gt.nzt) cycle
  if(knew.lt.1)   cycle
#endif
  do dj = -1,1
    jnew=j+dj
    if(jnew.gt.nyt)jnew=jnew-nxt
    if(jnew.lt.1)jnew=jnew+nxt
    do di = -1,1
#if defined(DIAGONALS_CONNECTED)
!drops are connected on the diagonals, so voids cannot be
      if(abs(di)+abs(dj)+abs(dk).gt.1) cycle
#endif
      inew=i+di
      if(inew.gt.nxt)inew=inew-nxt
      if(inew.lt.1)inew=inew+nxt
      if((top(inew,jnew,knew).eq.0).and.(s_drop(inew,jnew,knew).eq.0)) call carrier_fill(s_drop,inew,jnew,knew)
    enddo
  enddo
enddo
return
end subroutine carrier_fill
end module mod_flood_fill
