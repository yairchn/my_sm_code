    
subroutine upperbound

use vars
use params
implicit none
real coef
integer i,j,k
real tau_nudging
parameter ( tau_nudging = 3600.)

if(dolargescale) then

!
! if there is an "observed" sounding - nudge two highest levels to it
! to avoid problems with the upper boundary.
!

  coef = dtn / tau_nudging
  do k=nzm-1,nzm
    do j=1,ny
      do i=1,nx
         t(i,j,k)=t(i,j,k)-(t(i,j,k)-tg0(k)-gamaz(k))*coef
         q(i,j,k)=q(i,j,k)-(q(i,j,k)-qg0(k))*coef
      end do
    end do
  end do

else
!
!  otherwise, preserve the vertical gradients:
!
  coef = dz*adz(nzm)
  gamt0=(t0(nzm-1)-t0(nzm-2))/(z(nzm-1)-z(nzm-2))
  gamq0=(q0(nzm-1)-q0(nzm-2))/(z(nzm-1)-z(nzm-2))
  do j=1,ny
   do i=1,nx
     t(i,j,nzm)=t(i,j,nzm-1)+gamt0*coef
     q(i,j,nzm)=q(i,j,nzm-1)+gamq0*coef
   end do    
  end do 
           
end if

end   
