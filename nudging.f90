subroutine nudging
	
use vars
use params
implicit none

real coef
integer i,j,k
	
coef = 1./tauls

if(donudging_uv) then

    coef = 1. / tauls
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))*coef
           dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))*coef
        end do
      end do
    end do

endif


if(donudging_tq) then

    coef = dtn / tauls
    do k=1,nzm
      do j=1,ny
        do i=1,nx
           t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*coef
           q(i,j,k)=q(i,j,k)-(q0(k)-qg0(k))*coef
        end do
      end do
    end do

endif

end subroutine nudging
