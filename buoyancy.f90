
subroutine buoyancy()

use vars
use params
implicit none
	
integer i,j,k,kb
real du(nx,ny,nz,3)

if(docolumn) return

if(dostatis) then

  do k=1,nzm
    do j=1,ny
      do i=1,nx
         du(i,j,k,3)=dwdt(i,j,k,na)
      end do
    end do
  end do

endif

do k=2,nzm	
 kb=k-1
 do j=1,ny
  do i=1,nx
   dwdt(i,j,k,na)=dwdt(i,j,k,na)+0.5*(bet(k)* &
              (tabs0(k)*(0.61*(q(i,j,k)-q0(k))-1.61*qn(i,j,k)-qp(i,j,k)) &
             +(tabs(i,j,k)-tabs0(k))*(1.+0.61*q0(k))) + bet(kb)* &
              (tabs0(kb)*(0.61*(q(i,j,kb)-q0(kb))-1.61*qn(i,j,kb)-qp(i,j,kb)) &
             +(tabs(i,j,kb)-tabs0(kb))*(1.+0.61*q0(kb))) & 
             - betp(k)*p(i,j,k) - betp(kb)*p(i,j,kb))
  end do ! i
 end do ! j
end do ! k

	 
if(dostatis) then

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        du(i,j,k,1)=0.
        du(i,j,k,2)=0.
        du(i,j,k,3)=dwdt(i,j,k,na)-du(i,j,k,3)
      end do
    end do
  end do

  call stat_tke(du,tkelebuoy)
  call stat_mom(du,momlebuoy)
  call setvalue(twlebuoy,nzm,0.)
  call setvalue(qwlebuoy,nzm,0.)
  call setvalue(swlebuoy,nzm,0.)
  call stat_sw1(du,twlebuoy,qwlebuoy,swlebuoy)

endif

 
end subroutine buoyancy


