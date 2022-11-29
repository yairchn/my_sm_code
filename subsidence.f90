
subroutine subsidence()
	
use vars
implicit none

integer i,j,k,k1,k2
real rdz
real t_vtend, q_vtend
real t_tend(nx,ny,nzm), q_tend(nx,ny,nzm)

do k=2,nzm-1
  if(wsub(k).ge.0) then
     rdz=wsub(k)/(dz*adzw(k))	
     k1 = k
     k2 = k-1 
  else
     rdz=wsub(k)/(dz*adzw(k+1))       
     k1 = k+1
     k2 = k
  end if
  do j=1,ny
    do i=1,nx
      dudt(i,j,k,na) = dudt(i,j,k,na) - rdz*(u(i,j,k1)-u(i,j,k2)) 
      dvdt(i,j,k,na) = dvdt(i,j,k,na) - rdz*(v(i,j,k1)-v(i,j,k2)) 
      t_tend(i,j,k) =  - rdz * (t(i,j,k1)-t(i,j,k2))
      q_tend(i,j,k) =  - rdz * (q(i,j,k1)-q(i,j,k2))
    end do
  end do
end do
do k=2,nzm-1
  t_vtend = 0.
  q_vtend = 0.
  do j=1,ny
    do i=1,nx
      t(i,j,k) = t(i,j,k) + dtn * t_tend(i,j,k)
      q(i,j,k) = q(i,j,k) + dtn * q_tend(i,j,k)
      q(i,j,k) = max(0.,q(i,j,k))
      t_vtend = t_vtend + t_tend(i,j,k)
      q_vtend = q_vtend + q_tend(i,j,k)
    end do
  end do
  t_vtend = t_vtend / float(nx*ny) 
  q_vtend = q_vtend / float(nx*ny) 
  ttend(k) = ttend(k) + t_vtend
  qtend(k) = qtend(k) + q_vtend
end do 
	
end
