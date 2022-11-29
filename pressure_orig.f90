! Non-blocking receives before blocking sends

subroutine pressure_orig
	
!       Original pressure solver based on horizontal slabs
!       (C) 1998, 2002 Marat Khairoutdinov
!       Works only when the number of slabs is equal to the number of processors.
!       Therefore, the number of processors shouldn't exceed the number of levels nzm
!       Also, used for a 2D version 
!       For more processors for the given number of levels and 3D, use pressure_big

use vars
implicit none
	

integer, parameter :: npressureslabs = nsubdomains
integer, parameter :: nzslab = max(1,nzm / npressureslabs)
integer, parameter :: nx2=nx_gl+2, ny2=ny_gl+2*YES3D
integer, parameter :: n3i=3*nx_gl/2+1,n3j=3*ny_gl/2+1

real f(nx2,ny2,nzslab) ! global rhs and array for FTP coefficeients
real ff(nx+1,ny+2*YES3D,nzm)	! local (subdomain's) version of f
real buff_slabs(nxp1,nyp2,nzslab,npressureslabs)
real buff_subs(nxp1,nyp2,nzslab,nsubdomains) 
real bufp_slabs(0:nx,1-YES3D:ny,nzslab,npressureslabs)  
real bufp_subs(0:nx,1-YES3D:ny,nzslab,nsubdomains)  
common/tmpstack/f,ff,buff_slabs,buff_subs
equivalence (buff_slabs,bufp_slabs)
equivalence (buff_subs,bufp_subs)

real work(nx2,ny2),trigxi(n3i),trigxj(n3j) ! FFT stuff
integer ifaxj(100),ifaxi(100)

double precision a(nzm),b(nxp1,nyp2),c(nzm),e	
double precision xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty,eign(nxp1,nyp2)
double precision alfa(nxp1,nyp2,nzm-1),beta(nxp1,nyp2,nzm-1)

integer reqs_in(nsubdomains)
integer i, j, k, id, jd, m, n, it, jt, ii, jj, tag, rf
integer nyp22, n_in, count
integer iii(0:nx_gl),jjj(0:ny_gl)
logical flag(nsubdomains)
integer iwall,jwall

! check if the grid size allows the computation:

if(nsubdomains.gt.nzm) then
  if(masterproc) print*,'pressure_orig: nzm < nsubdomains. STOP'
  call task_abort
endif

if(mod(nzm,npressureslabs).ne.0) then
  if(masterproc) print*,'pressure_orig: nzm/npressureslabs is not round number. STOP'
  call task_abort
endif

!-----------------------------------------------------------------

if(docolumn) return

if(dowallx) then
  iwall=1
else
  iwall=0
end if
if(RUN2D) then  
  nyp22=1
  jwall=0
else
  nyp22=nyp2
  if(dowally) then
    jwall=2
  else
    jwall=0
  end if
endif
	
!-----------------------------------------------------------------
!  Compute the r.h.s. of the Poisson equation for pressure

call press_rhs()


!-----------------------------------------------------------------	 
!   Form the horizontal slabs of right-hand-sides of Poisson equation 
!   for the global domain. Request sending and receiving tasks.

! iNon-blocking receive first:

n_in = 0
do m = 0,nsubdomains-1

  if(rank.lt.npressureslabs.and.m.ne.nsubdomains-1) then

    n_in = n_in + 1
    call task_receive_float(bufp_subs(0,1-YES3D,1,n_in), &
                           nzslab*nxp1*nyp1,reqs_in(n_in))
    flag(n_in) = .false.
 
  endif

  if(rank.lt.npressureslabs.and.m.eq.nsubdomains-1) then

    call task_rank_to_index(rank,it,jt)	  
    n = rank*nzslab
    do k = 1,nzslab
     do j = 1,ny
       do i = 1,nx
         f(i+it,j+jt,k) = p(i,j,k+n)
       end do
     end do
    end do
  endif

end do ! m


! Blocking send now:


do m = 0,nsubdomains-1

  if(m.lt.npressureslabs.and.m.ne.rank) then

    n = m*nzslab + 1
    call task_bsend_float(m,p(0,1-YES3D,n),nzslab*nxp1*nyp1, 33)
  endif

end do ! m


! Fill slabs when receive buffers are full:

count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           call task_rank_to_index(rf,it,jt)	  
           do k = 1,nzslab
            do j = 1,ny
             do i = 1,nx
               f(i+it,j+jt,k) = bufp_subs(i,j,k,m)
             end do
            end do
           end do
	endif   
   endif
  end do
end do


!-------------------------------------------------
! Perform Fourier transformation for a slab:

if(rank.lt.npressureslabs) then

 call fftfax(nx_gl,ifaxi,trigxi)
 if(RUN3D) call fftfax(ny_gl,ifaxj,trigxj)

 do k=1,nzslab

  if(dowallx) then
   call cosft(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)
  else
   call fft991(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)
  end if

  if(RUN3D) then
   if(dowally) then
     call cosft(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
   else
     call fft991(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
   end if
  end if

 end do 

endif


! Synchronize all slabs:

call task_barrier()

!-------------------------------------------------
!   Send Fourier coeffiecients back to subdomains:

! Non-blocking receive first:

n_in = 0
do m = 0, nsubdomains-1
		
   call task_rank_to_index(m,it,jt)

   if(rank.lt.npressureslabs.and.m.eq.rank) then

     n = rank*nzslab
     do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          ff(i,j,k+n) = f(i+it,j+jt,k) 
        end do
      end do
     end do 

   end if

   if(m.lt.npressureslabs-1.or.m.eq.npressureslabs-1 &
                            .and.rank.ge.npressureslabs) then

     n_in = n_in + 1
     call task_receive_float(buff_slabs(1,1,1,n_in), &
                                nzslab*nxp1*nyp22,reqs_in(n_in))
     flag(n_in) = .false.	    
   endif

end do ! m

! Blocking send now:

do m = 0, nsubdomains-1

   call task_rank_to_index(m,it,jt)

   if(rank.lt.npressureslabs.and.m.ne.rank) then

     do k = 1,nzslab
      do j = 1,nyp22
       do i = 1,nxp1
         buff_subs(i,j,k,1) = f(i+it,j+jt,k)
       end do
      end do
     end do

     call task_bsend_float(m, buff_subs(1,1,1,1),nzslab*nxp1*nyp22,44)

   endif

end do ! m



! Fill slabs when receive buffers are complete:


count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           n = rf*nzslab           
           do k = 1,nzslab
             do j=1,nyp22
               do i=1,nxp1
                 ff(i,j,k+n) = buff_slabs(i,j,k,m)
               end do
             end do
           end do
	endif   
   endif
  end do
end do

!-------------------------------------------------
!   Solve the tri-diagonal system for Fourier coeffiecients 
!   in the vertical for each subdomain:

do k=1,nzm
    a(k)=rhow(k)/(adz(k)*adzw(k)*dz*dz)
    c(k)=rhow(k+1)/(adz(k)*adzw(k+1)*dz*dz)	 
end do 

call task_rank_to_index(rank,it,jt)
	
ddx2=dx*dx
ddy2=dy*dy
pii = dacos(-1.d0)
xny=ny_gl 	 
xnx=nx_gl
do j=1,nyp22-jwall
 if(dowally) then
    jd=j+jt-1
    facty = 1.d0
 else
    jd=(j+jt-0.1)/2.
    facty = 2.d0
 end if
 xj=jd
 do i=1,nxp1-iwall
  if(dowallx) then
    id=i+it-1
    factx = 1.d0
 else
    id=(i+it-0.1)/2.
    factx = 2.d0
 end if
  xi=id
  eign(i,j)=(2.d0*cos(factx*pii/xnx*xi)-2.d0)/ddx2+ & 
            (2.d0*cos(facty*pii/xny*xj)-2.d0)/ddy2
  if(id+jd.eq.0) then               
     b(i,j)=eign(i,j)*rho(1)-a(1)-c(1)
     alfa(i,j,1)=-c(1)/b(i,j)
     beta(i,j,1)=ff(i,j,1)/b(i,j)
  else
     b(i,j)=eign(i,j)*rho(1)-c(1)
     alfa(i,j,1)=-c(1)/b(i,j)
     beta(i,j,1)=ff(i,j,1)/b(i,j)
  end if
 end do 
end do 

do k=2,nzm-1
 do j=1,nyp22-jwall
  do i=1,nxp1-iwall 
    e=eign(i,j)*rho(k)-a(k)-c(k)+a(k)*alfa(i,j,k-1)
    alfa(i,j,k)=-c(k)/e
    beta(i,j,k)=(ff(i,j,k)-a(k)*beta(i,j,k-1))/e
  end do
 end do
end do

do j=1,nyp22-jwall
  do i=1,nxp1-iwall 
     ff(i,j,nzm)=(ff(i,j,nzm)-a(nzm)*beta(i,j,nzm-1))/ &
	        (eign(i,j)*rho(nzm)-a(nzm)+a(nzm)*alfa(i,j,nzm-1))
  end do
end do
	  
do k=nzm-1,1,-1
  do j=1,nyp22-jwall
    do i=1,nxp1-iwall	 
       ff(i,j,k)=alfa(i,j,k)*ff(i,j,k+1)+beta(i,j,k)
    end do
  end do  
end do 

call task_barrier()

!-----------------------------------------------------------------	 
!   Send the Fourier coefficient to the tasks performing
!   the inverse Fourier transformation:

! Non-blocking receive first:

n_in = 0
do m = 0,nsubdomains-1

  if(rank.lt.npressureslabs.and.m.ne.nsubdomains-1) then
    n_in = n_in + 1
    call task_receive_float(buff_subs(1,1,1,n_in), &
                              nzslab*nxp1*nyp22, reqs_in(n_in))
    flag(n_in) = .false.	    
  endif

  if(rank.lt.npressureslabs.and.m.eq.nsubdomains-1) then

    call task_rank_to_index(rank,it,jt)	  
    n = rank*nzslab
    do k = 1,nzslab
     do j = 1,nyp22-jwall
       do i = 1,nxp1-iwall
         f(i+it,j+jt,k) = ff(i,j,k+n)
       end do
     end do
    end do

  endif

end do ! m

! Blocking send now:

do m = 0,nsubdomains-1

  if(m.lt.npressureslabs.and.m.ne.rank) then
    n = m*nzslab+1
    call task_bsend_float(m,ff(1,1,n),nzslab*nxp1*nyp22, 33)
  endif

end do ! m


! Fill slabs when receive buffers are full:


count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           call task_rank_to_index(rf,it,jt)	  
           do k = 1,nzslab
            do j = 1,nyp22-jwall
             do i = 1,nxp1-iwall
                f(i+it,j+jt,k) = buff_subs(i,j,k,m)
             end do
            end do
           end do
	endif   
   endif
  end do
end do

!-------------------------------------------------
!   Perform inverse Fourier transformation:

if(rank.lt.npressureslabs) then

 do k=1,nzslab

  if(RUN3D) then
   if(dowally) then
     call cosft(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
   else
     call fft991(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
   end if
  end if
	 
  if(dowallx) then
   call cosft(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)
  else
   call fft991(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)
  end if

 end do 

endif

call task_barrier()

!-----------------------------------------------------------------	 
!   Fill the pressure field for each subdomain: 

do i=1,nx_gl
 iii(i)=i
end do
iii(0)=nx_gl
do j=1,ny_gl
 jjj(j)=j
end do
jjj(0)=ny_gl

! Non-blocking receive first:

n_in = 0
do m = 0, nsubdomains-1
		
   call task_rank_to_index(m,it,jt)

   if(m.lt.npressureslabs-1.or.  &
		m.eq.npressureslabs-1.and.rank.ge.npressureslabs) then

     n_in = n_in + 1
     call task_receive_float(bufp_slabs(0,1-YES3D,1,n_in), &
                                nzslab*nxp1*nyp1, reqs_in(n_in))
     flag(n_in) = .false.    

   endif

   if(rank.lt.npressureslabs.and.m.eq.rank) then

     n = rank*nzslab
     do k = 1,nzslab
      do j = 1-YES3D,ny
       jj=jjj(j+jt)
        do i = 0,nx
	 ii=iii(i+it)
          p(i,j,k+n) = f(ii,jj,k) 
        end do
      end do
     end do 

   end if

end do ! m


! Blocking send now:

do m = 0, nsubdomains-1

   call task_rank_to_index(m,it,jt)

   if(rank.lt.npressureslabs.and.m.ne.rank) then

     do k = 1,nzslab
      do j = 1-YES3D,ny
       jj=jjj(j+jt)
       do i = 0,nx
         ii=iii(i+it)
         bufp_subs(i,j,k,1) = f(ii,jj,k)
       end do
      end do
     end do

     call task_bsend_float(m, bufp_subs(0,1-YES3D,1,1), nzslab*nxp1*nyp1,44)

   endif

end do ! m

! Fill the receive buffers:

count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
	call task_test(reqs_in(m), flag(m), rf, tag)
        if(flag(m)) then 
	   count=count-1
           n = rf*nzslab           
           do k = 1,nzslab
            do j=1-YES3D,ny
             do i=0,nx
               p(i,j,k+n) = bufp_slabs(i,j,k,m)
             end do
            end do
           end do
        endif   
   endif
  end do
end do


call task_barrier()

!  Add pressure gradient term to the rhs of the momentum equation:

call press_grad()

end subroutine pressure_orig



