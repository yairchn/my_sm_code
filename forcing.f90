
subroutine forcing
	
use vars
use params
use simple_ocean, only: sst_evolve

implicit none

integer i,j,k,k1,k2
real coef, radtend

 do k=1,nzm
   total_water_ls = total_water_ls - &
              sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
 end do


! Large-scale sounding:

   k1=1
   k2=1
   do i=1,nsnd-1
    if(day.gt.daysnd(i)) then
      k2=i+1
      k1=i
    endif
   end do
   coef=(day-daysnd(k1))/(daysnd(k2)-daysnd(k1))
   do k=1,nzm
    tg0(k)=tsnd(k,k1)+(tsnd(k,k2)-tsnd(k,k1))*coef
    qg0(k)=qsnd(k,k1)+(qsnd(k,k2)-qsnd(k,k1))*coef
    qg0(k)=qg0(k)*1.e-3
! Note that ug0 and vg0 maybe reset if dolargescale is true)
    ug0(k)=usnd(k,k1)+(usnd(k,k2)-usnd(k,k1))*coef - ug
    vg0(k)=vsnd(k,k1)+(vsnd(k,k2)-vsnd(k,k1))*coef - vg
   end do

! Initialize tendencies:


do k=1,nzm
 ttend(k)=0.
 qtend(k)=0.
end do


! Large-Scale Advection Forcing:


if(dolargescale.and.time.gt.timelargescale) then

   k1=1
   k2=1
   do i=1,nlsf-1
    if(day.gt.dayls(i)) then
      k2=i+1
      k1=i
    endif
   end do
   coef=(day-dayls(k1))/(dayls(k2)-dayls(k1))

   do k=1,nzm
    ttend(k)=dtls(k,k1)+(dtls(k,k2)-dtls(k,k1))*coef
    qtend(k)=dqls(k,k1)+(dqls(k,k2)-dqls(k,k1))*coef
    ug0(k)=ugls(k,k1)+(ugls(k,k2)-ugls(k,k1))*coef - ug
    vg0(k)=vgls(k,k1)+(vgls(k,k2)-vgls(k,k1))*coef - vg
    wsub(k)=wgls(k,k1)+(wgls(k,k2)-wgls(k,k1))*coef
    do j=1,ny
     do i=1,nx
      t(i,j,k)=t(i,j,k)+ttend(k) * dtn
      q(i,j,k)=q(i,j,k)+qtend(k) * dtn
     end do
    end do
   end do 

   if(dosubsidence) call subsidence()

end if 


! Prescribed Radiation Forcing:


if(doradforcing.and.time.gt.timelargescale) then

   k1=1
   k2=1
   do i=1,nrfc-1
    if(day.gt.dayrfc(i)) then
      k2=i+1
      k1=i
    endif
   end do
   coef=(day-dayrfc(k1))/(dayrfc(k2)-dayrfc(k1))
   do k=1,nzm
    radtend=dtrfc(k,k1)+(dtrfc(k,k2)-dtrfc(k,k1))*coef
    radqrlw(k)=radtend*float(nx*ny)
    radqrsw(k)=0.
    do j=1,ny
     do i=1,nx
       t(i,j,k)=t(i,j,k)+radtend*dtn 
     end do
    end do
   end do

endif



! Surface flux forcing:

if(dosfcforcing.and.time.gt.timelargescale) then

   k1=1
   k2=1
   do i=1,nsfc-1
     if(day.gt.daysfc(i)) then
       k2=i+1
       k1=i
     endif
   end do
   if(doxy) then 
     if(masterproc)  print*,'doxy=',doxy,': no sfc forcing is allowed'
     call task_abort()
   end if
   coef=(day-daysfc(k1))/(daysfc(k2)-daysfc(k1))
   tabs_s=sstsfc(k1)+(sstsfc(k2)-sstsfc(k1))*coef
   fluxt0=(hsfc(k1)+(hsfc(k2)-hsfc(k1))*coef)/(rho(1)*cp)
   fluxq0=(lesfc(k1)+(lesfc(k2)-lesfc(k1))*coef)/(rho(1)*lcond)
   tau0=tausfc(k1)+(tausfc(k2)-tausfc(k1))*coef
   do j=1,ny
    do i=1,nx
      sstxy(i,j) = tabs_s
    end do
   end do

   if(dostatis) then
     sstobs = tabs_s  ! sst is not averaged over the sampling period
     lhobs = lhobs + fluxq0 * rho(1)*lcond
     shobs = shobs + fluxt0 * rho(1)*cp
   end if

endif

if(.not.dosfcforcing.and.dodynamicocean) then

   call sst_evolve()

endif


  do k=1,nzm
     total_water_ls = total_water_ls + &
                     sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
  end do


end subroutine forcing
