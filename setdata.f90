subroutine setdata()
	
use vars
use params
use simple_ocean, only: set_sst
implicit none
	
integer ndmax,n,i,j,k,iz,it,jt
real tmp, qsat0(nzm), dqsat0(nzm), om(nzm), dtabs, rel(nzm)
real presr(nz)	
parameter (ndmax = 1000)
real zz(ndmax),tt(ndmax),qq(ndmax),uu(ndmax),vv(ndmax) 
real zz1(ndmax),tt1(ndmax),qq1(ndmax),uu1(ndmax),vv1(ndmax) 
real rrr1,rrr2, pres1, pp(ndmax),ta(ndmax)
real pp1(ndmax)
real ratio_t1,ratio_t2,ratio_p1,ratio_p2
real tpert0(ndmax), qpert0(ndmax)
real latit,long
logical zgrid
!-------------------------------------------------------------
!	read subensemble perturbation file first:

if(doensemble) then
	
open(76,file='./'//trim(case)//'/tqpert',status='old',form='formatted')  
read(76,*)
  do j=0,nensemble
    read(76,*) i,n
    do i=1,n
      read(76,end=766,fmt=*) pp(i),tpert0(i),qpert0(i)
      tpert0(i)=tpert0(i)*(1000./pp(i))**(rgas/cp)
    end do
  end do
  close(76)
  if(masterproc) then
    print*,'Subensemble run. nensemble=',nensemble
    print*,'tpert:',(tpert0(i),i=1,n)
    print*,'qpert:',(qpert0(i),i=1,n)
  end if
  goto 767
766  print*,'Error: nensemble is too large.'  
  call task_abort()
767  continue
else
  do i=1,ndmax
    tpert0(i)=0.
    qpert0(i)=0.
  end do
end if


!**************************************************************
!	Read Initial Sounding

open(77,file='./'//trim(case)//'/snd',status='old',form='formatted')
read(77,*)

do while(.true.)

  read(77,err=55,end=55,fmt=*) rrr1, n, pres0
  do i=1,n
      read(77,*) zz(i),pp(i),tt(i),qq(i),uu(i),vv(i)
  end do	      
  read(77,err=55,end=55,fmt=*) rrr2, n, pres1
  do i=1,n
      read(77,*) zz1(i),pp1(i),tt1(i),qq1(i),uu1(i),vv1(i)
  end do	      

  if(day.ge.rrr1.and.day.le.rrr2) then
    if(zz(2).gt.zz(1)) then
      zgrid = .true.
      do i=1,n
      zz(i)=zz(i)+(zz1(i)-zz(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      tt(i)=tt(i)+(tt1(i)-tt(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      qq(i)=qq(i)+(qq1(i)-qq(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      uu(i)=uu(i)+(uu1(i)-uu(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      vv(i)=vv(i)+(vv1(i)-vv(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      tt(i)=tt(i)+tpert0(i)
      qq(i)=qq(i)+qpert0(i) 
      end do
    else if(pp(2).lt.pp(1)) then
      zgrid = .false.
      do i=1,n
      pp(i)=pp(i)+(pp1(i)-pp(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      tt(i)=tt(i)+(tt1(i)-tt(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      qq(i)=qq(i)+(qq1(i)-qq(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      uu(i)=uu(i)+(uu1(i)-uu(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      vv(i)=vv(i)+(vv1(i)-vv(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
      tt(i)=tt(i)+tpert0(i)
      qq(i)=qq(i)+qpert0(i) 
      ta(i)=tt(i)*(pp(i)/1000.)**(rgas/cp)
      end do
    else  
      if(masterproc) print*,'vertical grid is undefined...'
      stop
    end if
    pres0=pres0+(pres1-pres0)/(rrr2-rrr1+1.e-5)*(day-rrr1)      
    goto 56
  endif
  do i=1,n+1
    backspace(77)
  end do	      

end do

55 continue
if(masterproc) then
  print*,'Error: day is beyond the sounding time range'
  print*,day,rrr1,rrr2
end if
call task_abort()
56 continue	

close (77)

if(masterproc) then
  print *	
  print *,'surface pressure: ',pres0
endif  

! compute heights from pressure:

if(.not.zgrid) then
 zz(1) = rgas/ggr*ta(1)*log(pres0/pp(1))
 do i=2,n
  zz(i)=zz(i-1)+0.5*rgas/ggr*(ta(i)+ta(i-1))*log(pp(i-1)/pp(i))
 end do
end if  	
!-----------------------------------------------------------
!       Interpolate sounding into vertical grid:

presr(1)=(pres0/1000.)**(rgas/cp)
presi(1)=pres0
do k= 1,nzm
 do iz = 2,n
  if(z(k).le.zz(iz)) then
    t0(k)=tt(iz-1)+(tt(iz)-tt(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1))	
    q0(k)=qq(iz-1)+(qq(iz)-qq(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1))
    u0(k)=uu(iz-1)+(uu(iz)-uu(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1))  
    v0(k)=vv(iz-1)+(vv(iz)-vv(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1)) 
    goto 12
  endif
 end do

!  Utilize 1976 standard atmosphere for points above sounding:


 call atmosphere(z(k-1)/1000.,ratio_p1,rrr1,ratio_t1)
 call atmosphere(z(k)/1000.,ratio_p2,rrr1,ratio_t2)

 tabs0(k)=ratio_t2/ratio_t1*tabs0(k-1)
 presi(k+1)=presi(k)*exp(-ggr/rgas/tabs0(k)*(zi(k+1)-zi(k)))
 pres(k) = 0.5*(presi(k)+presi(k+1))
 prespot(k)=(1000./pres(k))**(rgas/cp)
 q0(k)=max(0.,2.*q0(k-1)-q0(k-2))
 u0(k)=u0(k-1)
 v0(k)=v0(k-1)
 goto 13
12 continue
 q0(k)=q0(k)*1.e-3
 tv0(k)=t0(k)*(1.+0.61*q0(k))
 presr(k+1)=presr(k)-ggr/cp/tv0(k)*(zi(k+1)-zi(k))
 presi(k+1)=1000.*presr(k+1)**(cp/rgas)
 pres(k) = exp(log(presi(k))+log(presi(k+1)/presi(k))* &
                             (z(k)-zi(k))/(zi(k+1)-zi(k)))
 prespot(k)=(1000./pres(k))**(rgas/cp)
 tabs0(k)=t0(k)/prespot(k)
13 continue
 ug0(k)=u0(k)
 vg0(k)=v0(k)
end do


! recompute pressure levels (for consistancy):

!	call pressz()
        
!-------------------------------------------------------------
!       Initial thernodynamic profiles: 
	
do k=1,nzm

  gamaz(k)=ggr/cp*z(k)
  t0(k) = tabs0(k)+gamaz(k) 
  ttend(k) = 0.
  qtend(k) = 0.
  wsub(k)=0.

  om(k) = (tabs0(k) - tbgmin)/(tbgmax - tbgmin)
  if(om(k).le.0.) then
     qsat0(k) = qsati(tabs0(k),pres(k))
     dqsat0(k) = dtqsati(tabs0(k),pres(k))
     tmp = lsub/cp
     om(k) = 0.
  else if(om(k).ge.1.) then
     qsat0(k) = qsatw(tabs0(k),pres(k))
     dqsat0(k) = dtqsatw(tabs0(k),pres(k))
     tmp = lcond/cp               
     om(k) = 1.
  else
     qsat0(k) = om(k)*qsatw(tabs0(k),pres(k))+(1.-om(k))*qsati(tabs0(k),pres(k))
     dqsat0(k) = om(k)*dtqsatw(tabs0(k),pres(k)) + &
                           (1.- om(k))*dtqsati(tabs0(k),pres(k))
     tmp = (om(k)*lcond + (1-om(k))*lsub)/cp
  endif
  if(qsat0(k).lt.q0(k)) then
     dtabs = tmp*(q0(k)-qsat0(k))/(1.+tmp*dqsat0(k))
     tabs0(k) = tabs0(k) + dtabs
     qsat0(k) = qsat0(k) + dqsat0(k) * dtabs
     qv0(k) = qsat0(k)
     qc0(k) = om(k)*(q0(k) - qv0(k))
     qi0(k) = (1.-om(k))*(q0(k) - qv0(k))
  else
     qv0(k) = q0(k)
     qc0(k) = 0.
     qi0(k) = 0.
  endif 

  rel(k) = qv0(k)/qsat0(k)*100.
  rho(k) = (presi(k)-presi(k+1))/(zi(k+1)-zi(k))/ggr*100.
  bet(k) = ggr/tabs0(k)
  betp(k)= ggr/(pres(k)*100.)
 
  u0(k) = u0(k) - ug
  v0(k) = v0(k) - vg
  ug0(k) = ug0(k) - ug
  vg0(k) = vg0(k) - vg
end do
do k=2,nzm
  rhow(k) =  (pres(k-1)-pres(k))/(z(k)-z(k-1))/ggr*100.
end do
rhow(1) = 2*rhow(2) - rhow(3)
rhow(nz)= 2*rhow(nzm) - rhow(nzm-1)

wsub(nz)=0.

if(masterproc) then
 print *,'Initial Sounding:'
 print *,'k height rho rhoi s h h* qt u v adz'
 do k=nzm,1,-1
   write(6,'(f7.1,2f7.3,f7.2,f7.2,4f7.2,4g11.4)')z(k),rho(k),rhow(k),t0(k), &
          t0(k)+lcond/cp*qv0(k), t0(k)+lcond/cp*qsatw(tabs0(k),pres(k)), &
		q0(k)*1.e3,u0(k),v0(k), adz(k)
 end do  

 print *  
 print *,'z  dz pres presi Tabs  s  tp  qt Qc Qi  REL' 
 do k = nzm,1,-1
  write(6,'(7f8.2,3f8.4,f8.2)')   z(k),zi(k+1)-zi(k),pres(k),presi(k),tabs0(k), &
     t0(k),tabs0(k)*prespot(k),q0(k)*1.e3, qc0(k)*1.e3,qi0(k)*1.e3,rel(k) 
 end do
endif

do k=1,nzm
 do j=1,ny
  do i=1,nx
   u(i,j,k)= u0(k)
   v(i,j,k)= v0(k)
   w(i,j,k)= 0.
   t(i,j,k)= t0(k)
   q(i,j,k)=q0(k)
   tabs(i,j,k)=tabs0(k)
   qn(i,j,k)=0.
   qp(i,j,k)=0.
   tke(i,j,k)=0.
   tk(i,j,k)=0.
   tkh(i,j,k)=0.
   p(i,j,k)=0.
   w(i,j,nz)=0.
   fluxbu(i,j)=0.
   fluxbv(i,j)=0.
   fluxbt(i,j)=0.
   fluxbq(i,j)=0.
   fluxtu(i,j)=0.
   fluxtv(i,j)=0.
   fluxtt(i,j)=0.
   fluxtq(i,j)=0.
   fzero(i,j) =0.
   precsfc(i,j)=0.
   prec_xy(i,j)=0.
   shf_xy(i,j)=0.
   lhf_xy(i,j)=0.
   lwns_xy(i,j)=0.
   swns_xy(i,j)=0.
   lwnsc_xy(i,j)=0.
   swnsc_xy(i,j)=0.
   lwnt_xy(i,j)=0.
   swnt_xy(i,j)=0.
   lwntc_xy(i,j)=0.
   swntc_xy(i,j)=0.
   solin_xy(i,j)=0.
   qocean_xy(i,j)=0.
   u200_xy(i,j) = 0.
   v200_xy(i,j) = 0.
   usfc_xy(i,j) = 0.
   vsfc_xy(i,j) = 0.
   w500_xy(i,j) = 0.
   pw_xy(i,j)=0.
   cw_xy(i,j)=0.
   iw_xy(i,j)=0.
   sstxy(i,j)=0.
   fcory(j) = fcor
   fcorzy(j) = fcorz
   ttend_cup(i,j,k) = 0.
   qtend_cup(i,j,k) = 0.
   utend_cup(i,j,k) = 0.
   vtend_cup(i,j,k) = 0.
  end do 
 end do 
end do 

dudt = 0.
dvdt = 0.
dwdt = 0.	   
tke = 0.
	
call setperturb()
if(doprecip) call precip_init()
if(docloud) call cloud()

if(doxy) then

 call task_rank_to_index(rank,it,jt) 

 do j=1,ny
   latit=latitude0+dy*(j+jt-(ny_gl+YES3D-1)/2-1)*2.5e-8*360. 
   fcory(j)= 4.*pi/86400.*sin(latit*pi/180.)
   fcorzy(j) = sqrt(4.*(2*pi/(3600.*24.))**2-fcory(j)**2)
   do i=1,nx
    long = longitude0+dx*(i+it-nx_gl/2-1)*2.5e-8*360. 
    latitude(i,j) = latit
    longitude(i,j) = long
   end do
 end do

 call set_sst(ocean_type)

else

 do j=1,ny
  do i=1,nx
    latitude(i,j) = latitude0
    longitude(i,j) = longitude0
  end do
 end do 

end if ! doxy

end

