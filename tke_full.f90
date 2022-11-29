
subroutine tke_full

!	this subroutine solves the TKE equation

use vars
use params
implicit none

real def2(nx,ny,nzm)	
real grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Pr1,Cee,Cs,Cs1
real buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
real lstarn, lstarp, bbb, omn, omp
real qsat,dqsat
integer i,j,k,kc,kb

Cs = 0.1944
Cs1 = 0.14
Pr = 3.0
Ck=0.1
Ce=Ck**3/Cs**4
Ces=Ce/0.7*3.0	

if(RUN3D) then
  call shear_prod3D(def2)
else
  call shear_prod2D(def2)
endif

do k=1,nzm      
  kb=k-1
  kc=k+1

  grd=dz*adz(k)

  betdz=bet(k)/dz/(adzw(kc)+adzw(k))
  Ce1=Ce/0.7*0.19
  Ce2=Ce/0.7*0.51
  if(k.eq.1) then
    kb=1
    kc=2
    betdz=bet(k)/dz/adzw(kc)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  if(k.eq.nzm) then
    kb=nzm-1
    kc=nzm
    betdz=bet(k)/dz/adzw(k)
    Ce1=Ces/0.7*0.19
    Ce2=Ces/0.7*0.51
  end if
  tkelediss(k) = 0.
  tkesbdiss(k) = 0.
  tkesbshear(k)= 0.
  tkesbbuoy(k) = 0.
  do j=1,ny
  do i=1,nx
!  SGS buoyancy flux

   omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg)) 
   omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr)) 
   lstarn = fac_cond+(1.-omn)*fac_fus 
   lstarp = fac_cond+(1.-omp)*fac_fus 

!   if(qn(i,j,kb).gt.0. .and. qn(i,j,k).gt.0. .and. qn(i,j,kc).gt.0.) then
   if(qn(i,j,k) .gt. 0.) then
     
      dqsat = omn*dtqsatw(tabs(i,j,k),pres(k))+ &
                             (1.-omn)*dtqsati(tabs(i,j,k),pres(k))
      qsat = omn*qsatw(tabs(i,j,k),pres(k))+(1.-omn)*qsati(tabs(i,j,k),pres(k))
      bbb = 1. + 0.61*qsat-qn(i,j,k) -qp(i,j,k)+1.61*tabs(i,j,k)*dqsat
      bbb = bbb / (1.+lstarn*dqsat)
      buoy_sgs=betdz*(bbb*(t(i,j,kc)-t(i,j,kb)) &
	+(bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))*(q(i,j,kc)-q(i,j,kb)) &
	+(bbb*lstarp - (1.+lstarp*dqsat)*tabs(i,j,k))*(qp(i,j,kc)-qp(i,j,kb)) )
   else

      bbb = 1.+0.61*q(i,j,k)-qp(i,j,k)
      buoy_sgs=betdz*( bbb*(t(i,j,kc)-t(i,j,kb)) &
        +0.61*tabs(i,j,k)*(q(i,j,kc)-q(i,j,kb)) &
	+(bbb*lstarp-tabs(i,j,k))*(qp(i,j,kc)-qp(i,j,kb)) )    	     
   end if


   tke(i,j,k)=max(0.,tke(i,j,k))

   if(buoy_sgs.le.0.) then
     smix=grd
   else
     smix=min(grd,max(0.1*grd,0.76*sqrt(tke(i,j,k)/buoy_sgs+1.e-10)))
   end if

   ratio=smix/grd
   Pr1=1.+2.*ratio
   Cee=Ce1+Ce2*ratio

   if(dosmagor) then

     tk(i,j,k)=sqrt(Ck**3/Cee*max(0.,def2(i,j,k)-Pr1*buoy_sgs))*smix**2
     if(.not.doscalar) tke(i,j,k) = (tk(i,j,k)/(Ck*smix))**2
     a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
     a_prod_bu=-(tk(i,j,k)+0.001)*Pr1*buoy_sgs
     a_diss=a_prod_sh+a_prod_bu

   else

     a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
     a_prod_bu=-(tk(i,j,k)+0.001)*Pr1*buoy_sgs
     a_diss=Cee/smix*tke(i,j,k)**1.5	   
     tke(i,j,k)=max(0.,tke(i,j,k)+dtn*(max(0.,a_prod_sh+a_prod_bu)-a_diss))
     tk(i,j,k)=Ck*smix*sqrt(tke(i,j,k))

   end if
	
   tkh(i,j,k)=Pr1*tk(i,j,k)

   tkelediss(k) = tkelediss(k) - a_prod_sh
   tkesbdiss(k) = tkesbdiss(k) + a_diss
   tkesbshear(k)= tkesbshear(k)+ a_prod_sh
   tkesbbuoy(k) = tkesbbuoy(k) + a_prod_bu

  end do ! i
  end do ! j

  tkelediss(k) = tkelediss(k)/float(nx*ny)

end do ! k

end


