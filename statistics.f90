
subroutine statistics()

use vars
use tracers
use params
use hbuffer
use isccp, only: isccp_get
implicit none	
	
	real mse(nzm)
	real dse(nzm)
	real sse(nzm)
	real tpz(nzm)
	real tlz(nzm)
	real tvz(nzm)
	real qcz(nzm)
	real qiz(nzm)
	real tez(nzm)
	real qvz(nzm)
	real qrz(nzm)
	real qsz(nzm)
	real qgz(nzm)
	real relhz(nzm)

	real u2z(nzm)
	real v2z(nzm)
	real w2z(nzm)
	real w22(nzm)
	real w3z(nzm)
	real skw(nzm)
	real t2z(nzm)
	real tqz(nzm)
	real q2z(nzm)
	real s2z(nzm)
	real qc2z(nzm)
	real qi2z(nzm)
	real qs2z(nzm)
	real tkez(nzm)
	real fadv(nz)
	real shear(nz)
	real shearx(nzm)
	real sheary(nzm)
	real presx(nzm)
	real presy(nzm)
	real twgrad(nzm)
	real qwgrad(nzm)
	real swgrad(nzm)
	
	
	real tvwle(nzm)
	real qcwle(nzm)
	real qiwle(nzm)
	real aup(nzm)
	real wcl(nzm)
	real ucl(nzm)
	real vcl(nzm)
	real tcl(nzm)
	real tacl(nzm)
	real tvcl(nzm)
	real qcl(nzm)
	real qccl(nzm)
	real qicl(nzm)
	real qpcl(nzm)
	real twcl(nzm)
	real swcl(nzm)
	real qwcl(nzm)
	real tvwcl(nzm)
	real qcwcl(nzm)
	real qiwcl(nzm)
	real scl(nzm)
	real wacl(nzm)	
	real cld(nzm)
 	real cldd(nzm)
	real hydro(nzm)
	real qsatwz(nzm)

	real tvirt(nx,ny,nzm)
	
	integer i,j,k,n,ntr
	real qcc,qii,qrr,qss,qgg,lstarn,lstarp,coef,coef1
	real factor_xy, factor_n, tmp(4)
        real buffer(nzm,6),buffer1(nzm,6)
	real prof1(nzm),prof2(nzm),prof3(nzm),prof4(nzm)	
	real cwp(nx,ny),cwpl(nx,ny),cwpm(nx,ny),cwph(nx,ny)
        real cwpmax,prec_conv, prec_strat
	real omn, omn1, omp, omg, qsat
	logical condition, condition_cl
	real zero(nzm)

	integer topind(nx,ny)	
	factor_xy = 1./float(nx*ny)
	factor_n = 1./float(nsubdomains)
	
!-----------------------------------------------
!	Mean thermodynamics profiles:
!-----------------------------------------------	
		
	do k=1,nzm
	 dse(k)=0.
	 mse(k)=0.
	 sse(k)=0.
	 tpz(k) = 0.
	 tlz(k) = 0.
	 tvz(k) = 0.
	 tez(k) = 0.
	 qvz(k) = 0.
	 qcz(k) = 0.
	 qiz(k) = 0.
	 qrz(k) = 0.
	 qsz(k) = 0.
	 qgz(k) = 0.
	 qsatwz(k)=0.
	 relhz(k)=0.
	 prof1(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
	 zero(k)=0.
	 do j=1,ny
	  do i=1,nx
	   omn = omegan(tabs(i,j,k))
	   omp = omegap(tabs(i,j,k))
	   omg = omegag(tabs(i,j,k))
	   lstarn=fac_cond+(1.-omn)*fac_fus
	   lstarp=fac_cond+(1.-omp)*fac_fus
	   qcc=qn(i,j,k)*omn
	   qii=qn(i,j,k)*(1.-omn)
	   qrr=qp(i,j,k)*omp
	   qss=qp(i,j,k)*(1.-omp)*(1.-omg)
	   qgg=qp(i,j,k)*(1.-omp)*omg
           qrz(k)=qrz(k)+qrr
           qsz(k)=qsz(k)+qss
           qgz(k)=qgz(k)+qgg
           qcz(k)=qcz(k)+qcc
           qiz(k)=qiz(k)+qii
	   prof1(k)=prof1(k)+qcc+qii
	   prof2(k)=prof2(k)+qrr+qss+qgg
	   prof3(k)=prof3(k)+qcc+qii+qrr+qss+qgg
	   tmp(1)=tabs(i,j,k)*prespot(k)
           tpz(k)=tpz(k)+tmp(1)
	   tlz(k)=tlz(k)+tmp(1)*(1.-fac_cond*qn(i,j,k)/tabs(i,j,k))
           tvirt(i,j,k)=tmp(1)*(1.+0.61*q(i,j,k)-1.61*(qn(i,j,k))-qp(i,j,k))
	   tvz(k)=tvz(k)+tvirt(i,j,k)
           tez(k)=tez(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k) &
                     +fac_cond*(q(i,j,k)-qcc-qii)-fac_fus*(qii+qss+qgg)
	   qvz(k) =qvz(k)+q(i,j,k)-qn(i,j,k)
     	   dse(k)=dse(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	   mse(k)=mse(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k) &
                	+fac_cond*(q(i,j,k)-qn(i,j,k))
	   qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
             (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	   sse(k)=sse(k)+t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k) &
                        +fac_cond*qsat
	   qsatwz(k) = qsatwz(k)+qsatw(tabs(i,j,k),pres(k))
	   relhz(k)=relhz(k)+(q(i,j,k)-qn(i,j,k))/qsatw(tabs(i,j,k),pres(k))
	  end do
	 end do
	end do	
	

	call hbuf_avg_put('TL',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	call hbuf_avg_put('TABS',tabs,1,nx, 1,ny, nzm,1.)
	call hbuf_avg_put('U',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1.)
	call hbuf_avg_put('V',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1.)
	call hbuf_avg_put('QT',q,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.e3)
	call hbuf_avg_put('TK',tk,0,nxp1, (1-YES3D),nyp1, nzm,1.)
	call hbuf_avg_put('TKH',tkh,0,nxp1, (1-YES3D),nyp1, nzm,1.)
	call hbuf_put('TABSOBS',tg0,1.)
	call hbuf_put('QVOBS',qg0,1.e3)
	call hbuf_put('UOBS',ug0,1.)
	call hbuf_put('VOBS',vg0,1.)
        call hbuf_put('WOBS',wsub,1.)
	call hbuf_put('TTEND',ttend,86400.)
	call hbuf_put('QTEND',qtend,86400.*1.e3)

	call hbuf_put('DSE',dse,factor_xy)
	call hbuf_put('MSE',mse,factor_xy)
	call hbuf_put('SSE',sse,factor_xy)
	call hbuf_put('THETA',tpz,factor_xy)
	call hbuf_put('THETAL',tlz,factor_xy)
	call hbuf_put('THETAV',tvz,factor_xy)
	call hbuf_put('THETAE',tez,factor_xy)

	call hbuf_put('PRES',pres,1.)
	call hbuf_put('RHO',rho,1.)
	call hbuf_put('QV',qvz,1.e3*factor_xy)
	call hbuf_put('QC',qcz,1.e3*factor_xy)
	call hbuf_put('QI',qiz,1.e3*factor_xy)
	call hbuf_put('QR',qrz,1.e3*factor_xy)
	call hbuf_put('QS',qsz,1.e3*factor_xy)
	call hbuf_put('QG',qgz,1.e3*factor_xy)
	call hbuf_put('QN',prof1,1.e3*factor_xy)
	call hbuf_put('QP',prof2,1.e3*factor_xy)
	call hbuf_put('QCOND',prof3,1.e3*factor_xy)
	call hbuf_put('QSAT',qsatwz,1.e3*factor_xy)
	call hbuf_put('RELH',relhz,100.*factor_xy)

!-------------------------------------------------------------
!	Fluxes:
!-------------------------------------------------------------

	do k=1,nzm
	  tmp(1) = dz/rhow(k)
	  tmp(2) = tmp(1) / dtn
	  uwsb(k) = uwsb(k) * tmp(1)
	  vwsb(k) = vwsb(k) * tmp(1)
	  twsb(k) = twsb(k) * tmp(1) * rhow(k) * cp
	  qwsb(k) = qwsb(k) * tmp(1) * rhow(k) * lcond
	  uwle(k) = uwle(k)*tmp(1) + uwsb(k)
	  vwle(k) = vwle(k)*tmp(1) + vwsb(k)
	  twle(k) = twle(k)*tmp(2)*rhow(k)*cp + twsb(k)
	  qwle(k) = qwle(k)*tmp(2)*rhow(k)*lcond + qwsb(k)
	  if(docloud.and.doprecip) then
	    qpwsb(k) = qpwsb(k) * tmp(1)*rhow(k)*lcond
     	    qpwle(k) = qpwle(k)*tmp(2)*rhow(k)*lcond + qpwsb(k)
	  else
	    qpwsb(k) = 0.
	    qpwle(k) = 0.
	    precflux(k) = 0.
	  endif
	  if(dotracers) then
           do ntr=1,ntracers
	    trwsb(k,ntr) = trwsb(k,ntr) * tmp(1)*rhow(k)
     	    trwle(k,ntr) = trwle(k,ntr) * tmp(2)*rhow(k) + trwsb(k,ntr)
           end do
          end if
	end do
	uwle(nz) = 0.
	vwle(nz) = 0.
	uwsb(nz) = 0.
	vwsb(nz) = 0.

	call hbuf_put('UW',uwle,factor_xy)
	call hbuf_put('VW',vwle,factor_xy)
	call hbuf_put('UWSB',uwsb,factor_xy)
	call hbuf_put('VWSB',vwsb,factor_xy)
	call hbuf_put('TLFLUX',twle,factor_xy)
	call hbuf_put('QTFLUX',qwle,factor_xy)
	call hbuf_put('TLFLUXS',twsb,factor_xy)
	call hbuf_put('QTFLUXS',qwsb,factor_xy)
	call hbuf_put('QPFLUXS',qpwsb,factor_xy)
	call hbuf_put('QPFLUX',qpwle,factor_xy)
	call hbuf_put('PRECIP',precflux,factor_xy/dt*dz*86400./(nstatis+1.e-5))
	
	do j=1,ny
	 do i=1,nx
	  precsfc(i,j)=precsfc(i,j)*dz/dt*3600./(nstatis+1.e-5)
	 end do
	end do

	call stat_precip(precsfc, prec_conv, prec_strat)
	p_conv = p_conv + prec_conv
	p_strat = p_strat + prec_strat

	do k=1,nzm
	 tvz(k) = 0.
	 qcz(k) = 0.
	 qiz(k) = 0.
	 qsatwz(k) = 0.
	 prof1(k)=0.
	 prof2(k)=0.
	 do j=1,ny
	  do i=1,nx
	    omn = omegan(tabs(i,j,k))
	    tvz(k) = tvz(k) + tvirt(i,j,k)
	    qcz(k) = qcz(k) + qn(i,j,k)*omn
	    qiz(k) = qiz(k) + qn(i,j,k)*(1.-omn)
	    qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                       (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	    qsatwz(k) = qsatwz(k)+qsat
	  end do
	 end do
	 tvz(k) = tvz(k)*factor_xy
	 qcz(k) = qcz(k)*factor_xy
	 qiz(k) = qiz(k)*factor_xy	 
	 qsatwz(k) = qsatwz(k)*factor_xy
	end do
	if(dompi) then
	  coef1 = 1./float(nsubdomains)
	  do k=1,nzm
	    buffer(k,1) = tvz(k)
	    buffer(k,2) = qcz(k)
	    buffer(k,3) = qiz(k)
	    buffer(k,4) = qsatwz(k)
	  end do
	  call task_sum_real(buffer,buffer1,nzm*4)
	  do k=1,nzm
	    tvz(k) = buffer1(k,1) * coef1
	    qcz(k) = buffer1(k,2) * coef1
	    qiz(k) = buffer1(k,3) * coef1
	    qsatwz(k) = buffer1(k,4) * coef1
	  end do
	end if ! dompi

	tvwle(1) = 0.
	qcwle(1) = 0.
	qiwle(1) = 0.
	do k=2,nzm
	 tvwle(k) = 0.
	 qcwle(k) = 0.
	 qiwle(k) = 0.
	 do j=1,ny
	  do i=1,nx
	    omn = omegan(tabs(i,j,k))
	    omn1 = omegan(tabs(i,j,k-1))
	    tvwle(k) = tvwle(k) + 0.5*w(i,j,k)* &
		(tvirt(i,j,k-1)-tvz(k-1)+tvirt(i,j,k)-tvz(k))
	    qcwle(k) = qcwle(k) + 0.5*w(i,j,k)* &
		(qn(i,j,k-1)*omn1-qcz(k-1)+ qn(i,j,k)*omn-qcz(k))	  
	    qiwle(k) = qiwle(k) + 0.5*w(i,j,k)* &
                (qn(i,j,k-1)*(1.-omn1)-qiz(k-1)+qn(i,j,k)*(1.-omn)-qiz(k))
	    prof1(k)=prof1(k)+rho(k)*0.5* &
                (w(i,j,k)**2+w(i,j,k+1)**2)*(t(i,j,k)-t0(k))
          end do
	 end do
	 tvwle(k) = tvwle(k)*rhow(k)*cp
	 qcwle(k) = qcwle(k)*rhow(k)*lcond
	 qiwle(k) = qiwle(k)*rhow(k)*lcond
	end do	

	call hbuf_put('TVFLUX',tvwle,factor_xy)
	call hbuf_put('QCFLUX',qcwle,factor_xy)
	call hbuf_put('QIFLUX',qiwle,factor_xy)

!---------------------------------------------------------
!	Mean turbulence related profiles:
!-----------------------------------------------------------


	do k=1,nzm

	 u2z(k) = 0.
	 v2z(k) = 0.
	 w2z(k) = 0.
	 w22(k) = 0.
	 w3z(k) = 0.
	 aup(k) = 0.
	 t2z(k) = 0.
	 tqz(k) = 0.
	 q2z(k) = 0.
	 s2z(k) = 0.
	 qc2z(k) = 0.
	 qi2z(k) = 0.
	 qs2z(k) = 0.
	 do j=1,ny
	  do i=1,nx
	    u2z(k) = u2z(k)+(u(i,j,k)-u0(k))**2	  
	    v2z(k) = v2z(k)+(v(i,j,k)-v0(k))**2
	    w2z(k) = w2z(k)+0.5*(w(i,j,k+1)**2+w(i,j,k)**2)
	    w22(k) = w22(k)+w(i,j,k)**2
	    w3z(k) = w3z(k)+0.5*(w(i,j,k+1)**3+w(i,j,k)**3)	  
	    t2z(k) = t2z(k)+(t(i,j,k)-t0(k))**2	  
	    tqz(k) = tqz(k)+(t(i,j,k)-t0(k))*(q(i,j,k)-q0(k))	  
	    q2z(k) = q2z(k)+(q(i,j,k)-q0(k))**2
	    s2z(k) = s2z(k)+(tke(i,j,k)-tke0(k))**2
	    if(w(i,j,k)+w(i,j,k+1).gt.0) aup(k) = aup(k) + 1
	  end do
	 end do
	 skw(k) = w3z(k)/(w2z(k)*factor_xy+1.e-5)**1.5
	 tkez(k)= 0.5*(u2z(k)+v2z(k)*YES3D+w2z(k))
	 tvwle(k) = tvwle(k) * bet(k) /(rho(k)*cp)
	 do j=1,ny
	  do i=1,nx
	    omn = omegan(tabs(i,j,k))
	    qc2z(k) = qc2z(k)+(qn(i,j,k)*omn-qcz(k))**2
	    qi2z(k) = qi2z(k)+(qn(i,j,k)*(1.-omn)-qiz(k))**2
	    qsat = omn*qsatw(tabs(i,j,k),pres(k))+ &
                    (1.-omn)*qsati(tabs(i,j,k),pres(k)) 
	    qs2z(k) = qs2z(k)+(qsat-qsatwz(k))**2
	  end do
	 end do
	 
	end do

	call hbuf_put('U2',u2z,factor_xy)
	call hbuf_put('V2',v2z,factor_xy)
	call hbuf_put('W2',w2z,factor_xy)
	call hbuf_put('W3',w3z,factor_xy)
	call hbuf_put('WSKEW',skw,factor_xy)
	call hbuf_put('AUP',aup,factor_xy)

	call hbuf_put('TL2',t2z,factor_xy)
	call hbuf_put('TQ',tqz,factor_xy)
	call hbuf_put('QT2',q2z,1.e6*factor_xy)
	call hbuf_put('QC2',qc2z,1.e6*factor_xy)
	call hbuf_put('QI2',qi2z,1.e6*factor_xy)
	call hbuf_put('QS2',qs2z,1.e6*factor_xy)
	
	call hbuf_put('TKE',tkez,factor_xy)
	
	fadv(1)=0.
	fadv(nz)=0.

!-----------------------------------------------------------------
!  TKE balance:

	shear(1)=0.
	shear(nz)=0.
	do k=2,nzm
          shear(k)=-( (uwle(k)-uwsb(k))*(u0(k)-u0(k-1)) &
            +(vwle(k)-vwsb(k))*(v0(k)-v0(k-1))*YES3D )*factor_xy /(dz*adzw(k))
	end do	
	do k=1,nzm
	  shear(k)=0.5*(shear(k)+shear(k+1)) 
	  tkeleadv(k)=tkeleadv(k)-shear(k)
	  tkelediff(k)=tkelediff(k)-tkelediss(k)
	end do

	call hbuf_put('ADVTR',tkeleadv,1.)
	call hbuf_put('PRESSTR',tkelepress,1.)
	call hbuf_put('BUOYA',tkelebuoy,1.)
	call hbuf_put('SHEAR',shear,1.)
	call hbuf_put('DISSIP',tkelediss,1.)
	call hbuf_put('DIFTR',tkelediff,1.)

!-----------------------------------------------------------------
!  Momentum flux balance:

!  UW advection d(w'w'u')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                         ( u(i,j,k-1)-u0(k-1)+u(i,j,k)-u0(k))
	  end do
	 end do
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 shearx(k)=momleadv(k,1)-coef
	 momleadv(k,1)=coef
	end do	


!  VW advection d(w'w'v')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &	
                         ( v(i,j,k-1)-v0(k-1)+v(i,j,k)-v0(k))
	  end do
	 end do
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 sheary(k)=momleadv(k,2)-coef
	 momleadv(k,2)=coef
	end do	


!  UW advection d(p'u')/dz:

	do k=1,nz
	 fadv(k)=0.
	 if(k.eq.1) then
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+(1.5*(u(i,j,k)-u0(k))*p(i,j,k)*rho(k)- &
                     0.5*(u(i,j,k+1)-u0(k+1))*p(i,j,k+1)*rho(k+1))
	    end do
	   end do
	 else if(k.eq.nz) then
	   do j=1,ny
	    do i=1,nx
	     fadv(k)=fadv(k)+(1.5*(u(i,j,k-1)-u0(k-1))*p(i,j,k-1)*rho(k-1)- &
                  0.5*(u(i,j,k-2)-u0(k-2))*p(i,j,k-2)*rho(k-2))
	    end do
	   end do
	 else
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+0.5*((u(i,j,k)-u0(k))*p(i,j,k)*rho(k)+ &
                       (u(i,j,k-1)-u0(k-1))*p(i,j,k-1)*rho(k-1))
	    end do
	   end do
	 end if
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 presx(k)=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	end do	


!  VW advection d(p'v')/dz:

	do k=1,nz
	 fadv(k)=0.
	 if(k.eq.1) then
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+(1.5*(v(i,j,k)-v0(k))*p(i,j,k)*rho(k)- & 
                    0.5*(v(i,j,k+1)-v0(k+1))*p(i,j,k+1)*rho(k+1))
	    end do
	   end do
	 else if(k.eq.nz) then
	   do j=1,ny
	    do i=1,nx
             fadv(k)=fadv(k)+(1.5*(v(i,j,k-1)-v0(k-1))*p(i,j,k-1)*rho(k-1)- &
                     0.5*(v(i,j,k-2)-v0(k-2))*p(i,j,k-2)*rho(k-2))
	    end do
	   end do
	 else
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+0.5*((v(i,j,k)-v0(k))*p(i,j,k)*rho(k)+ &
                       (v(i,j,k-1)-v0(k-1))*p(i,j,k-1)*rho(k-1))
	    end do
	   end do
	 end if
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 presy(k)=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	end do	

	do k=1,nzm
	  momlepress(k,1)=momlepress(k,1)-presx(k)
	  momlepress(k,2)=momlepress(k,2)-presy(k)
	  momlepress(k,3)=momlepress(k,3)-tkelepress(k)
	end do

	call hbuf_put('WUADV',momleadv(1,1),1.)
	call hbuf_put('WUANIZ',momlepress(1,1),1.)
	call hbuf_put('WUBUOY',momlebuoy(1,1),1.)
	call hbuf_put('WUSHEAR',shearx,1.)
	call hbuf_put('WUPRES',presx,1.)
	call hbuf_put('WUDIFF',momlediff(1,1),1.)

	call hbuf_put('WVADV',momleadv(1,2),1.)
	call hbuf_put('WVANIZ',momlepress(1,2),1.)
	call hbuf_put('WVBUOY',momlebuoy(1,2),1.)
	call hbuf_put('WVSHEAR',sheary,1.)
	call hbuf_put('WVPRES',presy,1.)
	call hbuf_put('WVDIFF',momlediff(1,2),1.)

	call hbuf_put('W2BUOY',momlebuoy(1,3),2.)
	call hbuf_put('W2ADV',momleadv(1,3),2.)
	call hbuf_put('W2REDIS',momlepress(1,3),2.)
	call hbuf_put('W2PRES',tkelepress,2.)
	call hbuf_put('W2DIFF',momlediff(1,3),2.)

!-----------------------------------------------------------
! T2 and Q2 variance budget:


	do k=1,nzm
	  q2lediff(k)=q2lediff(k)-q2lediss(k)
	  t2lediff(k)=t2lediff(k)-t2lediss(k)
	end do	
	

	call hbuf_put('T2ADVTR',t2leadv,1.)
	call hbuf_put('T2GRAD',t2legrad,1.)
	call hbuf_put('T2DISSIP',t2lediss,1.)
	call hbuf_put('T2DIFTR',t2lediff,1.)
	call hbuf_put('T2PREC',t2leprec,1.)

	call hbuf_put('Q2ADVTR',q2leadv,1.)
	call hbuf_put('Q2GRAD',q2legrad,1.)
	call hbuf_put('Q2DISSIP',q2lediss,1.)
	call hbuf_put('Q2DIFTR',q2lediff,1.)
	call hbuf_put('Q2PREC',q2leprec,1.)

!------------------------------------------------------------------
! HW and QW budgets:


	fadv(1)=0.
	fadv(nz)=0.


!  HW advection d(w'w'h')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                 	      ( t(i,j,k-1)-t0(k-1)+t(i,j,k)-t0(k))
	  end do
	 end do
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 twgrad(k)=twleadv(k)-coef
	 twleadv(k)=coef
	end do	


!  QW advection d(w'w'h')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                	      ( q(i,j,k-1)-q0(k-1)+q(i,j,k)-q0(k))
	  end do
	 end do
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 qwgrad(k)=qwleadv(k)-coef
	 qwleadv(k)=coef
	end do	



	call hbuf_put('TWADV',twleadv,factor_xy)
	call hbuf_put('TWDIFF',twlediff,factor_xy)
	call hbuf_put('TWGRAD',twgrad,factor_xy)
	call hbuf_put('TWBUOY',twlebuoy,factor_xy)
	call hbuf_put('TWPRES',twlepres,factor_xy)
	call hbuf_put('TWPREC',twleprec,factor_xy)

	call hbuf_put('QWADV',qwleadv,factor_xy)
	call hbuf_put('QWDIFF',qwlediff,factor_xy)
	call hbuf_put('QWGRAD',qwgrad,factor_xy)
	call hbuf_put('QWBUOY',qwlebuoy,factor_xy)
	call hbuf_put('QWPRES',qwlepres,factor_xy)
	call hbuf_put('QWPREC',qwleprec,factor_xy)


!  SW advection d(w'w's')/dz:

	if(doscalar) then

	  do k=2,nzm
	   fadv(k)=0.
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
          	        ( tke(i,j,k-1)-tke0(k-1)+tke(i,j,k)-tke0(k))
	    end do
	   end do
	  end do
	  do k=1,nzm
	   coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	   swgrad(k)=swleadv(k)-coef
	   swleadv(k)=coef
	  end do	

	  call hbuf_put('SWADV',swleadv,factor_xy)
	  call hbuf_put('SWDIFF',swlediff,factor_xy)
	  call hbuf_put('SWGRAD',swgrad,factor_xy)
	  call hbuf_put('SWBUOY',swlebuoy,factor_xy)
	  call hbuf_put('SWPRES',swlepres,factor_xy)
	  call hbuf_put('S2ADVTR',s2leadv,1.)
	  call hbuf_put('S2GRAD',s2legrad,1.)
	  call hbuf_put('S2DISSIP',s2lediss,1.)
	  call hbuf_put('S2DIFTR',s2lediff,1.)

	  do k=1,nzm
	   tmp(1) = dz/rhow(k)
	   tmp(2) = tmp(1) / dtn
	   tkewsb(k) = tkewsb(k) * tmp(1) * rhow(k)
	   tkewle(k) = tkewle(k)*tmp(2)*rhow(k) + tkewsb(k)
	  end do
	  call hbuf_put('SW',tkewle,factor_xy)
	  call hbuf_put('SWSB',tkewsb,factor_xy)
	  call hbuf_put('S2',s2z,1.e6*factor_xy)
	  call hbuf_avg_put('S',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

	 else

	  call hbuf_put('SWADV',zero,factor_xy)
	  call hbuf_put('SWDIFF',zero,factor_xy)
	  call hbuf_put('SWGRAD',zero,factor_xy)
	  call hbuf_put('SWBUOY',zero,factor_xy)
	  call hbuf_put('SWPRES',zero,factor_xy)
	  call hbuf_put('S2ADVTR',zero,1.)
	  call hbuf_put('S2GRAD',zero,1.)
	  call hbuf_put('S2DISSIP',zero,1.)
	  call hbuf_put('S2DIFTR',zero,1.)
	  call hbuf_put('SW',zero,factor_xy)
	  call hbuf_put('SWSB',zero,factor_xy)
	  call hbuf_put('S',zero,factor_xy)
	  call hbuf_put('S2',zero,factor_xy)

	end if


!---------------------------------------------------------
! SGS TKE Budget:

	if(.not.dosmagor) then
 	 call hbuf_put('ADVTRS',tkewle,factor_xy)
	 call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	else
	 if(doscalar) then
 	   call hbuf_put('TKES',zero,factor_xy)
	 else
	   call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	 end if
 	 call hbuf_put('ADVTRS',zero,factor_xy)
	end if
	 call hbuf_put('BUOYAS',tkesbbuoy,factor_xy)
	 call hbuf_put('SHEARS',tkesbshear,factor_xy)
	 call hbuf_put('DISSIPS',tkesbdiss,factor_xy)
	
	
!-------------------------------------------------------------
!	Cloud statistics:
!-------------------------------------------------------------

	do k=1,nzm
	 cld(k) = 0.
	 hydro(k) = 0.
	 wcl(k) = 0.
	 ucl(k) = 0.
	 vcl(k) = 0.
	 wacl(k) = 0.
	 tcl(k) = 0.
	 tacl(k) = 0.
	 tvcl(k)= 0.
	 qcl(k) = 0.
	 qccl(k)= 0.
	 qicl(k)= 0.
	 qpcl(k)= 0.
	 tvwcl(k)= 0.
	 twcl(k)= 0.
	 swcl(k)= 0.
	 qwcl(k)= 0.
	 qcwcl(k)= 0.
	 qiwcl(k)= 0.
	 scl(k) = 0.
	 prof1(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
	 prof4(k)=0.
	 dse(k)=0.
	 mse(k)=0.
	 sse(k)=0.
	 if(LES) then
	  coef=0.
	 else
	  coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
	 endif
	 do j=1,ny
	  do i=1,nx
	    if(qn(i,j,k).gt.coef) then
	      omn = omegan(tabs(i,j,k))
	      omp = omegap(tabs(i,j,k))
	      lstarn=fac_cond+(1.-omn)*fac_fus
	      lstarp=fac_cond+(1.-omp)*fac_fus
	      cld(k)=cld(k) + 1
	      hydro(k) = hydro(k) + 1
	      tmp(1)=0.5*(w(i,j,k+1)+w(i,j,k))
	      wcl(k) = wcl(k) + tmp(1)
	      ucl(k) = ucl(k) + u(i,j,k) 
	      vcl(k) = vcl(k) + v(i,j,k) 
	      if(tmp(1).gt.0.) then
		prof1(k)=prof1(k)+rho(k)*tmp(1)
	      else
	        prof2(k)=prof2(k)+rho(k)*tmp(1)
	      endif	
	      qcc=qn(i,j,k)*omn
	      qii=qn(i,j,k)*(1.-omn)
	      tmp(1)=t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	      dse(k)=dse(k)+tmp(1)	
	      mse(k)=mse(k)+tmp(1)+fac_cond*(q(i,j,k)-qn(i,j,k))	
	      tcl(k) = tcl(k) + t(i,j,k)
	      qcl(k) = qcl(k) + q(i,j,k)  
	      scl(k) = scl(k) + tke(i,j,k)
	      qccl(k) = qccl(k) + qcc
	      qicl(k) = qicl(k) + qii
	      qpcl(k) = qpcl(k) + qp(i,j,k)
	      tvcl(k) = tvcl(k) + tvirt(i,j,k)	 
	      tacl(k) = tacl(k) + tabs(i,j,k)	 
	      twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      swcl(k) = swcl(k) + tke(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qwcl(k) = qwcl(k) + q(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
	      qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))
	    elseif(qp(i,j,k).gt.1.e-4) then
	      hydro(k) = hydro(k) + 1
	      if(w(i,j,k)+w(i,j,k+1).lt.0.) &
 	         prof3(k)=prof3(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	    endif
	  end do
	 end do
	 if(cld(k).gt.0.) then
	   wacl(k) = wcl(k)
	   wcl(k) = wcl(k)/cld(k)
	   ucl(k) = ucl(k)/cld(k)
	   vcl(k) = vcl(k)/cld(k)
	   dse(k)=dse(k)/cld(k)
	   mse(k)=mse(k)/cld(k)
	   tcl(k) = tcl(k)/cld(k)
	   tvcl(k) = tvcl(k)/cld(k)
	   tacl(k) = tacl(k)/cld(k)
	   qcl(k) = qcl(k)/cld(k)
	   qccl(k) = qccl(k)/cld(k)
	   qicl(k) = qicl(k)/cld(k)
	   qpcl(k) = qpcl(k)/cld(k)
	   scl(k) = scl(k)/cld(k)
	   cloud_factor(k) = cloud_factor(k)+1
	 endif

	 prof4(k)=prof1(k)+prof2(k)+prof3(k)

	end do
		
	call hbuf_put('CLD',cld,factor_xy)
	call hbuf_put('HYDRO',hydro,factor_xy)
	call hbuf_put('WCLD',wcl,1.)
	call hbuf_put('UCLD',ucl,1.)
	call hbuf_put('VCLD',vcl,1.)
	call hbuf_put('DSECLD',dse,1.)
	call hbuf_put('MSECLD',mse,1.)
	call hbuf_put('TLCLD',tcl,1.)
	call hbuf_put('TVCLD',tvcl,1.)
	call hbuf_put('TACLD',tacl,1.)
	call hbuf_put('QTCLD',qcl,1.e3)
	call hbuf_put('QCCLD',qccl,1.e3)
	call hbuf_put('QICLD',qicl,1.e3)
	call hbuf_put('QPCLD',qpcl,1.e3)
	call hbuf_put('WCLDA',wacl,factor_xy)
	call hbuf_put('TLWCLD',twcl,factor_xy)
	call hbuf_put('TVWCLD',tvwcl,factor_xy)
	call hbuf_put('QTWCLD',qwcl,factor_xy*1.e3)
	call hbuf_put('QCWCLD',qcwcl,factor_xy*1.e3)
	call hbuf_put('QIWCLD',qiwcl,factor_xy*1.e3)
	call hbuf_put('MCUP',prof1,factor_xy)
	call hbuf_put('MCDNS',prof2,factor_xy)
	call hbuf_put('MCDNU',prof3,factor_xy)
	call hbuf_put('MC',prof4,factor_xy)

	if(doscalar) then
	 call hbuf_put('SCLD',tcl,1.)
	 call hbuf_put('SWCLD',swcl,factor_xy)
	else
	 call hbuf_put('SCLD',zero,1.)
	 call hbuf_put('SWCLD',zero,factor_xy)
	end if

!-------------------------------------------------------------
!	Updraft Core statistics:
!-------------------------------------------------------------


	do k=1,nzm
	 cld(k) = 0.
	 cldd(k) = 0.
	 wcl(k) = 0.
	 ucl(k) = 0.
	 vcl(k) = 0.
	 wacl(k) = 0.
	 tcl(k) = 0.
	 scl(k) = 0.
	 tvcl(k)= 0.
	 tacl(k)= 0.
	 qcl(k) = 0.
	 qccl(k)= 0.
	 qicl(k)= 0.
	 qpcl(k)= 0.
	 tvwcl(k)= 0.
	 twcl(k)= 0.
	 swcl(k)= 0.
	 qwcl(k)= 0.
	 qcwcl(k)= 0.
	 qiwcl(k)= 0.
	 dse(k)=0.
	 mse(k)=0.
	 prof1(k)=0.
         if(LES) then
          coef=0.
         else
          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
         endif
	 do j=1,ny
	  do i=1,nx
	    condition_cl = qn(i,j,k).gt.coef
     	    condition = tvirt(i,j,k).gt.tvz(k) 
	    if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).gt.2.
	    if(condition) then
	      omn = omegan(tabs(i,j,k))
	      omp = omegap(tabs(i,j,k))
	      cld(k)=cld(k) + 1
	      wcl(k) = wcl(k) + 0.5*(w(i,j,k+1)+w(i,j,k))
	      ucl(k) = ucl(k) + u(i,j,k) 
	      vcl(k) = vcl(k) + v(i,j,k) 
	      lstarn=fac_cond+(1.-omn)*fac_fus
	      lstarp=fac_cond+(1.-omp)*fac_fus
	      tmp(1)=t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	      dse(k)=dse(k)+tmp(1)	
	      mse(k)=mse(k)+tmp(1)+fac_cond*(q(i,j,k)-qn(i,j,k))	
	      tcl(k) = tcl(k) + t(i,j,k)
	      scl(k) = scl(k) + tke(i,j,k)
	      qcc=qn(i,j,k)*omn
	      qii=qn(i,j,k)*(1.-omn)
	      qcl(k) = qcl(k) + q(i,j,k)  
	      qccl(k) = qccl(k) + qcc
	      qicl(k) = qicl(k) + qii  
	      qpcl(k) = qpcl(k) + qp(i,j,k)  
	      tvcl(k) = tvcl(k) + tvirt(i,j,k)
	      tacl(k) = tacl(k) + tabs(i,j,k)
	      twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      swcl(k) = swcl(k) + tke(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qwcl(k) = qwcl(k) + q(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
	      qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))
	      prof1(k)=prof1(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	      if(condition_cl) then
	        cldd(k)=cldd(k)+1
	      end if
	    endif
	  end do
	 end do
	 if(cld(k).gt.0.) then
	   wacl(k) = wcl(k)
	   wcl(k) = wcl(k)/cld(k)
	   ucl(k) = ucl(k)/cld(k)
	   vcl(k) = vcl(k)/cld(k)
	   dse(k)=dse(k)/cld(k)
	   mse(k)=mse(k)/cld(k)
	   tcl(k) = tcl(k)/cld(k)
	   scl(k) = scl(k)/cld(k)
	   qcl(k) = qcl(k)/cld(k)
	   qccl(k) = qccl(k)/cld(k)
	   qicl(k) = qicl(k)/cld(k)
	   qpcl(k) = qpcl(k)/cld(k)
	   tvcl(k) = tvcl(k)/cld(k)
	   tacl(k) = tacl(k)/cld(k)
	   core_factor(k) = core_factor(k)+1
	 endif
	end do
		
	call hbuf_put('CORE',cld,factor_xy)
	call hbuf_put('CORECL',cldd,factor_xy)
	call hbuf_put('WCOR',wcl,1.)
	call hbuf_put('UCOR',ucl,1.)
	call hbuf_put('VCOR',vcl,1.)
	call hbuf_put('DSECOR',dse,1.)
	call hbuf_put('MSECOR',mse,1.)
	call hbuf_put('TLCOR',tcl,1.)
	call hbuf_put('TVCOR',tvcl,1.)
	call hbuf_put('TACOR',tacl,1.)
	call hbuf_put('QTCOR',qcl,1.e3)
	call hbuf_put('QCCOR',qccl,1.e3)
	call hbuf_put('QICOR',qicl,1.e3)
	call hbuf_put('QPCOR',qpcl,1.e3)
	call hbuf_put('WCORA',wacl,factor_xy)
	call hbuf_put('TLWCOR',twcl,factor_xy)
	call hbuf_put('TVWCOR',tvwcl,factor_xy)
	call hbuf_put('QTWCOR',qwcl,factor_xy*1.e3)
	call hbuf_put('QCWCOR',qcwcl,factor_xy*1.e3)
	call hbuf_put('QIWCOR',qiwcl,factor_xy*1.e3)

	if(doscalar) then
	 call hbuf_put('SCOR',tcl,1.)
	 call hbuf_put('SWCOR',swcl,factor_xy)
	else
	 call hbuf_put('SCOR',zero,1.)
	 call hbuf_put('SWCOR',zero,factor_xy)
	end if
!-------------------------------------------------------------
!	Cloud Downdraft Core statistics:
!-------------------------------------------------------------


	do k=1,nzm
	 cld(k) = 0.
	 cldd(k) = 0.
	 wcl(k) = 0.
	 ucl(k) = 0.
	 vcl(k) = 0.
	 wacl(k) = 0.
	 tcl(k) = 0.
	 scl(k) = 0.
	 tvcl(k)= 0.
	 tacl(k)= 0.
	 qcl(k) = 0.
	 qccl(k)= 0.
	 qicl(k)= 0.
	 qpcl(k)= 0.
	 tvwcl(k)= 0.
	 twcl(k)= 0.
	 swcl(k)= 0.
	 qwcl(k)= 0.
	 qcwcl(k)= 0.
	 qiwcl(k)= 0.
	 dse(k)=0.
	 mse(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
         if(LES) then
          coef=0.
         else
          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
         endif
	 do j=1,ny
	  do i=1,nx
	    condition_cl = qn(i,j,k).gt.coef .or. qp(i,j,k).gt.1.e-4 
     	    condition = tvirt(i,j,k).lt.tvz(k) 
	    if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).lt.-2.
	    if(condition) then
	      omn = omegan(tabs(i,j,k))
	      omp = omegap(tabs(i,j,k))
	      cld(k)=cld(k) + 1
	      wcl(k) = wcl(k) + 0.5*(w(i,j,k+1)+w(i,j,k))
	      ucl(k) = ucl(k) + u(i,j,k) 
	      vcl(k) = vcl(k) + v(i,j,k) 
	      lstarn=fac_cond+(1.-omn)*fac_fus
	      lstarp=fac_cond+(1.-omp)*fac_fus
	      tmp(1)=t(i,j,k)+lstarn*qn(i,j,k)+lstarp*qp(i,j,k)
	      dse(k)=dse(k)+tmp(1)	
	      mse(k)=mse(k)+tmp(1)+fac_cond*(q(i,j,k)-qn(i,j,k))	
	      tcl(k) = tcl(k) + t(i,j,k)
	      scl(k) = scl(k) + tke(i,j,k)
	      qcc=qn(i,j,k)*omn
	      qii=qn(i,j,k)*(1.-omn)
	      qcl(k) = qcl(k) + q(i,j,k)  
	      qccl(k) = qccl(k) + qcc
	      qicl(k) = qicl(k) + qii  
	      qpcl(k) = qpcl(k) + qp(i,j,k)  
	      tvcl(k) = tvcl(k) + tvirt(i,j,k)
	      tacl(k) = tacl(k) + tabs(i,j,k)
	      twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      swcl(k) = swcl(k) + tke(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qwcl(k) = qwcl(k) + q(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
	      qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
	      qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))
	      if(condition_cl) then
	        prof2(k)=prof2(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	        cldd(k)=cldd(k) + 1
	      else
	        prof3(k)=prof3(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	      end if
	    endif
	  end do
	 end do
	 if(cld(k).gt.0.) then
	   wacl(k) = wcl(k)
	   wcl(k) = wcl(k)/cld(k)
	   ucl(k) = ucl(k)/cld(k)
	   vcl(k) = vcl(k)/cld(k)
	   dse(k)=dse(k)/cld(k)
	   mse(k)=mse(k)/cld(k)
	   tcl(k) = tcl(k)/cld(k)
	   scl(k) = scl(k)/cld(k)
	   qcl(k) = qcl(k)/cld(k)
	   qccl(k) = qccl(k)/cld(k)
	   qicl(k) = qicl(k)/cld(k)
	   qpcl(k) = qpcl(k)/cld(k)
	   tvcl(k) = tvcl(k)/cld(k)
	   tacl(k) = tacl(k)/cld(k)
	   coredn_factor(k) = coredn_factor(k)+1
	 endif

	 prof4(k)=prof1(k)+prof2(k)+prof3(k)

	end do
		
	call hbuf_put('COREDN',cld,factor_xy)
	call hbuf_put('COREDNCL',cldd,factor_xy)
	call hbuf_put('UCORDN',ucl,1.)
	call hbuf_put('VCORDN',vcl,1.)
	call hbuf_put('WCORDN',wcl,1.)
	call hbuf_put('DSECORDN',dse,1.)
	call hbuf_put('MSECORDN',mse,1.)
	call hbuf_put('TLCORDN',tcl,1.)
	call hbuf_put('TVCORDN',tvcl,1.)
	call hbuf_put('TACORDN',tacl,1.)
	call hbuf_put('QTCORDN',qcl,1.e3)
	call hbuf_put('QCCORDN',qccl,1.e3)
	call hbuf_put('QICORDN',qicl,1.e3)
	call hbuf_put('QPCORDN',qpcl,1.e3)
	call hbuf_put('WCORDNA',wacl,factor_xy)
	call hbuf_put('TLWCORDN',twcl,factor_xy)
	call hbuf_put('TVWCORDN',tvwcl,factor_xy)
	call hbuf_put('QTWCORDN',qwcl,factor_xy*1.e3)
	call hbuf_put('QCWCORDN',qcwcl,factor_xy*1.e3)
	call hbuf_put('QIWCORDN',qiwcl,factor_xy*1.e3)
	call hbuf_put('MCRUP',prof1,factor_xy)
	call hbuf_put('MCRDNS',prof2,factor_xy)
	call hbuf_put('MCRDNU',prof3,factor_xy)
	call hbuf_put('MCR',prof4,factor_xy)

	if(doscalar) then
	 call hbuf_put('SCORDN',tcl,1.)
	 call hbuf_put('SWCORDN',swcl,factor_xy)
	else
	 call hbuf_put('SCORDN',zero,1.)
	 call hbuf_put('SWCORDN',zero,factor_xy)
	end if
!---------------------------------------------------------
!  Radiation and other stuff

	do j=1,ny
	 do i=1,nx
	   cwp(i,j)=0.
	   cwpl(i,j)=0.
	   cwpm(i,j)=0.
	   cwph(i,j)=0.
	   topind(i,j)=1
	 end do
	end do
	
	if(CEM) then
	  cwpmax=0.02
	else
	  cwpmax=0.0
	endif

	do k=nzm,1,-1
	 prof1(k)=(radqrlw(k)+radqrsw(k))*factor_xy
	 tmp(1)=rho(k)*adzw(k)*dz
	 do j=1,ny
	  do i=1,nx
	    cwp(i,j)=cwp(i,j)+tmp(1)*qn(i,j,k)
            if(pres(k).ge.700.) then
	      cwpl(i,j)=cwpl(i,j)+tmp(1)*qn(i,j,k)
            else if(pres(k).le.400.) then
	      cwph(i,j)=cwph(i,j)+tmp(1)*qn(i,j,k)
            else
	      cwpm(i,j)=cwpm(i,j)+tmp(1)*qn(i,j,k)
	    end if
	    if(cwp(i,j).gt.cwpmax.and.topind(i,j).eq.1)topind(i,j)=k
	  end do
	 end do
	end do

	do j=1,ny
	 do i=1,nx
	   if(cwp(i,j).gt.cwpmax) s_acld=s_acld+1.
	   if(cwpl(i,j).gt.cwpmax) s_acldl=s_acldl+1.
	   if(cwpm(i,j).gt.cwpmax) s_acldm=s_acldm+1.
	   if(cwph(i,j).gt.cwpmax) s_acldh=s_acldh+1.
	   if(tabs(i,j,topind(i,j)).lt.245.) s_acldcold=s_acldcold+1
	 end do
	end do

	do k=1,nzm
	 prof2(k)=0.
	 prof3(k)=0.	 
	 n=0
	 if(dolongwave.or.doshortwave) then
	   do j=1,ny
	     do i=1,nx
	       if(cwp(i,j).gt.cwpmax) then
	         n=n+1
	         prof2(k)=prof2(k)+misc(i,j,k)
	       else
	         prof3(k)=prof3(k)+misc(i,j,k)
	       endif
	     end do
	   end do 
	 end if
	end do


	call hbuf_put('HLADV',tadv,factor_xy*86400./dtn)
	call hbuf_put('HLDIFF',tdiff,factor_xy*86400./dtn)
	call hbuf_put('HLLAT',tlat+tlatqi,factor_xy*86400./dtn)
	call hbuf_put('HLRAD',prof1,86400.)
	call hbuf_put('QTADV',qadv+qifall,factor_xy*86400000./dtn)
	call hbuf_put('QTDIFF',qdiff,factor_xy*86400000./dtn)
	call hbuf_put('QTSINK',qpsrc,-factor_xy*86400000./dtn)
	call hbuf_put('QTSRC',qpevp,-factor_xy*86400000./dtn)
	call hbuf_put('QPADV',qpadv,factor_xy*86400000./dtn)
	call hbuf_put('QPDIFF',qpdiff,factor_xy*86400000./dtn)
	call hbuf_put('QPFALL',qpfall,factor_xy*86400000./dtn)
	call hbuf_put('QPSRC',qpsrc,factor_xy*86400000./dtn)
	call hbuf_put('QPEVP',qpevp,factor_xy*86400000./dtn)


        if(dotracers) then
          do ntr=1,ntracers
           call hbuf_put(trim(tracername(ntr))//'FLX',trwle(:,ntr),factor_xy)
           call hbuf_put(trim(tracername(ntr))//'FLXS',trwsb(:,ntr),factor_xy)
           call hbuf_put(trim(tracername(ntr))//'ADV',tradv(:,ntr),factor_xy*86400./dtn)
           call hbuf_put(trim(tracername(ntr))//'DIFF',trdiff(:,ntr),factor_xy*86400./dtn)
           call hbuf_put(trim(tracername(ntr))//'PHYS',trphys(:,ntr),factor_xy*86400./dtn)
          end do
        end if


	call hbuf_put('RADLWUP',radlwup,factor_xy)
	call hbuf_put('RADLWDN',radlwdn,factor_xy)
	call hbuf_put('RADSWUP',radswup,factor_xy)
	call hbuf_put('RADSWDN',radswdn,factor_xy)
	call hbuf_put('RADQRLW',radqrlw,factor_xy*86400.) 	
	call hbuf_put('RADQRSW',radqrsw,factor_xy*86400.) 	
	call hbuf_put('RADQR',prof1,86400.) 	
	call hbuf_put('RADQRC',prof2,86400./(n+1.e-5)) 	
	call hbuf_put('RADQRS',prof3,86400./(nx*ny-n+1.e-5))

!---------------------------------------------------------
!  Apparent heat/moisture sources/sinks

	tmp(1)=1./dtn

	do k=1,nzm
	 prof2(k)=0.
	 prof3(k)=0.	 
	 n=0
	 do j=1,ny
	  do i=1,nx
	    prof2(k)=prof2(k)+(tabs(i,j,k)-t01(k))*tmp(1)-ttend(k)-prof1(k)
	    prof3(k)=prof3(k)-fac_cond*((q(i,j,k)-q01(k))*tmp(1)-qtend(k))
	  end do
	 end do 
	end do

	call hbuf_put('Q1C',prof2,factor_xy*86400.) 	
	call hbuf_put('Q2',prof3,factor_xy*86400.) 	

        call isccp_get()

end
	
