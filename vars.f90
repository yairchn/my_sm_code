module vars

use grid

implicit none
!--------------------------------------------------------------------
! prognostic variables:

real u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
real v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
real w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
real t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! moist static energy
real q   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! total water
real qp	 (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! rain+snow
real tke (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! SGS TKE

!--------------------------------------------------------------------
! diagnostic variables:

real p      (0:nx, (1-YES3D):ny, nzm)     ! pressure
real tabs   (nx, ny, nzm)                 ! temperature
real qn     (nx, ny, nzm)                 ! cloud water+cloud ice 
real tk     (0:nxp1, (1-YES3D):nyp1, nzm) ! SGS eddyviscosity
real tkh    (0:nxp1, (1-YES3D):nyp1, nzm) ! SGS eddy conductivity
        
!--------------------------------------------------------------------
! time-tendencies for prognostic variables

real dudt   (nxp1, ny, nzm, 3)
real dvdt   (nx, nyp1, nzm, 3)
real dwdt   (nx, ny, nz,  3)

!----------------------------------------------------------------
! Temporary storage array:

	real misc(nx, ny, nz)
!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:

real fluxbu (nx, ny), fluxbv (nx, ny), fluxbt (nx, ny)
real fluxbq (nx, ny), fluxtu (nx, ny), fluxtv (nx, ny)
real fluxtt (nx, ny), fluxtq (nx, ny), fzero  (nx, ny)
real precsfc(nx,ny) ! surface precip. rate
                
!-----------------------------------------------------------------
! profiles 

real   t0(nzm), q0(nzm), qc0(nzm), qv0(nzm),tabs0(nzm), &
       qi0(nzm),tv0(nzm), rel0(nzm), u0(nzm), v0(nzm), &
       tg0(nzm), qg0(nzm), ug0(nzm), vg0(nzm), p0(nzm), &
       tke0(nzm), t01(nzm), q01(nzm)
!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics 

real  sstobs, lhobs, shobs
!----------------------------------------------------------------
!  Domain top stuff:

real   gamt0    ! gradient of t() at the top,K/m
real   gamq0    ! gradient of q() at the top,g/g/m

!-----------------------------------------------------------------
! reference vertical profiles:
 
real   prespot(nzm)  ! (1000./pres)**R/cp
real   rho(nzm)	  ! air density at pressure levels,kg/m3 
real   rhow(nz)   ! air density at vertical velocity levels,kg/m3
real   bet(nzm)	  ! = ggr/tv0
real   betp(nzm)  ! = ggr/p0
real   gamaz(nzm) ! ggr/cp*z
real   wsub(nz)   ! Large-scale subsidence velocity,m/s
real   qtend(nzm) ! Large-scale tendency for total water
real   ttend(nzm) ! Large-scale tendency for temp.

!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nmaxlsf ! dimension of lsf forcing arrays (time)
integer nmaxrfc ! dimension of rad forcing arrays (time) 
integer nmaxsfc ! dimension of sfc forcing arrays (time) 
integer nmaxsnd ! dimension of observed sounding arrays (time) 
parameter (nmaxrfc=1000,nmaxlsf=1000,nmaxsfc=1000, nmaxsnd=1000)
integer nlsf	! number of large-scale forcing profiles
integer nrfc	! number of radiative forcing profiles
integer nsfc	! number of surface forcing profiles
integer nsnd	! number of observed soundings

real   dqls(nzm,nmaxlsf) ! Large-scale tendency for total water
real   dtls(nzm,nmaxlsf) ! Large-scale tendency for temp.
real   ugls(nzm,nmaxlsf) ! Large-scale wind in X-direction
real   vgls(nzm,nmaxlsf) ! Large-scale wind in Y-direction
real   wgls(nzm,nmaxlsf) ! Large-scale subsidence velocity,m/s
real   pres0ls(nmaxlsf) ! Surface pressure, mb
real   dayls(nmaxlsf)    ! Large-scale forcing arrays time (days) 
real   dtrfc(nzm,nmaxrfc)! Radiative tendency for pot. temp.
real   dayrfc(nmaxrfc)   ! Radiative forcing arrays time (days) 
real   sstsfc(nmaxsfc)   ! SSTs
real   hsfc(nmaxsfc) 	 ! Sensible heat flux,W/m2
real   lesfc(nmaxsfc) 	 ! Latent heat flux,W/m2
real   tausfc(nmaxsfc) 	 ! Surface drag,m2/s2
real   daysfc(nmaxsfc)   ! Surface forcing arrays time (days) 
real   usnd(nzm,nmaxsnd) ! Observed zonal wind
real   vsnd(nzm,nmaxsnd) ! Observed meriod wind
real   tsnd(nzm,nmaxsnd) ! Observed Abs. temperature
real   qsnd(nzm,nmaxsnd) ! Observed Moisture
real   daysnd(nmaxsnd)
 
!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!
real sstxy(nx,ny)	!  surface temperature xy-distribution
real fcory(ny)      !  Coriolis parameter xy-distribution
real fcorzy(ny)      !  z-Coriolis parameter xy-distribution
real latitude(nx,ny)	     ! latitude (degrees)
real longitude(nx,ny)	     ! longitude(degrees)
real prec_xy(nx,ny) ! mean precip. rate for outout
real shf_xy(nx,ny) ! mean precip. rate for outout
real lhf_xy(nx,ny) ! mean precip. rate for outout
real lwns_xy(nx,ny) ! mean net lw at SFC
real swns_xy(nx,ny) ! mean net sw at SFC
real lwnsc_xy(nx,ny) ! clear-sky mean net lw at SFC
real swnsc_xy(nx,ny) ! clear-sky mean net sw at SFC
real lwnt_xy(nx,ny) ! mean net lw at TOA
real swnt_xy(nx,ny) ! mean net sw at TOA
real lwntc_xy(nx,ny) ! clear-sky mean net lw at TOA
real swntc_xy(nx,ny) ! clear-sky mean net sw at TOA
real solin_xy(nx,ny) ! solar TOA insolation
real pw_xy(nx,ny)   ! precipitable water
real cw_xy(nx,ny)   ! cloud water path
real iw_xy(nx,ny)   ! ice water path
real u200_xy(nx,ny) ! u-wind at 200 mb
real usfc_xy(nx,ny) ! u-wind at at the surface
real v200_xy(nx,ny) ! v-wind at 200 mb
real vsfc_xy(nx,ny) ! v-wind at the surface
real w500_xy(nx,ny) ! w at 500 mb
real qocean_xy(nx,ny) ! ocean cooling in W/m2

!----------------------------------------------------------------------
!  Tendencies due to convective parameterization:
!

real   qtend_cup(nx,ny,nzm) ! CUP tendency for total water
real   ttend_cup(nx,ny,nzm) ! CUP tendency for temp.
real   utend_cup(nx,ny,nzm) ! CUP tendency for u
real   vtend_cup(nx,ny,nzm) ! CUP tendency for v

!----------------------------------------------------------------------
!	Vertical profiles of quantities sampled for statitistics purposes:

real &
    twle(nz), twsb(nz), qwle(nz), qwsb(nz), tkewle(nz), &
    tkewsb(nz), qpwle(nz), qpwsb(nz), precflux(nz), &
    uwle(nz), uwsb(nz), vwle(nz), vwsb(nz), &
    cloud_factor(nz), core_factor(nz), coredn_factor(nz), &
    tkeleadv(nz), tkelepress(nz), tkelediss(nz), tkelediff(nz), &
    tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz), &
    tkelebuoy(nz), radlwup(nz), radlwdn(nz), radswup(nz), radswdn(nz), &
    radqrlw(nz), radqrsw(nz), w_max, s_acld, s_acldcold, s_ar, p_conv, p_strat,&
    s_acldl, s_acldm, s_acldh, s_acldisccp, &
    s_acldlisccp, s_acldmisccp, s_acldhisccp, &
    s_flns,s_flnt,s_flnsc,s_flntc,s_flds,s_fsns, &
    s_fsnt,s_fsnsc,s_fsntc,s_fsds,s_solin, & 
    t2leadv(nz),t2legrad(nz),t2lediff(nz),t2leprec(nz),t2lediss(nz), &
    q2leadv(nz),q2legrad(nz),q2lediff(nz),q2leprec(nz),q2lediss(nz), &
    s2leadv(nz),s2legrad(nz),s2lediff(nz),s2lediss(nz), &
    twleadv(nz),twlediff(nz),twlepres(nz),twlebuoy(nz),twleprec(nz), &
    qwleadv(nz),qwlediff(nz),qwlepres(nz),qwlebuoy(nz),qwleprec(nz), &
    swleadv(nz),swlediff(nz),swlepres(nz),swlebuoy(nz), &
    momleadv(nz,3),momlepress(nz,3),momlebuoy(nz,3), &
    momlediff(nz,3),tadv(nz),tdiff(nz),tlat(nz), tlatqi(nz), qifall(nz), &
    qadv(nz),qdiff(nz),qpadv(nz),qpdiff(nz),qpsrc(nz),qpfall(nz),qpevp(nz)


! register functions:


        real esatw,esati,dtesatw,dtesati
        real qsatw,qsati,dtqsatw,dtqsati
        external esatw,esati,dtesatw,dtesati,qsatw,qsati,dtqsatw,dtqsati

	real omegan, omegap, omegag
	external omegan, omegap, omegag


        integer lenstr
        external lenstr

! energy conservation diagnostics:
 
  double precision total_water_before, total_water_after
  double precision total_water_evap, total_water_prec, total_water_ls

end module vars