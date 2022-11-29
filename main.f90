program crm

!       Main module.

use vars
use hbuffer
use tracers
implicit none

integer k, icyc, nn, nstatsteps
real dummy(nzm)
!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
   call tracers_init() ! initialize tracers
   call setforcing()
elseif(nrestart.eq.1) then
   call read_all()
elseif(nrestart.eq.2) then
   call read_all()
   call setparm() ! overwrite the parameters
   call setforcing()
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif
if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()
	
!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + time/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

!----------------------------------------------------------
!       Nadging:

     call nudging()

!----------------------------------------------------------
!   	suppress turbulence near the upper boundary (spange):

     if(dodamping) call damping()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!----------------------------------------------------------
!      Update the subdomain's boundaries for velocity

     call boundaries(0)

!---------------------------------------------------------
!	SGS TKE equation:     	
	   
     if(dosgs) call tke_full()

!---------------------------------------------------------
!        Update boundaries for the SGS exchange coefficients:

     call boundaries(4)

!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!-----------------------------------------------
!   	surface fluxes:

     if(dosurface) then

       call surface()

     end if
!----------------------------------------------------------
!	SGS diffusion of momentum:

     if(dosgs) call diffuse_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()


!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for velocity fields to use for advection of scalars:

     call boundaries(1)

!        Update boundaries for scalars:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
     

     call advect_scalar(q,qadv,qwle,q2leadv,q2legrad,qwleadv,.true.)

     if(dosgs.and..not.dosmagor) then
      call advect_scalar(tke,dummy,tkewle,dummy,dummy,dummy,.false.)
     else if(doscalar) then
      call advect_scalar(tke,dummy,tkewle,s2leadv,s2legrad,swleadv,.true.)
     end if

     if(docloud.and.doprecip) then

       call advect_scalar(qp,qpadv,qpwle,dummy,dummy,dummy,.false.)

       do k=1,nzm
         total_water_prec = total_water_prec + &
                        sum(qp(1:nx,1:ny,k))*adz(k)*dz *rho(k)
       end do
 
       call precip_fall()

       do k=1,nzm
         total_water_prec = total_water_prec - &
                        sum(qp(1:nx,1:ny,k))*adz(k)*dz *rho(k)
       end do

     endif	

 ! advection of tracers:

     if(dotracers) then

        do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
        end do

     end if

!---------------------------------------------------------
!      diffusion of scalars :

!        Update boundaries for scalars:

      if(dosgs) call boundaries(3)

      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.)
     
      do k=1,nzm
         total_water_evap = total_water_evap - &
                         sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
      end do

      call diffuse_scalar(q,fluxbq,fluxtq,qdiff,qwsb, &
                          q2lediff,q2lediss,qwlediff,.true.)

      do k=1,nzm
         total_water_evap = total_water_evap + &
                         sum(q(1:nx,1:ny,k))*adz(k)*dz *rho(k)
      end do

      if(.not.dosmagor) then
          call diffuse_scalar(tke,fzero,fzero,dummy,tkewsb, &
                                    dummy,dummy,dummy,.false.)
      else if(doscalar) then
          call diffuse_scalar(tke,fluxbq,fluxtq,dummy,tkewsb, &
                           s2lediff,s2lediss,swlediff,.true.)
      end if

      if(docloud.and.doprecip) then
           call diffuse_scalar(qp,fzero,fzero,qpdiff,qpwsb, &
                           dummy,dummy,dummy,.false.)

      endif

 ! diffusion of tracers:

      if(dotracers) then

        call tracers_flux()

        do k = 1,ntracers

          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
                           dummy,dummy,dummy,.false.)
 
        end do

      end if

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:


      if(docloud) then
        call cloud()
        call ice_fall()
        if(doprecip) call precip_proc()
      end if

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) then 
	call radiation()     
      end if

!-----------------------------------------------------------
!       Compute field diagnostics and update the velocity field:

      call diagnose()

!----------------------------------------------------------
! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle	
          
!----------------------------------------------------------

   if(mod(nstep,nstatis).eq.0) call statistics()
	
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)
  
!----------------------------------------------------------
!----------------------------------------------------------

end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------
 	
call task_abort()

end program crm
