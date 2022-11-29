module simple_ocean

!------------------------------------------------------------
! Purpose:
!
! A collection of routines used to specify fixed 
! or compute interactive SSTs, like slab-ocean model, etc.
!
! Author: Marat Khairoutdinov
! Based on dynamic ocean impelemntation from the UW version of SAM.
!------------------------------------------------------------
use grid, only: nx,ny,nx_gl,dx,nsubdomains_x,rank,masterproc,dtn

implicit none

private 

! parameters of the sinusoidal SST destribution 
! along the X for Walker-type simulatons( ocean-type = 1):
real :: sst0 = 300.  ! maximum temperature (at the center of the domain)
real :: delta_sst0 = -1.5 ! amplitude of SST change

! Public methods:

public set_sst     ! set SST 
public sst_evolve ! evolve SST according to a model set by the ocean_type

CONTAINS


SUBROUTINE set_sst(ocean_type)

 use vars, only: sstxy
 use params, only: tabs_s

 integer, intent(in) :: ocean_type  

 real tmpx(nx), pii, lx
 integer i,j

 select case (ocean_type)

   case(0) ! fixed constant SST

      sstxy = tabs_s

   case(1) ! Sinusoidal distribution along the x-direction:

     lx = float(nx_gl)*dx
     do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
     end do
     pii = atan2(0.d0,-1.d0)
     do j=1,ny
       do i=1,nx
         sstxy(i,j) = sst0+delta_sst0*cos(2.*pii*tmpx(i)/lx)
       end do
     end do
     tabs_s = sst0

   case default

     if(masterproc) then
         print*, 'unknown ocean type in set_sst. Exitting...'
         call task_abort
     end if

 end select

end subroutine set_sst



SUBROUTINE sst_evolve
 use vars, only: sstxy, fluxbt, fluxbq, rhow,qocean_xy
 use params, only: rhor, cp, lcond
 use rad, only: swnsxy, lwnsxy

 real, parameter :: cpr = 4000.   ! Liquid Water Cp = 4000 J/kg/K
 real, parameter :: hml = 2.   ! Depth of the mixed layer (m)
 real, parameter :: Szero = -25.
 real, parameter :: deltaS = -50.
 real factor_cp, factor_lc, qoceanxy
 real tmpx(nx), lx
 integer i,j

      lx = float(nx_gl)*dx
      do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
      end do

      ! Define weight factors for the mixed layer heating due to
      ! the model's sensible and latent heat flux.
      factor_cp = rhow(1)*cp
      factor_lc = rhow(1)*lcond

      ! Use forward Euler to integrate the differential equation
      ! for the ocean mixed layer temperature: dT/dt = S - E.
      ! The source: CPT?GCSS WG4 idealized Walker-circulation 
      ! RCE Intercomparison proposed by C. Bretherton.
      do j=1,ny
         do i=1,nx
           qoceanxy = Szero + deltaS*abs(2.*tmpx(i)/lx - 1)
	   qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy

            sstxy(i,j) = sstxy(i,j) &
                 + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
                 - lwnsxy(i,j)               & ! LW Radiative Heating
                 - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                 - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                 + qoceanxy)            & ! Ocean Heating
                 /(rhor*cpr*hml)        ! Convert W/m^2 Heating to K/s
         end do
      end do

end subroutine sst_evolve


end module simple_ocean
