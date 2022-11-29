     
subroutine write_fields2D
	
use vars
use params
implicit none
character *80 filename,long_name
character *8 name
character *10 timechar
character *4 rankchar
character *6 filetype
character *10 units
character *10 c_dx, c_dy, c_time
integer i,j,nfields,nfields1
real tmp(nx,ny,nzm)
character*3 filestatus

nfields= 22
if(.not.doxy) nfields=nfields-1
if(.not.dolongwave) nfields = nfields-4
if(.not.doshortwave) nfields = nfields-5
if(.not.dodynamicocean) nfields=nfields-1
if(.not.doprecip) nfields = nfields-1
if(.not.docloud) nfields = nfields-2
if(SFC_FLX_FXD) nfields = nfields-2
nfields1=0

if(masterproc) then

  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do


! Make sure that the new run doesn't overwrite the file from the old run 

    filestatus='old'
    if(notopened2D.and.(nrestart.eq.0.or.nrestart.eq.3)) then
      filestatus='new'
    end if

    if(save2Dbin) then
      filetype = '.2Dbin'
    else
      filetype = '.2Dcom'
    end if
    if(save2Dsep) then
       filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
          rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype 
          open(46,file=filename,status='unknown',form='unformatted')
    else
       filename='./DATA3D/'//trim(case)//'_'//trim(caseid)//'_'// &
           rankchar(5-lenstr(rankchar):4)//filetype
	      
       if(nrestart.eq.0.and.notopened2D) then
          open(46,file=filename,status=filestatus,form='unformatted')	
       else
          open(46,file=filename,status=filestatus,form='unformatted', &
                                                     position='append')
       end if
       notopened2D=.false.
    end if

     if(save2Dbin) then

       write(46) nx,ny,nzm,nsubdomains, nsubdomains_x,nsubdomains_y,nfields
       write(46) dx
       write(46) dy
       write(46) nstep*dt/(3600.*24.)+day0

     else

       write(long_name,'(8i4)') nx,ny,nzm,nsubdomains,  &
                                        nsubdomains_x,nsubdomains_y,nfields
       write(c_dx,'(f10.0)') dx
       write(c_dy,'(f10.0)') dy
       write(c_time,'(f10.2)') nstep*dt/(3600.*24.)+day0
       write(46) long_name(1:32)
       write(46) c_time,c_dx,c_dy

     end if ! save2Dbin

end if! masterproc


! 2D fields:


if(doprecip) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=prec_xy(i,j)*dz/dt*86400./float(nsave2D)
      prec_xy(i,j) = 0.
    end do
   end do
  name='Prec'
  long_name='Surface Precip. Rate'
  units='mm/day'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(.not.SFC_FLX_FXD) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=shf_xy(i,j)*rhow(1)*cp/float(nsave2D)
      shf_xy(i,j) = 0.
    end do
   end do
  name='SHF'
  long_name='Sensible Heat Flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=lhf_xy(i,j)*rhow(1)*lcond/float(nsave2D)
      lhf_xy(i,j) = 0.
    end do
   end do
  name='LHF'
  long_name='Latent Heat Flux'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
        
end if

if(dolongwave) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwns_xy(i,j)/float(nsave2D)
       lwns_xy(i,j) = 0.
     end do
   end do
  name='LWNS'
  long_name='Net LW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwnsc_xy(i,j)/float(nsave2D)
       lwnsc_xy(i,j) = 0.
     end do
   end do
  name='LWNSC'
  long_name='Net clear-sky LW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwnt_xy(i,j)/float(nsave2D)
       lwnt_xy(i,j) = 0.
     end do
   end do
  name='LWNT'
  long_name='Net LW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=lwntc_xy(i,j)/float(nsave2D)
       lwntc_xy(i,j) = 0.
     end do
   end do
  name='LWNTC'
  long_name='Clear-Sky Net LW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(doshortwave) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=solin_xy(i,j)/float(nsave2D)
       solin_xy(i,j) = 0.
     end do
   end do
  name='SOLIN'
  long_name='Solar TOA insolation'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swns_xy(i,j)/float(nsave2D)
       swns_xy(i,j) = 0.
     end do
   end do
  name='SWNS'
  long_name='Net SW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swnsc_xy(i,j)/float(nsave2D)
       swnsc_xy(i,j) = 0.
     end do
   end do
  name='SWNSC'
  long_name='Net Clear-sky SW at the surface'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swnt_xy(i,j)/float(nsave2D)
       swnt_xy(i,j) = 0.
     end do
   end do
  name='SWNT'
  long_name='Net SW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=swntc_xy(i,j)/float(nsave2D)
       swntc_xy(i,j) = 0.
     end do
   end do
  name='SWNTC'
  long_name='Net Clear-Sky SW at TOA'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(docloud) then
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=cw_xy(i,j)/float(nsave2D)
       cw_xy(i,j) = 0.
     end do
   end do
  name='CWP'
  long_name='Cloud Water Path'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=iw_xy(i,j)/float(nsave2D)
       iw_xy(i,j) = 0.
     end do
   end do
  name='IWP'
  long_name='Ice Path'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=pw_xy(i,j)/float(nsave2D)
     pw_xy(i,j) = 0.
     end do
   end do
  name='PW'
  long_name='Precipitable Water'
  units='mm'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=usfc_xy(i,j)/float(nsave2D)
       usfc_xy(i,j) = 0.
     end do
   end do
  name='USFC'
  long_name='U at the surface'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=u200_xy(i,j)/float(nsave2D)
       u200_xy(i,j) = 0.
     end do
   end do
  name='U200'
  long_name='U at 200 mb'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
      tmp(i,j,1)=vsfc_xy(i,j)/float(nsave2D)
      vsfc_xy(i,j) = 0.
     end do
   end do
  name='VSFC'
  long_name='V at the surface'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=v200_xy(i,j)/float(nsave2D)
       v200_xy(i,j) = 0.
     end do
   end do
  name='V200'
  long_name='V at 200 mb'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j,1)=w500_xy(i,j)/float(nsave2D)
       w500_xy(i,j) = 0.
     end do
   end do
  name='W500'
  long_name='W at 500 mb'
  units='m/s'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

if(dodynamicocean) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=qocean_xy(i,j)/float(nsave2D)
      qocean_xy(i,j) = 0.
    end do
   end do
  name='QOCN'
  long_name='Deep Ocean Cooling'
  units='W/m2'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(doxy) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j,1)=sstxy(i,j)
    end do
   end do
  name='SST'
  long_name='Sea Surface Temperature'
  units='K'
  call compress3D(tmp,nx,ny,1,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

call task_barrier()


if(nfields.ne.nfields1) then
  if(masterproc) print*,'write_fields2D: error in nfields!!'
  call task_abort()
end if

if(masterproc) then
     close(46)
     if(save2Dsep.and.dogzip2D) call systemf('gzip -f '//filename)
     print*, 'Appending 2D data. file:'//filename
endif


end
