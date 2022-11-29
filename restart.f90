	subroutine write_all()
	
	use vars
	use params
	implicit none
	character *4 rankchar
	character *256 filename
	integer irank


	if(restart_sep) then

          write(rankchar,'(i4)') rank

          filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'


          open(65,file=trim(filename), status='unknown',form='unformatted')
          write(65) nsubdomains

	  call write_statement


	else
	  write(rankchar,'(i4)') nsubdomains
	  filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       if(masterproc) then
	      
	        open(65,file=trim(filename), status='unknown',form='unformatted')
	        write(65) nsubdomains

	       else

                open(65,file=trim(filename), status='unknown',form='unformatted',&
                   position='append')

	       end if

               call write_statement

             end if
	  end do

	end if ! restart_sep

	call task_barrier()

        return
        end
 
 
 
 
     
	subroutine read_all()
	
	use vars
	use params
	implicit none
	character *4 rankchar
	character *256 filename
	integer irank, ii


	if(restart_sep) then

           write(rankchar,'(i4)') rank

           filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'


           open(65,file=trim(filename), status='unknown',form='unformatted')
           read(65)

	   call read_statement


	else

	  write(rankchar,'(i4)') nsubdomains

	  filename='./RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          open(65,file=trim(filename), status='unknown',form='unformatted')

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       read (65)
 
               do ii=0,irank-1 ! skip records
                 read(65)
	       end do

	       call read_statement

             end if

	  end do

	end if ! restart_sep

	call task_barrier()
	

	prec_xy(:,:) = 0.
	lhf_xy(:,:) = 0.
	shf_xy(:,:) = 0.
	lwns_xy(:,:) = 0.
	swns_xy(:,:) = 0.
	lwnsc_xy(:,:) = 0.
	swnsc_xy(:,:) = 0.
	lwnt_xy(:,:) = 0.
	swnt_xy(:,:) = 0.
	lwntc_xy(:,:) = 0.
	swntc_xy(:,:) = 0.
	pw_xy(:,:) = 0.
	cw_xy(:,:) = 0.
	iw_xy(:,:) = 0.
	u200_xy(:,:) = 0.
	v200_xy(:,:) = 0.
	usfc_xy(:,:) = 0.
	vsfc_xy(:,:) = 0.
        w500_xy(:,:) = 0.
        solin_xy(:,:) = 0.
        qocean_xy(:,:) = 0.

        return
        end
 
 
 
        subroutine write_statement()

        use vars
        use tracers
        use params
        implicit none

             write(65)  &
               rank, dx, dy, dz, adz, adzw, at, bt, ct, dt, dtn, dt3, time, &
               day, day0, nstep, na, nb, nc, nadams, rank, rank, caseid, case, &
               u, v, w, t, q, qp, tke, tk, tkh, p, tabs, qn, dudt, dvdt, dwdt,&
               fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,&
               fzero, precsfc, t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0,&
               p0, tg0, qg0, ug0, vg0, z, pres, rho, rhow, bet, gamaz, &
               t01, q01, wsub, qtend, ttend, prespot, tke0, betp,tauls, &
               dqls,dtls,ugls,vgls,wgls,pres0ls,dayls, zi, presi, &
               dtrfc,dayrfc,sstsfc,hsfc,lesfc,tausfc,daysfc, &
               usnd,vsnd,tsnd,qsnd,daysnd,nlsf,nrfc,nsfc,nsnd, &
               dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
               dolongwave, doshortwave, dosgs, dosubsidence, &
               docoriolis, dosurface, dolargescale,doradforcing, &
               dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
               dosmagor, doscalar,doxy,dowallx,dowally, doperpetual, doseasons, &
               docup, docolumn, dt_cup, soil_wetness, dodynamicocean, ocean_type,&
               ttend_cup, qtend_cup, utend_cup, vtend_cup, &
               pres0, ug,vg,fcor,fcorz,tabs_s,z0,sstxy,fcory,fcorzy, &
               longitude,latitude,fluxt0,fluxq0,gamt0,gamq0, &
               tau0,timelargescale,gam3,gams1,gams2,gams3,gamr1,gamr2,gamr3,&
               gamg1, gamg2, gamg3,accrsc,accrsi,accrgc,accrgi,accrrc,coefice,&
               evaps1,evaps2,evapg1,evapg2,evapr1,evapr2, a_bg, a_pr, a_gr, &
               CEM, LES, OCEAN, LAND, SFC_FLX_FXD, SFC_FLX_FXD, SFC_TAU_FXD, &
               nx, ny, nz, grdf_x, grdf_y, grdf_z, &
               ntracers, tracername, tracerunits, tracer, fluxbtr, fluxttr 
               close(65)
               if(rank.eq.nsubdomains-1) then
                  print *,'Restart file was saved. nstep=',nstep
               endif


        return
        end




        subroutine read_statement()

        use vars
        use tracers
        use params
        implicit none
        integer  nx1, ny1, nz1, rank1, ntr


             read(65) &
               rank1, dx, dy, dz, adz, adzw, at, bt, ct, dt, dtn, dt3, time, &
               day, day0, nstep, na, nb, nc, nadams, rank, rank, caseid, case, &
               u, v, w, t, q, qp, tke, tk, tkh, p, tabs, qn, dudt, dvdt, dwdt,&
               fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq,&
               fzero, precsfc, t0, q0, qc0, qv0, qi0,tabs0, tv0, rel0, u0, v0,&
               p0, tg0, qg0, ug0, vg0, z, pres, rho, rhow, bet, gamaz, &
               t01, q01, wsub, qtend, ttend, prespot, tke0, betp,tauls, &
               dqls,dtls,ugls,vgls,wgls,pres0ls,dayls, zi, presi, &
               dtrfc,dayrfc,sstsfc,hsfc,lesfc,tausfc,daysfc, &
               usnd,vsnd,tsnd,qsnd,daysnd,nlsf,nrfc,nsfc,nsnd, &
               dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo, &
               dolongwave, doshortwave, dosgs, dosubsidence, &
               docoriolis, dosurface, dolargescale,doradforcing, &
               dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
               dosmagor, doscalar,doxy,dowallx,dowally, doperpetual, doseasons, &
               docup, docolumn, dt_cup, soil_wetness, dodynamicocean, ocean_type, &
               ttend_cup, qtend_cup, utend_cup, vtend_cup, &
               pres0, ug,vg,fcor,fcorz,tabs_s,z0,sstxy,fcory,fcorzy, &
               longitude,latitude,fluxt0,fluxq0,gamt0,gamq0, &
               tau0,timelargescale,gam3,gams1,gams2,gams3,gamr1,gamr2,gamr3,&
               gamg1, gamg2, gamg3,accrsc,accrsi,accrgc,accrgi,accrrc,coefice,&
               evaps1,evaps2,evapg1,evapg2,evapr1,evapr2, a_bg, a_pr, a_gr, &
               CEM, LES, OCEAN, LAND, SFC_FLX_FXD, SFC_FLX_FXD, SFC_TAU_FXD, &
               nx1, ny1, nz1, grdf_x, grdf_y, grdf_z, &
               ntr, tracername, tracerunits, tracer, fluxbtr, fluxttr 
             if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then
              print *,'Error: domain dims (nx,ny,nz) set by grid.f'
              print *,'do not correspond to ones in the restart file.'
              print *,'in executable:   nx, ny, nz:',nx,ny,nz
              print *,'in restart file: nx, ny, nz:',nx1,ny1,nz1
              print *,'Exiting...'
              call task_abort()
             endif
             if(ntr.ne.ntracers) then
               print*,'Error: number of tracers in restart file is not the same as ntracers.'
               print*,'ntracers=',ntracers,'   ntracers(in file)=',ntr
               print*,'Exiting...'
             end if
             close(65)
             if(rank.eq.nsubdomains-1) then
                 print *,'Case:',caseid
                 print *,'Restarting at step:',nstep
                 print *,'Time:',nstep*dt
             endif


        return
        end
 
