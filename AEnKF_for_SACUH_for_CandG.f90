!**********************************************************************
! AEnKF_for_SACUH.f90 - Program performs adaptive conditional bias-   *
!                       penalized ensemble Kalman filter for ensemble *
!                       streamflow prediction using the Sacramento    *
!                       soil moisture accounting (SAC, Burnash et     *
!                       al. 1973, Burnash 1995) and unit hydrograph   *
!                       routing (UH, Chow et al. 1988) models. For    *
!                       comparison the program also runs ensesmble    *
!                       Kalman filter (EnKF) following Lorentzen and  *
!                       Naeval (2011).                                *
!                       Version Aug 08, 2021, by D.-J. Seo at UTA/    *
!                       Hydrology and Water Resources Lab             *
!**********************************************************************
 
module pass0
real, dimension (:), allocatable :: pxv0,edmnd0,qobs0,ro_old_keep
integer :: kbeg,kend
save
end module pass0

module pass1
real :: varq,varp,vare,varw
real, dimension (:), allocatable :: zvar
save
end module pass1

module pass4
integer :: nhrp,nhre,nhrw,nsp,nse,nsw,nchp,ncht
character (len=132) :: path
character (len=5) :: segment
save
end module pass4

module pass5
integer, dimension (:), allocatable :: itime
save
end module pass5

module pass7
real, dimension (:), allocatable :: state1step,state1stepc
integer :: nn,ns,nm,nw,nwin,m
save
end module pass7

program AEnKF
use pass0
use pass1
use pass4
use pass5
use pass7
!
! Specify istep, the input timestep in hrs, and ltx, the maximum lead 
! time in number of timesteps.
!
  parameter (istep=6,ltx=120/istep)
 
integer, dimension (:), allocatable :: jtime,ktime,istore
real :: lztwm,lzfpm,lzfsm,lzpk,lzsk
real, dimension (:), allocatable :: pp,frac,frad,po,state,state1,pstate
real, dimension (:), allocatable :: qstate,temp,p,qsimu,qsimu0
real, dimension (:), allocatable :: qsimu0_dum,qsimu1_dum,ybar,ymean
real, dimension (:), allocatable :: ybar1,ymean1,zobs,zobs1
real, dimension (:), allocatable :: qsimf1_ens_mean,qsimf0_ens_mean
real, dimension (:), allocatable :: qsimf1_ens_crps,qsimf0_ens_crps
real, dimension (:), allocatable :: qsimf1_opt,qtrace
real, dimension (:,:), allocatable :: qsimf,qsimf0,qsimf0_dum
real, dimension (:,:), allocatable :: qsimf1_dum,y,c1,hhat1t,w1
real, dimension (:,:), allocatable :: amatrixi,p2state_opt,p2state
real, dimension (:,:), allocatable :: sigf,sigfs,rr,sigfi,sigu,ri
real, dimension (:,:), allocatable :: sigu1,htri,sigu1a,sigu_s,sigfid
real, dimension (:,:), allocatable :: sigfid2,htrih,temp1,temp4,gaint
real, dimension (:,:), allocatable :: h,zobsp,ht,yb,ym,y_save,yy
real, dimension (:,:), allocatable :: qsimf0_ens,qsimf1_ens,aic1tri
real, dimension (:,:), allocatable :: aic1trih,aic1trihai,hamatrixi
real, dimension (:,:), allocatable :: tempmn

  dimension ett(24),randnu(1)
  character basin(14)*5,mach*20,code*2,code1*1,code2*1,nwsid*5,&
            nwsidsac*20,nwsiduhg*20,map_name*132,qin1_name*132,&
            c_beg*2,c_end*2,cindex*1,c_strong*1,rfc*2
!
! Specify the 5-character NWS IDs of the basins 1 through 14 (see Seo 
! et al. 2021) for details.
!
  data basin/'NFDC1','ABRN1','ICLI4','GTBM3','COLI2','DAVI3',&
             'DLTC1','MONN7','NIMM5','GAXV2','TRYM7','MCNM6',&
             'ORAI3','GRDN6'/
!
! Specify the beginning basin number for the DA run, [1,14].
! 
  call getarg(1,c_beg)
  read(c_beg,'(i2)') iiibeg
  write(6,*) 'beginning basin number ',iiibeg
!
! Specify the ending basin number for the DA run, [1,14], iiend >= 
! iibeg.
!
  call getarg(2,c_end)
  read(c_end,'(i2)') iiiend
  write(6,*) 'ending basin number ',iiiend
!
! Use heteroscedastic modeling of observational uncertainty by default.
!
  ihetero=1
  write(code1,'(i1)') ihetero
!
! if homoscedastic observation error variance is desired (primarily for
! testing purposes), set ihetero to 0 above and specify the observation
! error variance values for MAP (varp), MAPE (vare), additive error to 
! total channel inflow (varw) and streamflow (varq) below.
!
  if(ihetero.eq.0) then
  varp=100.          !in mm^2  (for 6hr depth)
  vare=100.          !in mm^2  (for 6hr depth)
  varw=100.          !in mm^2  (for 6hr depth)
  varq=625.          !in cms^2 (for instantaneous flow)
  endif
!
! Specify 1 for strongly-constrained DA and 0 for weakly-constrained; 
! see Lee et al. (2019) and Shen et al. (2021) for explanation. It is
! recommended that the user run both and compare the results. They 
! often reveal model parametric and/or structural errors.
!
  call getarg(3,c_strong)
  read(c_strong,'(i1)') istrong
  if(istrong.eq.1) write(6,*) 'strongly-constrained DA'
  if(istrong.eq.0) write(6,*) 'weakly-constrained DA'
  write(code2,'(i1)') istrong
!
! Build a character string for heteroscedastic/homoscedsastic
! observation uncertainty (code1) and strongly- vs. weakly-constrained 
! DA (code2).
!
  code=code1//code2
!
! Specify the last character of the subdireectory for output files.
!
  if(istrong.eq.0) cindex='0'          !subdirectory name is '/output0'
  if(istrong.eq.1) cindex='1'          !subdirectory name is '/output1'
!
! run EnKF and AEnKF when iaenkf=0 and iaenkf=1, respectively.
!
  do iaenkf=0,1
!
! Set system clock to track time elapsed.
! 
  call system_clock(icount1, icount_rate1, icount_max1)
!
! Run DA for the range of basins selected, [iiibeg,iiiend]
!
  do iii=iiibeg,iiiend
!
! Specify the ensemble size.
!
  ns=200   
!
! Spcify the number of streamflow obs to be assimilated. The default is
! 1, i.e, use the streamflow obs valid at prediction time (see Fig 1 of
! Lee et al. 2019 for schematic of the assimilation window). for SAC-
! UH, nf > 1 is not recommended as it is very likely to present
! numerical singularity issues due to the highly correlated nature of 
! the convolution operation of UH.
!
  nf=1
!
! Specify the sampling frequency of streamflow obs to be assimilated.
! For nf=1, this parameter is irrelevant. If, e.g., nf=2, the stream-
! flow obs valid at timesteps k and k-ifrq are used where k is
! associated with the end of the assimilationn window or the prediction
! time.
!
  ifrq=1  
!
! Specify the random number seed.
!
  jseed=11111
!
! Get the 5-character NWS ID for this basin.
!
  segment=basin(iii)
  ncht=len_trim(segment)
  nwsid=segment(1:ncht)
  write(6,*) 'This segment is ',nwsid
!
! Specify the cutoff flow (qcut) for screening of significant events.
!
  call get_basin_info(nwsid,rfc,ncs,nwsidsac,ncu,nwsiduhg,qcut)
!
! Specify the directory path for output files
!
  if(iaenkf.eq.0) then
!
! for EnKF output, and
!
  mach='run_'//segment(1:ncht)
  nchm=len_trim(mach)
  path='./output'//cindex//'/'//segment(1:ncht)//'/'//&
  mach(1:nchm)//'_enkf_'//code//'/'
  nchp=len_trim(path)

  else if(iaenkf.eq.1) then
!
! for AEnKF output.
!
  mach='run_'//segment(1:ncht)
  nchm=len_trim(mach)
  path='./output'//cindex//'/'//segment(1:ncht)//'/'//&
  mach(1:nchm)//'_aenkf_'//code//'/'
  nchp=len_trim(path)

  endif
!
! Create the directory for output files.
!
  write(6,*) 'output directory path ', path(1:nchp)
  call system('mkdir -p '//path(1:nchp))
!
! Read UH.
!
  call get_uhg1(nwsid,rfc,ncu,nwsiduhg,nv,c_strong)
allocate(po(nv))
  po=0.
  call get_uhg2(nwsid,nv,po,idtr,area,cbflow,c_strong)
!
! Evaluate the double integral, integral integral u(time-t) u(time-s) 
! ds dt, for modeling uncertainty for additive error to TCI (see 
! Rafieeinasab et al. 2014, Shen et al. 2021). The resulting sumd is 
! used in subroutine get_obs_uncertainty below.
!
  sumd=0.
  do i=1,nv
  do j=1,nv
  sumd=sumd+po(nv-i+1)*po(nv-j+1)
  enddo
  enddo
!
! Specify the size of the assimilation window. The window ends at the
! prediction time, i.e., the current time (see Fig 1 of Lee et al. 2019
! for schematic of the assimilation window).
!
  nwin=nv-1+nf*ifrq
!
! Specify the time scale of adjustment for MAP, MAPE and additive error 
! to total channel inflow (TCL).
!
! Adjust at the scale of the entire assimilation window (i.e., uniform
! adjustment); this is the default.
!
  nhrp=nwin
  nhre=nwin
  nhrw=nwin
!
! Adjust at the scale of 2*timestep.
!
! nhrp=2
! nhre=2
! nhrw=2
!
! Adjust at the scale of the timestep.
!
! nhrp=1
! nhre=1
! nhrw=1
!
! If necessary, modify the assimilation window so that it is an integer
! multiple of the adjustment scale.
!
  rsub=float(nwin)/nhrp
  nsub=int(rsub)
  if(rsub-nsub.ge.0.5) nsub=nsub+1
  nwin=nsub*nhrp
!
! Specify the number of state variables for adjustment of MAP, MAPE and
! additive error to TCI.
!
  nsp=nwin/nhrp
  nse=nwin/nhre
  nsw=nwin/nhrw
!
! Error-check.
!
  if(nsp.ne.float(nwin)/nhrp) then
  write(6,*) 'nwin not an integer multiple of nsp...stop'
  stop
  endif
  if(nse.ne.float(nwin)/nhre) then
  write(6,*) 'nwin not an integer multiple of nse...stop'
  stop
  endif
  if(nsw.ne.float(nwin)/nhrw) then
  write(6,*) 'nwin not an integer multiple of nsw...stop'
  stop
  endif
!
! Specify the number of non-augmented state variables: 6 SAC states, 
! nsp+nse+nsw states for biases in MAP and MAPE and additive error to
! TCI.
!
  if(istrong.eq.0) then
  nn=6+nsp+nse+nsw                         !for weakly-constrained DA
  nw=3
  else
  nn=6+nsp+nse                             !for strongly-constrained DA
  nw=2
  endif

allocate(frac(nn))
allocate(frad(nn))
allocate(state1(nn))
!
! Prescribe the model dynamical error expressed as fraction of the 
! state.
!
  do i=1,nn
  frac(i)=0.01     !for get_perturbed_aug_states0 (for cold DA)
  frad(i)=0.01     !for get_perturbed_aug_states1 (for warm DA) 
!
! For TRYM7 and MCNM6, use smaller dynamical errors based on 
! sensitivity analysis.
!
  if(iii.eq.11.or.iii.eq.12) then
  frac(i)=0.0001
  frad(i)=0.0001
  endif
  enddo
!
! Read the SAC parameters and initial conditions (IC).
!
  call get_sac_params(nwsid,rfc,ncs,nwsidsac,rexp,lzpk,lzfpm,&
       pxadj,idt,pfree,zperc,riva,peadj,lztwm,rserv,adimp,uzk,&
       side,lzfsm,lzsk,uztwm,uzfwm,pctim,efc,ett,saved,parea,&
       nn,state1,c_strong)
!
! Specify the number of observations. For each timestep, there are
! MAP, MAPE and faux obs for TCI error plus the streamflow obs at
! the end of the assimilation window (i.e., if nf=1).
!
! Note that the TCI error is not observed in reality. As such, We 
! assume a priori that TCI has no errors, and assign a faux error of 0 
! over the entire assimilation window.
!
! Specify the total number of observations or measurements.
!
  nm=nw*nwin+nf
!
! Specify the total number of state variables, including the augmented.
!
  m=nn+nf

allocate(ro_old_keep(nv))
allocate(state1step(nn))
allocate(state1stepc(nn))
allocate(state(nn))
allocate(pstate(nn))
allocate(qstate(nn))
allocate(temp(nn))
allocate(p(nn))
allocate(qsimu(nf))
allocate(qsimu0(nf))
allocate(qsimu0_dum(nf))
allocate(qsimu1_dum(nf))
allocate(istore(nf))
allocate(qsimf(ltx+nwin,nn+2))
allocate(qsimf0(ltx+nwin,nn+2))
allocate(qsimf0_dum(ltx+nwin,nn+2))
allocate(qsimf1_dum(ltx+nwin,nn+2))
allocate(qsimf1_ens_mean(ltx+nwin))
allocate(qsimf0_ens_mean(ltx+nwin))
allocate(qsimf1_ens_crps(ltx+nwin))
allocate(qsimf0_ens_crps(ltx+nwin))
allocate(qsimf1_opt(ltx+nwin))
allocate(p2state_opt(nn,ns))
allocate(p2state(nn,ns))
allocate(zobs(nm))
allocate(zobs1(nm))
allocate(zobsp(nm,ns))
allocate(qtrace(ns))
allocate(qsimf0_ens(ltx+nwin,ns))
allocate(qsimf1_ens(ltx+nwin,ns))
allocate(ybar(m))
allocate(ymean(m))
allocate(ybar1(m))
allocate(ymean1(m))
allocate(yb(m,ns))
allocate(ym(m,ns))
allocate(sigf(m,m))
allocate(sigfs(m,m))
allocate(rr(m,m))
allocate(sigfi(m,m))
allocate(sigu(m,m))
allocate(ri(nm,nm))
allocate(sigu1(m,m))
allocate(htri(m,nm))
allocate(sigu1a(m,m))
allocate(sigu_s(m,m))
allocate(sigfid(m,m))
allocate(sigfid2(m,m))
allocate(htrih(m,m))
allocate(temp1(m,m))
allocate(temp4(m,m))
allocate(pp(nm))
allocate(zvar(nm))
allocate(y(m,ns))
allocate(y_save(m,ns))
allocate(gaint(nm,m))
allocate(h(nm,m))
allocate(ht(m,nm))
allocate(yy(nn,ns*2))
allocate(c1(nm,m))
allocate(hhat1t(nm,m))
allocate(w1(m,nm))
allocate(amatrixi(m,m))
allocate(aic1tri(m,nm))
allocate(aic1trih(m,m))
allocate(aic1trihai(m,m))
allocate(hamatrixi(nm,m))
allocate(tempmn(m,nm))
!
! Open files for MAP and streamflow (QIN) time series.
!
  map_name='new_map06_'//segment(1:ncht)
  nchx=len_trim(map_name)
  qin1_name=segment(1:ncht)//'.qin'
  nchq=len_trim(qin1_name)
!
! Read MAP.
!
  call get_map_length(nchx,map_name,jdata1)
allocate(jtime(jdata1))
  call get_map(nchx,map_name,jdt,jdata,jtime,jdata1)
!
! Read QIN.
!
  call get_qin_length(nchq,qin1_name,kdata1)
allocate(ktime(kdata1))
  call get_qin(nchq,qin1_name,kdata,ktime,kdata1)
!
! Identify the common period of record.
!
  call get_period_of_record(jdata,jtime,kdata,ktime,jbeg1,kbeg1,&
       jend1,kend1)
!
! Write the lengths and beginning and ending indices of the MAP and QIN
! time series.
!
  write(6,*) 'jdata,kdata,jbeg1,kbeg1,jend1,kend1'
  write(6,*)  jdata,kdata,jbeg1,kbeg1,jend1,kend1
!
! Write the beginning and ending times of MAP and QIN.
!
  write(6,*) jtime(jbeg1),ktime(kbeg1)
  write(6,*) jtime(jend1),ktime(kend1)

deallocate(jtime)
deallocate(ktime)
!
! Get the length of the common period of record.
!
mdata=max(jend1-jbeg1+1,kend1-kbeg1+1)

allocate(pxv0(mdata))
allocate(edmnd0(mdata))
allocate(qobs0(mdata))
allocate(itime(mdata))
!
! Get the MAP and QIN time series within the common period of record.
!
  call get_map_mape_qin(nchx,map_name,nchq,qin1_name,idt,ett,&
       nwsid,ndata,jbeg1,kbeg1,mdata,pxv0,edmnd0,qobs0,itime,dt)
!
! Error-check.
!
  if(ndata.gt.mdata) then
  write(6,*) 'ndata gt mdata...stop ',ndata,mdata
  stop
  endif
!
! Run SAC-UH over the entire common period of record to generate and 
! store the base (i.e., DA-less) single-valued SAC states.
!
  kbeg=1
  kend=ndata
!
! Initialize the carry-over from the UH operation.
!
  ro_old_keep=-0.1
!
! Write state variables while forward-integrating.
!
  iseed=jseed
  call randd(iseed,1,randnu(1))
  iseed=0

  knd=1
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,state1,&
             knd,qsimu,qsimf,0,iseed,gset,kf,ifrq,istore,nf,sum)
!
! Open output files
!
! for SAC states from base single-valued simulation,
!
  open(11,file=path(1:nchp)//segment(1:ncht)//'.statevar_out',&
  status='unknown',access='direct',form='unformatted',recl=24)
!
! for verifying observed flow,
!
  open(13,file=path(1:nchp)//segment(1:ncht)//'.fcst_obs_ctl',&
  status='unknown')
!
! for single-valued streamflow prediction,
!
  open(14,file=path(1:nchp)//segment(1:ncht)//'.fcst_sngl_ctl',&
  status='unknown')
!
! for EnKF ensemble streamflow prediction,
!
  open(22,file=path(1:nchp)//segment(1:ncht)//'.fcst_enkf_ctl',&
  status='unknown')
!
! for base ensemble streamflow prediction,
!
  open(25,file=path(1:nchp)//segment(1:ncht)//'.fcst_base_ctl',&
  status='unknown')
!
! for AEnKF ensemble streamflow prediction,
!
  open(32,file=path(1:nchp)//segment(1:ncht)//'.fcst_aenkf_ctl',&
  status='unknown')
!
! for diagnostics,
!
  open(42,file=path(1:nchp)//segment(1:ncht)//'.dfs',& 
  status='unknown')
!
! for mean CRPS of base ensemble streamflow prediction,
!
  open(51,file=path(1:nchp)//segment(1:ncht)//'.crps_base',&
  status='unknown')
!
! for mean CRPS of EnKF ensemble streamflow prediction, and
! 
  open(52,file=path(1:nchp)//segment(1:ncht)//'.crps_enkf',&
  status='unknown')
!
! for mean CRPS of AEnKF ensemble streamflow prediction.
!
  open(53,file=path(1:nchp)//segment(1:ncht)//'.crps_aenkf',&
  status='unknown')
!
! Initialize valid time for the preceding assimilation cycle.
!
  itime_old=-9999
!
! Perform assimilation in each timestep.
!
c20 : do i=nwin,ndata-ltx
!
! Get the current time.
!
  call time1(itime(i),iyr,mon,iday,ihr)
!
! Skip DA if streamflow observation valid at this time step does not
! exist.
!
  if(qobs0(i).lt.0.) cycle c20
!
! Skip if there is a missing streamflow observation within the forecast
! horizon.
!
  do j=1,ltx
  if(qobs0(i+j).lt.0.) cycle c20
  enddo
!
! Get the beginning and ending times of the assimilation window.
!
  kbeg=i-nwin+1
  kbeg_save=kbeg
  if(kbeg.lt.1) cycle c20
  kend=i
!
! Get the beginning and ending times of the carry-over window for 
! convolution of UH.
!
  ibeg=i-nv+1
  if(ibeg.lt.1) cycle c20
  iend=kbeg-1
!
! Skip DA if there are missing input data in the union of windows.
!
  do k=min(ibeg,kbeg),max(iend,kend)
  if(pxv0(k).lt.0.) cycle c20
  if(edmnd0(k).lt.0.) cycle c20
  if(qobs0(k).lt.0.) cycle c20
  enddo
!
! Check if there exists a significant flow in the next 5 days.
!
  nhr_look_forward=120/istep
  isig=0
c22 : do j=1,nhr_look_forward
  if(i+j.gt.ndata) cycle c22
  if(qobs0(i+j).gt.qcut) isig=isig+1
  enddo c22
!
! Check if there exists a significant flow in the preceding 5 days.
!
  nhr_look_backward=120/istep
  jsig=0
c26 : do j=1,nhr_look_backward
  if(i-j.lt.1) cycle c26
  if(qobs0(i-j).gt.qcut) jsig=jsig+1
  enddo c26
!
! Skip DA if no significant flow exists within the 10-day window.
!
  if(isig+jsig.eq.0) cycle c20
!
! Determine rising vs. receding flow for conditional evaluation.
!
  irise=-9999
  drise=-9999.
  if(qobs0(i+1).gt.qobs0(i)) irise= 1
  if(qobs0(i+1).lt.qobs0(i)) irise=-1
  if(qobs0(i+1).eq.qobs0(i)) irise= 0
  drise=qobs0(i+1)-qobs0(i)
!
! Specify observational uncertainties.
!
  call get_obs_uncertainty(i,nwsid,nwin,nw,nf,nm,varq,varp,vare,&
       varw,sumd,zvar,ihetero,ifrq,ndata)
!
! Initialize carry-over for UH convolution.
! 
  ro_old_keep=-0.1
!
! Get time elapsed (in hrs) since the last DA.
!
  call diftime(itime_old,itime(i),lag)
!
! Specify the a priori SAC states from the direct access file.
!
  read(11,rec=kbeg) (state1(j),j=1,6)
  do j=1,nsp
  state1(6+j)=pxadj
  enddo
  do j=1,nse
  state1(6+nsp+j)=peadj
  enddo
  if(nw.gt.2) then
  do j=1,nsw
  state1(6+nsp+nse+j)=0.
  enddo
  endif
!
! Check feasible range before generating DA-less single-value forecast.
!
  call check_range(2,state1,uztwm,uzfwm,lztwm,lzfsm,lzfpm,nsp,nse,&
                   nsw,nn,nw) 
!
! Save the state variables before being overwritten in subroutine 
! sacuh.
!
  state=state1
!
! Predict flow based on the raw model ICs; this is the single-valued
! base simulation.
!
  knd=3
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,state,&
             knd,qsimu0,qsimf0,0,iseed,gset,kf,ifrq,istore,nf,sum)

  if(lag.ne.istep) then
!
! This is a cold assimilation cycle, i.e., DA was not run in the 
! preceding timestep.
!
! Perturb state1 to obtain pstate, and propagate it one timestep to
! obtain y; this y is used in the next assimilation cycle.
!
  call get_perturbed_aug_states0(uztwm,uzfwm,lztwm,lzfsm,lzfpm,&
       nsp,nse,nsw,state1,iseed,p,gset,frac,pstate,&
       uzk,pctim,adimp,riva,zperc,rexp,lzsk,lzpk,&
       pfree,side,saved,parea,dt,nv,po,pxadj,peadj,knd,&
       y,kf,ifrq,nf)

  else 
!
! This is a warm assimilation cycle, i.e., DA was run in the preceding
! timestep.
!
! Build the observation vector for MAP, MAPE, additive error to TCI and
! streamflow. Note that the prediction time is i and hence, if nf=1,
! only the streamflow obs valid at time step i is assimilated.
!
  do k=kbeg,kend
  i1=k-kbeg+1
  zobs(    i1)=  pxv0(k)
  zobs(   nwin+i1)=edmnd0(k)
  if(nw.gt.2) then
  zobs( 2*nwin+i1)=   0.
  endif
  enddo
  do j=1,nf
  if(i+(j-nf)*ifrq.lt.1.or.i+(j-nf)*ifrq.gt.ndata) then
  write(6,*) 'qobs0 out of range...stop 2'
  stop
  endif
  zobs(nw*nwin+ j)=qobs0(i+(j-nf)*ifrq)
  enddo
!
! Perturb observations.
!
  do ii=1,ns
  call get_perturbs_for_obs(nm,iseed,pp,zvar,gset)
  do j=1,nm
  zobsp(j,ii)=zobs(j)+pp(j)
!
! Check non-negativity of MAP, MAPE and streamflow obs.
!
  if(j.le.(nw-1)*nwin.or.j.gt.nw*nwin) then
  if(zobsp(j,ii).lt.0.) zobsp(j,ii)=0.
  endif
  enddo
  enddo
!
! Calculate sample mean and variance of perturbed obs.
!
  do ik=1,nm
  rmean_ens=0.
  std_ens=0.
!
! Calculate the ensemble mean of the ik-th obs.
!
  do ii=1,ns
  rmean_ens=rmean_ens+zobsp(ik,ii)
  enddo
  rmean_ens=rmean_ens/ns
!
! Calculate ensembles standard deviation of the ik-th obs.
!
  std_ens=0.
  do ii=1,ns
  std_ens=std_ens+(zobsp(ik,ii)-rmean_ens)**2
  enddo
  std_ens=sqrt(std_ens/(ns-1))
!
! Shift and scale the perturbed observations to match the mean of 
! oserved value and variance of zvar, i.e., the prescribed observation
! uncertainty.
!
  do ii=1,ns
  zobsp(ik,ii)=((zobsp(ik,ii)-rmean_ens)/std_ens)*sqrt(zvar(ik))+&
  zobs(ik)
  enddo
  enddo
!
! If perturbed MAP, MAPE or streamflow is negative, set to zero.
!
  do ii=1,ns
  do ik=1,nwin
  if (zobsp(ik    ,ii).lt.0) zobsp(ik    ,ii)=0.
  if (zobsp(ik+   nwin,ii).lt.0) zobsp(ik+   nwin,ii)=0.
  if(nw.gt.2) then
  if (zobsp(ik+ 2*nwin,ii).lt.0) zobsp(ik+ 2*nwin,ii)=0.
  endif
  enddo
  do ik=1,nf
  if (zobsp(ik+nw*nwin,ii).lt.0) zobsp(ik+nw*nwin,ii)=0.
  enddo
  enddo
!
! Save the ensemble state, y, from the preceding time step valid at the
! beginning of the assimilation window.
!
  do ii=1,ns 
  do ik=1,m
  y_save(ik,ii)=y(ik,ii)
  enddo
  enddo
!
! Initialize the optimal alpha and degrees of freedom for noise (dfn).
!
  ialpha_opt=-999
  dfn_opt=9999999.
!
! Specify the maximum number of iterations allowed.
!
  if(iaenkf.eq.0) ie=1                                           !EnKF
  if(iaenkf.eq.1) ie=100                                         !AEnKF

c666 : do ialpha=1,ie

  if(ialpha.eq.1) alpha=0.
!
! Prescribe the ensemble IC.
!
  do ii=1,ns
  do ik=1,m
  y(ik,ii)=y_save(ik,ii)
  enddo
  enddo
!   
! Generate model-simulated observations within the assimilation window
! to augment the state vector.
!
  do ii=1,ns
 
  do j=1,nn
  temp(j)=y(j,ii)
  enddo

  knd=2
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,temp,&
             knd,qsimu0_dum,qsimf0_dum,0,iseed,gset,kf,ifrq,&
             istore,nf,sum)

  do j=1,nf
  y(nn   + j,ii)=qsimu0_dum(j)
  enddo

  enddo
! 
! generate DA-less ensemble forecast.
!
! Get perturbations for the state variables. Note that, for every 
! assimilation cycle, the states from the original single-valued run
! are perturbed, rather than those propogated from the preceding 
! assimilation cycle. This explains the large mean CRPS for DA-less
! analysis and very short-term prediction.
!
  do ii=1,ns

  call get_perturbs_for_states0(uztwm,uzfwm,lztwm,lzfsm,&
                                lzfpm,nsp,nse,&
                                nsw,state1,iseed,p,gset,frac,nn,nw)
!
! Add the perturbations to the state variables.
!
  temp=state1+p
!
! Check the feasible range.
! 
  call check_range(2,temp,uztwm,uzfwm,lztwm,lzfsm,&
                   lzfpm,nsp,nse,nsw,nn,nw)
!
! Make a prediction run of SACUH for this ensemble member.
!
  knd=3
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,temp,&
             knd,qsimu0_dum,qsimf0_dum,1,iseed,gset,kf,ifrq,&
             istore,nf,sum)

  do j=1,ltx+nwin
  qsimf0_ens(j,ii)=qsimf0_dum(j,1)    !base ensemble prediction for
              !streamflow
  enddo

  enddo
!
! Calculate the ensemble mean for DA-less base ensemble prediction.
!
  do j=1,ltx+nwin
  qsimf0_ens_mean(j)=0.
  do ii=1,ns
  qsimf0_ens_mean(j)=qsimf0_ens_mean(j)+qsimf0_ens(j,ii)
  enddo
  qsimf0_ens_mean(j)=qsimf0_ens_mean(j)/ns
  enddo
!
! Calculate selected quantiles for plotting.
!
  do j=nwin,nwin
  do ii=1,ns
  qtrace(ii)=qsimf0_ens(j,ii)
  enddo
  call hpsort(ns,qtrace)
!
! Get the 99th, 90th, 10th and 1st percentile flows for base ensemble
! prediction.
!
  threshold_frac=0.99
  m1=int(threshold_frac*ns)
  b99=qtrace(m1)
  threshold_frac=0.90
  m1=int(threshold_frac*ns)
  b90=qtrace(m1)
  threshold_frac=0.10
  m1=int(threshold_frac*ns)
  b10=qtrace(m1)
  threshold_frac=0.01
  m1=int(threshold_frac*ns)
  b01=qtrace(m1)
  enddo
!
! Calculate mean CRPS for base ensemble prediction using Hersbach 
! (2000).
!
  do j=1,ltx+nwin
  do ii=1,ns
  qtrace(ii)=qsimf0_ens(j,ii)
  enddo
  call hpsort(ns,qtrace)
  call calc_crps(ns,qtrace,qobs0(i+j-nwin),crps)
  qsimf0_ens_crps(j)=crps
  enddo
!
! Write the base ensemble prediction.
!
  if(iaenkf.eq.0) then
  write(51,70) (qsimf0_ens_crps(j),j=nwin,nwin+ltx)
  endif
!
! Calculate ensemblemean forecast, ybar, and save the ensemble
! forecast in yb.
!
  do ik=1,m
  ybar(ik)=0.
  do ii=1,ns
  ybar(ik)=ybar(ik)+y(ik,ii)
  yb(ik,ii)=y(ik,ii)
  enddo
  ybar(ik)=ybar(ik)/ns
  enddo
!
! Calculate the forecast error covariance.
!
  do ik=1,m
  do ij=1,m
  sigf(ik,ij)=0.
  do ii=1,ns
  sigf(ik,ij)=sigf(ik,ij)+y(ik,ii)*y(ij,ii)
  enddo
  sigf(ik,ij)=sigf(ik,ij)/ns-ybar(ik)*ybar(ij)
  enddo
  enddo

  del=0.
  100 do ik=1,m
  do ij=1,m
  if(ik.eq.ij) then
  sigf(ik,ij)=sigf(ik,ij)+del
  rr(ik,ij)=1.
  else
  rr(ik,ij)=0.
  endif
  sigfs(ik,ij)=sigf(ik,ij)
  enddo
  enddo
!
! Check positive definiteness of sigf by inverting.
!
  call lsolve_cd(m,sigfs,m,rr,sigfi,ierr)

  if(ierr.ne.0) then
!
! If not positive definite, add a small positive number to the 
! diagonal.
!
  del=del+0.0001
  write(6,*) 'ierr ne 0 sigfs from lsolve_cd 0...add del ',ierr,del
  go to 100
  endif
!
! Calculate gain
!
  if(iaenkf.eq.1) then 
!
! for AEnKF, and
!
  call eval_aenkf_gain(nm,m,sigf,zvar,gaint,alpha,nwin,nf,sigu,&
                       sigu1,nw,ierr,h,ht,ri)

  if(alpha.gt.0.) then
!
! Error encountered...terminate the iteration.
!
  if(ierr.ne.0) exit c666
  endif

  else
!
! for EnKF.
!
  call eval_enkf_gain(nm,m,sigf,zvar,gaint,nwin,nf,nw,ierr,h,ht,&
                      ri,sigu)

  if(ierr.ne.0) then
  write(6,*) 'ierr ne 0 from eval_enkf_gain...stop ',ierr
  stop
  endif

  endif
!
! Update states
!
  do ij=1,nn
  do ii=1,ns

  sumg=0.
!
! due to MAP,
!
  do jj=1,nsp
  ibeg=(jj-1)*nhrp+1
  iend=jj*nhrp
  do ik=ibeg,iend
  resid=zobsp(ik,ii)-pxv0(kbeg+ik-1)*y(6+jj,ii)
  sumg=sumg+gaint(ik,ij)*resid
  enddo
  enddo
!
! due to MAPE,
!
  do jj=1,nse
  ibeg=nwin+(jj-1)*nhre+1
  iend=nwin+jj*nhre
  do ik=ibeg,iend
  resid=zobsp(ik,ii)-edmnd0(kbeg+ik-nwin-1)*y(6+nsp+jj,ii)
  sumg=sumg+gaint(ik,ij)*resid
  enddo
  enddo
!
! due to error to TCI,
!
  if(nw.eq.3) then

  do jj=1,nsw
  ibeg=2*nwin+(jj-1)*nhrw+1
  iend=2*nwin+jj*nhrw
  do ik=ibeg,iend
  resid=zobsp(ik,ii)-y(6+nsp+nse+jj,ii)
  sumg=sumg+gaint(ik,ij)*resid
  enddo
  enddo
!
! due to Qobs for the weakly-constrained, and, if applicable,
!
  do jj=1,nf
  do ik=3*nwin+1,3*nwin+nf
  resid=zobsp(ik,ii)-y(6+nsp+nse+nsw+jj,ii)
  sumg=sumg+gaint(ik,ij)*resid
  enddo
  enddo

  endif

  if(nw.eq.2) then
!
! due to Qobs for the strongly-constrained.
!
  do jj=1,nf
  do ik=2*nwin+1,2*nwin+nf
  resid=zobsp(ik,ii)-y(6+nsp+nse+jj,ii)
  sumg=sumg+gaint(ik,ij)*resid
  enddo
  enddo

  endif
 
  p2state(ij,ii)=y(ij,ii)+sumg

  enddo
  enddo
!
! State updating is complete. The updated ensemble states are held in
! p2state. Prepare to make DA-aided ensemble prediction.
!
  do ii=1,ns
!
! Copy a member of the updated ensemble state.
!
  do j1=1,nn
  qstate(j1)=p2state(j1,ii)
  enddo
!
! Check for feasible range.
!
  call check_range(2,qstate,uztwm,uzfwm,lztwm,lzfsm,&
                   lzfpm,nsp,nse,nsw,nn,nw)

  do j1=1,nn
  y(j1,ii)=qstate(j1)
  enddo

  enddo
!
! Calculate ensemble mean of updated states, ymean, and save the
! updated ensemble analysis in ym.
!
  do j1=1,nn
  ymean(j1)=0.
  do ii=1,ns
  ymean(j1)=ymean(j1)+y(j1,ii)
  ym(j1,ii)=y(j1,ii)
  enddo
  ymean(j1)=ymean(j1)/ns
  enddo
!
! Initialize.
!
  dfn_ave=0.

  do j=1,nf
  ymean(nn+j)=0.
  enddo

  do ii=1,ns
!
! Copy each member of the ensemble IC for prediction.
!
  do j1=1,nn
  qstate(j1)=y(j1,ii)
  enddo

  knd=3
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,qstate,&
             knd,qsimu1_dum,qsimf1_dum,1,iseed,gset,kf,ifrq,&
             istore,nf,sum)
!
! Save the updated ensemble streamflow analysis in ym.
!
  do j=1,nf
  ymean(nn+j)=ymean(nn+j)+qsimu1_dum(j)
  ym(nn+j,ii)=qsimu1_dum(j)
  enddo
!
! Save the DA-aided ensemble prediction in qsimf1_ens.
!
  do j=1,ltx+nwin
  qsimf1_ens(j,ii)=qsimf1_dum(j,1)     !DA-aided ens prediction
  enddo

  dfn_ave=dfn_ave+sum

  enddo
!
! Calculate the ensemble mean of updated augmented states. Below, ymean
! holds the updated ensemble mean states.
!
  do j=1,nf
  ymean(nn+j)=ymean(nn+j)/ns
  enddo
!
! Calculate the sample analysis error covariance.
! 
  do ik=1,m
  do ij=1,m
  sigu_s(ik,ij)=0.
  do ii=1,ns
  sigu_s(ik,ij)=sigu_s(ik,ij)+ym(ik,ii)*ym(ij,ii)
  enddo
  sigu_s(ik,ij)=sigu_s(ik,ij)/ns-ymean(ik)*ymean(ij)
  enddo
  enddo
!
! Check modeled updated error covariance, sigu, vs. the actual, sigu_s.
!
  sumu_t=0.
  sumu1t=0.
  sumu_s=0.
!     do j=1,m
  do j=m,m             !for streamflow only
  sumu_t=sumu_t+sigu(j,j)
  sumu1t=sumu1t+sigu1(j,j)
  sumu_s=sumu_s+sigu_s(j,j)
  enddo
!
! Calculate mean dfn.
!
  dfn_ave=dfn_ave/ns
!
! Calculate ensemble mean.
!
  do j=1,ltx+nwin
  qsimf1_ens_mean(j)=0.
  do ii=1,ns
  qsimf1_ens_mean(j)=qsimf1_ens_mean(j)+qsimf1_ens(j,ii)
  enddo
  qsimf1_ens_mean(j)=qsimf1_ens_mean(j)/ns
  enddo

  if(iaenkf.eq.1) then
!
! Optimize alpha for AEnKF.
!
  alpha_sum=0.
!
! Evaluate terms that are not ensemble member-specific.
!
  call eval_opt_alpha1(m,alpha,sigfi,h,ht,nm,ri,sigu1,&
       sigfid,sigfid2,htri,temp1,sigu1a,temp4,htrih)

  do ii=1,ns

  do j1=1,m
  ymean1(j1)=ym(j1,ii)       !updated ensemble analysis
  ybar1(j1)=yb(j1,ii)        !ensemble forecast
  enddo

  do i1=1,nm
  zobs1(i1)=zobsp(i1,ii)
  enddo
!
! Evaluate terms that are ensemble member-specific.
!
  alpha1=alpha
  call eval_opt_alpha2(m,alpha1,h,nm,ymean1,ybar1,&
       zobs1,deriv1m,deriv2m,htri,temp1,temp4,htrih,ierr)
  if(ierr.ne.0) exit c666

  alpha_sum=alpha_sum+alpha1

  enddo
!
! Calculate average alpha; this is the new alpha.
!
  alpha=alpha_sum/ns

  endif
!
! Calculate mean CRPS.
!
  do j=1,ltx+nwin
  do ii=1,ns
  qtrace(ii)=qsimf1_ens(j,ii)
  enddo

  call hpsort(ns,qtrace)

  call calc_crps(ns,qtrace,qobs0(i+j-nwin),crps)
  qsimf1_ens_crps(j)=crps
  enddo
!
! Write verifying obs, base single-valued prediction, base ensemble 
! predictionn, EnKF prediction, and mean CRPS for EnKF prediction.
!
  if(iaenkf.eq.0) then
!
! Write verifying observed flow.
!
  write(13,70) (qobs0(j),j=i  ,i+ltx)
   70 format(121f9.2)
!
! Write single-valued prediction.
!
  write(14,70) (qsimf0(j,1),j=nwin  ,nwin+ltx)
!
! Write base ensemble prediction.
!
  write(25,70) (qsimf0_ens_mean(j),j=nwin  ,nwin+ltx)
!
! Write EnKF prediction.
!
  write(22,70) (qsimf1_ens_mean(j),j=nwin  ,nwin+ltx)
!
! Write mean CRPS for EnKF prediction.
!
  write(52,70) (qsimf1_ens_crps(j),j=nwin,nwin+ltx)

  endif
!
! Calculate percent reduction in dfn relative to the current optimum.
!
  per_redu=100.*(dfn_opt-dfn_ave)/dfn_opt

  if(per_redu.gt.0.) then 
!
! AEnKF reduced dfn; update the optimum. Note that this condition is 
! always met initially (i.e., with alpha=0).
!
  sumu_topt=sumu_t
  sumu1topt=sumu1t
  sumu_sopt=sumu_s

  ialpha_opt=ialpha
  alpha_opt=alpha
  dfn_opt=dfn_ave

  do ii=1,ns
  do j1=1,nn
  p2state_opt(j1,ii)=p2state(j1,ii)
  enddo
  enddo

  do j=1,ltx+nwin
  qsimf1_opt(j)=qsimf1_ens_mean(j)
  enddo

!     do j=1,ltx+nwin
  do j=nwin,nwin             !analysis only

  do ii=1,ns
  qtrace(ii)=qsimf1_ens(j,ii)
  enddo

  call hpsort(ns,qtrace)
!
! calculate selected percentiles for plotting.
!
  threshold_frac=0.99
  m1=int(threshold_frac*ns)
  q99=qtrace(m1)
  threshold_frac=0.90
  m1=int(threshold_frac*ns)
  q90=qtrace(m1)
  threshold_frac=0.10
  m1=int(threshold_frac*ns)
  q10=qtrace(m1)
  threshold_frac=0.01
  m1=int(threshold_frac*ns)
  q01=qtrace(m1)

  enddo

  else
  write(6,*) 'no reduction in dfn...exit ',per_redu 
  exit c666 

  endif 

  enddo c666
!
! Write the AEnKF ensemble mean prediction.
!
  write(32,70) (qsimf1_opt(j),j=nwin  ,nwin+ltx)
!
! Write mean CRPS of AEnKF ensemble prediction.
!
  write(53,70) (qsimf1_ens_crps(j),j=nwin,nwin+ltx)
!
! Write miscellaneous results for diagnostics.
!
  write(42,*) itime(i),irise,ialpha_opt,alpha_opt,rdum1,&
  qobs0(i),rdum2,qsimf1_opt(nwin),q01,q10,q90,q99,qsimf0(nwin,1),&
  idum1,idum2,drise,b01,b10,b90,b99,sumu_topt,sumu1topt,sumu_sopt
!
! Copy the optimal ensemble state for the next time step; its
! feasibility is checked in subroutine get_perturbed_aug_states1 below.
!
  do ii=1,ns
  do j1=1,nn
  p2state(j1,ii)=p2state_opt(j1,ii)
  enddo
  enddo
!
! Propagate the ensemble state one timestep for the next assimilation
! cycle.
!
  call get_perturbed_aug_states1(uztwm,uzfwm,lztwm,lzfsm,lzfpm,&
       nsp,nse,nsw,iseed,gset,frad,p2state,&
       uzk,pctim,adimp,riva,zperc,rexp,lzsk,lzpk,&
       pfree,side,saved,parea,dt,nv,po,pxadj,peadj,knd,&
       y,kf,ifrq,nf)

  endif 
!
! At this stage, the ensemble state valid at the beginning of the 
! assimilation window, y, already exist.
!
  itime_old=itime(i)

  enddo c20

  close(11)
  close(13)
  close(14)
  close(22)
  close(25)
  close(42)
  close(32)
  close(51)
  close(52)
  close(53)
!
! deallocate arrays
!
deallocate(po) 
deallocate(frac) 
deallocate(frad) 
deallocate(state1) 
deallocate(ro_old_keep)  
deallocate(state1step)  
deallocate(state1stepc) 
deallocate(state) 
deallocate(pstate) 
deallocate(qstate) 
deallocate(temp)  
deallocate(p)  
deallocate(qsimu)
deallocate(qsimu0)  
deallocate(qsimu0_dum) 
deallocate(qsimu1_dum) 
deallocate(istore)  
deallocate(qsimf)  
deallocate(qsimf0) 
deallocate(qsimf0_dum) 
deallocate(qsimf1_dum)
deallocate(qsimf1_ens_mean)
deallocate(qsimf0_ens_mean)
deallocate(qsimf1_ens_crps)
deallocate(qsimf0_ens_crps)
deallocate(qsimf1_opt)
deallocate(p2state_opt) 
deallocate(p2state)  
deallocate(zobs)  
deallocate(zobs1)  
deallocate(zobsp)  
deallocate(qtrace)
deallocate(qsimf0_ens) 
deallocate(qsimf1_ens) 
deallocate(ybar)
deallocate(ymean)
deallocate(ybar1)
deallocate(ymean1)
deallocate(yb)
deallocate(ym)
deallocate(sigf)
deallocate(sigfs)
deallocate(rr)
deallocate(sigfi)
deallocate(sigu)
deallocate(ri)
deallocate(sigu1)
deallocate(htri)
deallocate(sigu1a)
deallocate(sigu_s)
deallocate(sigfid)
deallocate(sigfid2)
deallocate(htrih)
deallocate(temp1)
deallocate(temp4)
deallocate(pp)
deallocate(zvar)
deallocate(y)
deallocate(y_save)
deallocate(gaint)
deallocate(h)
deallocate(ht)
deallocate(yy)
deallocate(pxv0)
deallocate(edmnd0)
deallocate(itime)
deallocate(qobs0)
deallocate(c1)
deallocate(hhat1t)
deallocate(w1)
deallocate(amatrixi)
deallocate(aic1tri)
deallocate(aic1trih)
deallocate(aic1trihai)
deallocate(hamatrixi)
deallocate(tempmn)

  enddo !do iii loop
 
  call system_clock(icount2, icount_rate2, icount_max2)
!
! Calculate time elapsed in sec for EnKF.
!
  if(iaenkf.eq.0) then
  enkf_time=float(icount2)/icount_rate1&
       -float(icount1)/icount_rate2
  write(6,*) '  enkf_time ',enkf_time
  endif
!
! Calculate time elapsed in sec for AEnKF.
!
  if(iaenkf.eq.1) then
  aenkf_time=float(icount2)/icount_rate1&
         -float(icount1)/icount_rate2
  write(6,*) '  aenkf_time ',aenkf_time
  endif

  enddo !do iaenkf loop
!
! Calculate the ratio of elapsed time for AEnKF to EnKF.
!
  write(6,*) 'AEnKF-to-EnKF computing time ratio ',aenkf_time/&
  enkf_time

  stop
  end program AEnKF
!**********************************************************************
  function gasdev(iseed,gset)
!**********************************************************************
! function generates standard normal random deviates
! from Numerical Recipes (Press et al. 1998)
!
  dimension randnu(1)
  data iset/0/
  if (iset.eq.0) then
1     continue
    call randd(iseed,1,randnu(1))
    v1=2.*randnu(1)-1.
    call randd(iseed,1,randnu(1))
    v2=2.*randnu(1)-1.
    r=v1**2+v2**2
    if(r.ge.1.)go to 1
    fac=sqrt(-2.*log(r)/r)
    gset=v1*fac
    gasdev=v2*fac
    iset=1
  else
    gasdev=gset
    iset=0
  endif
  return
end function gasdev
!**********************************************************************
    subroutine randd(seed,n,vector)
!**********************************************************************
! this random number generator generates random numbers in ]0,1[
! note that if the seed value is zero on the first call, a default
! value of 1369 will be used  in a linear congruential generator to
! generate 55 odd integers for the array 'itab()'. these values are
! preserved by a common statement, so that they may be used in sub-
! sequent calls by setting the seed to zero.if the value of 'seed'
! is greater than zero in a call to the subroutine, then the array
! 'itab' will be initialized and a new seed value will be returned
! by the subroutine. best results are obtained by making the initial
! call with a seed of your choice and then setting the seed to '0'
! for all subsequent calls.
! from Deutsch and Journel (1992)
!
  dimension vector(*)
  common /unusual/itab(55),n1,n2,nseed
  integer rn1,seed
!
! test to see if 55 odd integers must be generated.
!
  if((seed.gt.0).or.(nseed.lt.1)) then
    nseed = seed
    if(seed.le.0) nseed  = 7931
    do i=1,55
      rn1=mod(nseed*9069,32768)
      if(mod(rn1,2).eq.0) rn1 = rn1-1
      itab(i) = rn1
      nseed = rn1
    enddo
    n1 = 0
    n2 = 24
  endif
!
! generate "n" random components for the vector "vector"
!
  do i=1,n
    itab(55-n1) = mod(itab(55-n2)*itab(55-n1),32768)
    vector(i) = abs(float(itab(55-n1))/float(32768))
    n1 = mod(n1+1,55)
    n2 = mod(n2+1,55)
  enddo
  if(seed.gt.0) seed=nseed

  return
  end subroutine randd
!**********************************************************************
  subroutine diftime(itime1,itime2,lag)
!**********************************************************************
! subroutine computes the time difference in hours between the current 
! (itime1) and the previous time (itime2)
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  if(itime1.lt.0) then
!
! this is the very first calculation
!
  lag=0
  return
  endif

  if(itime2.lt.itime1) then
  write(6,*) 'diftime - itime2 lt itime1...stop ',itime1,itime2
  stop
  endif

  call time1(itime1,iyr1,mon1,iday1,ihr1)
  call time1(itime2,iyr2,mon2,iday2,ihr2)

  jday1=julian(iyr1,mon1,iday1)
  jhr1=jday1*24+ihr1
  jday2=julian(iyr2,mon2,iday2)
  jhr2=jday2*24+ihr2
  if(iyr1.eq.iyr2) then
  lag=jhr2-jhr1
  else
  jday=julian(iyr1,12,31)
  jhr=jday*24
  lag=jhr2+jhr-jhr1
  endif

  return
  end subroutine diftime
!**********************************************************************
  function julian(iyr,mon,iday)
!**********************************************************************
! given year, month, and day, function returns Julian day
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension isumn(12),isuml(12)
  data isumn/0,31,59,90,120,151,181,212,243,273,304,334/
  data isuml/0,31,60,91,121,152,182,213,244,274,305,335/
  yr=iyr
  if(iyr/4-yr/4..eq.0.) then
  julian=isuml(mon)+iday
  else
  julian=isumn(mon)+iday
  endif

  return
  end function julian
!**********************************************************************
  subroutine time1(itime,iyr,mon,iday,ihr)
!**********************************************************************
! given yyyymmddhh, subroutine returns year, month, day and hour
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  iyr=itime/1000000
  mon=itime-iyr*1000000
  mon=mon/10000
  iday=itime-iyr*1000000-mon*10000
  iday=iday/100
  ihr=itime-iyr*1000000-mon*10000-iday*100

  return
  end subroutine time1
!**********************************************************************
  subroutine sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,state,&
             knd,qsimu,qsimf,ipert,iseed,gset,kf,ifrq,istore,nf,sum)
!**********************************************************************
! subroutine runs SAC and UH over the time horizon of [kbeg,jend]
! version Aug 08, 2021 by D.-J. Seo at UTA/HWRL
!
use pass0
use pass1
use pass4
use pass5
use pass7
real, dimension (:), allocatable :: ro_old
real :: lztwm,lzfpm,lzfsm,lzpk,lzsk
  parameter(istep=6,ltx=120/istep)
  dimension po(nv),qsimf(ltx+nwin,nn+2)
  dimension state(nn),qsimu(nf),istore(nf)

  tiny=0.00001

allocate(ro_old(nv))

  jend=0
  if(knd.eq.1) then
  open(11,file=path(1:nchp)//segment(1:ncht)//'.statevar_out',&
  status='unknown',access='direct',form='unformatted',recl=24)
  endif
!
! initialize
!
  qsim=-0.1
!
! specify the end of the simulation horizon
!
  if(knd.le.2) jend=kend          !assimilation
  if(knd.eq.3) jend=kend+ltx          !forecast

  open(9,file=path(1:nchp)//'qo_vs_qs.out_'//segment(1:ncht),&
  status='unknown')

  sump=0.
  sume=0.
  sumw=0.
  sumq=0.

  knt=0
  pbiasln=0.
  ebiasln=0.
  wnor=0.

  kf=0

c10 : do k=kbeg,jend

  if(k.eq.kbeg+1) then
!
! save the 1 step-ahead prediction of the states
!
  do j=1,nn
  state1step(j)=state(j)
  enddo
  endif
!
! run SAC-SMA for one time step 
!
  if(knd.eq.1) then
!
! save the SAC states into a direct-access file and specify MAP and 
! MAPE
!
  write(11,rec=k) (state(i),i=1,6)
  pxv=pxv0(k)*pxadj
  edmnd=edmnd0(k)*peadj
  if(nw.gt.2) w=0.

  else
 
    if(k.le.kend) then

    i=k-kbeg+1
    jp=(i-1)/nhrp+1
    je=(i-1)/nhre+1
    if(nw.gt.2) jw=(i-1)/nhrw+1

    pxv=pxv0(k)*state(6+jp)
    edmnd=edmnd0(k)*state(6+nsp+je)
    if(nw.gt.2) w=state(6+nsp+nse+jw)

    sump=sump+(1./zvar(i))*(pxv0(k)**2)*(pxadj-state(6+jp))**2
    sume=sume+&
    (1./zvar(nwin+i))*(edmnd0(k)**2)*(peadj-state(6+nsp+je))**2
    if(nw.gt.2) then
    sumw=sumw+&
    (1./zvar(2*nwin+i))*(1.**2)*(0.-state(6+nsp+nse+jw))**2
    endif

    else
 
  knt=knt+1

  if(knd.eq.2) then
  write(6,*) 'k gt kend and knd eq 2 in subroutine sacuh...stop'
  stop
  endif

  if(ipert.eq.1) then
  call get_perturbs_for_pxpetci(iseed,pxadj,peadj,per,eer,wer,&
                                pbiasln,ebiasln,wnor,knt,gset)
  pxv=pxv0(k)*per
  edmnd=edmnd0(k)*eer
  if(nw.gt.2) w=0.+wer
  else
  pxv=pxv0(k)*pxadj
  edmnd=edmnd0(k)*peadj
  if(nw.gt.2) w=0.
  endif

  if(k-kbeg+1.gt.ltx+nwin) then
  write(6,*)  'qsimf out of range in subroutine sacuh'
  stop
  endif

  qsimf(k-kbeg+1,1)=-0.1

  endif

  endif
!
! specify the observed flow and data quality flag
!
  qobs=qobs0(k)
!
! run SAC for one timestep
! 
  call sac(state,edmnd,uztwm,uzfwm,lztwm,lzfpm,lzfsm,saved,adimp,&
       pxv,pctim,uzk,lzpk,lzsk,dt,zperc,rexp,pfree,side,riva,parea,&
       tci,nn)
!
! adjust TCI
!
  if(nw.gt.2) then
  tci=tci+w
  if(tci.lt.0.) tci=tiny
  endif
!
! SAC-SMA has been run for one time step: update the runoff entries in
! the UH convolution
!
  do it=1,(nv-1)
  i=(nv-1)-it+1
  ro_old(i+1)=ro_old(i)
  enddo
  ro_old(1)=tci
!
! check if all runoff entries are filled before the convolution
! operation
!
  null=0
  do i=1,nv
  if(ro_old(i).lt.0.) null=null+1
  enddo
!
! calculate flow in cms via UH
!
  if(null.eq.0) then
!
! perform convolution
!
  rox=0.
  q=0.
  do i=1,nv
  q=q+ro_old(i)*po(i)
  if(ro_old(i).gt.rox) rox=ro_old(i)
  enddo

  if(knd.eq.1) write(9,23)itime(k),qobs,q,pxv,edmnd,tci
   23 format(i15,3f12.2,2f9.5)
!
! evaluate the objective function value
!
  if(k.gt.kend-nf*ifrq.and.k.le.kend) then
   
  if(mod(k-kend,ifrq).eq.0) then
  kf=kf+1
  istore(kf)=k-kend
  qsimu(kf)=q
  if(knd.eq.3) sumq=sumq+(q-qobs)**2/zvar(nw*nwin+kf)
  endif

  endif

  if(knd.eq.3) qsimf(k-kbeg+1,1)=q

  endif

  enddo c10

  if(knd.eq.3) then
  if(nw.eq.3) sum=sump+sume+sumw+sumq
  if(nw.eq.2) sum=sump+sume+sumq
!     sum=sumq
  endif

  close(9)

  if(knd.eq.1) close(11)
  if(knd.eq.2) close(24)

  if(kf.ne.nf) then
  write(6,*) 'kf ne nf in subroutine sacuh...stop ',kf,nf
  stop
  endif

deallocate(ro_old)

  return
  end subroutine sacuh
!**********************************************************************
  subroutine get_perturbs_for_pxpetci(iseed,pxadj,peadj,per,eer,&
             wer,pbiasln,ebiasln,wnor,knt,gset)
!**********************************************************************
! subroutine perturbs multiplicative biases in MAP aand MAPE, pxadj and
! peadj, respectively, and addtive error to TCI
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension vary(3),p(3),rho(3)

  vary(1)=0.1 
  vary(2)=0.01 
  vary(3)=0.001 
!
! specify lag-1 correlation for the AR(1) model for pxadj, peadj, and
! additive error to TCI
!
  rho(1)=0.96
  rho(2)=0.96
  rho(3)=0.96
!
! assume pxadj and peadj ~ LN(avey,vary)
!
  avey=0.
  do i=1,2
  if(i.eq.1) avey=pxadj
  if(i.eq.2) avey=peadj
  varx=-2.*alog(avey)+alog(vary(i)+avey**2)
  avex=2.*alog(avey)-0.5*alog(vary(i)+avey**2)
!
! ln(pxadj) and ln(peadj) ~ N(avex,varx)
!
  if(knt.eq.1) then
  p(i)=sqrt(varx)*gasdev(iseed,gset)+avex
  else
  p(i)=avex+rho(i)*(pbiasln-avex)&
      +sqrt(1.-rho(i)**2)*sqrt(varx)*gasdev(iseed,gset)
  endif
  if(i.eq.1) pbiasln=p(i)
  if(i.eq.2) ebiasln=p(i)
  p(i)=exp(p(i))
  enddo
!
! perturb additive error to TCI, w, w ~ N(0,vary(3))
!
  if(knt.eq.1) then
  p(3)=sqrt(vary(3))*gasdev(iseed,gset)+0.
  else
  p(3)=rho(3)*wnor&
  +sqrt(1.-rho(3)**2)*sqrt(vary(3))*gasdev(iseed,gset)
  endif
  wnor=p(3)
!
! check if they fall within p hysically-realistic range
!
  smin=0.2
  smax=5.
  do j=1,2
  if(p(j).lt.smin) p(j)=smin
  if(p(j).gt.smax) p(j)=smax
  enddo

  per=p(1)
  eer=p(2)
  wer=p(3)

  return
  end subroutine get_perturbs_for_pxpetci
!**********************************************************************
  subroutine hpsort(nn,zd)
!**********************************************************************
! subroutine performs heap-sorting of zd in the ascending order
! from Press et al. (1998)
!
  dimension zd(nn)
  l=nn/2+1
  ir=nn
   15 if(l.gt.1) go to 10
  go to 25
   10 l=l-1
  rra=zd(l)
  go to 30
   25 rra=zd(ir)
  zd(ir)=zd(1)
  ir=ir-1
  if(ir.eq.1) go to 40
  go to 30
   40 zd(1)=rra
  go to 5
   30 i=l
  j=l+l
   20 if(j.le.ir) go to 50
  go to 60
   50 if(j.lt.ir) go to 70
  go to 80
   70 if(zd(j).lt.zd(j+1)) j=j+1
   80 if(rra.lt.zd(j)) go to 90
  go to 100
   90 zd(i)=zd(j)
  i=j
  j=j+j
  go to 20
  100 j=ir+1
  go to 20
   60 zd(i)=rra
  go to 15
    5 return
  end subroutine hpsort
!**********************************************************************
  subroutine check_range(itype,state,uztwm,uzfwm,lztwm,&
             lzfsm,lzfpm,nsp,nse,nsw,nn,nw)
!**********************************************************************
! subroutine checks feasible range of state variables
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
real :: lztwm,lzfpm,lzfsm
  dimension state(nn)

  if(state(1).lt.  0. ) state(1)=0.
  if(state(1).gt.uztwm) state(1)=uztwm

  if(state(2).lt.  0. ) state(2)=0.
  if(state(2).gt.uzfwm) state(2)=uzfwm

  if(state(3).lt.  0. ) state(3)=0.
  if(state(3).gt.lztwm) state(3)=lztwm

  if(state(4).lt.  0. ) state(4)=0.
  if(state(4).gt.lzfsm) state(4)=lzfsm

  if(state(5).lt.  0. ) state(5)=0.
  if(state(5).gt.lzfpm) state(5)=lzfpm

  if (state(6).lt.state(1)) state(6)=state(1)

  if (state(6).gt.(uztwm+lztwm)) state(6)=uztwm+lztwm

  if(itype.eq.1) return

  smin=0.1
  smax=10.
  do j=1,nsp
  if(state(6+j).lt.smin) state(6+j)=smin
  if(state(6+j).gt.smax) state(6+j)=smax
  enddo

  do j=1,nse
  if(state(6+nsp+j).lt.smin) state(6+nsp+j)=smin
  if(state(6+nsp+j).gt.smax) state(6+nsp+j)=smax
  enddo

  if(nw.gt.2) then
  smin=-10.
  smax= 10.
  do j=1,nsw
  if(state(6+nsp+nse+j).lt.smin) state(6+nsp+nse+j)=smin
  if(state(6+nsp+nse+j).gt.smax) state(6+nsp+nse+j)=smax
  enddo
  endif

  return
  end subroutine check_range
!**********************************************************************
  subroutine get_perturbs_for_obs(nm,iseed,pp,zvar,gset)
!**********************************************************************
! subroutine returns perturbation for observation
! A. Rafieeinasab Jul 24, 2016 at UTA/HWRL
!
  dimension pp(nm),zvar(nm)
 
  do i=1,nm
  pp(i)=sqrt(zvar(i))*gasdev(iseed,gset)
  enddo

  return
  end subroutine get_perturbs_for_obs
!*********************************************************************
  subroutine get_perturbed_aug_states0(uztwm,uzfwm,lztwm,lzfsm,&
             lzfpm,nsp,nse,nsw,state1,iseed,p,gset,frac,pstate,&
             uzk,pctim,adimp,riva,zperc,rexp,lzsk,lzpk,&
             pfree,side,saved,parea,dt,nv,po,pxadj,peadj,knd,&
             y,kf,ifrq,nf)
!*********************************************************************
! subroutine perturbs states for cold DA
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
use pass7
real :: lztwm,lzfpm,lzfsm,lzpk,lzsk
integer, dimension (:), allocatable :: istore
real, dimension (:), allocatable :: qsimu0_dum
real, dimension (:,:), allocatable :: qsimf0_dum

  parameter (istep=6,ltx=120/istep)
  dimension y(m,ns),po(nv),state1(nn),pstate(nn),frac(nn),p(nn)

allocate(istore(nf))
allocate(qsimu0_dum(nf))
allocate(qsimf0_dum(ltx+nwin,nn+2))
!
! one step-predict from time step k-1 to time step k
!
  do ii=1,ns
!
! perturb the states at the (k-1)st timestep, state1, before
! propagating them to the k-th timestep
!
  call get_perturbs_for_states0(uztwm,uzfwm,lztwm,lzfsm,&
       lzfpm,nsp,nse,&
       nsw,state1,iseed,p,gset,frac,nn,nw)
!
! get perturbed states (no need to range-check; checked in
! subroutine state_perturbations)
!
  do j=1,nn
  pstate(j)=state1(j)+p(j)
  enddo
!
! propagate the perturbed states at the (k-1)st timestep, pstate,
! one-step forward to the k-th timestep to get M(Xk-1+pi)
!
  call check_range(2,pstate,uztwm,uzfwm,lztwm,lzfsm,&
                   lzfpm,nsp,nse,nsw,nn,nw)

  knd=2
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,pstate,&
             knd,qsimu0_dum,qsimf0_dum,0,iseed,gset,kf,ifrq,&
             istore,nf,sum)
!
! save the one step-propagated state variables for the next time step
!
  do j=1,nn
  y(j,ii)=state1step(j)
  enddo

  enddo

deallocate(istore)
deallocate(qsimu0_dum)
deallocate(qsimf0_dum)

  return
  end subroutine get_perturbed_aug_states0
!**********************************************************************
  subroutine get_perturbed_aug_states1(uztwm,uzfwm,lztwm,lzfsm,&
             lzfpm,nsp,nse,nsw,iseed,gset,frad,p2state,&
             uzk,pctim,adimp,riva,zperc,rexp,lzsk,lzpk,&
             pfree,side,saved,parea,dt,nv,po,pxadj,peadj,knd,&
             y,kf,ifrq,nf)
!**********************************************************************
! subroutine perturb states for warm DA
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
use pass7
integer, dimension (:), allocatable :: istore
real, dimension (:), allocatable :: qstate2
real, dimension (:), allocatable :: err
real, dimension (:), allocatable :: tempp
real, dimension (:), allocatable :: qsimu0_dum
real, dimension (:,:), allocatable :: qsimf0_dum
real :: lztwm,lzfpm,lzfsm,lzpk,lzsk

  parameter (istep=6,ltx=120/istep)
  dimension y(m,ns),po(nv),p2state(nn,ns),frad(nn)
allocate(qstate2(nn))
allocate(err(nn))
allocate(tempp(nn))
allocate(istore(nf))
allocate(qsimu0_dum(nf))
allocate(qsimf0_dum(ltx+nwin,nn+2))
!
! ensemble-forecast from time step k to time step k+1
!
  do ii=1,ns 
!
! propagate the perturbed states at the k-th timestep, p, one-step
! forward to the (k+1)st timestep
!
  do j1=1,nn
  qstate2(j1)=p2state(j1,ii)
  enddo
!
! check range of the state variables
!
  call check_range(2,qstate2,uztwm,uzfwm,lztwm,lzfsm,&
                   lzfpm,nsp,nse,nsw,nn,nw)

  knd=2
  call sacuh(uztwm,uzfwm,uzk,pctim,adimp,&
             riva,zperc,rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
             side,saved,parea,dt,nv,po,pxadj,peadj,qstate2,&
             knd,qsimu0_dum,qsimf0_dum,0,iseed,gset,kf,ifrq,&
             istore,nf,sum)
!
! save the forecast states valid at the (k+1)st time step
!
  do j=1,nn
  y(j,ii)=state1step(j)
  enddo
!
! generate model error; model error is a function of state variables
!
  call get_perturbs_for_states1(nn,err,iseed,state1step,nsp,nse,&
                                gset,frad)
!
! add model error
!
  do j=1,nn
  tempp(j)=y(j,ii)+err(j)
  enddo
!
! check feasible range
!
  call check_range(2,tempp,uztwm,uzfwm,lztwm,lzfsm,&
                   lzfpm,nsp,nse,nsw,nn,nw)

  do j=1,nn
  y(j,ii)=tempp(j)
  enddo

  enddo 

deallocate(qstate2)
deallocate(err)
deallocate(tempp)
deallocate(istore)
deallocate(qsimu0_dum)
deallocate(qsimf0_dum)

  return
  end subroutine get_perturbed_aug_states1
!**********************************************************************
  subroutine get_perturbs_for_states1(nn,err,iseed,state1step,&
             nsp,nse,gset,frad)
!**********************************************************************
! subroutine generate model errors
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension err(nn),state1step(nn),frad(nn)
!
! for SAC states
!
  do i=1,6
  err(i)=frad(i)*state1step(i)
  rnormal=gasdev(iseed,gset)
  err(i)=sqrt(2.)*err(i)*rnormal
  enddo
!
! for bias for MAP
!
  do i=6+1,6+nsp
  rnormal=gasdev(iseed,gset)
  err(i)=(frad(i)**1)*rnormal
  enddo
!
! for bias for MAPE
!
  do i=6+nsp+1,6+nsp+nse
  rnormal=gasdev(iseed,gset)
  err(i)=(frad(i)**1)*rnormal
  enddo
!
! for additive error to TCI
!
  do i=6+nsp+nse+1,nn
  rnormal=gasdev(iseed,gset)
  err(i)=(frad(i)**1)*rnormal
  enddo

  return
  end subroutine get_perturbs_for_states1
!**********************************************************************
  subroutine get_perturbs_for_states0(uztwm,uzfwm,lztwm,&
             lzfsm,lzfpm,nsp,nse,nsw,state3,iseed,err,gset,frac,&
             nn,nw)
!**********************************************************************
! subroutine perturbs IC with normal errors
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
real :: lztwm,lzfpm,lzfsm
  dimension err(nn),state3(nn),frac(nn)
!
! for SAC states
!
  do i=1,6
  err(i)=frac(i)*state3(i)
  rnormal=gasdev(iseed,gset)
  err(i)=sqrt(2.)*err(i)*rnormal
  enddo
!
! for bias for MAP
!
  do i=6+1,6+nsp
  rnormal=gasdev(iseed,gset)
  err(i)=(frac(i)**1)*rnormal
  enddo
!
! for bias for MAPE
!
  do i=6+nsp+1,6+nsp+nse
  rnormal=gasdev(iseed,gset)
  err(i)=(frac(i)**1)*rnormal
  enddo
!
! for additive error to TCI
!
  do i=6+nsp+nse+1,nn
  rnormal=gasdev(iseed,gset)
  err(i)=(frac(i)**1)*rnormal
  enddo

  do i=1,nn
  err(i)=state3(i)+err(i)
  enddo
!
! check for feasible range
!
  call check_range(2,err,uztwm,uzfwm,lztwm,lzfsm,&
                   lzfpm,nsp,nse,nsw,nn,nw)
!
! return only the perturbations
!
  do i=1,nn
  err(i)=err(i)-state3(i)
  enddo

  return
  end subroutine get_perturbs_for_states0
!**********************************************************************
  subroutine lsolve_cd(n,as,m,bbs,xxs,ierr)
!**********************************************************************
! solves linear system via Cholesky decomposition
! from Press et al. (1998)
!
  dimension as(n,n),bbs(n,m),xxs(n,m)

real*8, dimension (:,:), allocatable :: a,xx,bb
real*8, dimension (:), allocatable :: b,p,x

allocate(a(n,n))
allocate(b(n))
allocate(p(n))
allocate(x(n))
allocate(bb(n,m))
allocate(xx(n,m))

  do i=1,n
  do j=1,n
  a(i,j)=dble(as(i,j))
  enddo
  enddo

  do i=1,n
  do j=1,m
  bb(i,j)=dble(bbs(i,j))
  enddo
  enddo
!
! feed the upper matrix only
!
  call choldc(a,n,p,ierr)

  if(ierr.ne.0) return
!
! solve the linear system column by column for inversion
!
  do j=1,m
!
! copy entries in the j-th column of the identity matrix
!
  do i=1,n
  b(i)=bb(i,j)
  enddo
!
! solve the linear system via Cholesky decomposition
!
  call cholsl(a,n,p,b,x)
!
! copy the column by column solution to matrix X
!
  do i=1,n
  xx(i,j)=x(i)
  enddo

  enddo

  do i=1,n
  do j=1,m
  xxs(i,j)=sngl(xx(i,j))
  enddo
  enddo

deallocate(a)
deallocate(b)
deallocate(p)
deallocate(x)
deallocate(bb)
deallocate(xx)

  return
  end subroutine lsolve_cd
!**********************************************************************
  SUBROUTINE choldc(a,n,p,ierr)
!**********************************************************************
! Given a positive-definite symmetric matrix a(1:n,1:n), with physical
! dimension np, this routine constructs its Cholesky decomposition,
! A = L·LT . On input, only the upper triangle of a need be given; 
! it is not modified. The Cholesky factor L is returned in the lower 
! triangle of a, except for its diagonal elements which are returned 
! in p(1:n).
! from Press et al. (1998)
!
  REAL*8 a(n,n),p(n)
  REAL*8 sum

  ierr=0
  do i=1,n
  do j=i,n
  sum=a(i,j)
  do k=i-1,1,-1
  sum=sum-a(i,k)*a(j,k)
  enddo
  if(i.eq.j)then
  if(sum.le.0.d0) then
  write(6,*) 'choldc failed ',sum
  ierr=1
  return
  endif
  p(i)=dsqrt(sum)
  else
  a(j,i)=sum/p(i)
  endif
  enddo
  enddo

  return
  END SUBROUTINE choldc
!**********************************************************************
  SUBROUTINE cholsl(a,n,p,b,x)
!**********************************************************************
! Solves the set of n linear equations A · x = b, where a is a 
! positive-definite symmetric matrix with physical dimension np. 
! a and p are input as the output of the routine choldc.
! Only the lower triangle of a is accessed. b(1:n) is input as the
! right-hand side vector. The solution vector is returned in x(1:n).
! a, n, np, and p are not modified and can be left in place for 
! successive calls with different right-hand sides b. b is not modified
! unless you identify b and x in the calling sequence, which is 
! allowed.
! from Press et al. (1998)
!
  REAL*8 a(n,n),p(n),b(n),x(n)
  REAL*8 sum

  do i=1,n
  sum=b(i)
  do k=i-1,1,-1
  sum=sum-a(i,k)*x(k)
  enddo
  x(i)=sum/p(i)
  enddo
  do i=n,1,-1
  sum=x(i)
  do k=i+1,n
  sum=sum-a(k,i)*x(k)
  enddo
  x(i)=sum/p(i)
  enddo

  return
  END SUBROUTINE cholsl
!**********************************************************************
  subroutine calc_crps(n,fcst,obs,crps)
!**********************************************************************
! subroutine calculates mean CRPS based on Hersbach (2000)
! Version Aug 08, 2021 by D.-J. Seo at UTA/HWRL
!
  dimension fcst(0:n-1),alpha(0:n),beta(0:n)

  do j=0,n
  alpha(j)=0.
  beta(j)=0.
  enddo

     if(obs.le.fcst(0)) then
!
! verifying obs is less than fcst(0)
!
     beta(0)=beta(0)+(fcst(0)-obs)

     do k=1,n-1
     beta(k)=beta(k)+(fcst(k)-fcst(k-1))
     enddo

     else if(obs.ge.fcst(n-1)) then
!
! verifying obs is greater than fcst(n-1)
!
     do k=1,n-1
     alpha(k)=alpha(k)+(fcst(k)-fcst(k-1))
     enddo

     alpha(n)=alpha(n)+(obs-fcst(n-1))

     else
!
! verifying observation falls between fcst(0) and fcst(n-1)
!
       klow=0
       do k=1,n-1
       if(obs.ge.fcst(k-1).and.obs.le.fcst(k)) then
       klow=k
       exit
       endif
       enddo
!
! verifying obs is between fcst(klow-1) and fcst(klow)
!
     if(klow.ge.2) then
     do k=1,klow-1
     alpha(k)=alpha(k)+(fcst(k)-fcst(k-1))
     enddo
     endif

     alpha(klow)=alpha(klow)+(obs-fcst(klow-1))
     beta(klow)=beta(klow)+(fcst(klow)-obs)

     if(klow+1.le.n-1) then
     do k=klow+1,n-1
     beta(k)=beta(k)+(fcst(k)-fcst(k-1))
     enddo
     endif

  endif

  crps=0.
  do k=0,n
  if(k.eq.0) then
  crps=crps+beta(k)
  else if(k.eq.n) then
  crps=crps+alpha(k)
  else
  prob=float(k)/(n+1.)
  crps=crps+alpha(k)*prob**2+beta(k)*(1.-prob)**2
  endif
  enddo

  return
  end subroutine calc_crps
!**********************************************************************
  subroutine get_map(nch,infile,jdt,k,jtime,jdata1)
!**********************************************************************
! subroutine reads MAP data
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension jtime(jdata1)
  character infile*132,string*132,site1*5,cdum2*8

  open(10,file=infile(1:nch),status='old')

   10 read(10,40) string
   40 format(a132)
  if(string(1:1).eq.'$') go to 10
  if(string(15:18).ne.'MAPX') then
  endif
  if(string(25:26).eq.'MM') then
  fac1=1
  else if(string(25:26).eq.'IN') then
  fac1=25.4
  else
  write(6,*) 'unknown unit for MAPX...stop ',string(25:26)
  stop
  endif
!
! get the timestep
!
  read(string(29:30),'(i2)') jdt

  site1=string(35:39)
  read(10,40) string
!
! read MAP data
!
  k=0
  jday_old=-999
   20 read(10,50,end=30) cdum2,jmon,jyr,jday,pxv
   50 format(a8,4x,i2,i2,2x,i2,f10.3)
  k=k+1
  if(jday.ne.jday_old) knt2=0
  knt2=knt2+1
  jday_old=jday
  if(jyr.gt.21) then
  jyr=jyr+1900
  else
  jyr=jyr+2000
  endif
  jtime(k)=jyr*1000000+jmon*10000+jday*100+knt2*6
  go to 20
   30 close(10)

  return
  end subroutine get_map
!**********************************************************************
  subroutine get_qin(nch,infile,k,ktime,kdata1)
!**********************************************************************
! subroutine reads streamflow data
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension ktime(kdata1)
  character infile*132,string*132,usgsid*8,cdum1*8

  open(10,file=infile(1:nch),status='old')

   10 read(10,40) string
   40 format(a132)
  if(string(1:1).eq.'$') go to 10
  if(string(15:17).ne.'QIN') then
  write(6,*) 'not QIN time series...stop'
  stop
  endif
  if(string(25:27).eq.'CMS') then
  fac2=1
  else if(string(25:27).eq.'CFS') then
  write(6,*) 'in CFS'
  fac2=0.3048**3
  else
  write(6,*) 'unknown unit for QIN...stop ',string(25:27)
  stop
  endif
  if(string(30:30).ne.'1') then
  write(6,*) 'not hourly data...stop'
  stop
  endif
  usgsid=string(35:42)
  read(10,40) string
!
! read streamflow data
!
  k=0
  iday_old=-999
   20 read(10,50,end=30) cdum1,imon,iyr,iday,qobs
   50 format(a8,4x,i2,i2,2x,i2,f10.3)
  k=k+1
  if(iday.ne.iday_old) knt2=0
  knt2=knt2+1
  iday_old=iday
  if(iyr.gt.21) then
  iyr=iyr+1900
  else
  iyr=iyr+2000
  endif
  ktime(k)=iyr*1000000+imon*10000+iday*100+knt2*6
  go to 20
   30 close(10)

  return
  end subroutine get_qin
!**********************************************************************
  subroutine get_period_of_record(ndata,jtime,kdata,ktime,jbeg,&
                                  kbeg,jend,kend)
!**********************************************************************
! subroutine finds the beginning and ending times of the common period
! Version Aug 08, 2021 by D.-J. Seo at UTA/HWRL
!
  dimension jtime(ndata),ktime(kdata)
!
! get the common beginning date
!
  if(jtime(1).gt.ktime(1)) then
  do i=2,kdata
  if(jtime(1).eq.ktime(i)) go to 10
  enddo
  write(6,*) 'no matching beginning date...stop'
  stop
   10 jbeg=1
  kbeg=i
  else if(jtime(1).lt.ktime(1)) then
  do i=2,ndata
  if(jtime(i).eq.ktime(1)) go to 20
  enddo
  write(6,*) 'no matching beginning date...stop'
  stop
   20 jbeg=i
  kbeg=1
  else
  jbeg=1
  kbeg=1
  endif

  if(jtime(jbeg).ne.ktime(kbeg)) then
  write(6,*) 'jtime(jbeg).ne.ktime(kbeg)...stop ',jbeg,&
      jtime(jbeg),kbeg,ktime(kbeg)
  stop
  endif
!
! get the common ending date
!
  if(jtime(ndata).lt.ktime(kdata)) then
  do i=kdata-1,1,-1
  if(jtime(ndata).eq.ktime(i)) go to 30
  enddo
  write(6,*) 'no matching ending date...stop'
  stop
   30 jend=ndata
  kend=i
  else if(jtime(ndata).gt.ktime(kdata)) then
  do i=ndata-1,1,-1
  if(jtime(i).eq.ktime(kdata)) go to 40
  enddo
  write(6,*) 'no matching ending date...stop'
  stop
   40 jend=i
  kend=kdata
  else
  jend=ndata
  kend=kdata
  endif

  if(jtime(jbeg).ne.ktime(kbeg)) then
  write(6,*) 'jtime(jend).ne.ktime(kend)...stop ',jend,&
      jtime(jend),kend,ktime(kend)
  stop
  endif

  return
  end subroutine get_period_of_record
!**********************************************************************
  subroutine sac(state,edmnd,uztwm,uzfwm,lztwm,lzfpm,lzfsm,saved,&
             adimp,pxv,pctim,uzk,lzpk,lzsk,dt,zperc,rexp,pfree,&
             side,riva,parea,tci,nn)
!**********************************************************************
! subroutine executes the Sacramento model for one timestep
! from Eric Anderson (1979) of NWS/HRL
! version Aug 08, 2021, adapted by D.-J. Seo at UTA/HWRL
!
real :: lztwm,lzfpm,lzfsm,lzpk,lzsk
  dimension state(nn)

  simpvt=0.
  sintft=0.
  sgwfp=0.
  sgwfs=0.
  srecht=0.
  srost=0.
  srodt=0.
  srot=0.
  sett=0.
  se1=0.
  se3=0.
  se4=0.
  se5=0.
!
! SAC starts here
!
!     compute et from upper zone.
  e1=edmnd*(state(1)/uztwm)
  red=edmnd-e1
!     red is residual evap demand
  state(1)=state(1)-e1
  e2=0.0
  if(state(1).lt.0.) then
!     e1 can not exceed state(1)
  e1=e1+state(1)
  state(1)=0.0
  red=edmnd-e1
    if(state(2).lt.red) then
!   e2 is evap from state(2).
    e2=state(2)
    state(2)=0.0
    red=red-e2
    else
    e2=red
    state(2)=state(2)-e2
    red=0.0
      if((state(1)/uztwm).lt.(state(2)/uzfwm)) then
!
!     upper zone free water ratio exceeds upper zone
!     tension water ratio, thus transfer free water to tension
!
      uzrat=(state(1)+state(2))/(uztwm+uzfwm)
      state(1)=uztwm*uzrat
      state(2)=uzfwm*uzrat
      endif
   endif
  else
    if((state(1)/uztwm).lt.(state(2)/uzfwm)) then
!   upper zone free water ratio exceeds upper zone
!   tension water ratio, thus transfer free water to tension
    uzrat=(state(1)+state(2))/(uztwm+uzfwm)
    state(1)=uztwm*uzrat
    state(2)=uzfwm*uzrat
    endif
  endif
  if (state(1).lt.0.00001) then
  state(1)=0.0
  endif
  if (state(2).lt.0.00001) then
  state(2)=0.0
  endif
!
!     compute et from the lower zone.
!     compute et from state(3) (e3)
  e3=red*(state(3)/(uztwm+lztwm))
  state(3)=state(3)-e3
  if(state(3).lt.0.0) then
!     e3 can not exceed state(3)
  e3=e3+state(3)
  state(3)=0.0
  endif
  ratlzt=state(3)/lztwm
  ratlz=(state(3)+state(5)+state(4)-saved)/(lztwm+lzfpm+lzfsm-saved)
  if(ratlzt.lt.ratlz) then
!     resupply lower zone tension water from lower
!     zone free water if more water available there.
  del=(ratlz-ratlzt)*lztwm
!     transfer from state(4) to state(3).
  state(3)=state(3)+del
  state(4)=state(4)-del
    if(state(4).lt.0.0) then
!   if transfer exceeds state(4) then remainder comes from state(5)
    state(5)=state(5)+state(4)
    state(4)=0.0
    endif
  endif
  if (state(3).lt.0.00001) then
  state(3)=0.0
  endif
!
!     compute et from adimp area.-e5
  e5=e1+(red+e2)*((state(6)-e1-state(1))/(uztwm+lztwm))
!     adjust state(6),additional impervious area storage, for
!     evaporation
  state(6)=state(6)-e5
  if(state(6).lt.0.0) then
!     e5 can not exceed state(6).
  e5=e5+state(6)
  state(6)=0.0
  endif
  e5=e5*adimp
!     e5 is et from the area adimp.
!.......................................
!     compute percolation and runoff amounts.
  twx=pxv+state(1)-uztwm
!     twx is the time interval available moisture in excess
!     of uztw requirements.
  if(twx.lt.0.0) then
!     all moisture held in uztw--no excess.
  state(1)=state(1)+pxv
  twx=0.0
  else
!     moisture available in excess of uztw storage.
  state(1)=uztwm
  endif
  state(6)=state(6)+pxv-twx
!
!     compute impervious area runoff.
  roimp=pxv*pctim
!     roimp is runoff from the minimum impervious area.
  simpvt=simpvt+roimp
!
!     initialize time interval sums.
  sbf=0.0
  ssur=0.0
  sif=0.0
  sperc=0.0
  sdro=0.0
  spbf=0.0
!
!     determine computational time increments for the basic time
!     interval
  rninc=1.0+0.2*(state(2)+twx)
  if(rninc.lt.0.) then
  write(6,*) 'rninc lt 0....stop ',rninc
  stop
  endif
  ninc=int(rninc)

! 
!  ninc=1
!

!     ninc=number of time increments that the time interval
!     is divided into for further
!     soil-moisture accounting.  no one increment
!     will exceed 5.0 millimeters of state(2)+pav
  dinc=(1.0/ninc)*dt
!     dinc=length of each increment in days.
  pinc=twx/ninc
!     pinc=amount of available moisture for each increment.
!     compute free water depletion fractions for
!     the time increment being used-basic depletions
!     are for one day
  duz=1.0-((1.0-uzk)**dinc)
  dlzp=1.0-((1.0-lzpk)**dinc)
  dlzs=1.0-((1.0-lzsk)**dinc)
!
!     start incremental do loop for the time interval.
!
  do 240 i=1,ninc
  adsur=0.0
!     compute direct runoff (from adimp area).
  ratio=(state(6)-state(1))/lztwm
  if (ratio.lt.0.0) then
  ratio=0.0
  endif
  addro=pinc*(ratio**2)
!     addro is the amount of direct runoff from the area adimp.
!
!     compute baseflow and keep track of time interval sum.
  bf=state(5)*dlzp
  state(5)=state(5)-bf
  if (state(5).le.0.0001) then
  bf=bf+state(5)
  state(5)=0.0
  endif
  sbf=sbf+bf
  spbf=spbf+bf
  bf=state(4)*dlzs
  state(4)=state(4)-bf
  if(state(4).le.0.0001) then
  bf=bf+state(4)
  state(4)=0.0
  endif
  sbf=sbf+bf
!
!     compute percolation-if no water available then skip
  if((pinc+state(2)).le.0.01) then
  state(2)=state(2)+pinc
  else
  percm=lzfpm*dlzp+lzfsm*dlzs
  perc=percm*(state(2)/uzfwm)
  defr=1.0-((state(3)+state(5)+state(4))/(lztwm+lzfpm+lzfsm))
!     defr is the lower zone moisture deficiency ratio
  fr=1.0
!     fr is the change in percolation withdrawal due to frozen ground
  fi=1.0
!     fi is the change in interflow withdrawal due to frozen ground
  ifrze=0
  if (ifrze.ne.0) then
  uzdefr=1.0-((state(1)+state(2))/(uztwm+uzfwm))
!     call fgfr1(defr,fr,uzdefr,fi)
  endif
  perc=perc*(1.0+zperc*(defr**rexp))*fr
!     note...percolation occurs from state(2) before pav is added.
  if(perc.ge.state(2)) then
!     percolation rate exceeds state(2).
  perc=state(2)
!     percolation rate is less than state(2).
  endif
  state(2)=state(2)-perc
!     check to see if percolation exceeds lower zone deficiency.
  check=state(3)+state(5)+state(4)+perc-lztwm-lzfpm-lzfsm
  if(check.gt.0.0) then
  perc=perc-check
  state(2)=state(2)+check
  endif
  sperc=sperc+perc
!     sperc is the time interval summation of perc
!
!     compute interflow and keep track of time interval sum.
!     note...pinc has not yet been added
  del=state(2)*duz*fi
  sif=sif+del
  state(2)=state(2)-del
!     distribe percolated water into the lower zones
!     tension water must be filled first except for the pfree area
!     perct is percolation to tension water and percf is percolation
!     going to free water.
  perct=perc*(1.0-pfree)
  if ((perct+state(3)).le.lztwm) then
  state(3)=state(3)+perct
  percf=0.0
  else
  percf=perct+state(3)-lztwm
  state(3)=lztwm
!
!  distribute percolation in excess of tension
!  requirements among the free water storages.
  endif
  percf=percf+perc*pfree
  if(percf.ne.0.0) then
  hpl=lzfpm/(lzfpm+lzfsm)
!     hpl is the relative size of the primary storage
!     as compared with total lower zone free water storage.
  ratlp=state(5)/lzfpm
  ratls=state(4)/lzfsm
!     ratlp and ratls are content to capacity ratios, or
!     in other words, the relative fullness of each storage
  fracp=(hpl*2.0*(1.0-ratlp))/((1.0-ratlp)+(1.0-ratls))
!     fracp is the fraction going to primary.
    if (fracp.gt.1.0) then
    fracp=1.0
    endif
  percp=percf*fracp
  percs=percf-percp
!     percp and percs are the amount of the excess
!     percolation going to primary and supplemental
!  storges,respectively.
  state(4)=state(4)+percs
    if(state(4).gt.lzfsm) then
    percs=percs-state(4)+lzfsm
    state(4)=lzfsm
    endif
  state(5)=state(5)+(percf-percs)
!     check to make sure state(5) does not exceed lzfpm.
    if (state(5).gt.lzfpm) then
    excess=state(5)-lzfpm
    state(3)=state(3)+excess
    state(5)=lzfpm
    endif
!
!     distribute pinc between state(2) and surface runoff.
  endif
    if(pinc.ne.0.0) then 
!   check if pinc exceeds uzfwm
      if((pinc+state(2)).le.uzfwm) then
!     no surface runoff
      state(2)=state(2)+pinc
      else
!
!     compute surface runoff (sur) and keep track of time interval
!     sum
      sur=pinc+state(2)-uzfwm
      state(2)=uzfwm
      ssur=ssur+sur*parea
      adsur=sur*(1.0-addro/pinc)
!     adsur is the amount of surface runoff which comes
!     from that portion of adimp which is not
!     currently generating direct runoff.  addro/pinc
!     is the fraction of adimp currently generating
!     direct runoff.
      ssur=ssur+adsur*adimp
!
!     adimp area water balance -- sdro is the idt sum of
!      direct runoff.
      endif
    endif
  endif
  state(6)=state(6)+pinc-addro-adsur
  if (state(6).gt.(uztwm+lztwm)) then
  addro=addro+state(6)-(uztwm+lztwm)
  state(6)=uztwm+lztwm
  endif
  sdro=sdro+addro*adimp
  if (state(6).lt.0.00001) then
  state(6)=0.0
  endif
  240 continue
!
!     end of incremental do loop.
!
!     compute sums and adjust runoff amounts by the area over
!     which they are generated.
  eused=e1+e2+e3
!     eused is the et from parea which is 1.0-adimp-pctim
  sif=sif*parea
!
!     separate channel component of baseflow
!     from the non-channel component
  tbf=sbf*parea
!     tbf is total baseflow
  bfcc=tbf*(1.0/(1.0+side))
!     bfcc is baseflow, channel component
  bfp=spbf*parea/(1.0+side)
  bfs=bfcc-bfp
  if(bfs.lt.0.0) then
  bfs=0.0
  endif
  bfncc=tbf-bfcc
!     bfncc is baseflow,non-channel component
!
!     add to monthly sums.
  sintft=sintft+sif
  sgwfp=sgwfp+bfp
  sgwfs=sgwfs+bfs
  srecht=srecht+bfncc
  srost=srost+ssur
  srodt=srodt+sdro
!
!     compute total channel inflow for the time interval.
  tci=roimp+sdro+ssur+sif+bfcc
!
!     compute e4-et from riparian vegetation.
  e4=(edmnd-eused)*riva
!
!     subtract e4 from channel inflow
  tci=tci-e4
  if(tci.lt.0.0) then
  e4=e4+tci
  tci=0.0
  endif

  srot=srot+tci
!
!     compute total evapotranspiration-tet
  eused=eused*parea
  tet=eused+e5+e4
  sett=sett+tet
  se1=se1+e1*parea
  se3=se3+e3*parea
  se4=se4+e4
  se5=se5+e5
!     check that state(6).ge.state(1)
  if (state(6).lt.state(1)) then
  state(6)=state(1)
  endif

  return
  end subroutine sac
!**********************************************************************
  subroutine get_sac_params(nwsid,rfc,ncs,nwsidsac,rexp,lzpk,&
             lzfpm,pxadj,idt,pfree,zperc,riva,peadj,lztwm,rserv,&
             adimp,uzk,side,lzfsm,lzsk,uztwm,uzfwm,pctim,efc,ett,&
             saved,parea,nn,state,c_strong)
!**********************************************************************
! subroutine reads SAC parameters
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
real :: lztwm,lzfpm,lzfsm,lzpk,lzsk
  dimension ett(24),state(nn)
  character nwsid*5,nwsidsac*20,string*132,rfc*2,c_strong*1
!
! extract SAC parameters and ET demand curve
!
  open(10,file='SACSMA_'//nwsidsac(1:ncs)//'_UpdateStates.xml',&
  status='old')

  open(20,file=nwsid//'.deck1_'//c_strong,status='unknown')

c70 : do i=1,10000

  read(10,30,end=20) string
   30 format(a132)
!
! extract SAC parameters
!
c40 : do j=1,132
  if(j+7.gt.132) cycle c40
 
  if(string(j:j+7).eq.'dblValue') then

c50 : do k=j+8,132
  if(k+7.gt.132) cycle c50

  if(string(k:k+7).eq.'dblValue') then
  write(20,*) string(j+9:k-3)
  endif

  enddo c50

  endif

  enddo c40
!
! extract RUNOFF_COMPONENT_INTERVAL and SASC_INPUT_OPTION
!
c160 : do j=1,132
  if(j+7.gt.132) cycle c160

  if(string(j:j+7).eq.'intValue') then

c170 : do k=j+8,132
  if(k+7.gt.132) cycle c170

  if(string(k:k+7).eq.'intValue') then
  write(20,*) string(j+9:k-3)
  cycle c70
  endif

  enddo c170

  endif

  enddo c160
!
! extract ET demand curve
!
c60 : do j=1,132
  if(string(j:j+5).eq.'row A=') then
  write(20,*) string(j+7:j+9)
  cycle c70
  endif
  enddo c60

  enddo c70
   20 close(10)
  close(20)

  if(rfc.eq.'cn'.or.&
     rfc.eq.'lm'.or.&
     rfc.eq.'mb'.or.&
     rfc.eq.'nc') then
  open(20,file=nwsid//'.deck1_'//c_strong,status='old')
  idt=6
  read(20,*) rexp
  read(20,*) lzpk
  read(20,*) uzk
  read(20,*) side
  read(20,*) lzfpm
  read(20,*) pxadj
  read(20,*) pfree
  read(20,*) lzsk
  read(20,*) lzfsm
  read(20,*) zperc
  read(20,*) riva
  read(20,*) peadj
  read(20,*) uztwm
  read(20,*) uzfwm
  read(20,*) lztwm
  read(20,*) rserv
  read(20,*) pctim
  read(20,*) adimp
  read(20,*) efc
  do i=1,12
  read(20,*) ett(i)
  enddo
  close(20)
  endif

  if(rfc.eq.'ne') then
  open(20,file=nwsid//'.deck1_'//c_strong,status='old')
  read(20,*) rexp
  read(20,*) lzpk
  read(20,*) lzfpm
  read(20,*) pxadj
  read(20,*) idt
  read(20,*) pfree
  read(20,*) zperc
  read(20,*) riva
  read(20,*) peadj
  read(20,*) lztwm
  read(20,*) rserv
  read(20,*) adimp
  read(20,*) uzk
  read(20,*) side
  read(20,*) lzfsm
  read(20,*) lzsk
  read(20,*) ismzc_interval
  if(ismzc_interval.ne.idt) then
  write(6,*) 'ismzc_interval ne idt...stop ',nwsid
  write(6,*) ismzc_interval,idt
  stop
  endif
  read(20,*) uztwm
  read(20,*) uzfwm
  read(20,*) pctim
  read(20,*) efc
  do i=1,12
  read(20,*) ett(i)
  enddo
  close(20)
  endif

  if(rfc.eq.'se'.or.&
     rfc.eq.'oh') then
  open(20,file=nwsid//'.deck1_'//c_strong,status='old')
  idt=6
  read(20,*) rexp
  read(20,*) lzpk
  read(20,*) uzk
  read(20,*) side
  read(20,*) lzfpm
  read(20,*) idt
  read(20,*) pxadj
  read(20,*) pfree
  read(20,*) lzsk
  read(20,*) lzfsm
  read(20,*) ismzc_interval
  if(ismzc_interval.ne.idt) then
  write(6,*) 'ismzc_interval ne idt...stop ',nwsid
  write(6,*) ismzc_interval,idt
  stop
  endif
  read(20,*) zperc
  read(20,*) riva
  read(20,*) peadj
  read(20,*) uztwm
  read(20,*) lztwm
  read(20,*) uzfwm
  read(20,*) rserv
  read(20,*) pctim
  read(20,*) adimp
  read(20,*) efc
  do i=1,12
  read(20,*) ett(i)
  enddo
  close(20)
  endif

  saved=rserv*(lzfpm+lzfsm)
  parea=1.0-pctim-adimp

  state(1)=0.1*uztwm
  state(2)=0.1*uzfwm
  state(3)=0.1*lztwm
  state(4)=0.1*lzfsm
  state(5)=0.1*lzfpm
  state(6)=0.1*(uztwm+lztwm)

  if (state(1).gt.uztwm) state(1)=uztwm
  if (state(2).gt.uzfwm) state(2)=uzfwm
  if (state(3).gt.lztwm) state(3)=lztwm
  if (state(4).gt.lzfsm) state(4)=lzfsm
  if (state(5).gt.lzfpm) state(5)=lzfpm
  if (state(6).gt.(uztwm+lztwm)) state(6)=uztwm+lztwm
  if (state(6).lt.state(1)) state(6)=state(1)

  return
  end subroutine get_sac_params
!**********************************************************************
  subroutine get_uhg1(nwsid,rfc,ncu,nwsiduhg,nv,c_strong)
!**********************************************************************
! subroutine reads part 1 of UH deck
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  character nwsid*5,nwsiduhg*20,string*132,rfc*2,c_strong*1
!
! extract UH
!
  open(10,file='/home/djseo/6hrly/R/'//rfc//'/'//&
  'UNITHG_'//nwsiduhg(1:ncu)//'_UpdateStates.xml',status='old')

  open(20,file=nwsid//'.deck2_'//c_strong,status='unknown')

  nv=0
c100 : do i=1,10000

  read(10,30,end=80) string
   30 format(a132)
!
! extract UHG_DURATION and UHG_INTERVAL
!
c130 : do j=1,132
  if(j+7.gt.132) cycle c130

  if(string(j:j+7).eq.'intValue') then

c110 : do k=j+8,132
  if(k+7.gt.132) cycle c110

  if(string(k:k+7).eq.'intValue') then
  write(20,*) string(j+9:k-3)
  endif

  enddo c110

  endif

  enddo c130
!
! extract DRAINAGE_AREA and CONSTANT_BASE_FLOW
!
c90 : do j=1,132
  if(j+7.gt.132) cycle c90

  if(string(j:j+7).eq.'dblValue') then

c140 : do k=j+8,132
  if(k+7.gt.132) cycle c140

  if(string(k:k+7).eq.'dblValue') then
  write(20,*) string(j+9:k-3)
  cycle c100
  endif

  enddo c140

  endif

  enddo c90
!
! extract UHG_ORDINATES
!
c120 : do j=1,132
  if(string(j:j+5).eq.'row A=') then

c150: do k=j+7,132

  if(string(k:k).eq.'"') then
  write(20,*) string(j+7:k-1)
  nv=nv+1
  cycle c100
  endif

  enddo c150

  endif

  enddo c120

  enddo c100
   80 close(10)
  close(20)

  return
  end subroutine get_uhg1
!**********************************************************************
  subroutine get_uhg2(nwsid,nv,po,idtr,area,cbflow,c_strong)
!**********************************************************************
! subroutine reads part 2 of UH deck
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension po(nv)
  character nwsid*5,c_strong*1
  open(20,file=nwsid//'.deck2_'//c_strong,status='old')
  read(20,*) idtr
  read(20,*) idtq
  if(idtq.ne.idtr) then
  write(6,*) 'idtq ne idtr...stop ',idtq,idtr
  stop
  endif
  read(20,*) area
  read(20,*) cbflow     !CONSTANT_BASE_FLOW
  do i=1,nv
  read(20,*) po(i)
  po(i)=po(i)*0.3048**3/25.4
  enddo
  close(20)

  return
  end subroutine get_uhg2
!**********************************************************************
  subroutine get_map_length(nch,infile,k)
!**********************************************************************
! subroutine finds the length of the MAP time series
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  character infile*132,string*132,site1*5,cdum2*8
  open(10,file=infile(1:nch),status='old')

   12 read(10,11) string
   11 format(a132)
  if(string(1:1).eq.'$') go to 12
  if(string(25:26).eq.'MM') then
  fac1=1
  else if(string(25:26).eq.'IN') then
  fac1=25.4
  else
  write(6,*) 'unknown unit for MAPX...stop ',string(25:26)
  stop
  endif
!
! get the timestep
!
  read(string(29:30),'(i2)') jdt

  site1=string(35:39)
  read(10,11) string
!
! read precipitation data
!
  k=0
  jday_old=-999
   17 read(10,14,end=16) cdum2,jmon,jyr,jday,pxv
   14 format(a8,4x,i2,i2,2x,i2,f10.3)
  k=k+1
  if(jday.ne.jday_old) knt2=0
  knt2=knt2+1
  jday_old=jday
  if(jyr.gt.21) then
  jyr=jyr+1900
  else
  jyr=jyr+2000
  endif
  go to 17
   16 close(10)

  return
  end subroutine get_map_length
!**********************************************************************
  subroutine get_qin_length(nch,infile,k)
!**********************************************************************
! subroutine finds the length of the QIN time series
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  character infile*132,string*132,usgsid*8,cdum1*8

  open(40,file=infile(1:nch),status='old')
   15 read(40,11) string
   11 format(a132)
  if(string(1:1).eq.'$') go to 15
  if(string(15:17).ne.'QIN') then
  write(6,*) 'not QIN time series...stop'
  stop
  endif
  if(string(25:27).eq.'CMS') then
  fac2=1
  else if(string(25:27).eq.'CFS') then
  write(6,*) 'in CFS'
  fac2=0.3048**3
  else
  write(6,*) 'unknown unit for QIN...stop ',string(25:27)
  stop
  endif
  if(string(30:30).ne.'1') then
  write(6,*) 'not hourly data...stop'
  stop
  endif
  usgsid=string(35:42)
  read(40,11) string
!
! read streamflow data
!
  k=0
  iday_old=-999
   17 read(40,14,end=16) cdum1,imon,iyr,iday,qobs
   14 format(a8,4x,i2,i2,2x,i2,f10.3)
  k=k+1
  if(iday.ne.iday_old) knt2=0
  knt2=knt2+1
  iday_old=iday
  if(iyr.gt.21) then
  iyr=iyr+1900
  else
  iyr=iyr+2000
  endif
  go to 17
   16 close(40)

  return
  end subroutine get_qin_length
!**********************************************************************
  subroutine get_basin_info(nwsid,rfc,ncs,nwsidsac,ncu,nwsiduhg,&
                            qcut)
!**********************************************************************
! subroutine specifies the NWS ID for SAC and UH decks and cutoff flow
! flow identification of significant events
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  character rfc*2,nwsid*5,nwsidsac*20,nwsiduhg*20
!
! Basin 1
!
  if(nwsid.eq.'NFDC1') then                             !multiple zones
  rfc='cn'
  nwsidsac='NFDC1H_NFDC1HLF'
  ncs=len_trim(nwsidsac)
  nwsiduhg='NFDC1H_NFDC1'
  ncu=len_trim(nwsiduhg)
  qcut=130.
!
! Basin 2
!
  else if(nwsid.eq.'ABRN1') then
  rfc='mb'
  nwsidsac='ABRN1S_1610'
  ncs=len_trim(nwsidsac)
  nwsiduhg='ABRN1S_1610'
  ncu=len_trim(nwsiduhg)
  qcut=10.
!
! Basin 3
!
  else if(nwsid.eq.'ICLI4') then
  rfc='mb'
  nwsidsac='ICLI4S_1618'
  ncs=len_trim(nwsidsac)
  nwsiduhg='ICLI4S_1618'
  ncu=len_trim(nwsiduhg)
  qcut=75.
!
! Basin 4
!
  else if(nwsid.eq.'GTBM3') then
  rfc='ne'
  nwsidsac='GTBM3SNE_GTBM3SNE'
  ncs=len_trim(nwsidsac)
  nwsiduhg='GTBM3SNE_GTBM3SNE'
  ncu=len_trim(nwsiduhg)
  qcut=60.
!
! Basin 5
!
  else if(nwsid.eq.'COLI2') then
  rfc='nc'
  nwsidsac='COLI2_COLI2'
  ncs=len_trim(nwsidsac)
  nwsiduhg='COLI2_COLI2'
  ncu=len_trim(nwsiduhg)
  qcut=50.
!
! Basin 6
!
  else if(nwsid.eq.'DAVI3') then
  rfc='nc'
  nwsidsac='DAVI3_DAVI3'
  ncs=len_trim(nwsidsac)
  nwsiduhg='DAVI3_DAVI3'
  ncu=len_trim(nwsiduhg)
  qcut=25.
!
! Basin 7
!
  else if(nwsid.eq.'DLTC1') then
  rfc='cn'
  nwsidsac='DLTC1H_DLTC1HLF'
  ncs=len_trim(nwsidsac)
  nwsiduhg='DLTC1H_DLTC1'
  ncu=len_trim(nwsiduhg)
  qcut=175.
!
! Basin 8
!
  else if(nwsid.eq.'MONN7') then
  rfc='se'
  nwsidsac='MONN7_MONN7'
  ncs=len_trim(nwsidsac)
  nwsiduhg='MONN7_MONN7'
  ncu=len_trim(nwsiduhg)
  qcut=300.
!
! Basin 9
!
  else if(nwsid.eq.'NIMM5') then
  rfc='nc'
  nwsidsac='NIMM5_NIMM5'
  ncs=len_trim(nwsidsac)
  nwsiduhg='NIMM5_NIMM5'
  ncu=len_trim(nwsiduhg)
  qcut=25.
!
! Basin 10
!
  else if(nwsid.eq.'GAXV2') then
  rfc='oh'
  nwsidsac='GAXV2_GAXV2'
  ncs=len_trim(nwsidsac)
  nwsiduhg='GAXV2_GAXV2'
  ncu=len_trim(nwsiduhg)
  qcut=175.
!
! Basin 11
!
  else if(nwsid.eq.'TRYM7') then
  rfc='nc'
  nwsidsac='TRYM7_TRYM7'
  ncs=len_trim(nwsidsac)
  nwsiduhg='TRYM7_TRYM7'
  ncu=len_trim(nwsiduhg)
  qcut=120.
!
! Basin 12
!
  else if(nwsid.eq.'MCNM6') then
  rfc='se'
  nwsidsac='MCNM6_MCNM6'
  ncs=len_trim(nwsidsac)
  nwsiduhg='MCNM6_MCNM6'
  ncu=len_trim(nwsiduhg)
  qcut=175.
!
! Basin 13
!
  else if(nwsid.eq.'ORAI3') then
  rfc='oh'
  nwsidsac='ORAI3_ORAI3'
  ncs=len_trim(nwsidsac)
  nwsiduhg='ORAI3_ORAI3'
  ncu=len_trim(nwsiduhg)
  qcut=75.
!
! Basin 14
!
  else if(nwsid.eq.'GRDN6') then
  rfc='ne'
  nwsidsac='GRDN6HUD_GRDN6HUD'
  ncs=len_trim(nwsidsac)
  nwsiduhg='GRDN6HUD_GRDN6HUD'
  ncu=len_trim(nwsiduhg)
  qcut=160.

  else
  write(6,*) 'unknown location...stop ',nwsid
  stop
  endif

  return
end subroutine get_basin_info
!**********************************************************************
  subroutine get_coefr_consr(nwsid,coefr,consr)
!**********************************************************************
! subroutine specifies uncertainty parameters for additive error to TCI
! see Shen et al. (2021) for details
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  character nwsid*5

  if(nwsid.eq.'NFDC1') then
  coefr=0.49
  consr=1.68
  else if(nwsid.eq.'ABRN1') then
  coefr=0.44
  consr=1.73
  else if(nwsid.eq.'ICLI4') then
  coefr=0.42
  consr=1.35
  else if(nwsid.eq.'GTBM3') then
  coefr=0.37
  consr=0
  else if(nwsid.eq.'COLI2') then
  coefr=0.74
  consr=4.95
  else if(nwsid.eq.'DAVI3') then
  coefr=0.20
  consr=0.
  else if(nwsid.eq.'DLTC1') then
  coefr=0.39
  consr=5.95
  else if(nwsid.eq.'MONN7') then
  coefr=0.43
  consr=0
  else if(nwsid.eq.'NIMM5') then
  coefr=0.35
  consr=0
  else if(nwsid.eq.'GAXV2') then
  coefr=0.32
  consr=0
  else if(nwsid.eq.'TRYM7') then
  coefr=0.38
  consr=3.81
  else if(nwsid.eq.'MCNM6') then
  coefr=0.21
  consr=6.65
  else if(nwsid.eq.'ORAI3') then
  coefr=0.31 
  consr=0.
  else if(nwsid.eq.'GRDN6') then
  coefr=0.35
  consr=0.
  else
  write(6,*) 'unknown location...stop ',nwsid
  stop
  endif

  return
  end subroutine get_coefr_consr
!*********************************************************************
  subroutine get_obs_uncertainty(i,nwsid,nwin,nw,nf,nm,varq,&
             varp,vare,varw,sumd,zvar,ihetero,ifrq,ndata)
!*********************************************************************
! subroutine specifies observational uncertainties
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
use pass0
  dimension zvar(nm)
  character nwsid*5

  do ik=1,nwin
!
! for precip
!
  if(ihetero.eq.1) then
!
! model the error standard deviation of observed MAP
!
  coefp=0.39
  consp=0.25
  zstdp=coefp*pxv0(kbeg+ik-1)+consp
  zvar(ik)=zstdp**2
  else
  zvar(ik)=varp
  endif
!
! for PE
!
  if(ihetero.eq.1) then
  zvar(nwin+ik)=1.
  else
  zvar(nwin+ik)=vare
  endif
!
! for TCI error
!
  if(nw.eq.3) then

  if(ihetero.eq.1) then
!
! model the error std of runoff in terms of observed flow
!
  call get_coefr_consr(nwsid,coefr,consr)
  zstdr=coefr*qobs0(i)+consr
  zvar(2*nwin+ik)=zstdr**2/sumd
  else
  zvar(2*nwin+ik)=varw
  endif

  endif

  enddo
!
! for streamflow
!
  do ik=1,nf

  if(i+(ik-nf)*ifrq.lt.1.or.i+(ik-nf)*ifrq.gt.ndata) then
  write(6,*) 'qobs0 out of range...stop 1'
  stop
  endif

  if(ihetero.eq.1) then
!
! model the error standard deviation of observed streamflow
!
  coefq=0.15
  consq=0.25
  zstdq=coefq*qobs0(i+(ik-nf)*ifrq)+consq
  zvar(nw*nwin+ik)=zstdq**2
  else
  zvar(nw*nwin+ik)=varq
  endif
  enddo

  return
  end subroutine get_obs_uncertainty
!**********************************************************************
  subroutine get_map_mape_qin(nchx,map_name,nchq,qin1_name,idt,&
             ett,nwsid,ndata,jbeg1,kbeg1,mdata,pxv0,edmnd0,qobs0,&
             itime,dt)
!**********************************************************************
! subroutine reads MAP, MAPE and QIN data
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
  dimension ett(24),pxv0(mdata),edmnd0(mdata),qobs0(mdata),&
            itime(mdata)
  character map_name*132,string*132,site1*5,cdum2*8,qin1_name*132,&
            cdum1*8,usgsid*8,nwsid*5
!
! read any preceding precipitation data before the common period starts
!
  open(10,file=map_name(1:nchx),status='old')
   12 read(10,11) string
   11 format(a132)
  if(string(1:1).eq.'$') go to 12
  if(string(15:18).ne.'MAPX') then
  endif
  if(string(25:26).eq.'MM') then
  fac1=1
  else if(string(25:26).eq.'IN') then
  fac1=25.4
  else
  write(6,*) 'unknown unit for MAPX...stop'
  stop
  endif
!
! get the timestep
!
  read(string(29:30),'(i2)') jdt
  write(6,*) jdt

  site1=string(35:39)
  read(10,11) string
!
! read precipitation data before the common period starts
!
  do i=1,jbeg1-1
  read(10,14) cdum2,jmon,jyr,jday,pxv
  enddo
!
! read any preceding streamflow data before the common period starts
!
  open(20,file=qin1_name(1:nchq),status='old')
   15 read(20,11) string
  if(string(1:1).eq.'$') go to 15
  if(string(15:17).ne.'QIN') then
  write(6,*) 'not QIN time series...stop'
  stop
  endif
  if(string(25:27).eq.'CMS') then
  fac2=1
  else if(string(25:27).eq.'CFS') then
  write(6,*) 'in CFS'
  fac2=0.3048**3
  else
  write(6,*) 'unknown unit for QIN...stop ',string(25:27)
  stop
  endif
  usgsid=string(35:42)
  read(20,11) string
  write(6,*) string
  write(6,*) 'header read'
  write(6,*) 'kbeg1 ',kbeg1
!
! read streamflow data
!
  do i=1,kbeg1-1
  read(20,14) cdum1,imon,iyr,iday,qobs
  enddo
!
! initialize
!
  istart=0
!
! get precipitation, evapotranspiration demand, flow, and flow data
! quality flag for each hour
!
  ndata=0
  jday_old=-999

  do k=1,1000000
!
! read precipitation
!
  read(10,14,end=28) cdum2,jmon,jyr,jday,pxv
  if(jday.ne.jday_old) knt2=0
  knt2=knt2+1
  ihour=knt2*6
  if(nwsid.eq.'CREC1') ihour=knt2*1
  if(ihour.eq.0) ihour=24
  jday_old=jday
!
! cannot allow missing precipitation data 
!
  if(pxv.lt.0.) pxv=0.
!
! read observed flow in cms and quality flag
!
  read(20,14,end=28) cdum1,mon,iyr,iday,qobs
   14 format(a8,4x,i2,i2,2x,i2,f10.3)

  dt=idt/24.0
  edmnd=ett(jmon)*dt
!
! store precipitation, evapotranspiration demand, observed flow, and
! flow data quality flag
!
  pxv0(k)=pxv
  edmnd0(k)=edmnd
  qobs0(k)=qobs*fac2
  if(jyr.gt.21) jyr=jyr+1900
  if(jyr.le.21) jyr=jyr+2000
  itime(k)=jyr*1000000+jmon*10000+jday*100+ihour
  ndata=ndata+1

  enddo
!
! close the precipitation and flow files
!
   28 close(10)
  close(20)

  return
  end subroutine get_map_mape_qin
!**********************************************************************
  subroutine eval_aenkf_gain(nm,m,sigf,zvar,gaint,alpha,nwin,nf,&
                             sigu,sigu1,nw,ierr,h,ht,ri)
!**********************************************************************
! Subroutine evaluates AEnKF gain and AKF error covariance. See Seo et
! al. (2021) for details.
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
! nm     - number of observations (includes the faux obs)
! m      - number of state variables (includes the augmented)
! sigf   - forecast error covariance matrix
! zvar   - observation error variance vector
! gaint  - AKF gain matrix
! alpha  - alpha
! nwin   - length of the assimilation window in units of timesteps
! nf     - number of streamflow obs to be assimilated (1 recommended 
!          for SAC-UH)
! sigu   - AKF-updated error covariance matrix
! sigu1  - KF-updated apparent error covariance matrix associated with
!          sigf1
! nw     - 3 if weakly-constrained, 2 if strongly-constrained
! ierr   - error flag; a successful call returns ierr=0
! h      - structure matrix in the observation equation
! ht     - transpose of h
! ri     - inverse of the observation error covariance matrix
!
use pass0
use pass4

real, dimension (:,:), allocatable :: r
real, dimension (:,:), allocatable :: rj
real, dimension (:,:), allocatable :: sigusigu2i
real, dimension (:,:), allocatable :: sigf1
real, dimension (:,:), allocatable :: sigf2
real, dimension (:,:), allocatable :: sigu2
real, dimension (:,:), allocatable :: sigu2i
real, dimension (:,:), allocatable :: gain2t

  dimension sigf(m,m),zvar(nm),gaint(nm,m),sigu(m,m),&
            sigu1(m,m),h(nm,m),ht(m,nm),ri(nm,nm)

allocate(r(nm,nm))
allocate(rj(m,m))
allocate(sigusigu2i(m,m))
allocate(sigf1(m,m))
allocate(sigf2(m,m))
allocate(sigu2(m,m))
allocate(sigu2i(m,m))
allocate(gain2t(nm,m))

  ierr=0
!
! initialize the arrays
!
  h=0.
  ht=0.
  sigf1=0.
  sigf2=0.
  sigu=0.
  sigu1=0.
  sigu2=0.
  gaint=0.
  gain2t=0.
  sigu2i=0.
  sigusigu2i=0.

  call get_h_and_r(nm,m,nf,nw,nwin,h,ht,zvar,r,ri)
!
! inflate error covariance matrix
!
  do i=1,m
  do j=1,m
  sigf1(i,j)=(1.+alpha)*sigf(i,j)
  sigf2(i,j)=((1.+alpha)**2)*sigf(i,j)
  enddo
  enddo

  call get_kf_solution(nm,m,h,sigf1,ht,r,gaint,sigu1,ierr)

  if(ierr.ne.0) then
  write(6,*) 'ierr ne 0 from get_kf_solution 1...stop ',alpha
  stop
  endif

  call get_kf_solution(nm,m,h,sigf2,ht,r,gain2t,sigu2,ierr)

  if(alpha.gt.0) then
  if(ierr.ne.0) then
!
! this can occur due to numerical errors
!
  write(6,*) 'ierr ne 0 from get_kf_solution 2 ',alpha
  go to 10
  return
  endif
  endif
!
! invert sigu2
!
  do i=1,m
  do j=1,m
  rj(i,j)=0.
  if(i.eq.j) then
  rj(i,j)=1.
  endif
  enddo
  enddo

  call lsolve_cd(m,sigu2,m,rj,sigu2i,ierr)

  if(alpha.gt.0.) then
  if(ierr.ne.0) go to 10
  endif
!
! evaluated  updated error covariance
!
  sigusigu2i=matmul(sigu1,sigu2i)
  sigu=matmul(sigusigu2i,sigu1)

  do i1=1,m
  if(sigu(i1,i1).lt.0.) then
  write(6,*) 'negative diagonal element in sigu',i1
  ierr=1
  go to 10
  return
  endif
  enddo

   10 continue

deallocate(r)
deallocate(rj)
deallocate(sigusigu2i)
deallocate(sigf1)
deallocate(sigf2)
deallocate(sigu2)
deallocate(sigu2i)
deallocate(gain2t)

  return
  end subroutine eval_aenkf_gain
!**********************************************************************
  subroutine get_kf_solution(nm,m,h,sigf,ht,r,gaint,sigu,ierr)
!**********************************************************************
! Subroutine evaluates KF solution.
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
! nm     - number of observations (includes the faux obs)
! m      - number of state variables (includes the augmented)
! h      - structure matrix in the observation equation
! sigf   - forecast error covariance matrix
! ht     - transpose of h
! r      - observation error covariance matrix
! gaint  - KF gain matrix
! sigu   - KF-updated error covariance matrix
! ierr   - error flag; nonzero means an unsuccessful call
!
real, dimension (:,:), allocatable :: ui
real, dimension (:,:), allocatable :: hsigf
real, dimension (:,:), allocatable :: hsigfht
real, dimension (:,:), allocatable :: gain
real, dimension (:,:), allocatable :: gainh
real, dimension (:,:), allocatable :: uimgainh
real, dimension (:,:), allocatable :: rphsigfht

  dimension h(nm,m),sigf(m,m),ht(m,nm),r(nm,nm),gaint(nm,m),&
            sigu(m,m)

allocate(ui(m,m))
allocate(hsigf(nm,m))
allocate(hsigfht(nm,nm))
allocate(gain(m,nm))
allocate(gainh(m,m))
allocate(uimgainh(m,m))
allocate(rphsigfht(nm,nm))
!
! evaluate r+h*sigf*ht
!
  hsigf=matmul(h,sigf)
  hsigfht=matmul(hsigf,ht)

  rphsigfht=r+hsigfht
!
! solve for nxm Kalman gain
!
  call lsolve_cd(nm,rphsigfht,m,hsigf,gaint,ierr)

  if(ierr.ne.0) then
  write(6,*) 'ierr from lsolve_cd in get_kf_solution ne 0...stop ',&
  stop
  endif

  gain=transpose(gaint)
!
! evaluate sigu
!
  gainh=matmul(gain,h)

  do i=1,m
  do j=1,m
  ui(i,j)=0.
  if(i.eq.j) ui(i,j)=1.
  enddo
  enddo

  uimgainh=ui-gainh
  sigu=matmul(uimgainh,sigf)

deallocate(ui)
deallocate(hsigf)
deallocate(hsigfht)
deallocate(gain)
deallocate(gainh)
deallocate(uimgainh)
deallocate(rphsigfht)

  return
  end subroutine get_kf_solution
!**********************************************************************
  subroutine eval_opt_alpha1(m,alpha,sigfi,h,ht,nm,ri,sigu1,&
             sigfid,sigfid2,htri,temp1,sigu1a,temp4,htrih)
!**********************************************************************
! Subroutine evalulates the first- and second-order derivatives of 
! degrees of freedom for noise with respect to alpha to determine dfn-
! minimizing alpha; this is part 1 of the calculations before the
! ensemble loop. See Seo et al. (2021) for details.
!
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
! m       - number of state variables (including the augmented)
! alpha   - alpha
! sigfi   - inverse of forecast error covariance matrix
! h       - structure matrix in the observation equation
! ht      - transpose of h
! nm      - number of observations (including the faux obs)
! ri      - inverse of observation error covariance matrix
! sigu1   - KF-updated error covariance associated with (1+alpha)*sigf
! sigfid  - inverse of (1+alpha)*sigf
! sigfid2 - inverse of (1+alpha**2)*sigf
! htri    - h*ri
! temp1   - sigu1*sigfid2
! sigu1a  - derivative of sigu1 with respect to alpha
! temp4   - array holding intermidiate result for use by eval_opt_alpha2
! htrih   - ht*ri*h
!
real, dimension (:,:), allocatable :: temp3

  dimension h(nm,m),ht(m,nm),ri(nm,nm),sigu1(m,m),sigfid(m,m),&
            sigfid2(m,m),htri(m,nm),temp1(m,m),sigu1a(m,m),&
            sigfi(m,m),temp4(m,m),htrih(m,m)

allocate(temp3(m,m))
!
! evaluate inverse of (1+alpha)*sigf and {(1+alpha)**2}*sigf
!
  do i1=1,m
  do j1=1,m
  sigfid(i1,j1)=sigfi(i1,j1)/(1.+alpha)
  sigfid2(i1,j1)=sigfi(i1,j1)/(1.+alpha)**2
  enddo
  enddo
!
! evaluate ht*ri 
!
  htri=matmul(ht,ri)
!
! evaluate sigu1*sigfid2
!
  temp1=matmul(sigu1,sigfid2)
!
! evaluate d(sigu1)/d(alpha)=sigu1*sigfid2*sigu1
!
  sigu1a=matmul(temp1,sigu1)
!
! evaluate (-2/(1+alpha))*sigu1+dsigu1/dalpha
!
  temp3=(-2./(1.+alpha))*sigu1+sigu1a
!
! evaluate {(-2/(1+alpha))*sigu1+dsigu1/dalpha}*sigfid2
!
  temp4=matmul(temp3,sigfid2)
!
! evaluate ht*ri*h
!
  htrih=matmul(htri,h)

deallocate(temp3)

  return
  end subroutine eval_opt_alpha1
!**********************************************************************
  subroutine eval_opt_alpha2(m,alpha,h,nm,ymean,ybar,zobs,deriv1m,&
                             deriv2m,htri,temp1,temp4,htrih,ierr)
!**********************************************************************
! Subroutine evalulates gradient and hessian of degrees of freedom for
! noise with respect to alpha to determine dfn-minimizing alpha; this
! is part 2 of the calculation for each ensemble member. See Seo et 
! al. (2021) for details.
!
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
! m       - number of state variables (including the augmented)
! alpha   - current alpha on input, updated alpha on output
! h       - structure matrix in the observatation equation
! nm      - number of observaations (including the faux obs)
! ymean   - AEnKF analysis
! ybar    - one step-ahead model prediction
! zobs    - vector of all observations
! deriv1m - first-order derivative of dfn with respect to alpha
! deriv2m - second-order derivative of dfn with respect to alpha
! htri    - ht*ri
! temp1   - sigu1*sigfid2
! temp4   - {(-2/(1+alpha))*sigu1+dsigu1/dalpha}*sigfid2
! htrih   - ht*ri*h
! ierr    - error flag; a successful call returns ierr=0
!
real, dimension (:), allocatable :: deriv1
real, dimension (:), allocatable :: deriv2
real, dimension (:), allocatable :: temp5
real, dimension (:), allocatable :: temp6
real, dimension (:), allocatable :: wxf2
real, dimension (:), allocatable :: term1
real, dimension (:), allocatable :: hxhat
real, dimension (:), allocatable :: residu
real, dimension (:), allocatable :: term2

  dimension h(nm,m),ymean(m),ybar(m),zobs(nm),htri(m,nm),&
            temp1(m,m),temp4(m,m),htrih(m,m)

allocate(deriv1(m))
allocate(deriv2(m))
allocate(temp5(m))
allocate(temp6(m))
allocate(wxf2(m))
allocate(term1(m))
allocate(hxhat(nm))
allocate(residu(nm))
allocate(term2(m))

  ierr=0

  dlow=-1000.
  dupp= 1000.
!
! evaluate (AEnKF analysis-ybar)
!
  wxf2=ymean-ybar
!
! evaluate d(AEnKF analysis)/d(alpha)
! =sigu1*sigfid2*(AEnKF analysis-ybar)
!
  deriv1=matmul(temp1,wxf2)
!
! evaluate {(-2/(1+alpha))*sigu1+dsigu1/dalpha}*sigfid2*(xhat-ybar)
!
  temp5=matmul(temp4,wxf2)
!
! evaluate sigu1*sigfid2*d(AEnKF analysis)/d(alpha)
!
  temp6=matmul(temp1,deriv1)
!
! evaluate (d^2AEnKF analysis/dalpha^2)=
! {(-2/(1+alpha))*sigu1+dsigu1/dalpha}*sigfid2*(AEnKF analysis-ybar)+
! sigu1*sigfid2*d(AEnKF analysis)/d(alpha)
!
  deriv2=temp5+temp6
!
! evaluate h*AEnKF analysis
!
  hxhat=matmul(h,ymean)
!
! evaluate z-h*AEnKF analysis
!
  residu=zobs-hxhat
!
! evaluate ht*ri*(z-h*AEnKF analysis)
!
  term2=matmul(htri,residu)
!
! evaluate d(dfn)/d(alpha)
!
  deriv1m=0.
  do i1=1,m
  deriv1m=deriv1m+deriv1(i1)*term2(i1)
  enddo
  deriv1m=-2.*deriv1m

  if(deriv1m.lt.dlow.or.deriv1m.gt.dupp) then
  alpha=0.
  return
  endif
!
! evaluate -2*(d^2 AEnKF analysis/d alpha^2)*ht*ri*(z-h*AEnKF analysis)
!
  deriv2m=0.
  do i1=1,m
  deriv2m=deriv2m+deriv2(i1)*term2(i1)
  enddo
  deriv2m=-2.*deriv2m
!
! evaluate 2*(dAEnKF analysis/dalpha)*htrih*dAEnKF analysis/dalpha
!
  term1=matmul(htrih,deriv1)
!
! evaluate d^2dfn/dalpha^2=
!-2*(d^2AEnKF analysis/dalpha^2)*htri*(z-h*xlat)
!+2*(dAEnKF analysis/dalpha)*htrih*dAEnKF analysis/dalpha
!
  deriv3m=0.
  do i1=1,m
  deriv3m=deriv3m+deriv1(i1)*term1(i1)
  enddo
  deriv2m=deriv2m+2.*deriv3m

  if(deriv2m.lt.dlow.or.deriv2m.gt.dupp) then
  alpha=0.
  return
  endif

  if(deriv2m.eq.0.) then
  write(6,*) 'deriv2m eq 0...set to 0.01 ',deriv2m
  deriv2m=0.01
  ierr=0
  return
  endif
!
! calculate alpha_new=alpha_old-deriv1m/deriv2m
!
  alpha=alpha-deriv1m/deriv2m

  if(alpha.lt.0.) alpha=0.

deallocate(deriv1)
deallocate(deriv2)
deallocate(temp5)
deallocate(temp6)
deallocate(wxf2)
deallocate(term1)
deallocate(hxhat)
deallocate(residu)
deallocate(term2)

  return
  end subroutine eval_opt_alpha2
!**********************************************************************
  subroutine get_h_and_r(nm,m,nf,nw,nwin,h,ht,zvar,r,ri)
!**********************************************************************
! subroutine specifies the observation structure matrix, h, its 
! transpose, observation error matrix, r, and its inverse, ri
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
use pass0
use pass4

  dimension h(nm,m),ht(m,nm),zvar(nm),r(nm,nm),ri(nm,nm)
!
! specify the observation structure matrix
!
  h=0.
!
! for MAP
!
  do jj=1,nsp
  ibeg=(jj-1)*nhrp+1
  iend=jj*nhrp
  do ii=ibeg,iend
  h(ii,6+jj)=pxv0(kbeg+ii-1)
  enddo
  enddo
!
! for MAPE
!
  do jj=1,nse
  ibeg=nwin+(jj-1)*nhre+1
  iend=nwin+jj*nhre
  do ii=ibeg,iend
  h(ii,6+nsp+jj)=edmnd0(kbeg+ii-(nwin+1))
  enddo
  enddo
!
! for TCI
!
  if(nw.eq.3) then
  do jj=1,nsw
  ibeg=2*nwin+(jj-1)*nhrw+1
  iend=2*nwin+jj*nhrw
  do ii=ibeg,iend
  h(ii,6+nsp+nse+jj)=1.
  enddo
  enddo
!
! for Qobs
!
  do ii=3*nwin+1,3*nwin+nf
  do jj=1,nf
  if(ii-3*nwin.eq.jj) h(ii,6+nsp+nse+nsw+jj)=1.
  enddo
  enddo
  endif

  if(nw.eq.2) then
  do ii=2*nwin+1,2*nwin+nf
  do jj=1,nf
  if(ii-2*nwin.eq.jj) h(ii,6+nsp+nse+jj)=1.
  enddo
  enddo
  endif

  ht=transpose(h)
!
! construct the measurement error covariance matrix
!
  do i=1,nm
  do j=1,nm
  r(i,j)=0.
  ri(i,j)=0.
  if(i.eq.j) then
  r(i,j)=zvar(i)
  ri(i,j)=1./zvar(i)
  endif
  enddo
  enddo

  return
  end subroutine get_h_and_r
!**********************************************************************
  subroutine eval_enkf_gain(nm,m,sigf,zvar,gaint,nwin,nf,nw,ierr,&
                            h,ht,ri,sigu)
!**********************************************************************
! subroutine evaluates enkf gain
! version Aug 08, 2021, by D.-J. Seo at UTA/HWRL
!
! nm     - number of observations (includes the faux obs)
! m      - number of state variables (includes the augmented)
! sigf   - forecast error covariance matrix
! zvar   - observation error variance vector
! gaint  - KF gain matrix
! alpha  - alpha
! nwin   - length of the assimilation window in units of timesteps
! nf     - number of streamflow obs to be assimilated (use 1 for SAC-UH)
! nw     - 3 if weakly constrained, 2 if strongly constrained
! ierr   - error flag; nonzero means an unsuccessful call
! h      - structure matrix in the observation equation
! ht     - transpose of h
! ri     - inverse of the observation error covariance matrix
! sigu   - KF-updated error covariance matrix
!
use pass0
use pass4

real, dimension (:,:), allocatable :: gain
real, dimension (:,:), allocatable :: gainh
real, dimension (:,:), allocatable :: sigfht
real, dimension (:,:), allocatable :: hsigfht
real, dimension (:,:), allocatable :: hsigfhtpr
real, dimension (:,:), allocatable :: sigfhtt
real, dimension (:,:), allocatable :: hsigf
real, dimension (:,:), allocatable :: gainhsigf
real, dimension (:,:), allocatable :: r

  dimension sigf(m,m),zvar(nm),gaint(nm,m),h(nm,m),ht(m,nm),&
            ri(nm,nm),sigu(m,m)

  ierr=0

allocate(gain(m,nm))
allocate(gainh(m,m))
allocate(sigfht(m,nm))
allocate(hsigfht(nm,nm))
allocate(hsigfhtpr(nm,nm))
allocate(sigfhtt(nm,m))
allocate(hsigf(nm,m))
allocate(gainhsigf(m,m))
allocate(r(nm,nm))
!
! specify the number of measurements
!
  call get_h_and_r(nm,m,nf,nw,nwin,h,ht,zvar,r,ri)
  sigfht=matmul(sigf,ht)
  hsigfht=matmul(h,sigfht)
  hsigfhtpr=hsigfht+r
  sigfhtt=transpose(sigfht)
  call lsolve_cd(nm,hsigfhtpr,m,sigfhtt,gaint,ierr) !for enkf

  if(ierr.ne.0) then
  write(6,*) 'ierr enkf failed: impossible event...stop ',ierr
  stop
  endif

  gain=transpose(gaint)
  gainh=matmul(gain,h)
!
! calculate the degrees of freedom for signal
!
  dfs=0.
  do ig=1,m
  dfs=dfs+gainh(ig,ig)
  enddo
!
! multiply gain and hsigf
!
  hsigf=transpose(sigfht)
!
! subtract gainhsigf from sigf
!
  sigu=sigf-gainhsigf

deallocate(gain)
deallocate(gainh)
deallocate(sigfht)
deallocate(hsigfht)
deallocate(hsigfhtpr)
deallocate(sigfhtt)
deallocate(hsigf)
deallocate(gainhsigf)
deallocate(r)

  return
  end subroutine eval_enkf_gain
