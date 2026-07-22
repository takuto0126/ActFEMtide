!# modified on 2017.07.13 for multiple sources
!# Generated on 2017.05.10
module param_inv
use modelpart
use param        ! added on May 13, 2017
implicit none

type obs_src          ! 2017.07.13
 integer(4) :: nobs_s ! nobs for one source 2017.07.13
 character(50),allocatable,dimension(:) :: obsfile ! (nobs) 2017.07.13
 integer(4),   allocatable,dimension(:) :: obsindex! (nobs) 2017.07.13
end type        ! 2017.07.13

type param_inversion
 integer(4)    :: nobs           ! same as g_param 2017.07.13
 integer(4)    :: nfreq          ! same as g_param 2017.07.13
 integer(4)    :: nsr_inv        ! 2017.07.13
 integer(4)    :: ndat           !defined in readdata 2017.07.14
 integer(4)    :: ialphaflag     ! 1 for L-curve, 2 for cooling strategy 2017.07.19
 !# for L-curve strategy                     ialphaflag = 1
 integer(4)    :: nalpha
 real(8),allocatable,dimension(:) :: alpha
 !# for Kordy et al. (2016) cooling strategy ialphaflag = 2
 real(8)       :: alpha_init     ! 2017.07.14
 !#
 character(50) :: g_faceinfofile ! added on March 10, 2017
 !# init cond or model
 integer(4)    :: icondflag_ini  ! 2018.03.18
 character(50) :: g_initmodelconn! 2018.03.18
 character(50) :: g_initcondfile ! added on May 13, 2017
 character(50) :: g_initmodelfile! 2018.03.18
 !# ref cond or model
 integer(4)    :: icondflag_ref  ! 2018.03.18
 character(50) :: g_refmodelconn ! 2018.03.18
 character(50) :: g_refcondfile  ! 2017.07.19
 character(50) :: g_refmodelfile ! 2018.03.18
 character(50) :: outputfolder   ! 2017.07.25
 integer(4),   allocatable,dimension(:) :: srcindex ! (nsr_inv) 2017.07.13
 type(obs_src),allocatable,dimension(:) :: obsinfo  ! (nsr_inv) 2017.07.13
 logical,allocatable,dimension(:,:,:) :: data_avail !(nfreq,nobs,nsr_inv) 2017.07.13
 real(8)       :: errorfloor     ! [%] added on May 13, 2017
 !# output level 0:default, 1;output Jacobian  2018.03.29
 integer(4)    :: ioutlevel      ! 2018.03.29
 real(8)       :: finalrms       ! optional , if less than 0.1, reset to 1.0 2018.06.14
end type

type data_vec  ! added on May 13, 2017
 integer(4)    :: nsr_inv   ! 2017.07.13
 integer(4)    :: nobs
 integer(4)    :: nfreq
 integer(4)    :: ndat ! depends on data_avail 2017.07.14
 logical,allocatable,dimension(:,:,:) :: data_avail !(nfreq,nobs,nsr_inv) 2017.07.13
 real(8),allocatable,dimension(:)     :: dvec
 real(8),allocatable,dimension(:)     :: error
end type

contains

!################################################## readparam
! modified on 2017.07.13 for multiple sources
subroutine readparainv(g_param_inv,g_modelpara,g_param,sparam,g_data) ! 2018.03.20
 implicit none
 type(param_forward),  intent(in)  :: g_param
 type(param_source),   intent(in)  :: sparam     ! 2018.03.20
 type(param_inversion),intent(out) :: g_param_inv
 type(modelpara),      intent(out) :: g_modelpara
 type(data_vec),       intent(out) :: g_data
 character(50)  :: a ! 2017.06.06
 integer(4) :: input=5, nobs_s,i,j,nsr_inv
 integer(4) :: icondflag ! 1: cond, 2: model

!#[1]## read files
 read(input,*)  ! for header
 write(*,*) "input path for faceinfo file"
 read(input,'(20x,a)') g_param_inv%g_faceinfofile ! face info
 write(*,*) "faceinfo file is",g_param_inv%g_faceinfofile
 !# ref cond/model file
 write(*,*) "[ref] Input 1 (condfile), or 2 (modelconnect and model file)" ! 2018.03.17
 read(input,'(20x,i10)') g_param_inv%icondflag_ref                 ! 2018.03.17
 if ( g_param_inv%icondflag_ref .eq. 1 ) then                      ! 2018.03.17
  write(*,*) "input path for reference conductivity structure file"! 2017.07.19
  read(input,'(20x,a)') g_param_inv%g_refcondfile ! ref cond file    2017.07.19
  write(*,*) "refcondfile=",g_param_inv%g_refcondfile              ! 2017.07.19
 elseif (g_param_inv%icondflag_ref .eq. 2 ) then                   ! 2018.03.18
  write(*,*) "Input modelconnect file and modelfile"               ! 2018.03.18
  read(input,'(20x,a)') g_param_inv%g_refmodelconn ! ref cond file   2018.03.18
  read(input,'(20x,a)') g_param_inv%g_refmodelfile ! ref cond file   2018.03.18
  write(*,*) "refmodelconnect",g_param_inv%g_refmodelconn          ! 2018.03.18
  write(*,*) "refmodelfile",   g_param_inv%g_refmodelfile          ! 2018.03.18
 end if
 !# init cond/model file
 write(*,*) "[ini] Input 1 (condfile), or 2 (modelconnect and model file)" ! 2018.03.18
 read(input,'(20x,i10)') g_param_inv%icondflag_ini                 ! 2018.03.18
 if ( g_param_inv%icondflag_ini .eq. 1 ) then                      ! 2018.03.18
  write(*,*) "input path for initial conductivity structure file"
  read(input,'(20x,a)') g_param_inv%g_initcondfile ! init cond file
  write(*,*) "initcondfile=",g_param_inv%g_initcondfile
 elseif (g_param_inv%icondflag_ini .eq. 2 ) then                   ! 2018.03.18
  write(*,*) "Input modelconnect file and modelfile"               ! 2018.03.18
  read(input,'(20x,a)') g_param_inv%g_initmodelconn ! ref cond file  2018.03.18
  read(input,'(20x,a)') g_param_inv%g_initmodelfile ! ref cond file  2018.03.18
  write(*,*) "inimodelconnect",g_param_inv%g_initmodelconn         ! 2018.03.18
  write(*,*) "inimodelfile",   g_param_inv%g_initmodelfile         ! 2018.03.18
 end if
 write(*,*) "input output folder for inversion result" ! 2017.07.25
 read(input,'(20x,a)') g_param_inv%outputfolder         ! 2017.07.25
 write(*,*) "output folder for inversion is",g_param_inv%outputfolder!2017.07.25
 write(*,*) "input error floor [0-1]"
 read(input,'(20x,g15.7)') g_param_inv%errorfloor ! error floor
 write(*,*) "errorfloor=",g_param_inv%errorfloor
 write(*,*) "1 for L-curve strategy; 2 for cooling strategy (ialphaflag)"
 read(input,'(20x,i10)') g_param_inv%ialphaflag
 !# L curve strategy # 2017.07.19
 if ( g_param_inv%ialphaflag .eq. 1) then! L-curve strategy
  write(*,*) " # of alpha"
  read(input,'(20x,i10)') g_param_inv%nalpha
  allocate(g_param_inv%alpha(g_param_inv%nalpha))
  do i=1,g_param_inv%nalpha
   read(input,'(20x,g15.7)')g_param_inv%alpha(i)
  end do
 !# cooling strategy # 2017.07.19
 else if (g_param_inv%ialphaflag .eq. 2) then!
  write(*,*) "initial alpha [real]"
  read(input,'(20x,g15.7)') g_param_inv%alpha_init !2017.07.14
  write(*,*) "initial alpha [real]=",g_param_inv%alpha_init
 else
  write(*,*) "g_param_inv%ialphaflag should be 1 or 2",g_param_inv%ialphaflag
  stop
 end if
!# alpha input end # 2017.07.19
 write(*,*) "input # of sources for inversion" ! 2017.07.13
 read(input,'(20x,i)') nsr_inv                ! 2017.07.13
 g_param_inv%nsr_inv = nsr_inv               ! 2017.07.13
 allocate(g_param_inv%srcindex(nsr_inv))      ! 2017.07.13
 allocate(g_param_inv%obsinfo(nsr_inv) )      ! 2017.07.13

!#[2]## read data files
 g_param_inv%nobs  = g_param%nobs
 g_param_inv%nfreq = g_param%nfreq

 do i=1,nsr_inv
  !# srcindex(i)
  write(*,*) "input source index specified in forward ctl file"
  read(input,'(20x,i)') g_param_inv%srcindex(i)
  !write(*,*) "srcindex=",g_param_inv%srcindex(i)
  !# nobs_s
  write(*,*) "input # of observations for inversion for current source"
  read(input,'(20x,i)') nobs_s
  if ( g_param_inv%nobs .lt. nobs_s) goto 100
  g_param_inv%obsinfo(i)%nobs_s = nobs_s           ! 2017.07.13
  allocate(g_param_inv%obsinfo(i)%obsfile(nobs_s)) ! 2017.07.13
  allocate(g_param_inv%obsinfo(i)%obsindex(nobs_s))! 2017.07.13
  do j=1,nobs_s
   read(input,'(20x,i5,a)') g_param_inv%obsinfo(i)%obsindex(j),&
				&   g_param_inv%obsinfo(i)%obsfile(j)
!   write(*,*) g_param_inv%obsinfo(i)%obsfile(j),&
!   &                       g_param_inv%obsinfo(i)%obsindex(j)
  end do ! nobs_s loop
 end do  ! nsr_inv loop
 call readdata(g_param_inv,g_data,g_param,sparam) ! 2018.03.20

!#[3]## read model partitions
 CALL readmodelpara(input,g_modelpara)

!#[4]## read output level  2018.03.29
 g_param_inv%ioutlevel=0
 write(*,*) "Input 0 for default or 1 for output Jacobian matrix"
 read(input,'(20x,i10)',err=90) g_param_inv%ioutlevel
 90 continue
 write(*,*) "enter final rms for convergence"        ! 2018.06.14
 read(input,'(20x,g15.7)',err=91) g_param_inv%finalrms ! 2018.06.14
 goto 92                   ! 2018.06.14
 91 continue               ! 2018.06.14
 g_param_inv%finalrms=1.d0 ! 2018.06.14
 92 continue
 write(*,*) "ioutlevel =",g_param_inv%ioutlevel
 write(*,*) "finalrms =",g_param_inv%finalrms    ! 2018.06.14

return

100 continue
 write(*,*) "GEGEGE # of nobs_s should be smaller than nobs specified for forward"
 write(*,*) "nobs",g_param_inv%nobs,"nobs_s",nobs_s
 stop

end

!#################################################### readdata
!# modified on 2018.03.20 for output obsfile with error adopted
!# modified on 2017.07.13 for multiple sources
!# Modified on 2017.05.18
!# Coded on May 13, 2017
subroutine readdata(g_param_inv,g_data,g_param,sparam) ! 2018.03.20
implicit none
type(param_inversion),    intent(inout)    :: g_param_inv ! 2017.07.14 data_avail
type(param_forward),      intent(in)       :: g_param     ! 2018.03.20
type(param_source),       intent(in)       :: sparam      ! 2018.03.20
type(data_vec),           intent(out)      :: g_data
character(50),allocatable,dimension(:)     :: obsfile
real(8),      allocatable,dimension(:)     :: dat,err
real(8),      allocatable,dimension(:,:,:) :: dd, ee      ! 2017.07.13
integer(4),   allocatable,dimension(:,:,:) :: iflag       ! 2017.07.13
logical,      allocatable,dimension(:,:,:) :: data_avail  ! 2017.07.13
real(8)        :: d1,dr,e1,er,f1,ed, errorfloor
integer(4)     :: i,j,k,ndat,nsr_inv ! 2017.07.13
integer(4)     :: nobs,nfreq
integer(4)     :: nobs_s,iobs         ! 2017.07.13
integer(4)     :: ndatmax,icount      ! 2017.07.13
character(50)  :: head,obsfile_inv    ! 2018.03.20
character(50)  :: site,sour           ! 2018.03.20
real(8),      allocatable,dimension(:)     :: freq              ! 2018.03.20
real(8),      allocatable,dimension(:,:,:) :: dat_obs           ! 2018.03.20
real(8),      allocatable,dimension(:,:,:) :: err_obs           ! 2018.03.20
integer(4)                                 :: nhead,nsite,nsour ! 2018.03.20
integer(4),   allocatable,dimension(:)     :: srcindex          ! 2018.03.20

!#[1]## set
 nobs       = g_param_inv%nobs        ! same as g_param 2017.07.13
 nfreq      = g_param_inv%nfreq       ! same as g_param 2017.07.13
 nsr_inv    = g_param_inv%nsr_inv     ! 2017.07.13
 errorfloor = g_param_inv%errorfloor  ! [%]
 allocate(srcindex(nsr_inv), freq(nfreq) ) ! 2017.03.20
 srcindex   = g_param_inv%srcindex           ! 2017.02.30
 head       = g_param_inv%outputfolder  ! 2018.03.20
 nhead      = len_trim(head)            ! 2018.03.20

!#[2]## read and gen dat,err
 allocate(      dd(nfreq,nobs,nsr_inv)) ! 2017.07.13
 allocate(      ee(nfreq,nobs,nsr_inv)) ! 2017.07.13
 allocate( dat_obs(nfreq,nobs,nsr_inv)) ! 2018.03.20
 allocate( err_obs(nfreq,nobs,nsr_inv)) ! 2018.03.20
 allocate(   iflag(nfreq,nobs,nsr_inv)) ! 2017.07.13

!#[2-1]## just read
do k=1,nsr_inv
 nobs_s = g_param_inv%obsinfo(k)%nobs_s
 do i=1, nobs_s            ! 2017.07.13
  open(1,file=g_param_inv%obsinfo(k)%obsfile(i)) ! 2017.07.13
  write(*,*) "file=",g_param_inv%obsinfo(k)%obsfile(i)
  iobs = g_param_inv%obsinfo(k)%obsindex(i)

  do j=1,nfreq
!   write(*,*) "j,i,k",j,i,k,"iobs=",iobs
   read(1,*) freq(j),iflag(j,iobs,k),dd(j,iobs,k),ee(j,iobs,k)
  end do
  close(1)
 end do ! nobs_s   loop
end do  ! nsr_inv loop

!#[2-2]## assemble data vec
!# count ndata
ndatmax = nsr_inv * nobs * (nfreq - 1) ! 2017.07.13
allocate(dat(ndatmax),err(ndatmax))
allocate(data_avail(nfreq,nobs,nsr_inv)) !2017.07.13
data_avail(:,:,:)=.false.                        ! 2017.07.13
icount = 0

do j=2,nfreq
 do k=1,nsr_inv
 nobs_s = g_param_inv%obsinfo(k)%nobs_s
 do i=1,nobs_s
  iobs = g_param_inv%obsinfo(k)%obsindex(i)
   if ( iflag(j,iobs,k) .eq. 1 ) then    ! 2017.07.13
    data_avail(1,iobs,k) = .true.        ! 2017.07.17
    data_avail(j,iobs,k) = .true.        ! 2017.07.13
    icount = icount + 1                  ! 2017.07.13
    dr = dd(1,iobs,k)
    d1 = dd(j,iobs,k)
    er = ee(1,iobs,k)
    e1 = ee(j,iobs,k)
    ed = sqrt((e1/dr)**2.d0 + (d1/(dr**2.)*er)**2.d0)   ! 2018.03.20
    if ( ed/abs(d1/dr) .lt. errorfloor ) ed=abs(d1/dr)*errorfloor
    dat( icount )     = d1/dr ! 2017.07.13
    err( icount )     = ed    ! 2017.07.13
    dat_obs(j,iobs,k) = d1/dr ! 2018.03.20
    err_obs(j,iobs,k) = ed    ! 2018.03.20
   end if                     ! 2017.07.13
  end do
  close(1)
 end do
end do
 ndat = icount           ! 2017.07.13
 write(*,*) "ndat=",ndat ! 2017.07.13

!# output obs.dat.inv files with error adopted 2018.03.20
do iobs=1,nobs
   site  = g_param%obsname(iobs)
   nsite = len_trim(site)
  do k=1,nsr_inv
   sour       = sparam%sourcename(srcindex(k))
   nsour      = len_trim(sour)
   obsfile_inv= head(1:nhead)//site(1:nsite)//"_"//sour(1:nsour)//"_obs.inv"

   open(1,file=obsfile_inv)
   do j=2,nfreq
    if ( data_avail(j,iobs,k) ) then
     write(1,*) freq(j),dat_obs(j,iobs,k),err_obs(j,iobs,k)
     end if
   end do
   close(1)

 end do
end do

!#[3]## output
! g_data
 g_data%nobs  = nobs
 g_data%nfreq = nfreq
 g_data%ndat  = nsr_inv ! 2017.07.14
 g_data%ndat  = ndat
 allocate(g_data%dvec(ndat),g_data%error(ndat))
 allocate(g_data%data_avail(nfreq,nobs,nsr_inv)) ! 2017.07.13
 g_data%dvec       = dat(1:ndat) ! 2017.07.13
 g_data%error      = err(1:ndat) ! 2017.07.13
 g_data%data_avail = data_avail  ! 2017.07.13
! g_param_inv
 allocate(g_param_inv%data_avail(nfreq,nobs,nsr_inv)) ! 2017.07.17
 g_param_inv%data_avail = data_avail  ! 2017.07.14
 g_param_inv%ndat       = ndat

return
end

!################################################# initializedatavec
!# 2017.07.14 modified for multisource inversion
subroutine initializedatavec(g_param_inv,g_data)
implicit none
type(param_inversion),intent(in)     :: g_param_inv ! 2017.07.14
type(data_vec),       intent(out)    :: g_data
!===============================================================
integer(4)                           :: nfreq,nobs,nsr_inv ! 2017.07.14
logical,allocatable,dimension(:,:,:) :: data_avail  ! 2017.07.14
integer(4)                           :: i,j,k,ndat

!#[0]## set
ndat    = g_param_inv%ndat    ! 2017.07.14
nobs    = g_param_inv%nobs    ! 2017.07.14
nfreq   = g_param_inv%nfreq   ! 2017.07.14
nsr_inv = g_param_inv%nsr_inv ! 2017.07.14
allocate(data_avail(nfreq,nobs,nsr_inv)) ! 2017.07.14
data_avail = g_param_inv%data_avail ! 2017.07.14

!#[1]## output
g_data%nobs    = nobs
g_data%nfreq   = nfreq
g_data%nsr_inv = nsr_inv     ! 2017.07.14
g_data%ndat    = ndat        ! 2017.07.14
allocate(g_data%dvec(ndat) ) ! 2017.07.14
allocate(g_data%error(ndat)) ! 2017.07.14
allocate(g_data%data_avail(nfreq,nobs,nsr_inv)) ! 2017.07.14
g_data%data_avail = data_avail

return
end

end module param_inv
