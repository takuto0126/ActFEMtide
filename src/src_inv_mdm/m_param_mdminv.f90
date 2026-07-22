!# 2017.10.15
module param_mdminv
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
 integer(4)    :: ndat1          ! former data in readdata 2017.10.15
 integer(4)    :: ndat2          ! later  data in readdata 2017.10.15
 integer(4)    :: ndat           ! total data 2017.10.15
 integer(4)    :: ialphaflag     ! 1 for L-curve, 2 for cooling strategy 2017.07.19
 !# for L-curve strategy                     ialphaflag = 1
 integer(4)    :: nalpha
 real(8),allocatable,dimension(:) :: alpha
 !# for Kordy et al. (2016) cooling strategy ialphaflag = 2
 real(8)       :: alpha_init1     ! 2017.11.01
 real(8)       :: alpha_init2     ! 2017.11.01
 !#
 character(50) :: g_faceinfofile ! added on March 10, 2017
 character(50) :: g_initcondfile ! added on May 13, 2017
 character(50) :: g_refcondfile  ! 2017.07.19
 character(50) :: outputfolder   ! 2017.07.25
 integer(4),   allocatable,dimension(:) :: srcindex ! (nsr_inv) 2017.07.13
 type(obs_src),allocatable,dimension(:) :: obsinfo1  ! (nsr_inv) 2017.10.15
 type(obs_src),allocatable,dimension(:) :: obsinfo2  ! (nsr_inv) 2017.10.15
 logical,allocatable,dimension(:,:,:) :: data_avail1 !(nfreq,nobs,nsr_inv)2017.07.13
 logical,allocatable,dimension(:,:,:) :: data_avail2 !(nfreq,nobs,nsr_inv)2017.07.13
 real(8)       :: errorfloor     ! [%] added on May 13, 2017

 !# change area for dmin 2017.10.15
 real(8)       :: xminc  ! c: change
 real(8)       :: xmaxc  ! c: change
 real(8)       :: yminc  ! c: change
 real(8)       :: ymaxc  ! c: change
 real(8)       :: zminc  ! c: change
 real(8)       :: zmaxc  ! c: change
end type

type data_vec  ! added on May 13, 2017
 integer(4)    :: nsr_inv   ! 2017.07.13
 integer(4)    :: nobs
 integer(4)    :: nfreq
 integer(4)    :: ndat1 ! depends on data_avail 2017.10.15
 integer(4)    :: ndat2 ! depends on data_avail 2017.10.15
 integer(4)    :: ndat  ! depends on data_avail  2017.10.15
 logical,allocatable,dimension(:,:,:) :: data_avail1 !(nfreq,nobs,nsr_inv)2017.07.13
 real(8),allocatable,dimension(:)     :: dvec1
 real(8),allocatable,dimension(:)     :: error1
 logical,allocatable,dimension(:,:,:) :: data_avail2 !(nfreq,nobs,nsr_inv)2017.07.13
 real(8),allocatable,dimension(:)     :: dvec2
 real(8),allocatable,dimension(:)     :: error2
 real(8),allocatable,dimension(:)     :: dvec        ! total
 real(8),allocatable,dimension(:)     :: error       ! total
end type

contains

!################################################## readparam
! modified on 2017.07.13 for multiple sources
subroutine readparainv(g_param_inv,g_modelpara,g_param,g_data)
 implicit none
 type(param_forward),  intent(in)  :: g_param
 type(param_inversion),intent(out) :: g_param_inv
 type(modelpara),      intent(out) :: g_modelpara
 type(data_vec),       intent(out) :: g_data
 character(50)  :: a ! 2017.06.06
 integer(4) :: input=5, nobs_s1,nobs_s2,i,j,nsr_inv

!#[1]## read files
 read(input,*)  ! for header
 write(*,*) "input path for faceinfo file"
 read(input,'(20x,a)') g_param_inv%g_faceinfofile ! face info
 write(*,*) "faceinfo file is",g_param_inv%g_faceinfofile
 write(*,*) "input path for reference conductivity structure file" ! 2017.07.19
 read(input,'(20x,a)') g_param_inv%g_refcondfile ! ref cond file     2017.07.19
 write(*,*) "refcondfile=",g_param_inv%g_refcondfile        !       2017.07.19
 write(*,*) "input path for initial conductivity structure file"
 read(input,'(20x,a)') g_param_inv%g_initcondfile ! init cond file
 write(*,*) "initcondfile=",g_param_inv%g_initcondfile
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
  read(input,'(20x,g15.7)') g_param_inv%alpha_init1 !2017.11.01 for m0
  read(input,'(20x,g15.7)') g_param_inv%alpha_init2 !2017.11.01 for dm
  write(*,*) "initial alpha1 [real]=",g_param_inv%alpha_init1 ! 2017.11.01
  write(*,*) "initial alpha2 [real]=",g_param_inv%alpha_init2 ! 2017.11.01
 else
  write(*,*) "g_param_inv%ialphaflag should be 1 or 2",g_param_inv%ialphaflag
  stop
 end if
!# alpha input end # 2017.07.19
 write(*,*) "input # of sources for inversion" ! 2017.07.13
 read(input,'(20x,i)') nsr_inv                 ! 2017.07.13
 g_param_inv%nsr_inv = nsr_inv                 ! 2017.07.13
 allocate(g_param_inv%srcindex(nsr_inv))       ! 2017.07.13
 allocate(g_param_inv%obsinfo1(nsr_inv) )      ! 2017.10.15 former data
 allocate(g_param_inv%obsinfo2(nsr_inv) )      ! 2017.10.15 latter data

!#[2]## read data files
 g_param_inv%nobs  = g_param%nobs
 g_param_inv%nfreq = g_param%nfreq

 do i=1,nsr_inv
  !# srcindex(i)
  write(*,*) "input source index specified in forward ctl file"
  read(input,'(20x,i)') g_param_inv%srcindex(i)
  !write(*,*) "srcindex=",g_param_inv%srcindex(i)

  !#[Former data] nobs_s1  2017.10.15
  write(*,*) "input # of observations for inversion for current source"
  read(input,'(20x,i)') nobs_s1 ! 2017.10.15
  write(*,*) "nobs for former is",nobs_s1 ! 2017.10.15
  if ( g_param_inv%nobs .lt. nobs_s1) goto 100
  g_param_inv%obsinfo1(i)%nobs_s = nobs_s1            ! 2017.10.15
  allocate(g_param_inv%obsinfo1(i)%obsfile( nobs_s1)) ! 2017.07.13
  allocate(g_param_inv%obsinfo1(i)%obsindex(nobs_s1)) ! 2017.07.13
  do j=1,nobs_s1
   read(input,'(20x,i5,a)') g_param_inv%obsinfo1(i)%obsindex(j),&
				&   g_param_inv%obsinfo1(i)%obsfile(j)
!   write(*,*) g_param_inv%obsinfo(i)%obsfile(j),&
!   &                       g_param_inv%obsinfo(i)%obsindex(j)
  end do ! nobs_s loop

  !#[Later data] nobs_s1  2017.10.15
  write(*,*) "input # of observations for inversion for current source"
  read(input,'(20x,i)') nobs_s2 ! 2017.10.15
  write(*,*) "nobs for later is",nobs_s2 ! 2017.10.15
  if ( g_param_inv%nobs .lt. nobs_s2) goto 100
  g_param_inv%obsinfo2(i)%nobs_s = nobs_s2            ! 2017.10.15
  allocate(g_param_inv%obsinfo2(i)%obsfile( nobs_s2)) ! 2017.07.13
  allocate(g_param_inv%obsinfo2(i)%obsindex(nobs_s2)) ! 2017.07.13
  do j=1,nobs_s2
   read(input,'(20x,i5,a)') g_param_inv%obsinfo2(i)%obsindex(j),&
				&   g_param_inv%obsinfo2(i)%obsfile(j)
!   write(*,*) g_param_inv%obsinfo(i)%obsfile(j),&
!   &                       g_param_inv%obsinfo(i)%obsindex(j)
  end do ! nobs_s loop

 end do  ! nsr_inv loop
 call readdata(g_param_inv,g_data) ! 2017.10.15

!#[3]## read area for dmin
  read(*,'(20x,g15.7)') g_param_inv%xminc
  read(*,'(20x,g15.7)') g_param_inv%xmaxc
  read(*,'(20x,g15.7)') g_param_inv%yminc
  read(*,'(20x,g15.7)') g_param_inv%ymaxc
  read(*,'(20x,g15.7)') g_param_inv%zminc
  read(*,'(20x,g15.7)') g_param_inv%zmaxc

!#[4]## read model partitions
 CALL readmodelpara(input,g_modelpara)

return

100 continue
 write(*,*) "GEGEGE # of nobs_s should be smaller than nobs specified for forward"
 write(*,*) "nobs",g_param_inv%nobs,"nobs_s1, s2",nobs_s1, nobs_s2
 stop

end

!#################################################### readdata
!# modified on 2017.07.13 for multiple sources
!# Modified on 2017.05.18
!# Coded on May 13, 2017
subroutine readdata(g_param_inv,g_data)
implicit none
type(param_inversion),    intent(inout)    :: g_param_inv ! 2017.07.14 data_avail
type(data_vec),           intent(out)      :: g_data
character(50),allocatable,dimension(:)     :: obsfile
real(8),      allocatable,dimension(:)     :: dat ,err      ! 2017.10.15
real(8),      allocatable,dimension(:)     :: dat1,err1     ! 2017.10.15
real(8),      allocatable,dimension(:)     :: dat2,err2     ! 2017.10.15
real(8),      allocatable,dimension(:,:,:) :: dd1, ee1      ! 2017.10.15
integer(4),   allocatable,dimension(:,:,:) :: iflag1        ! 2017.10.15
logical,      allocatable,dimension(:,:,:) :: data_avail1   ! 2017.10.15
real(8),      allocatable,dimension(:,:,:) :: dd2, ee2      ! 2017.10.15
integer(4),   allocatable,dimension(:,:,:) :: iflag2        ! 2017.10.15
logical,      allocatable,dimension(:,:,:) :: data_avail2   ! 2017.10.15
real(8)    :: d1,dr,e1,er,f1,ed, errorfloor
integer(4) :: i,j,k
integer(4) :: ndat,ndat1,ndat2,nsr_inv       ! 2017.11.01
integer(4) :: nobs,nfreq
integer(4) :: nobs_s1,nobs_s2,iobs           ! 2017.07.13
integer(4) :: ndatmax,icount1,icount2        ! 2017.07.13

write(*,*) "Start read data"
!#[1]## set
 nobs     = g_param_inv%nobs     ! same as g_param 2017.07.13
 nfreq    = g_param_inv%nfreq    ! same as g_param 2017.07.13
 nsr_inv  = g_param_inv%nsr_inv  ! 2017.07.13
 errorfloor = g_param_inv%errorfloor ! [%]

!#[2]## read and gen dat,err
 allocate(    dd1(nfreq,nobs,nsr_inv)) ! 2017.10.15
 allocate(    ee1(nfreq,nobs,nsr_inv)) ! 2017.10.15
 allocate( iflag1(nfreq,nobs,nsr_inv)) ! 2017.10.15
 allocate(    dd2(nfreq,nobs,nsr_inv)) ! 2017.10.15
 allocate(    ee2(nfreq,nobs,nsr_inv)) ! 2017.10.15
 allocate( iflag2(nfreq,nobs,nsr_inv)) ! 2017.10.15

!#[2-1]## just read
do k=1,nsr_inv

 !#[Former data]## 2017.10.15
 nobs_s1 = g_param_inv%obsinfo1(k)%nobs_s
 do i=1, nobs_s1            ! 2017.07.13
  open(1,file=g_param_inv%obsinfo1(k)%obsfile(i)) ! 2017.07.13
  write(*,*) "file=",g_param_inv%obsinfo1(k)%obsfile(i)
  iobs = g_param_inv%obsinfo1(k)%obsindex(i)
  do j=1,nfreq
!   write(*,*) "j,i,k",j,i,k,"iobs=",iobs
   read(1,*) f1,iflag1(j,iobs,k),dd1(j,iobs,k),ee1(j,iobs,k)
  end do
  close(1)
 end do ! nobs_s   loop

 !#[Later data]## 2017.10.15
 nobs_s2 = g_param_inv%obsinfo2(k)%nobs_s
 do i=1, nobs_s2            ! 2017.07.13
  open(1,file=g_param_inv%obsinfo2(k)%obsfile(i)) ! 2017.07.13
  write(*,*) "file=",g_param_inv%obsinfo2(k)%obsfile(i)
  iobs = g_param_inv%obsinfo2(k)%obsindex(i)
  do j=1,nfreq
!   write(*,*) "j,i,k",j,i,k,"iobs=",iobs
   read(1,*) f1,iflag2(j,iobs,k),dd2(j,iobs,k),ee2(j,iobs,k)
  end do
  close(1)
 end do ! nobs_s   loop

end do  ! nsr_inv loop

!#[2-2]## assemble data vec
!# count ndata
ndatmax = nsr_inv * nobs * (nfreq - 1) * 2 ! 2017.11.01
allocate(dat(ndatmax), err(ndatmax))
allocate(dat1(ndatmax),err1(ndatmax))
allocate(dat2(ndatmax),err2(ndatmax))
allocate(data_avail1(nfreq,nobs,nsr_inv))  ! 2017.10.15
allocate(data_avail2(nfreq,nobs,nsr_inv))  ! 2017.10.15
data_avail1(:,:,:)=.false.                 ! 2017.07.13
data_avail2(:,:,:)=.false.                 ! 2017.07.13
icount1 = 0
icount2 = 0

do j=2,nfreq
 do k=1,nsr_inv

 !#[Former data]## 2017.10.15
 nobs_s1 = g_param_inv%obsinfo1(k)%nobs_s
 do i=1,nobs_s1
  iobs = g_param_inv%obsinfo1(k)%obsindex(i)
   if ( iflag1(j,iobs,k) .eq. 1 ) then    ! 2017.07.13
    data_avail1(1,iobs,k) = .true.        ! 2017.07.17
    data_avail1(j,iobs,k) = .true.        ! 2017.07.13
    icount1 = icount1 + 1                  ! 2017.07.13
    dr = dd1(1,iobs,k)
    d1 = dd1(j,iobs,k)
    er = ee1(1,iobs,k)
    e1 = ee1(j,iobs,k)
    ed = sqrt((e1/dr)**2.d0 + (d1/dr*er)**2.d0)
    if ( ed/abs(d1/dr) .lt. errorfloor ) ed=abs(d1/dr)*errorfloor
    dat1( icount1 ) = d1/dr ! 2017.07.13
    err1( icount1 )  = ed   ! 2017.07.13
   end if                 ! 2017.07.13
  end do

 !#[Later data]## 2017.10.15
 nobs_s2 = g_param_inv%obsinfo2(k)%nobs_s
 do i=1,nobs_s2
  iobs = g_param_inv%obsinfo2(k)%obsindex(i)
   if ( iflag2(j,iobs,k) .eq. 1 ) then    ! 2017.07.13
    data_avail2(1,iobs,k) = .true.        ! 2017.07.17
    data_avail2(j,iobs,k) = .true.        ! 2017.07.13
    icount2 = icount2 + 1                  ! 2017.07.13
    dr = dd2(1,iobs,k)
    d1 = dd2(j,iobs,k)
    er = ee2(1,iobs,k)
    e1 = ee2(j,iobs,k)
    ed = sqrt((e1/dr)**2.d0 + (d1/dr*er)**2.d0)
    if ( ed/abs(d1/dr) .lt. errorfloor ) ed=abs(d1/dr)*errorfloor
    dat2( icount2 ) = d1/dr ! 2017.07.13
    err2( icount2 ) = ed   ! 2017.07.13
   end if                 ! 2017.07.13
  end do

 end do
end do

 ndat1 = icount1           ! 2017.10.15
 ndat2 = icount2           ! 2017.10.15
 ndat  = ndat1 + ndat2     ! 2017.10.31
 dat(1:ndat1)      = dat1(1:ndat1)      ! 2017.10.15
 dat(ndat1+1:ndat) = dat2(1:ndat2)      ! 2017.10.15
 err(1:ndat1)      = err1(1:ndat1)      ! 2017.10.15
 err(ndat1+1:ndat) = err2(1:ndat2)      ! 2017.10.15
 write(*,*) "ndat1",ndat,"ndat2",ndat2,"ndat",ndat ! 2017.10.15

do i=1,ndat
 write(*,*) "dat,err=",dat(i),err(i)
end do

!#[3]## output
! g_data
 g_data%nobs  = nobs
 g_data%nfreq = nfreq
 g_data%ndat  = nsr_inv ! 2017.07.14
 g_data%ndat1 = ndat1
 g_data%ndat2 = ndat2
 g_data%ndat  = ndat
 allocate(g_data%dvec1(ndat1),g_data%error1(ndat))
 allocate(g_data%dvec2(ndat2),g_data%error2(ndat))
 allocate(g_data%dvec( ndat ),g_data%error( ndat))
 allocate(g_data%data_avail1(nfreq,nobs,nsr_inv)) ! 2017.07.13
 allocate(g_data%data_avail2(nfreq,nobs,nsr_inv)) ! 2017.07.13
 g_data%dvec1       = dat1(1:ndat1) ! 2017.07.13
 g_data%error1      = err1(1:ndat1) ! 2017.07.13
 g_data%data_avail1 = data_avail1   ! 2017.07.13
 g_data%dvec2       = dat2(1:ndat2) ! 2017.07.13
 g_data%error2      = err2(1:ndat2) ! 2017.07.13
 g_data%data_avail2 = data_avail2   ! 2017.07.13
 !
 g_data%dvec        = dat(1:ndat)   ! 2017.07.13
 g_data%error       = err(1:ndat)   ! 2017.07.13
! g_param_inv
 allocate(g_param_inv%data_avail1(nfreq,nobs,nsr_inv)) ! 2017.07.17
 allocate(g_param_inv%data_avail2(nfreq,nobs,nsr_inv)) ! 2017.07.17
 g_param_inv%data_avail1 = data_avail1  ! 2017.07.14
 g_param_inv%ndat1       = ndat1
 g_param_inv%data_avail2 = data_avail2 ! 2017.07.14
 g_param_inv%ndat2       = ndat2
 g_param_inv%ndat        = ndat

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
logical,allocatable,dimension(:,:,:) :: data_avail1  ! 2017.07.14
logical,allocatable,dimension(:,:,:) :: data_avail2  ! 2017.07.14
integer(4)                           :: i,j,k,ndat
integer(4)                           :: ndat1,ndat2  ! 2017.11.01

!#[0]## set
ndat    = g_param_inv%ndat    ! 2017.07.14
ndat1   = g_param_inv%ndat1   ! 2017.11.01
ndat2   = g_param_inv%ndat2   ! 2017.11.01
nobs    = g_param_inv%nobs    ! 2017.07.14
nfreq   = g_param_inv%nfreq   ! 2017.07.14
nsr_inv = g_param_inv%nsr_inv ! 2017.07.14
allocate(data_avail1(nfreq,nobs,nsr_inv)) ! 2017.07.14
allocate(data_avail2(nfreq,nobs,nsr_inv)) ! 2017.07.14
data_avail1 = g_param_inv%data_avail1 ! 2017.07.14
data_avail2 = g_param_inv%data_avail2 ! 2017.07.14

!#[1]## output
g_data%nobs    = nobs
g_data%nfreq   = nfreq
g_data%nsr_inv = nsr_inv     ! 2017.07.14
g_data%ndat    = ndat        ! 2017.07.14
g_data%ndat1   = ndat1       ! 2017.11.01
g_data%ndat2   = ndat2       ! 2017.11.01
allocate(g_data%dvec(ndat) ) ! 2017.07.14
allocate(g_data%error(ndat)) ! 2017.07.14
allocate(g_data%data_avail1(nfreq,nobs,nsr_inv)) ! 2017.07.14
allocate(g_data%data_avail2(nfreq,nobs,nsr_inv)) ! 2017.07.14
g_data%data_avail1 = data_avail1 ! 2017.10.31
g_data%data_avail2 = data_avail2 ! 2017.10.31

return
end

end module param_mdminv
