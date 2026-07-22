!# modified for multisource amp phase inversion on 2017.08.31
!# Generated on 2017.06.08
module param_jointinv
use modelpart
use param        ! added on May 13, 2017
use param_mt     ! 2022.10.21
implicit none

 character(2),dimension(5) :: comp=(/"Bx","By","Bz","Ex","Ey"/) ! 2018.10.04

type obs_src          ! 2017.08.31
 integer(4) :: nobs_s ! nobs for one source 2017.07.13
 character(50),allocatable,dimension(:,:)   :: ampfile ! (5,nobs)   2018.10.04
 character(50),allocatable,dimension(:,:)   :: phafile ! (5,nobs)   2017.10.04
 integer(4),   allocatable,dimension(:,:,:) :: obsindex! (2,5,nobs) 2017.10.04
end type              ! 2017.08.31

type obs_mt ! 2021.12.27
integer(4) :: nobs_mt
character(50),allocatable,dimension(:) :: impfile    ! 2021.12.27
character(50),allocatable,dimension(:) :: imperrfile ! 2021.12.27
integer(4),   allocatable,dimension(:) :: impindex ! (nobs_mt)2021.12.27
integer(4),   allocatable,dimension(:) :: errindex ! (nobs_mt)2021.12.27
end type

type param_joint  ! 2021.12.25
 integer(4)    :: ijoint         ! 1:ACTIVE, 2:MT, 3:Joint 2022.10.14
 integer(4)    :: nobs           ! same as g_param 2017.07.13
 integer(4)    :: nobs_mt        ! 2022.01.04
 integer(4)    :: nfreq          ! same as g_param 2017.07.13
 integer(4)    :: nfreq_mt       ! 2022.01.04
 integer(4)    :: nfreq_tot      ! 2022.10.18
 integer(4),allocatable,dimension(:,:)    :: index ! (2,nfreq_tot)2022.10.18
 integer(4)    :: nsr_inv        ! 2017.07.13
 integer(4)    :: ndat           !defined in readdata_ap 2017.07.14
 integer(4)    :: ndat_mt        !defined in readdata_mt 2022.01.02
 !#--------------------------------------------------------------- type of roughness term
 integer(4)    :: itype_roughness
 !# itype_roughness = 1 : Smoothest Model (SM)            2017.12.11
 !# itype_roughness = 2 : Minimum Support (MS)            2017.12.11
 !# itype_roughness = 3 : Minimum Support Gradient (MSG)  2017.12.11
 !# itype_roughness = 4 : Minimum gradient support (MGS)  2017.12.11
 real(8)      :: beta   ! when itype_roughenss >= 2       2017.12.11
 !#--------------------------------------------------------------- choice of hyperparameter
 integer(4)    :: ialphaflag
 !# ialphaflag = 1 : L-curve  (Usui et al., 2017)          2017.12.11
 !# ialphaflag = 2 : Cooling strategy with given alpha0    2017.12.11
 !# ialphaflag = 3 : General cooling strategy for data-space method (Minami et al., 2018)
 !# ialphaflag = 4 : Modified cooling method of Grayver et al (2013)
 !#                  At each iteration, rambda = [ratio of spectral]/ite is determined
 !#                  ratio is calculated from data-space equation
 !#--------------------------------------------------------------- choice of hyperparameter
 integer(4)    :: nalpha
 real(8),allocatable,dimension(:) :: alpha
 !# for Kordy et al. (2016) cooling strategy ialphaflag = 2
 real(8)       :: alpha_init     ! 2017.07.14
 !# for ialphaflag = 3 (Minami cooling) and ialphaflag = 4 (Modified Grayver) ! 2018.01.09
 real(8)       :: gamma          ! init alpha = gamma * alpha_from_Lanczos 2017.12.22
 !#-----------------------------  Upper and lower bound information 2018.01.18
 integer(4)    :: iboundflag     !0: off, 1:on, 2:Transformed model variable 2018.01.22
 real(8)       :: p_param        ! only for iboundflag = 2 2018.01.22 see Grayver et al. (2013)
 real(8)       :: logrho_upper   ! upper bound 2018.01.18
 real(8)       :: logrho_lower   ! lower bound 2018.01.18
 integer(4)    :: BMIdenomflag   ! 1: model, 2:transfered model is used for denominator in BMI generation in MSG, MS
 !#-----------------------------
 character(50) :: g_faceinfofile ! added on March 10, 2017
 !# init cond or model
 !# icondflag_ini = 0 : homo geneous input, sigma for ini should be applied as sigma_ini 2018.10.04
 !# icondflag_ini = 1 : initial cond is supplied as cond file    2018.10.04
 !# icondflag_ini = 2 : initial cond is supplied as model file   2018.10.04
 integer(4)    :: icondflag_ini  ! 2018.06.21
 real(8)       :: sigmahomo_ini  ! 2018.10.04
 character(50) :: g_initmodelconn! 2018.06.21
 character(50) :: g_initcondfile ! added on May 13, 2017
 character(50) :: g_initmodelfile! 2018.06.21
 !# ref cond or model
 integer(4)    :: icondflag_ref  ! 2018.06.21
 real(8)       :: sigmahomo_ref  ! 2018.10.04
 character(50) :: g_refmodelconn ! 2018.06.21
 character(50) :: g_refcondfile  ! 2017.07.19
 character(50) :: g_refmodelfile ! 2018.06.21
!#
 character(50) :: outputfolder   ! 2017.07.25
 integer(4),   allocatable,dimension(:)   :: srcindex   ! (nsr_inv) 2017.07.13
 type(obs_src),allocatable,dimension(:)   :: obsinfo    ! (nsr_inv) 2017.07.13
 type(obs_mt)                             :: obsmtinfo ! 2021.12.27
 integer(4)         ,dimension(5)         :: iflag_comp ! 2018.10.04
 logical,allocatable,dimension(:,:,:,:,:) :: data_avail !(2,5,nfreq,nobs,nsr_inv) 2018.10.04
 logical,allocatable,dimension(:,:,:,:)   :: data_avail_mt !(2,4,nfreq,nobs_mt) 2022.01.02
 real(8)                                  :: errorfloor_act  ! [%] added on 2021.12.27
 real(8)                                  :: errorfloor_mt   ! 2021.12.27
 integer(4)                               :: mt_imp_ap_flag  ! 0:imp, 1:amp, phase 2021.12.27

 !# output level 0:default, 1;output Jacobian  2018.06.25
 integer(4)    :: ioutlevel      ! 2018.06.25
 real(8)       :: finalrms       ! optional , if less than 0.1, reset to 1.0 2018.06.25
 !# replace flag 0: nothing, 1: replace reference model when rms is larger than previous value
 integer(4)    :: iflag_replace  ! 2018.06.26
end type

type data_vec_ap  ! added on 2017.08.31
!  | bx amp_obs1 (src1,freq1)|
!  | bx pha_obs1 (src1,freq1)|
!  | by amp_obs1 (src1,freq1)|
!  | by pha_obs1 (src1,freq1)|
!  | bz amp_obs1 (src1,freq1)|
!  | bz pha_obs1 (src1,freq1)|
!  | ex amp_obs1 (src1,freq1)|
!  | ex pha_obs1 (src1,freq1)|
!  | ey amp_obs1 (src1,freq1)|
!  | ey pha_obs1 (src1,freq1)|
!d=| bx amp_obs2 (src1,freq1)|
!  | bx pha_obs2 (src1,freq1)|
!  | ...                     |
!  | ey amp_obs2 (src1,freq1)|
!  | ey pha_obs2 (src1,freq1)|
!  |  ....                   |
!  | bx amp_obs1 (src2,freq1)|
!  | bx pha_obs1 (src2,freq1)|
!  |  ....                   |
!  | ey amp_obs1 (src2,freq1)|
!  | ey pha_obs1 (src2,freq1)|
!  |  ....                   |
!  | bx amp_obs1 (src1,freq2)|
!  | bx pha_obs1 (src1,freq2)|
!  |  ....                   |
!  | ey amp_obs2 (src2,freq2)|
!  | ey pha_obs2 (src2,freq2)|
 integer(4) :: nsr_inv   ! 2017.08.31
 integer(4) :: nobs   = 0 ! 2022.10.14
 integer(4) :: nfreq  = 0 ! 2022.10.14
 integer(4) :: ndat   = 0   ! depends on data_avail 2022.10.14
 !data_avail(2,5,nfreq,nobs,nsr_inv) 1:amp, 2:phase added on 2017.08.31; 5 comp added on 2018.10.04
 integer(4),allocatable,dimension(:,:,:,:,:) :: idata! (2,5,nfreq,nobs,nsr_inv) 2022.11.01
 logical,allocatable,dimension(:,:,:,:,:) :: data_avail ! 2018.10.04
 integer(4),         dimension(5)         :: iflag_comp ! e.g. (/0,0,1,0,0/) for ACTIVE 2018.10.04
 real(8),allocatable,dimension(:)         :: dvec  ! dvec(ndat)
 real(8),allocatable,dimension(:)         :: error ! error(ndat)
end type

type data_vec_mt ! 2021.12.27
!  | zxx real_obs1 (freq1)|
!  | zxx imag_obs1 (freq1)|
!  | zxy real_obs1 (freq1)|
!  | zxy imag_obs1 (freq1)|
!  | zyx real_obs1 (freq1)|
!  | zyx imag_obs1 (freq1)|
!  | zyy real_obs1 (freq1)|
!d=| zyy imag_obs1 (freq1)|
!  | zxx real_obs2 (freq1)|
!  | zxx imag_obs2 (freq1)|
!  | ...                  |
!  | zxx real_obs1 (freq2)|
!  | zxx imag_obs1 (freq2)|
!  | ...                  |
!  | zyy real_obsM (freqN)|
!  | zyy imag_obsM (freqN)|
integer(4) :: nobs_mt = 0  ! 2022.10.14
integer(4) :: nfreq_mt = 0 ! 2022.10.14
integer(4) :: ndat_mt = 0     ! depends on data_avail 2022.10.14
!data_avail_mt(2,4,nobs,nfreq) 1:real/img, 2:zxx->zyy added on 2017.08.31; 5 comp added on 2018.10.04
integer(4),allocatable,dimension(:,:,:,:) :: idata_mt      !(2,4,nobs_mt,nfreq_mt) 2022.12.12
logical,   allocatable,dimension(:,:,:,:) :: data_avail_mt !(2,4,nobs_mt,nfreq_mt) 2021.12.27
real(8),allocatable,dimension(:)       :: dvec_mt  ! dvec(ndat_mt)
real(8),allocatable,dimension(:)       :: error_mt ! error(ndat_mt)
end type ! 2021.12.27

contains

!############################################# setnec 2022.10.14
subroutine setnec(ijoint,ACT,MT)
implicit none
 integer(4),intent(in) :: ijoint
 logical,   intent(out) :: ACT, MT
 ACT=.false.
 MT=.false.
 if ( ijoint == 1 ) ACT = .true. ! only ACTIVE
 if ( ijoint == 2 ) MT  = .true.  ! only MT
 if ( ijoint == 3 ) then
   ACT = .true. ! ACTIVE and MT
   MT  = .true. 
   end if
return
end
!################################################## readparam
!# 2017.08.31 modified for multisource inversion
!# 2017.06.08
subroutine readparaJOINTINV(ijoint,g_param_joint,g_modelpara,g_param,sparam,g_param_mt,g_data_ap,g_data_mt) !2022.10.14
 use caltime      ! 2017.12.22
 implicit none
 type(param_forward),    intent(in)  :: g_param
 type(param_forward_mt), intent(in)  :: g_param_mt     ! 2022.10.21
 type(param_source),     intent(in)  :: sparam         ! 2018.06.26
 type(param_joint),      intent(out) :: g_param_joint  ! 2017.08.31
 type(modelpara),        intent(out) :: g_modelpara
 type(data_vec_ap),      intent(out) :: g_data_ap
 type(data_vec_mt),      intent(out) :: g_data_mt      ! 2021.12.27 
 character(50)                       :: a              ! 2017.06.06
 integer(4),             intent(in)  :: ijoint  ! 1:ACTIVE, 2:MT, 3:ACTIVE + MT 2022.10.14
 integer(4)                          :: input=10, nobs_s,i,j,nsr_inv ! 2017.08.31
 integer(4)                          :: icomp          ! 2018.10.04
 integer(4)                          :: nobs_mt
 integer(4)                          :: ikeep          ! 2021.10.12
 character(100)              :: paramfile ! 2020.09.29
 integer(4),       parameter :: n = 1000  ! 2020.09.28
 character(200),dimension(n) :: lines     ! 2020.09.28
 type(watch) :: t_watch ! 2017.12.22

 call watchstart(t_watch) ! 2017.12.22

!#[-1]##
 g_param_joint%ijoint = ijoint ! 1;ACTIVE, 2:MT, 3:Joint

!#[0]## read inversion parameter file 2020.09.29
write(*,*) "" ! 2020.09.29
write(*,*) "<Please input the inversion parameter file>" ! 2020.09.29
read(*,'(a)') paramfile           ! 2020.09.29
ikeep = 1 !2021.10.04
call readcontrolfile(paramfile,n,lines,ikeep) ! 2021.10.04 see src/common/m_param.f90

open(input,file="tmp.ctl")

!#[1]## read files 2017.08.31
! read(input,*)  ! for header
 write(*,*) "" ! 2020.09.29
 write(*,*) "<Input path for faceinfo file>" ! 2020.09.29
 read(input,'(20x,a)') g_param_joint%g_faceinfofile       ! face info
 write(*,*) "faceinfo file : ",g_param_joint%g_faceinfofile         ! 2018.10.05
 !#
 !# ref cond/model file
 write(*,*) "" ! 2020.09.29
 write(*,*) "< reference conductivity model >"
 write(*,40) " < Input icondflag_ref: 1 (condfile), or 2 (modelconnect and model file) >" ! 2018.03.17
 read(input,11) g_param_joint%icondflag_ref                 ! 2018.06.21
 write(*,42) " icondflag_ref =",g_param_joint%icondflag_ref   ! 2020.09.29
 if     ( g_param_joint%icondflag_ref .eq. 0 ) then                  ! 2018.10.04
  write(*,*) "" ! 2020.09.29
  write(*,*) "< input homogeneous conductivity for reference model >"    ! 2018.10.04
  read(input,12) g_param_joint%sigmahomo_ref              ! 2018.10.04
  write(*,*) "homogeneous conductivity : ",g_param_joint%sigmahomo_ref," [S/m]"    ! 2018.10.04
 elseif ( g_param_joint%icondflag_ref .eq. 1 ) then                  ! 2018.10.04
  write(*,*) "" ! 2020.09.29
  write(*,*) "< input path for reference conductivity structure file >"  ! 2018.06.21
  read(input,10) g_param_joint%g_refcondfile ! ref cond file    2018.06.21
  write(*,*) "refcondfile : ",g_param_joint%g_refcondfile              ! 2018.06.21
 elseif (g_param_joint%icondflag_ref .eq. 2 ) then                   ! 2018.06.21
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input modelconnect file and modelfile >"                 ! 2018.06.21
  read(input,10) g_param_joint%g_refmodelconn ! ref cond file   2018.06.21
  read(input,10) g_param_joint%g_refmodelfile ! ref cond file   2018.06.21
  write(*,*) "refmodelconnect : ",g_param_joint%g_refmodelconn          ! 2018.06.21
  write(*,*) "refmodelfile    : ",g_param_joint%g_refmodelfile          ! 2018.06.21
 end if
 !# init cond/model file
 write(*,*) "" ! 2020.09.29
 write(*,*) "< initial conductivity model >"
 write(*,40) " < Input icondflag_ini: 1 (condfile), or 2 (modelconnect and model file) >" ! 2020.09.29
 read(input,11) g_param_joint%icondflag_ini                 ! 2018.06.21
 write(*,42) " icondflag_ini =",g_param_joint%icondflag_ini   ! 2020.09.29
 if     ( g_param_joint%icondflag_ini .eq. 0 ) then                  ! 2018.10.04
  write(*,*) "" ! 2020.09.29
  write(*,*) "< input homogeneous conductivity for initial model >"      ! 2018.10.04
  read(input,12) g_param_joint%sigmahomo_ini              ! 2018.10.04
  write(*,*) "sigmahomo_ini : ",g_param_joint%sigmahomo_ini,"[S/m]"    ! 2018.10.04
 elseif ( g_param_joint%icondflag_ini .eq. 1 ) then                  ! 2018.06.21
  write(*,*) "" ! 2020.09.29
  write(*,*) "< input path for initial conductivity structure file >"
  read(input,10) g_param_joint%g_initcondfile ! init cond file
  write(*,41) " initcondfile : ",g_param_joint%g_initcondfile
 elseif (g_param_joint%icondflag_ini .eq. 2 ) then                   ! 2018.06.21
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input modelconnect file and modelfile >"               ! 2018.06.21
  read(input,10) g_param_joint%g_initmodelconn ! ref cond file  2018.06.21
  read(input,10) g_param_joint%g_initmodelfile ! ref cond file  2018.06.21
  write(*,*) "inimodelconnect",g_param_joint%g_initmodelconn         ! 2018.06.21
  write(*,*) "inimodelfile",   g_param_joint%g_initmodelfile         ! 2018.06.21
 end if
!#
 write(*,*) "" ! 2020.09.29
 write(*,*) "< input output folder for inversion result >"    ! 2017.07.25
 read(input,10) g_param_joint%outputfolder         ! 2017.07.25
 write(*,41) " output folder for inversion: ",g_param_joint%outputfolder!2017.07.25

 write(*,*) "" ! 2020.09.29
 write(*,40) " < Input roughness type >"                              ! 2020.09.29
 write(*,40) " < 1:Smoothest,2:minimum support,3:Minimum support gradient,4:Minimum gradient support>"
 read(input,11) g_param_joint%itype_roughness    ! 2017.12.11
 write(*,'(a,i3)')   " itype_roughness =",g_param_joint%itype_roughness ! 2020.09.17
 if ( g_param_joint%itype_roughness .ge. 2 ) then
  write(*,*) "" ! 2020.09.29
  write(*,*) "<Input beta for roughenss>"
  read(input,12) g_param_joint%beta
  write(*,*) "Beta =",g_param_joint%beta
 end if
 write(*,*) "" ! 2020.09.29
 write(*,*) "< Input ialphaflag >"
 write(*,40) " < 1 for L-curve strategy; 2 for cooling w alpha, 3: cooling (Minami) (ialphaflag) >" ! 2017.12.11
 read(input,11) g_param_joint%ialphaflag
 write(*,'(a,i3)') " ialphaflag  =",g_param_joint%ialphaflag ! 2020.09.17

!# L curve strategy # 2017.07.19
 if ( g_param_joint%ialphaflag .eq. 1) then! L-curve strategy
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input # of alpha >"
  read(input,11) g_param_joint%nalpha
  allocate(g_param_joint%alpha(g_param_joint%nalpha))
  do i=1,g_param_joint%nalpha
   read(input,12)g_param_joint%alpha(i)
   write(*,*) i,"alpha",g_param_joint%alpha(i)
  end do

 !# cooling strategy # 2017.07.19
 else if (g_param_joint%ialphaflag .eq. 2 ) then     ! 2017.08.31 cooling strategy with given alpha
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input initial alpha [real] >"     ! 2018.06.25
  read(input,12) g_param_joint%alpha_init ! 2017.07.14
  write(*,'(a,f8.3)') " initial alpha [real] :",g_param_joint%alpha_init ! 2020.09.17
  write(*,*) "" ! 2020.09.29
  write(*,40) " < Enter iflag_replace 0:nothing, 1:replace ref model when rms is larger than previous> "! 2020.09.29
  read(input,11)  g_param_joint%iflag_replace    ! 2018.06.26
  write(*,42) " iflag_replace  = ", g_param_joint%iflag_replace ! 2018.06.26
 else if (g_param_joint%ialphaflag .eq. 3 ) then ! 2017.12.13 cooling with Lanczos
 ! gamma is parameter : init alpha=gamma*alpha_from_Lanczos
 write(*,*) "" ! 2020.09.29
 write(*,*) "< Input the parameter gamma (see Grayver et al., 2013)> "
 read(input,12) g_param_joint%gamma   ! 2017.12.22
  write(*,40) "Initial alpha is determined from ratio of spectral radii (Minami cooling)"
 else if (g_param_joint%ialphaflag .eq. 4 ) then ! 2017.12.22
  write(*,*) "Modified version of Grayver (2013) is adopted!"!2018.01.09
  write(*,*) "" ! 2020.09.29
  write(*,*)  "<Input the parameter gamma (see Grayver et al., 2013)>"
 read(input,12) g_param_joint%gamma   ! 2018.01.09
 else
  write(*,*) "g_param_joint%ialphaflag should be 1 or 2",g_param_joint%ialphaflag
  stop
 end if
!#--------------------------------------------iboundflag: model bound! 2018.01.18
write(*,*) "" ! 2020.09.29
write(*,*) "< Input iboundflag 0:off, 1:on, 2: Transformed variable >" ! 2018.01.18
read(input,11) g_param_joint%iboundflag        ! 2018.01.18
write(*,42) "iboundflag =",g_param_joint%iboundflag      ! 2018.06.25
if ( g_param_joint%iboundflag .eq. 1 ) then             ! 2018.01.18
 read(input,12) g_param_joint%logrho_upper   ! 2018.01.18 log10(Ohm.m)
 read(input,12) g_param_joint%logrho_lower   ! 2018.01.18 log10(Ohm.m)
end if
if ( g_param_joint%iboundflag .eq. 2 ) then             ! 2018.01.22 when Transform
 read(input,12)  g_param_joint%p_param       ! 2018.01.22
 read(input,12)  g_param_joint%logrho_upper  ! 2018.01.22 log10(Ohm.m)
 read(input,12)  g_param_joint%logrho_lower  ! 2018.01.22 log10(Ohm.m)
 write(*,*) "P_param=",     g_param_joint%p_param       ! 2018.01.23
 write(*,*) "logrho_upper=",g_param_joint%logrho_upper  ! 2018.01.23
 write(*,*) "logrho_lower=",g_param_joint%logrho_lower  ! 2018.01.23
 write(*,*) "" ! 2020.09.29
 write(*,40) "<Input BMIdenomflag = 1 for model and 2 transformed model for denominator of BMI>" ! 2018.02.05
 read(input,11)  g_param_joint%BMIdenomflag    ! 2018.02.05
 write(*,*) "BMIdenomflag : ",g_param_joint%BMIdenomflag! 2018.02.05
end if
!# alpha input end # 2017.08.31

!=================================================================== ACTIVE information
if ( ijoint == 1 .or. ijoint == 3 ) then ! ACTIVE or Joint inversion
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input # of sources for inversion >" ! 2017.08.31
 !################################################################ read nsr_inv
  read(input,11) nsr_inv                 ! 2017.08.31
  g_param_joint%nsr_inv = nsr_inv               ! 2017.08.31
  write(*,'(a,i3)') " nsr_inv =",g_param_joint%nsr_inv   ! 2018.06.25
  allocate(g_param_joint%srcindex(nsr_inv))     ! 2017.08.31
  allocate(g_param_joint%obsinfo(nsr_inv) )     ! 2017.08.13

 !################################################################ read iflag_comp(5)
  write(*,*) "" ! 2020.09.29
  write(*,40) " <Input 0 or 1 for use of Bx,By,Bz,Ex,Ey component as (5i2)>"   ! 2018.10.04
  write(*,40) " <When 1 is chosen for a component, you should supply data file for all src & site>"
  read(input,'(20x,5i2)') g_param_joint%iflag_comp(1:5)                        ! 2018.10.04

 !################################################################  read  errorfloor_act
  write(*,*) "" ! 2021.12.27
  write(*,*) "< input error floor for active data [0-1] >" ! 2021.12.27
  read(input,12) g_param_joint%errorfloor_act ! error floor
  write(*,'(a,f8.3)') " errorfloor_act      =",g_param_joint%errorfloor_act ! 2021.12.22

 !#[2]## read data files ! 2018.10.04
 g_param_joint%nobs     = g_param%nobs
 g_param_joint%nfreq    = g_param%nfreq
 g_param_joint%nfreq_mt = g_param_mt%nfreq ! 2022.10.21

 do i=1,nsr_inv ! 2017.08.31 #################################### nsr_inv loop start
  !# srcindex(i)
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input source index specified in forward ctl file >"  ! 2020.09.29
 !############################################################### read src index
  read(input,11) g_param_joint%srcindex(i)
  write(*,*) "srcindex=",g_param_joint%srcindex(i)
  !# nobs_s
  write(*,*) "" ! 2020.09.29
  write(*,*) "< Input # of observations for inversion for current source >" ! 2020.09.29
  read(input,11) nobs_s
  write(*,*) "Source",i,"nobs_s :",nobs_s                  ! 2018.06.25
  if ( g_param_joint%nobs .lt. nobs_s) goto 100
   g_param_joint%obsinfo(i)%nobs_s = nobs_s                 ! 2017.07.13
   allocate(g_param_joint%obsinfo(i)%ampfile(5,nobs_s))     ! 2018.10.04
   allocate(g_param_joint%obsinfo(i)%phafile(5,nobs_s))     ! 2018.10.04
   allocate(g_param_joint%obsinfo(i)%obsindex(2,5,nobs_s))  ! 2018.10.04

   do icomp = 1,5                                           ! 2018.10.04
     if ( g_param_joint%iflag_comp(icomp) .eq. 1 ) then     ! 2018.10.04 for comp to be used
       !# number of amp and phase files should be the same, i.e. nobs_s.
       do j=1,nobs_s !  amplitude files for i-th source 2017.08.31
         read(input,'(20x,i5,a)') g_param_joint%obsinfo(i)%obsindex(1,icomp,j),&
				 &   g_param_joint%obsinfo(i)%ampfile(icomp,j) ! 2018.10.04
       end do ! nobs_s loop
       do j=1,nobs_s !  phase     files for i-th source 2017.08.31
         read(input,'(20x,i5,a)') g_param_joint%obsinfo(i)%obsindex(2,icomp,j),&
				&   g_param_joint%obsinfo(i)%phafile(icomp,j) ! 2018.10.04
       end do ! nobs_s loop
     end if
   end do ! component loop 2018.10.04

 end do   ! nsr_inv loop
 call readdata_ap(g_param,sparam,g_param_joint,g_data_ap) ! 2018.10.04 read active data

end if ! read ACTIVE data end ! 2022.10.14 

!#[3]## read mt parameters and mt data 2021.12.27
if ( ijoint == 2 .or. ijoint == 3 ) then ! MT or Joint inversion 2022.10.14
!##[mt_imp_ap_flag]
 write(*,*) "" ! 2021.12.27
 write(*,*) "< input mt_imp_ap_flag for mt data [0:imp, 1:amp, phase] >" ! 2021.12.27
 read(input,12) g_param_joint%mt_imp_ap_flag ! error floor
 write(*,'(a,i10)') " mt_imp_ap_flag      =",g_param_joint%mt_imp_ap_flag ! 2021.12.22

 !##[errorfloor_mt]
 write(*,*) "" ! 2021.12.27
 write(*,*) "< input error floor for active data [0-1] >" ! 2021.12.27
 read(input,12) g_param_joint%errorfloor_mt ! error floor
 write(*,'(a,f8.3)') " errorfloor_mt      =",g_param_joint%errorfloor_mt ! 2021.12.22

 !##[nobs_mt]
 write(*,*) "" ! 2021.12.27
 write(*,*) "< input # of mt observatories for inv >" ! 2021.12.27
 read(input,12) g_param_joint%nobs_mt ! error floor
 write(*,'(a,i5)') " nobs_mt      =",g_param_joint%nobs_mt ! 2021.12.22
 nobs_mt = g_param_joint%nobs_mt
 g_param_joint%obsmtinfo%nobs_mt = nobs_mt    ! 2021.12.27
 allocate(g_param_joint%obsmtinfo%impfile(   nobs_mt)) ! 2021.12.27
 allocate(g_param_joint%obsmtinfo%imperrfile(nobs_mt)) ! 2021.12.27
 allocate(g_param_joint%obsmtinfo%impindex(  nobs_mt)) ! 2021.12.27
 allocate(g_param_joint%obsmtinfo%errindex(  nobs_mt)) ! 2021.12.27

 !##[ mt file names]
 do i=1,nobs_mt 
   read(input,'(20x,i5,a)') g_param_joint%obsmtinfo%impindex(i),&
   &              g_param_joint%obsmtinfo%impfile(i) 
   write(*,*) g_param_joint%obsmtinfo%impfile(i)
  end do
 do i=1,nobs_mt 
   read(input,'(20x,i5,a)') g_param_joint%obsmtinfo%errindex(i),&
   &              g_param_joint%obsmtinfo%imperrfile(i) 
   write(*,*) g_param_joint%obsmtinfo%imperrfile(i) 
 end do
 
 !#[4]## read mt data
   call readdata_mt(g_param_joint,g_param_mt,g_data_mt)! see below 2022.12.12

end if ! MT data read end ! 2022.10.14

!#[5]## read model partitions
 CALL readmodelpara(input,g_modelpara) !see ../common/m_param.f90

!#[4]## read output level  2018.06.25
 g_param_joint%ioutlevel=0
 write(*,*) "" ! 2020.09.29
 write(*,*) "< Input 0 for default or 1 for output Jacobian matrix >"
 read(input,11)  g_param_joint%ioutlevel
 write(*,'(a,i3)') " ioutlevel =", g_param_joint%ioutlevel

 write(*,*) "" ! 2020.09.29
 write(*,*) "< Enter final rms for convergence >"            ! 2018.06.25
 read(input,12) g_param_joint%finalrms        ! 2018.06.25
 write(*,'(a,f8.3)') " finalrms  =",  g_param_joint%finalrms        ! 2018.06.25

 !# check

 call watchstop(t_watch) ! 2017.12.22
 write(*,'(a,f9.4,a)') " ### READPARAINVAP END!! ### Time=",t_watch%time," [min]" ! 2017.12.22
 write(*,*) ! 2020.09.29
return

100 continue
 write(*,*) "GEGEGE # of nobs_s should be smaller than nobs specified for forward"
 write(*,*) "nobs",g_param_joint%nobs,"nobs_s",nobs_s
 stop

10 format(20x,a)
11 format(20x,i10)
12 format(20x,g15.7)
21 format(20x,2i10)
22 format(20x,2g12.5)
32 format(20x,3g12.5)
40 format(a)
41 format(a,a)
42 format(a,i3)

end subroutine
!#################################################### readdata_mt
subroutine readdata_mt(g_param_joint,g_param_mt,g_data_mt) ! 2022.12.12
use constants    ! dmu,pi,r2d 2022.12.12
implicit none
type(param_joint),     intent(inout)  :: g_param_joint ! set data_avail_mt,ndat_mt
type(data_vec_mt),     intent(out)    :: g_data_mt 
type(param_forward_mt),intent(in)     :: g_param_mt    ! 2022.12.12
integer(4) :: nfreq_mt,nobs_mt,i,j,k,l
integer(4) :: i1,i2   ! 2022.12.12
integer(4) :: iavail(4),ii,ndat_mt
real(8)    :: freq,a(8),e(8),ra(4),rp(4),ea(4),ep(4),coef,omega! 2022.12.12
real(8)    :: errorfloor_mt,assq     ! 2022.12.12
complex(8) :: zxx,zxy,zyx,zyy,zz(4),iunit=(0.d0,1.d0)  ! 2022.12.12
character(100)                            :: filename1       ! 2022.01.22
character(100)                            :: filename2       ! 2022.01.22
real(8),   allocatable,dimension(:)       :: dvec_mt,err_mt
logical,   allocatable,dimension(:,:,:,:) :: data_avail
integer(4),allocatable,dimension(:,:,:,:) :: idata_mt    ! 2022.12.12
real(8),   allocatable,dimension(:,:,:,:) :: data,err
character(3)                              :: z(4)
character(4)                              :: realimag(2)
character(50)                             :: head,site    ! 2022.12.12


!#[1]##
z(1:4)=(/"zxx","zxy","zyx","zyy"/)
realimag(1:2)=(/"real","imag"/)
nfreq_mt=g_param_joint%nfreq_mt             ! 2022.10.21
nobs_mt=g_param_joint%nobs_mt
errorfloor_mt = g_param_joint%errorfloor_mt ! 2022.12.12
allocate(data_avail(2,4,nobs_mt,nfreq_mt))
allocate(  idata_mt(2,4,nobs_mt,nfreq_mt))  ! 2022.12.12
allocate(      data(2,4,nobs_mt,nfreq_mt))
allocate(       err(2,4,nobs_mt,nfreq_mt))
data_avail=.false.

!#[2]## read and gen data vector
!#[2-1] ## read    2022.12.12
do i=1,nobs_mt

  open(1,file=trim(g_param_joint%obsmtinfo%impfile(   i)),status='old')
  open(2,file=trim(g_param_joint%obsmtinfo%imperrfile(i)),status='old')

  do j=1,nfreq_mt
    read(1,*) freq,a(1:8)
    read(2,*) freq,iavail(1:4),e(1:8)
    if ( iavail(1) ) data_avail(1:2,1,i,j)=.true. ! zxx available
    if ( iavail(2) ) data_avail(1:2,2,i,j)=.true. ! zxy available
    if ( iavail(3) ) data_avail(1:2,3,i,j)=.true. ! zyz available
    if ( iavail(4) ) data_avail(1:2,4,i,j)=.true. ! zyy available
    data(1:2,1,i,j)=a(1:2) ! zxx
    data(1:2,2,i,j)=a(3:4) ! zxy
    data(1:2,3,i,j)=a(5:6) ! zyx
    data(1:2,4,i,j)=a(7:8) ! zyy
    err(1:2,1,i,j)=e(1:2) ! zxx
    err(1:2,2,i,j)=e(3:4) ! zyz
    err(1:2,3,i,j)=e(5:6) ! zyx
    err(1:2,4,i,j)=e(7:8) ! zyy
    !# reflect errorfloor_mt  ! 2022.12.12
      zxx = a(1) + iunit*a(2) ! 2022.12.12
      zxy = a(3) + iunit*a(4) ! 2022.12.12
      zyx = a(5) + iunit*a(6) ! 2022.12.12
      zyy = a(7) + iunit*a(8) ! 2022.12.12
      assq = abs(sqrt((zxx**2. + zxy**2. + zyx**2. + zyy**2.)/4.)) ! 2022.12.12
      do k=1,4                                                      ! 2022.12.12
       if ( err(1,k,i,j) .lt. assq*errorfloor_mt ) err(1:2,k,i,j)= assq*errorfloor_mt ! 2022.12.12
      end do ! 2022.12.12
  end do
  close(1)
  close(2)
end do

!#[2-2]## assemble data vector  2022.12.12
allocate(dvec_mt(2*4*nobs_mt*nfreq_mt))
allocate( err_mt(2*4*nobs_mt*nfreq_mt))
ii       = 0
idata_mt = 0      ! 2022.12.12
coef = 1.d0          ! Case for input in [mV/km/nT] for Nakaya 2022.12.13
!coef = 1./dmu*1.e-3 ! case for input in [V/A]=[Ohm],coef is used for conversion from [Ohm]=[V/A] -> [mV/km/nT]
open(1,file=trim(g_param_joint%outputfolder)//"dvec_mt.dat") ! 2022.12.05
do i=1,nfreq_mt
 do j=1,nobs_mt
  do k=1,4 ! zxx -> zyy
   do l=1,2 ! real, imag
    if( data_avail(l,k,j,i) ) then
     ii=ii+1
     idata_mt(l,k,j,i) = ii           ! 2022.12.12
     dvec_mt(ii) = coef*data(l,k,j,i) ! Impedance [Ohm] -> [mV/km/nT] 2022.12.12
      err_mt(ii) = coef* err(l,k,j,i) ! Impedance [Ohm] -> [mV/km/nT] 2022.12.12
      ! write(*,'(i5,2x,a,2x,a,2f15.7,2(a,i3))') ii,z(k),realimag(l),dvec_mt(ii),err_mt(ii),&
      ! & " obs#",j," freq#",i !
      write(1,'(i5,2x,a,2x,a,2f15.7,2(a,i3))') ii,z(k),realimag(l),dvec_mt(ii),err_mt(ii),&
      & " obs#",j," freq#",i
    end if
    end do
   end do
  end do
 end do
close(1) ! 2022.12.05
ndat_mt=ii ! set ndat

!#[2-3]## out obsfiles with adopted errors 2022.12.12
head       = g_param_joint%outputfolder     ! 2017.07.25
do j=1,nobs_mt
  site = g_param_mt%obsname(j)
  filename1 = trim(head)//trim(site)//"_MTIN.dat"     ! 2022.01.02
  filename2 = trim(head)//trim(site)//"_MT_impIN.dat" ! 2022.01.02
  open(1,file=filename1,status='replace')
  open(2,file=filename2,status='replace')
  do i=1,nfreq_mt
    freq = g_param_mt%freq(i)  ! 2022.12.12
    zz(1:4) = (99999., 99999.) ! as dummy
    e(1:8)  = 999999.          ! as dummy
    ra(:)   = 999999. ; rp(:)=999999. ! as dummy
    ea(:)   = 0. ; ep(:)=0.           ! as dummy
    omega = 2.*pi*freq
    coef = dmu/omega*1.d+6
    do k=1,4
      if (data_avail(1,k,j,i)) then
       i1 = idata_mt(1,k,j,i)
       i2 = idata_mt(2,k,j,i)
       zz(k) = dvec_mt(i1) + iunit*dvec_mt(i2) ! = [mV/km]/[nT] 2022.12.12
       e(2*k-1) = err_mt(i1)
       e(2*k  ) = err_mt(i2)
       ra(k) = coef*(abs(zz(k)))**2.
       rp(k) = atan2(imag(zz(k)),real(zz(k)))*r2d
       ea(k) = 2.*coef*abs(zz(k))*max(e(2*k-1), e(2*k))
       ep(k) = 90.
       if ( ra(k) .ge. ea(k) ) ep(k)=asin(ea(k)/(ra(k)+1.e-5))*r2d
      end if   
    end do
    write(1,'(17g15.7)') freq,(ra(l),rp(l),l=1,4),(ea(l),ep(l),l=1,4)
    write(2,'(17g15.7)') freq,(real(zz(l)),imag(zz(l)),l=1,4),(e(2*k-1),e(2*k),k=1,4)
  end do
  close(1)
  close(2)
end do

!# set output
 g_data_mt%ndat_mt  = ndat_mt
 g_data_mt%nobs_mt  = nobs_mt
 g_data_mt%nfreq_mt = nfreq_mt
 allocate( g_data_mt%data_avail_mt(2,4,nobs_mt,nfreq_mt))
 g_data_mt%data_avail_mt = data_avail
 allocate(g_data_mt%dvec_mt(ndat_mt))
 allocate(g_data_mt%error_mt( ndat_mt))
 g_data_mt%dvec_mt  = dvec_mt(1:ndat_mt)
 g_data_mt%error_mt = err_mt(1:ndat_mt)

! g_param_joint 
 allocate(g_param_joint%data_avail_mt(2,4,nobs_mt,nfreq_mt)) ! 2022.01.02
 g_param_joint%data_avail_mt = data_avail  ! 2022.01.02
 g_param_joint%ndat_mt       = ndat_mt     ! 2022.01.02

return
end

  !#################################################### readdata_ap
!# Modified on 2018.10.04 for multisource amp phase inversion
subroutine readdata_ap(g_param,sparam,g_param_joint,g_data_ap)
implicit none
type(param_forward),      intent(in)           :: g_param         ! 2018.06.26
type(param_source),       intent(in)           :: sparam          ! 2018.06.26
type(param_joint),  intent(inout)        :: g_param_joint   ! 2017.08.31
type(data_vec_ap),        intent(out)          :: g_data_ap       ! 2017.08.31
character(50),allocatable,dimension(:)         :: ampfile, phafile! 2017.08.31
real(8),      allocatable,dimension(:)         :: dat,err
real(8),      allocatable,dimension(:,:,:,:,:) :: dd, ee      ! 2017.10.04
integer(4),   allocatable,dimension(:,:,:,:,:) :: iflag       ! 2018.10.04
logical,      allocatable,dimension(:,:,:,:,:) :: data_avail  ! 2018.10.04
integer(4),   allocatable,dimension(:,:,:,:,:) :: idata       ! 2018.10.04
real(8)       :: d1,dr,e1,er,f1,ed, errorfloor
integer(4)    :: i,j,k,l,ndat,nsr_inv,icomp ! 2018.10.04
integer(4)    :: nobs,nfreq
integer(4)    :: nobs_s,iobs,iobs1,iobs2    ! 2017.08.31
integer(4)    :: ndatmax,icount             ! 2017.07.13
real(8)       :: pi,r2d
!# for output
character(50) :: ampname, phaname           ! 2018.06.26
character(50) :: head, site, sour           ! 2018.06.26
integer(4)    :: nhead,nsite,nsour          ! 2018.06.26
character(50) :: apfile        ! 2018.06.26
character(3)  :: ap(2)         ! 2018.06.26
integer(4)    :: ii            ! 2018.06.26
character(2)  :: num           ! 2018.06.26
real(8),   allocatable,dimension(:)    :: freq      ! 2018.06.26
integer(4),allocatable,dimension(:)    :: srcindex  ! 2018.06.26
real(8)       :: epfloor       ! 2018.06.26
integer(4)    :: iflag_comp(5) ! 2018.10.04

!#[1]## set
 nobs       = g_param_joint%nobs     ! same as g_param 2017.08.31
 nfreq      = g_param_joint%nfreq    ! same as g_param 2017.08.31
 nsr_inv    = g_param_joint%nsr_inv  ! 2017.08.31
 errorfloor = g_param_joint%errorfloor_act ! [%] 2021.12.27
 pi         = 4.d0 *datan(1.d0)        ! 2017.08.31
 r2d        = 180.d0/pi                ! 2017.08.31
 allocate(freq(nfreq))                       ! 2018.06.26
 allocate(srcindex(nsr_inv))                 ! 2018.06.26
 freq       = g_param%freq                   ! 2018.06.26
 srcindex   = g_param_joint%srcindex         ! 2018.06.26
 iflag_comp = g_param_joint%iflag_comp       ! 2018.10.04

!#[2]## read and gen dat,err
 allocate(         dd(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.04
 allocate(         ee(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.04
 allocate(      iflag(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.04
 allocate( data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.04
 allocate(      idata(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.04
 ndatmax = 2 * 5 * nsr_inv * nobs * nfreq      ! 2018.10.04
 allocate( dat(ndatmax),err(ndatmax)       )   ! 2017.08.31
 data_avail(:,:,:,:,:)=.false.                 ! 2017.08.31

!#[2-1]## just read ! 2017.08.31
do k=1,nsr_inv
 nobs_s = g_param_joint%obsinfo(k)%nobs_s
 do icomp = 1,5            ! 2018.10.04
  if ( iflag_comp(icomp) .eq. 0 ) cycle        ! 2018.10.04 read only components with iflag_comp = 1
 do i=1, nobs_s            ! 2017.07.13
  open(1,file=g_param_joint%obsinfo(k)%ampfile(icomp,i)) ! 2018.10.04
  open(2,file=g_param_joint%obsinfo(k)%phafile(icomp,i)) ! 2018.10.04
!  write(*,*) "file=",g_param_joint%obsinfo(k)%ampfile(i)
!  write(*,*) "file=",g_param_joint%obsinfo(k)%phafile(i)
  iobs1 = g_param_joint%obsinfo(k)%obsindex(1,icomp,i)   ! 2018.10.04
  iobs2 = g_param_joint%obsinfo(k)%obsindex(2,icomp,i)   ! 2018.10.04
  do j=1,nfreq ! amp  ! 2018.06.26
!   write(*,*) "j,i,k",j,i,k,"iobs=",iobs
   read(1,*) f1,iflag(1,icomp,j,iobs1,k),dd(1,icomp,j,iobs1,k),ee(1,icomp,j,iobs1,k)!amp,2017.08.31
!   write(*,'(f15.7,i5,2f15.7)') f1,iflag(1,icomp,j,iobs1,k),dd(1,icomp,j,iobs1,k),ee(1,icomp,j,iobs1,k)!amp,2017.08.31
  end do
  do j=1,nfreq ! phase ! 2018.06.26
   read(2,*) f1,iflag(2,icomp,j,iobs2,k),dd(2,icomp,j,iobs2,k),ee(2,icomp,j,iobs2,k)!pha,2017.08.31
!   write(*,'(f15.7,i5,2f15.7)') f1,iflag(2,icomp,j,iobs2,k),dd(2,icomp,j,iobs2,k),ee(2,icomp,j,iobs2,k)!pha,2017.08.31
  end do
  close(1) ; close(2) ! 2017.08.31
 end do ! nobs_s   loop
 end do ! comp loop     2018.10.04
end do  ! nsr_inv loop

!#[2-2]## assemble data vec 2018.10.04
!# count ndata
icount = 0
idata  = 0 ! 2018.06.26
do j=1,nfreq ! 2017.08.31
 do k=1,nsr_inv
 nobs_s = g_param_joint%obsinfo(k)%nobs_s
 do i=1,nobs_s
  do icomp = 1,5                          ! 2018.10.04
   if ( iflag_comp(icomp) .eq. 0 ) cycle  ! 2018.10.04
   do l=1,2 ! 1 for amp, 2 for phase      ! 2017.08.31
   iobs = g_param_joint%obsinfo(k)%obsindex(l,icomp,i) ! 2018.10.04
   if ( iflag(l,icomp,j,iobs,k) .eq. 1 ) then   ! 2018.10.04
    data_avail(l,icomp,j,iobs,k) = .true.       ! 2018.10.04
    icount = icount + 1                   ! 2017.07.13
    idata(l,icomp,j,iobs,k)  = icount     ! 2018.10.04
    d1 = dd(l,icomp,j,iobs,k)             ! 2018.10.04
    e1 = ee(l,icomp,j,iobs,k)             ! 2018.10.04
    !# reflect errorfloor
    ! log10(a+/-b) = log10(a) +/- b/a/log10 2017.08.31
    if ( l .eq. 1) then               ! amp 2017.08.31
     dat( icount ) = log10(d1)        !     2017.08.31
     err( icount ) = e1/d1/log(10.)   ! 2017.08.31
     ! e1 /d1 < errorfloor then err = d1*errorfloor/d1/log(10)
     if ( e1/d1 .lt. errorfloor ) err(icount) = errorfloor/log(10.)
    else if ( l .eq. 2) then ! phase
     ! This program uses only error for amp ! 2017.08.31
     dat( icount ) = d1      ! 2017.08.31
     !# error for phase is asin(err_log10(amp)*amp*log(10))*d2r because err_log10(amp)=err/amp/log(10)
     !# epfloor = asin(d1*errorfloor/d1)*r2d
!     err( icount ) = asin(err(icount-1)*log(10.))*r2d ! 2017.08.31 commented out on 2018.06.26
     epfloor       = asin(errorfloor)*r2d  ! 2018.06.26
     err( icount ) = max(epfloor,e1)       ! 2018.06.26
    else
     goto 99
    end if ! if for l=1,2      2017.08.31
   end if  ! if for iflag!     2017.08.31
   end do  ! phase loop l=1,2  2017.08.31
   end do  ! component loop    2018.10.04
  end do   ! site loop
 end do    ! src loop
end do     ! freq loop
 ndat = icount           ! 2017.07.31
 write(*,'(a,i5)') " ndat=",ndat ! 2020.09.17

ap(1)="amp"
ap(2)="pha"

if ( .true. ) then ! 2017.12.22
write(*,*) "Assembled data vector for inversion : "   ! 2020.09.29
write(*,'(32x,a)') "log10(amp [nT/A]) / phase [deg]"  ! 2020.09.29
 icount = 0
 do j=1,nfreq ! 2018.10.08
  do k=1,nsr_inv
   nobs_s = g_param_joint%obsinfo(k)%nobs_s
  do i=1,nobs_s
   do icomp = 1,5                          ! 2018.10.04
    if ( iflag_comp(icomp) .eq. 0 ) cycle  ! 2018.10.04
    do l=1,2 ! 1 for amp, 2 for phase ! 2017.08.31
    iobs = g_param_joint%obsinfo(k)%obsindex(l,icomp,i) ! 2018.10.04
    if (data_avail(l,icomp,j,iobs,k)) then       ! 2018.10.04
    icount = icount + 1                   ! 2017.07.13 do i=1,ndat
    write(*,10) icount,"obs#",i,trim(g_param%obsname(i)),comp(icomp),ap(l),"dat,err=",dat(icount),err(icount),"src",k,"freq",j
    end if
    end do
   end do
  end do
 end do
end do
end if
10 format(    i3,1x, a,  i3,1x,a,1x,a,    1x,a,1x,    a,      2f15.7,       1x,a,i3,1x,a,i3) ! 20200929

!#[2-3]## out obs files with adopted errors 2018.06.26
 head       = g_param_joint%outputfolder     ! 2018.06.26
 nhead      = len_trim(head)
 num(1:2)   = "IN"                           ! IN means "input data with adopted error"
 do k=1,nsr_inv
  sour     = sparam%sourcename(srcindex(k))
  nsour    = len_trim(sour)
  nobs_s   = g_param_joint%obsinfo(k)%nobs_s
 do i=1,nobs_s
 do icomp = 1,5 ! 2018.10.04
  if ( iflag_comp(icomp) .eq. 0 ) cycle ! 2018.10.04
 do l=1,2
   iobs  = g_param_joint%obsinfo(k)%obsindex(l,icomp,i) ! 2018.10.04
   site  = g_param%obsname(iobs)
   nsite = len_trim(site)
   apfile=head(1:nhead)//site(1:nsite)//"_"//sour(1:nsour)//"_"//comp(icomp)//ap(l)//num(1:2)//".dat"
   open(1,file=apfile) ! create newfile
   do j=1,nfreq
     !# amp
     if ( data_avail(l,icomp,j,iobs,k) ) then ! 2018.10.04
      ii = idata(l,icomp,j,iobs,k)            ! 2018.10.04
      d1 = dat(ii)
      e1 = err(ii)
      if ( l .eq. 1 ) d1 = 10**dat(ii)
      if ( l .eq. 1 ) e1 = (10**dat(ii))*log(10.)*err(ii)
      write(1,'(g15.7,a,2g15.7)') freq(j)," 1 ", d1, e1 ! amp of bz 2017.09.04
     else
      write(1,'(g15.7,a)') freq(j)," 1 0 0 "     ! no data
     end if
   end do ! freq loop      2018.10.04
   close(1)
 end do   ! amp phase loop 2018.06.26
 end do   ! component loop 2018.10.04
 end do   ! site loop      2018.06.26
 end do   ! src loop       2018.06.26

!#[3]## output
! g_data
 g_data_ap%nobs       = nobs
 g_data_ap%nfreq      = nfreq
 g_data_ap%nsr_inv    = nsr_inv ! 2017.09.06
 g_data_ap%ndat       = ndat
 allocate(g_data_ap%dvec(ndat),g_data_ap%error(ndat))
 allocate(g_data_ap%data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.04
 g_data_ap%dvec       = dat(1:ndat) ! 2017.07.13
 g_data_ap%error      = err(1:ndat) ! 2017.07.13
 g_data_ap%data_avail = data_avail  ! 2017.07.13
 g_data_ap%iflag_comp = iflag_comp  ! 2018.10.04
 g_data_ap%idata      = idata       ! (l,icomp,j,iobs,k) 2018.10.04

! g_param_joint
 allocate(g_param_joint%data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2018,10,04
 g_param_joint%data_avail = data_avail  ! 2017.07.14
 g_param_joint%ndat       = ndat


return
99 continue
write(*,*) "GEGEGE l should be 1 or 2, l=",l
stop
end subroutine

!################################################# initializedatavec
!# 2017.09.04 modified for multisource inversion
subroutine initializedatavec(g_param_joint,g_data_ap)
implicit none
type(param_joint),intent(in)       :: g_param_joint ! 2017.09.04
type(data_vec_ap),intent(out)      :: g_data_ap     ! 2017.08.31
!===============================================================
integer(4)                               :: nfreq,nobs,nsr_inv ! 2017.07.14
logical,allocatable,dimension(:,:,:,:,:) :: data_avail  ! 2018.10.05
integer(4)                               :: i,j,k,ndat

!#[0]## set
ndat    = g_param_joint%ndat    ! 2017.08.31
nobs    = g_param_joint%nobs    ! 2017.08.31
nfreq   = g_param_joint%nfreq   ! 2017.08.31
nsr_inv = g_param_joint%nsr_inv ! 2017.08.31
allocate(data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.05
data_avail = g_param_joint%data_avail ! 2017.08.31

!#[1]## output
g_data_ap%nobs    = nobs     ! 2017.08.31
g_data_ap%nfreq   = nfreq    ! 2017.08.31
g_data_ap%nsr_inv = nsr_inv  ! 2017.08.31
g_data_ap%ndat    = ndat     ! 2017.08.31
allocate(g_data_ap%dvec(ndat) ) ! 2017.08.31
allocate(g_data_ap%error(ndat)) ! 2017.08.31
allocate(g_data_ap%data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.05
g_data_ap%data_avail = data_avail

return
end subroutine

!################################################# initializedatavec
!# 2017.09.04 modified for multisource inversion
subroutine initializedatavecmt(g_param_joint,h_data_mt)
  implicit none
  type(param_joint),intent(in)       :: g_param_joint ! 2017.09.04
  type(data_vec_mt),intent(out)      :: h_data_mt     ! 2017.08.31
  !===============================================================
  integer(4)                         :: nfreq_mt,nobs_mt ! 2017.07.14
  logical,allocatable,dimension(:,:,:,:) :: data_avail_mt  ! 2018.10.05
  integer(4)                               :: i,j,k,ndat_mt
  
  !#[0]## set
  ndat_mt  = g_param_joint%ndat_mt    ! 2017.08.31
  nobs_mt  = g_param_joint%nobs_mt    ! 2017.08.31
  nfreq_mt = g_param_joint%nfreq_mt   ! 2017.08.31
  allocate(data_avail_mt(2,4,nobs_mt,nfreq_mt)) ! 2018.10.05
  data_avail_mt = g_param_joint%data_avail_mt ! 2017.08.31
  
  !#[1]## output
  h_data_mt%nobs_mt  = nobs_mt     ! 2022.01.04
  h_data_mt%nfreq_mt = nfreq_mt    ! 2022.01.04
  h_data_mt%ndat_mt  = ndat_mt     ! 2022.01.04
  allocate(h_data_mt%dvec_mt(ndat_mt) ) ! 2022.01.04
  allocate(h_data_mt%error_mt(ndat_mt)) ! 2022.01.04
  allocate(h_data_mt%data_avail_mt(2,4,nobs_mt,nfreq_mt)) ! 2018.10.05
  h_data_mt%data_avail_mt = data_avail_mt
  
  return
  end subroutine

end module param_jointinv ! 2021.12.25
