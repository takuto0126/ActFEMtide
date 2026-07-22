! Modified on 2017.08.31 for multiple sources
! Coded on 2017.06.08
! to use amp and phase of bz
program inversion_joint
!
use param
use mesh_type
use line_type
use face_type
use matrix
use modelpart
use constants        ! added on 2017.05.14
use iccg_var_takuto  ! added on 2017.05.14
use outresp
use jacobian_joint   ! added on 2022.01.05
use freq_mpi_joint   ! added on 2017.06.07FF
use caltime          ! added on 2017.09.11
use spectral         ! added on 2017.12.12
use modelroughness   ! added on 2017.12.12
!# Modification for joint inversion from inv_ap 2021.12.25
use shareformpi_joint! added on 2021.12.25
use surface_type     ! 2021.12.25
use param_mt         ! 2021.12.25
use param_jointinv   ! 2021.12.25
!
implicit none
type(param_forward)    :: g_param      ! see m_param.f90
type(param_source)     :: sparam       ! see m_param.f90
type(param_cond)       :: g_cond       ! see m_param.f90 : goal conductivity
type(param_cond)       :: h_cond       ! see m_param.f90 : initial structure
type(param_cond)       :: r_cond       ! see m_param.f90 : ref 2017.07.19
type(param_joint)      :: g_param_joint! see m_param_jointinv.f90 ! 2021.12.25
type(mesh)             :: g_mesh       ! see m_mesh_type.f90
type(mesh)             :: h_mesh       ! z file; see m_mesh_type.f90
type(line_info)        :: g_line       ! see m_line_type.f90
type(face_info)        :: g_face       ! see m_line_type.f90
type(modelpara)        :: g_modelpara  ! see m_modelpart.f90
type(model)            :: g_model_ref  ! see m_modelpart.f90
type(model)            :: g_model_ini  ! see m_modelpart.f90 2017.08.31
type(model)            :: h_model      ! see m_modelpart.f90 2017.05.17
type(model)            :: pre_model    ! see m_modelpart.f90, 2017.08.31
type(data_vec_ap)      :: g_data       ! observed data  ; see m_param_jointinv.f90
type(data_vec_ap)      :: h_data       ! calculated data; see m_param_jointinv.f90
type(data_vec_mt)      :: g_data_mt    ! observed   data 2021.12.27
type(data_vec_mt)      :: h_data_mt    ! calculated data 2021.12.27
type(freq_info_joint)  :: g_freq_joint ! see m_freq_mpi_joint.f90  ! 2022.10.20
type(real_crs_matrix)  :: RTR,CD       ! see m_matrix.f90 20
type(real_crs_matrix)  :: BMI          ! see m_matrix.f90 2017.12.13
type(real_crs_matrix)  :: BM           ! see m_matrix.f90 2017.12.25
type(real_crs_matrix)  :: PT(5)        ! (1:5)[nobs,nlin] ! 2018.10.04, Bx,By,Bz,Ex,Ey
type(real_crs_matrix)  :: R,RI         ! see m_matrix.f90 2017.12.18
type(global_matrix)    :: A            ! see m_iccg_var_takuto.f90
type(real_crs_matrix)  :: coeffobs(2,3)! see m_matrix.f90 ; 1 for edge, 2 for face
integer(4)                                  :: ijoint ! 1:ACTIVE, 2: MT, 3: Joint 2022.10.14
logical                                     :: MT, ACT ! 1: if ACTIVE, 2:if MT necessary 2022.10.14
integer(4)                                  :: ite,i,j,k,errno,node
integer(4)                                  :: nline, ntet, nobs_act, ndat ! 2017.08.31
integer(4)                                  :: ndat_mt                 ! 2022.01.04
integer(4)                                  :: nsr_inv, nmodel         ! 2017.09.04
integer(4)                                  :: nmodelactive            ! 2018.06.25
real(8)                                     :: omega, freq_tot_ip
integer(4)                                  :: nfreq_tot, nfreq_act, nfreq_tot_ip ! 2022.10.20
integer(4)                                  :: nfreq_act_ip,nfreq_mt_ip ! 2022.10.20
real(8)                                     :: nrms, nrms0, alpha      ! 2017.09.08
real(8)                                     :: misfit       ! 2017.12.22
complex(8),    allocatable,dimension(:,:)   :: fp,fs        ! (nline,nsr_inv) 2017.08.31
type(obsfiles)                              :: files        ! see m_outresp.f90
type(respdata),allocatable,dimension(:,:,:) :: resp5        ! 2017.08.31resp5(5,nsr,nfreq_ip)
type(respdata),allocatable,dimension(:,:,:) :: tresp        ! 2017.08.31 tresp(5,nsr,nfreq)
integer                                     :: access       ! 2017.05.15
type(complex_crs_matrix)                    :: ut(5)        ! 2018.10.04, Bx,By,Bz,Ex,Ey
type(amp_phase_dm),allocatable,dimension(:) :: g_apdm       ! 2017.06.07
type(amp_phase_dm),allocatable,dimension(:) :: gt_apdm      ! 2017.06.07
type(real_crs_matrix)                       :: JJ           ! Jacobian matrix, 2017.05.17
real(8)                                     :: rough1,rough2! 2017.12.13
integer(4)                                  :: ialpha       ! 2017.07.19
real(8)                                     :: nrms_init    ! 2017.07.19
character(50)                               :: head         ! 2017.07.25
character(1)                                :: num2         ! 2017.09.11
type(watch)                                 :: t_watch      ! 2017.09.11
type(watch)                                 :: t_watch0     ! 2018.03.02
real(8)                                     :: nrms_ini     ! 2018.06.25
real(8)                                     :: frms         ! 2018.06.25
!## modification from inv_ap to inv_joint  2021.12.25 ======================================
complex(8),    allocatable,dimension(:,:)   :: fs_mt        ! (nline,2) 2021.12.30
integer(4)                                  :: nfreq_mt     ! 2022.10.20
type(param_forward_mt)                      :: g_param_mt   ! 2021.12.25
type(surface)                               :: g_surface(6) ! 2021.12.25
type(param_cond)                            :: i_cond       ! see m_param.f90 2021.12.25
integer(4),    allocatable,dimension(:,:)   :: n4           ! 2021.12.27
type(real_crs_matrix)                       :: CD_mt        ! 2021.12.29
type(complex_crs_matrix)                    :: ut_mt(4)     ! [nobs_mt,nline]*4 2021.12.30
type(real_crs_matrix)                       :: PT_mt(4)     ! [nobs_mt,nline]*4 Bx,By,Ex,Ey
type(real_crs_matrix)                       :: coeffobs_mt(2,3)! see m_matrix.f90 2021.12.30
type(respdata),allocatable,dimension(:,:,:) :: resp5_mt     ! 2021.12.30
type(respdata),allocatable,dimension(:,:,:) :: tresp_mt     ! 2022.01.02
type(respmt),  allocatable,dimension(:)     ::  imp_mt      ! 2021.12.30
type(respmt),  allocatable,dimension(:)     :: timp_mt      ! 2022.01.02
integer(4)                                  :: nsr_mt=2     ! 2021.12.30
character(1) :: num
real(8)      :: nrms_mt   ! 2022.01.04
real(8)      :: misfit_mt ! 2022.01.04
real(8)      :: nrms_mt_ini     ! 2018.06.25
type(mt_dm),   allocatable,dimension(:)     :: g_mtdm  ! 2022.01.05
type(mt_dm),   allocatable,dimension(:)     :: gt_mtdm ! 2022.01.05
integer(4)                                  :: nobs_mt ! 2022.01.05
type(real_crs_matrix)                       :: JJ_mt   ! Jacobian matrix, 2022.01.05
!##===========================================================================
integer(4) :: ip,np, itemax = 20, iflag, ierr=0, i_act,i_mt
integer(4) :: kmax ! maximum lanczos procedure 2017.12.13
integer(4) :: nalpha, ialphaflag  ! 2017.09.08
integer(4) :: itype_roughness     ! 2017.12.13
!#[Explanation]---------------------------------------------------------
 !# Algorithm
 !# Phi(m) = F(m) + alpha*R(m)
 !# F(m)     = 1/2*(d(m)-d_obs)^T Cd^-1 (d(m)-d_obs)
 !# dF/dm    = J(m)^T Cd^-1 (d(m)-d_obs)
 !# d2F/dm2  = J(m)^T Cd^-1 J(m)
 !# R(m)     = 1/2*(m - m_ref)^T[BM](m-m_ref)
 !# dR/dm    = [BM](m-m_ref)
 !# d2R/dm2  = [BM]
 !#
 !# dPhi(m)/dm = 0
 !# <=> dPhi(m_k)/dm + d2Phi(m_k)/dm2*dm_k = 0
 !#  v
 !# m_k+1 - m_ref
 !# = Cm J beta
 !# ((DQ)^T(DQ))Cm = I
 !# [J Cm J^T + alpha*Cd ] beta = X      (eq. 1)
 !# X = J*(m_k - m_ref) - (d(m_k) - d_obs)
 !#
 !# Cd : [Nd,Nd] : diagonal
 !# Cd(i,i) = error(i)**2
 !#
 !#
 !# J = d(d(m))/dm : [Nd,Nm]
 !# J = [ J_f1[Nd_f1,Nm] ]
 !#     [ J_f2[Nd_f2,Nm] ]
 !#         ;
 !#     [ J_fN[Nd_fN,Nm] ]
 !#
 !# Nd_f = N_obs (for ACTIVE observation)
 !#
 !#---------------------------------------------------------
 !# Program flow
 !# [1] read data and parameters
 !# [2] gen Cm from RTR
 !# [3] gen Cd from errors
 !# [4] --- forward modelling ( Nfreq ) (by PARDISO)
 !# [5] gen Jacobian          ( Nfreq ) (by PARDISO)
 !# [6] cal data misfit F(m)
 !# [7] solve beta in eq. 1 (by PARDISO)
 !# [8] m_k+1 = m_ref + Cm*J*beta

!#[-1]## START MPI on 2017.05.29
 CALL MPI_INIT(errno)
 CALL MPI_COMM_RANK(mpi_comm_world, ip, errno)  ! ip starts with 0 in the following
 CALL MPI_COMM_SIZE(mpi_comm_world, np, errno)
 write(*,'(a,i3,a,i3)') " ip =",ip," /",np ! 2020.09.17

if ( ip .eq. 0) then !################################################# ip = 0
   ! ip=0 start

!#[-1]## select inversion type
  read(*,'(i5)') ijoint ! 1: only active, 2: only MT, 3: both active and MT data 2022.10.14
  call declareinversiontype(ijoint,ierr) ! 2022.10.14
  if ( ierr .ne. 0 ) goto 998
  call setnec(ijoint,ACT,MT) ! 2022.10.14

!#[0]## read parameters
  if(ACT) CALL READPARAM(g_param,sparam,g_cond) ! include READCOND for initial model 2022.10.14
  if(MT ) CALL READPARAM_MT(g_param_mt,i_cond) ! read MT param 2022.10.14 (i_cond is not used)
  CALL READPARAJOINTINV(ijoint,g_param_joint,g_modelpara,g_param,sparam,g_param_mt,g_data,g_data_mt)
  g_param_joint%nobs_mt  = g_param_mt%nobs  ! 2022.01.04
  g_param_joint%nfreq_mt = g_param_mt%nfreq ! 2022.01.04
  !m_param_jointinv.f90 2022.10.22

!#[1]## read mesh
  CALL READMESH_TOTAL(g_mesh,g_param%g_meshfile)
  if ( access( g_param%z_meshfile, " ") .eq. 0 ) then! if exist, 2017.05.15
    CALL READMESH_TOTAL(h_mesh,g_param%z_meshfile)
    CALL PREPZSRCOBS(h_mesh,g_param,sparam)          ! see below, include kill h_mesh
  end if
  CALL GENXYZMINMAX(g_mesh,g_param) ! see below 2021.12.27
       g_param_mt%xyzminmax = g_param%xyzminmax ! 2021.12.27
  CALL READLINE(g_param%g_lineinfofile,g_line)       ! see m_line_type.f90

  !=================================================================
  nline   = g_line%nline      ! 2021.09.14
  ntet    = g_mesh%ntet       ! 2107.09.08
  node    = g_mesh%node       ! 2017.09.08
  n4      = g_mesh%n4         ! allocate n4(ntet,4) here 2020.10.04

!#[2]## Make face of 3D mesh
  CALL MKFACE(  g_face,node,ntet,4, n4) ! make g_line
  CALL MKN4FACE(g_face,node,ntet,   n4) ! make g_line%n6line
  !  CALL READFACE(g_param_joint%g_faceinfofile,g_face) ! see m_face_type.f90
  CALL FACE2ELEMENT(g_face)                          ! cal g_face%face2ele(1:2,1:nface)

!#[3]## Make g_surface for 3DMT calculation only when MT is required
 !Extract Faceinfo and fine boundary line of 3D mesh for TM mode calculation
 !        Face2        ^ y
 !        ----         |
 ! Face3 |    | Face5  |
 !        ----         ----> x
 !        Face4
 if(MT ) CALL EXTRACT6SURFACES(g_mesh,g_line,g_face,g_surface) ! ../src_2D/m_surface_type.f90
 if(MT ) CALL FINDBOUNDARYLINE(g_mesh,g_surface)  ! Find boudnary line m_surface_type.f90
!=================================================================

!#[4]## prepare initial and ref cond, note sigma_air of ref and init is from g_cond
  CALL DUPLICATESCALARCOND(g_cond,r_cond)            ! see m_param.f90,  2017.08.31
  CALL DUPLICATESCALARCOND(g_cond,h_cond)            ! see m_param.f90,  2017.08.31
  CALL deallocatecond(g_cond)                        ! see m_param.f90,  2017.08.31
  CALL READREFINITCOND(r_cond,h_cond,g_param_joint,g_mesh)  ! 2018.10.04 homo init is added

!#[5]## generate model space see m_modelpart.f90
  CALL genmodelspace(g_mesh,g_modelpara,g_model_ref,g_param,h_cond)! m_modelpart.f90
  CALL modelparam(g_model_ref)            ! generate random model m_modelpart.f90
  CALL OUTMODEL(g_param_joint,g_model_ref,g_mesh,0,0)! to show model space 
  CALL assignmodelrho(r_cond,g_model_ref)            ! 2017.08.31 renew g_model_ref
  g_model_ini = g_model_ref                          ! 2017.08.31
  CALL assignmodelrho(h_cond,g_model_ini)            ! 2017.08.31

!#[6]## cal BMI for SM
  itype_roughness = g_param_joint%itype_roughness  ! 2017.12.25 see m_param_jointinv.f90
  if     ( itype_roughness == 1 )  then  ! SM: smoothest model
    CALL GENBMI_SM(BM,BMI,g_face,g_mesh,g_model_ref)  ! 2017.06.14 m_modelroughenss.f90
  elseif ( itype_roughness == 3 )  then  ! MSG: minimum support gradient, 2017.12.18
    CALL GENRI_MSG(R,RI,g_face,g_mesh,g_model_ref) ! 2017.12.25 m_modelroughenss.f90
  end if

!#[7]## cal Cd : ACTIVE -> CD, MT -> CD_mt
 CALL GENCD(g_data,g_data_mt,CD,CD_mt)                          ! see below

998 continue
end if !################################################################# ip = 0 end
CALL MPI_BCAST(ierr,1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! share ierr 2022.10.14
if ( ierr .ne. 0 ) goto 999 ! 2022.10.14

!#[8]## share params, mesh, and line (see m_shareformpi.f90)
  CALL shareapinv(g_param,sparam,h_cond,g_mesh,g_line,g_param_joint,g_model_ini,ip) ! 2018.10.04
  CALL sharemt(g_param_mt,g_surface,ip) ! 2021.12.30

!#[8.5]## linkglobalmodel2surface
  call linkglobalmodel2surface(g_model_ini,g_surface(2:6)) ! 2022.01.16

!#[9]## prepare A of surface for 2D TM calculation
  CALL PREPAOFSURFACE(g_surface(2:5),4,ip) ! allocate A and table_dof for 2DMT

!#[10]## cal Pt: matrix for get data from simulation results
  CALL PREPOBSCOEFF(g_param,g_mesh,g_line,coeffobs) ! for x,y,z component, see below
  CALL PREPOBSCOEFF_MT(g_param_mt,g_mesh,g_line,coeffobs_mt,ip)! for bx,by,bz,ex,ey component
  CALL PREPAREPT_JOINT(coeffobs,coeffobs_mt,g_param_joint,PT,PT_mt)! 2018.10.04 see below

!#[11]## Preparation for global stiff matrix
  ijoint       = g_param_joint%ijoint           ! 2022.10.14
  call setnec(ijoint,ACT,MT)                    ! 2022.10.14 see m_param_joint.f90
  if (ip == 0) write(*,'(a)') "  ip  |  ACT   MT" ! 2022.12.05
  CALL MPI_BARRIER(mpi_comm_world, errno)       ! 2022.12.05
  write(*,'(i4,1x,2l6)') ip,ACT,MT                          ! 2022.12.05
  nline        = g_line%nline
  nobs_act     = g_param%nobs                   ! 2022.10.20
  nobs_mt      = g_param_mt%nobs                ! 2022.01.05
  ntet         = g_mesh%ntet
  nmodel       = g_model_ini%nmodel             ! 2018.06.25
  nmodelactive = g_model_ini%nmodelactive       ! 2018.06.25
  nsr_inv      = g_param_joint%nsr_inv          ! 2017.08.31
  CALL SET_ICCG_VAR(ntet,nline,g_line,A,ip)     ! see below, 2017.06.05
  ! write(*,'(a,2i6,a,i3)') "nmodel,nmodelactive",nmodel,nmodelactive,"ip",ip ! 2018.06.26

!#[12]## set frequency
  call SETFREQIPJOINT(g_param,g_param_mt,ip,np,g_freq_joint) ! see m_freq_mpi_joint.f90 2022.10.20
  nfreq_tot    = g_freq_joint%nfreq_tot    ! ACT + MT 2022.10.20
  nfreq_tot_ip = g_freq_joint%nfreq_tot_ip ! ACT + MT for ip
  nfreq_act    = g_freq_joint%nfreq_act    ! 2022.10.20
  nfreq_act_ip = g_freq_joint%nfreq_act_ip(ip) ! 2022.10.20
  nfreq_mt     = g_freq_joint%nfreq_mt         ! 2022.10.20
  nfreq_mt_ip  = g_freq_joint%nfreq_mt_ip(ip)  ! 2022.10.20

  if (ip == 0 ) write(*,*) "-----------------------------------------------"
  CALL MPI_BARRIER(mpi_comm_world, errno)                      ! 2022.12.05
  write(*,'(2(a,i3))') "  nfreq_mt_ip  =",nfreq_mt_ip, " | ip =",ip ! 2022.12.05
  CALL MPI_BARRIER(mpi_comm_world, errno)                      ! 2022.12.05
  if (ip == 0 ) write(*,*) "-----------------------------------------------"
  CALL MPI_BARRIER(mpi_comm_world, errno)                      ! 2022.12.05
  write(*,'(2(a,i3))') "  nfreq_act_ip =",nfreq_act_ip," | ip =",ip ! 2022.12.05
  CALL MPI_BARRIER(mpi_comm_world, errno)                      ! 2022.12.05

!#[13]## allocate respdata and open output files for each observatory
  if (ACT) then ! 2022.12.05
        allocate( resp5(5,nsr_inv,nfreq_act_ip)  )  ! 2020.10.20
        allocate( tresp(5,nsr_inv,nfreq_act   )  )  ! 2020.10.20
        CALL ALLOCATERESP(g_param,nsr_inv,resp5,tresp,ip,nfreq_act,nfreq_act_ip) ! 2020.10.20
  end if

!#[14]## allocate respdata for MT
  if (MT) then ! 2022.12.05
         allocate(  resp5_mt(5,nsr_mt,nfreq_mt_ip) ) ! 2022.10.20
         allocate(  tresp_mt(5,nsr_mt,nfreq_mt   ) ) ! 2022.10.20
         allocate(    imp_mt(         nfreq_mt_ip) ) ! 2022.10.20
         allocate(   timp_mt(         nfreq_mt   ) ) ! 2022.10.20
         CALL ALLOCATERESP_MT(nobs_mt,nsr_mt,resp5_mt, imp_mt,ip,nfreq_mt_ip) !2022.10.20
         CALL ALLOCATERESP_MT(nobs_mt,nsr_mt,tresp_mt,timp_mt,ip,nfreq_mt)    !2022.10.20
  end if

 !  CALL PREPRESPFILES(g_param,files,resp5,nfreq) ! 2017.05.18

!#[15]## initialize data vector
  if (ACT) CALL initializedatavec(g_param_joint,h_data)      ! m_param_joint.f90, 2017.08.31
  if (MT)  CALL initializedatavecmt(g_param_joint,h_data_mt) ! m_param_joint.f90 2022.01.04
  ndat     = g_data%ndat                                     ! 2017.08.31
  ndat_mt  = g_data_mt%ndat_mt                               ! 2022.01.04

  if (ACT) then ! 2022.12.05
         allocate(g_apdm(nfreq_act_ip))
         call allocateapdm(nobs_act,nfreq_act_ip,nmodelactive,nsr_inv,g_apdm ) ! 2018.06.25
         allocate(gt_apdm(nfreq_act)  )
         call allocateapdm(nobs_act,nfreq_act,   nmodelactive,nsr_inv,gt_apdm) ! 2018.06.25
         allocate( fp(nline,nsr_inv),fs(nline,nsr_inv) )      ! 2017.08.31
  end if
  if (MT) then  ! 2022.12.05
         allocate( gt_mtdm(nfreq_mt), g_mtdm(nfreq_mt_ip))    ! 2022.01.05
         call allocatemtdm(nobs_mt,nfreq_mt_ip,nmodelactive, g_mtdm) ! 2022.01.05
         call allocatemtdm(nobs_mt,nfreq_mt,   nmodelactive,gt_mtdm) ! 2022.01.05
         allocate( fs_mt(nline,2)) ! 2021.12.30
  end if

!#[16]## alpha loop start ======================================= alpha loop start
 ialphaflag = g_param_joint%ialphaflag ! 2021.12.25
 ialpha     = 1
 nalpha     = 1
 if ( ialphaflag .eq. 1 ) nalpha = g_param_joint%nalpha ! L-curve 2021.12.25
 ! write(*,*) "ip",ip,"nalpha=",nalpha

 do ialpha = 1,nalpha ! nalpha > 1 when L-curve 2017.09.08

 h_model     = g_model_ini                        ! initial model, 2017.08.31
 iflag       = 0                                  ! 2017.09.11,  1: end of iteration
 alpha       = 0.d0                               ! 2017.12.13

!#[17]## open rms file (only ip = 0)
 if ( ip .eq. 0 ) then

  call watchstart(t_watch)                                    ! 2017.09.11
  !  write(*,*) " watch start alpha = ",alpha                   ! 2017.09.11

  head = g_param_joint%outputfolder                           ! 2021.12.25
  if ( ialphaflag .eq. 1 ) then                      !  L-curve 2017.09.11
   write(num2,'(i1)') ialpha                                  ! 2017.09.11
   head = head(1:len_trim(head))//"a"//num2(1:1)//"/"         ! 2017.08.31
  end if
  if (ACT) open(21,file=head(1:len_trim(head))//"rms.dat")    ! 2022.10.14
  if (MT)  open(22,file=head(1:len_trim(head))//"rms_mt.dat") ! 2022.10.14

 end if ! ip = 0 end

!#[18]## iteration loop start ========================================== iteration loop start
 do ite = 1,10 !itemax
!#[19]## output cond and model file, convert h_model to h_cond 
  if ( ip .eq. 0) then
     call watchstart(t_watch0)                   ! 2017.03.02
     !CALL OUTCOND(g_param_joint,h_cond,g_mesh,ite,ialpha)!commented out 2021.12.25
     write(*,*) ""
     write(*,'(a,i3,a)') " ======= iteration step =",ite," START!! ======" ! 2020.09.29
     CALL OUTMODEL(g_param_joint,h_model,g_mesh,ite,ialpha) ! 2017.12.25
  end if

 CALL MPI_BARRIER(mpi_comm_world, errno) ! 2020.09.29
 call model2cond(h_model,h_cond)  ! see m_modelpart.f90 2017.12.22
 !call showcond(h_cond,1)

!#[20]## nfreq_ip loop start ========================================== freq loop start
 do i=1,nfreq_tot_ip ! forward and genjacobian1 loop  2022.10.20
   freq_tot_ip = g_freq_joint%freq_tot_ip(i,ip) ! 2022.10.20
   i_act       = g_freq_joint%i_act_ip(i,ip)    ! 2022.10.20
   i_mt        = g_freq_joint%i_mt_ip(i,ip)     ! 2022.10.20
   ACT = .false. ;  MT = .false.                ! 2022.10.20
   if ( i_act > 0 ) ACT=.true.                  ! 2022.10.20
   if ( i_mt  > 0 ) MT=.true.                   ! 2022.10.20
   if ( freq_tot_ip .lt. 0. ) cycle             ! 2022.10.20
   write(*,'(a,i2)') "--------------------------------------------------- ite =",ite ! 2022.12.05
   write(*,10) "  ip =",ip," /",np," freq =",freq_tot_ip," [Hz] start!!"
   write(*,'(3(a,i3))')"  i     =",i,    " /",nfreq_tot_ip," | ip =",ip ! 2022.12.05    
   write(*,'(3(a,i3))')"  i_mt  =",i_mt, " /",g_freq_joint%nfreq_mt_ip(ip)," | ip =",ip!2022.12.05    
   write(*,'(3(a,i3))')"  i_act =",i_act," /",g_freq_joint%nfreq_act_ip(ip)," | ip =",ip!2022.12.05    
   write(*,'(a)') "-----------------------------------------------------------" ! 2022.12.05
  
  omega=2.d0*pi*freq_tot_ip ! 2022.01.04

!#[21]## COND3DTO2D : set 2D conductivity at four side surfaces
   if (MT) then ! 2022.10.14 only for MT case
     CALL COND3DTO2D(g_mesh,g_surface,h_cond) ! see m_surface_type.f90 2021.12.30
     do j=2,5 ! surface 2DTM loop, 2:north,3:west,4:south,5:east surface
       CALL forward_2DTM(g_surface(j),freq_tot_ip,g_param_mt,h_cond,ip,j)!## solve 2DTM for BC
       if (ip .eq. 0 .and. .false. ) then ! 2D*.dat
       write(num,'(i1)') j
       open(1,file="2D"//num//".dat")
       write(1,'(i4,2g15.7)') (k,g_surface(j)%bs(k),k=1,g_surface(j)%nline)
       close(1)
       end if
     end do ! end surface loop
   end if ! 2022.10.14

!#[22]## 3D ACTIVE and MT forward and obtain [fs], [fs_mt], [ut] and [ut_mt] 
   CALL forward_joint(ACT,MT,A,g_mesh,g_line,g_surface,nline,nsr_inv,fs,fs_mt,&
      & omega,sparam,g_param,g_param_joint,h_cond,PT,PT_mt,ut,ut_mt,ip,np)!2020.12.30
   if (ip .eq. 0 .and. i .eq. 1 .and. .false.) then ! fs_mt01.dat
     open(1,file="fs_mt01.dat")
     write(1,'(i6,4f15.7)') (j,fs_mt(j,1:2),j=1,nline)
     close(1)
     end if
   fp=0.d0 ! 2018.06.25

!#[23]## cal ACTIVE resp  [resp5]
  if(ACT) CALL CALOBSEBCOMP(fp,fs,nline,nsr_inv,omega,coeffobs,resp5(:,:,i_act),g_param_joint)

!#[24]## cal MT response  [resp5_mt] and [imp_mt]
  !# cal E and B at obs
  if(MT ) CALL CALOBSEBCOMP_MT(fs_mt,nline,2,omega,coeffobs_mt,resp5_mt(:,:,i_mt),ip) !see below
  !# cal MT impedance
  if(MT ) CALL CALRESPMT( resp5_mt(:,1:2,i_mt),imp_mt(i_mt),omega,ip) ! 2022.12.05

!#[25]## generate d|amp|dm and d(pha)dm for jacobian
  if(ACT) then
     CALL genjacobian1(nobs_act,nline,nsr_inv,ut,fs,PT,h_model,g_mesh,g_line,&
                  &  omega,g_apdm(i_act),g_param_joint,ip,np)     ! 2022.10.20
  end if
  if(MT ) then
     write(*,'(a,i2)') " ### genjacobian1_mt st main ###  ip =",ip
     if ( ip == 4 .or. ip == 1 ) then
     write(*,'(2(a,i8))')     "    nobs_mt",nobs_mt,    " | ip =",ip
     write(*,'(2(a,i8))')     "      nline",nline,      " | ip =",ip
     write(*,'(2(a,i8))')     "size(fs_mt)",size(fs_mt)," | ip =",ip
     write(*,'(2(a,i8))')     "size(ut_mt)",size(ut_mt)," | ip =",ip
     write(*,'(2(a,i8))')     "size(PT_mt)",size(PT_mt)," | ip =",ip
     write(*,'(a,f8.3,a,i8)') "      omega",omega,      " | ip =",ip
     write(*,'(2(a,i8))')     "       i_mt",i_mt,       " | ip =",ip
     write(*,'(2(a,i2),a,i8)')     "ip",ip,"np",np," | ip =",ip
     end if
     CALL genjacobian1_mt(nobs_mt,nline,ut_mt,fs_mt,PT_mt,h_model,g_mesh,&
                  &g_line,omega,g_mtdm(i_mt),g_param_joint,ip,np) !2022.10.20
     write(*,'(a,i2)') " ### genjacobian1_mt en main ###  ip =",ip
  end if
  write(*,'(a,i2)') " ### genjacobian1 / genjacobian1_mt / both end!! ### ip =",ip! 2022.12.05
  end do ! nfreq_tot_ip loop end

if (ip .eq. 0 .and. MT  ) then ! 2022.10.14
  open(1,file="dzxxdm1_before.dat")
  call complexcrsout(g_mtdm(1)%dzxxdm,1) !m_matrix.f90
  close(1)
  open(1,file="dzxydm1_before.dat")
  call complexcrsout(g_mtdm(1)%dzxydm,1) !m_matrix.f90
  close(1)
  end if

!============================================================= freq loop end
  call setnec(ijoint,ACT,MT) ! ip where freq is not assigned is also necessary for MPI share

!#[26]## GATHER ACTIVE and MT results to ip = 0
  if(ACT)CALL SENDRECVRESULT(resp5,tresp,ip,np,nfreq_act,nfreq_act_ip,g_freq_joint,nsr_inv) 
  if(ACT)CALL SENDBZAPRESULT_AP(g_apdm,gt_apdm,nobs_act,nfreq_act,nfreq_act_ip,g_freq_joint,nsr_inv,ip,np,g_param_joint)
  if(MT) CALL SENDRECVIMP(imp_mt,timp_mt,ip,np,nfreq_mt,nfreq_mt_ip,g_freq_joint)! 2022.01.02
  if(MT) CALL SENDRESULTINV_MT(g_mtdm,gt_mtdm,nobs_mt,nfreq_mt,nfreq_mt_ip,g_freq_joint,ip,np,g_param_joint) ! 2022.01.05

  CALL MPI_BARRIER(mpi_comm_world, errno) ! 2017.09.03

!#[27]## OUTPUT ACTIVE and MT responsed for ite=============================  ip=0 start
 if ( ip .eq. 0 ) then
   if(ACT)CALL OUTOBSFILESINV(g_param,g_param_joint,nsr_inv,sparam,tresp,nfreq_act,ite,ialpha)!2017.09.11
   if(MT )CALL OUTOBSFILESINV_MT(g_param_mt,g_param_joint,timp_mt,nfreq_mt,ite,ialpha)

   !#[28]## gen dvec
     if(ACT)CALL GENDVEC_AP(g_param_joint,nsr_inv,tresp,nfreq_act,h_data,g_data) ! 2018.10.08 [h_data]
     if(MT )CALL GENDVEC_MT(g_param_joint,timp_mt,nfreq_mt,h_data_mt) ! 2022.01.04 [h_data_mt]

   !#[29]## cal rms and misfit term ! 2022.01.02
     if(ACT)CALL CALRMS_AP(g_data,h_data,Cd,misfit,nrms)   ! cal rms, 2017.09.08
     if(MT )CALL CALRMS_MT(g_data_mt,h_data_mt,Cd_mt,misfit_mt,nrms_mt)   !2022.01.04
     if ( ACT .and. ite .eq. 1 ) nrms_ini    = nrms    ! 2022.10.14
     if ( MT  .and. ite .eq. 1 ) nrms_mt_ini = nrms_mt ! 2022.10.14

   !#[30]## cal roughness 2018.01.22
     call CALROUGHNESS(g_param_joint,h_model,g_model_ref,rough1,rough2,BM,R) !m_modelroughenss

   !#[31]## check RMS and set alpha
    !#[1]## output misfit and nrms
    frms = g_param_joint%finalrms  ! 2018.06.25
    if ( frms .lt. 0.1  ) frms=1.0 ! 2018.06.25
    if ( nrms .lt. frms ) then     ! 2018.06.25
     if (ACT) call OUTRMS(21,ite,nrms,misfit,alpha,rough1,rough2,1) !Converged  2017.09.08
     if (MT)  call OUTRMS(22,ite,nrms_mt,misfit_mt,alpha,rough1,rough2,1)!converged  2017.09.08
     iflag = 1 ;  goto 80                           ! iflag = 1 means "end" 2017.12.20
    else
     if (ACT) call OUTRMS(21,ite,nrms,misfit,alpha,rough1,rough2,0) ! Not converged  2017.09.08
     if (MT ) call OUTRMS(22,ite,nrms_mt,misfit_mt,alpha,rough1,rough2,0)!Not converged22.01.04
    end if
    if ( ite .ge. 2 .and. nrms_ini .lt. nrms ) then ! 2018.01.25 stop when nrms is larger
     write(*,*) "GEGEGE, rms exceeded the initial rms under the current inversion setting "!2022.10.31
     write(*,*) "Consider of changing the current inversion setting " ! 2022.10.31
     iflag = 1 ; goto 80                            ! 2018.01.25
    end if                                          ! 2018.01.25

   !#[32]## [Option]## if nrms increased, replace the reference model
     if ( g_param_joint%iflag_replace .eq. 1 ) then ! 2018.06.26
       if ( ite .eq. 1 ) nrms0 = nrms/0.9             ! initial     2018.06.26
       if ( nrms/nrms0 .gt. 1.0  ) then               !             2018.06.26
         if (g_param_joint%ialphaflag    .eq. 2 ) then ! cooling strategy 2017.07.19
           g_model_ref = pre_model
           h_model     = pre_model
           goto 80
         end if
       end if
     end if

   !-------------------------------------------------------- followings are for next model
    ! ialphaflag = 1 : L-curve method
    ! ialphaflag = 2 : cooling with given alpha0
    ! ialphaflag = 3 : cooling with automatic alpha0 (Minami et al. 2018)

   !#[33]## Generate BMI for next model : see m_modelroughenss.f90
    if (      itype_roughness .eq. 2 ) then ! MS:  Minimum support
     call GENBMI_MS( h_model,g_model_ref,g_param_joint,BM,BMI,ite,ialphaflag)! 2017.12.25
    else if ( itype_roughness .eq. 3 ) then ! MSG: Minimum support gradient
     call GENBMI_MSG(h_model,g_model_ref,g_param_joint,R,RI,BM,BMI,ite,ialphaflag) ! 2017.12.25
    end if

   !#[34]## generate JJ (jacobian)
    if(ACT)call genjacobian2(g_param_joint,nfreq_act,gt_apdm,JJ) !m_jacobian_joint.f90 2017.12.11
    if(MT )call genjacobian2_MT(g_param_joint,nfreq_mt,gt_mtdm,JJ_mt) ! m_jacobian_joint.f90 2022.01.05

    if ( ACT .and. g_param_joint%ioutlevel .eq. 1 ) then  ! 2018.06.25
     CALL GENMODELVOLUME(h_model,g_mesh) !see ../src_inv/m_modelpart.f90 2022.11.01
     CALL OUTJACOB(g_param_joint,g_data,JJ,ite,h_model,g_param,sparam)!m_jacobian_joint.f90 2018.06.25
    end if                                      ! 2018.06.25

    if ( ACT .and. g_param_joint%iboundflag .eq. 2) then ! 2018.01.22 J -> J'
     call TRANSJACOB(g_param_joint,JJ,h_model) ! 2018.01.22 m_jacobian_joint.f90
    end if

   !#[35]# Initial alpha
   if ( ite == 1 ) then   ! 2017.12.13
     if ( ialphaflag .eq. 1 ) alpha = g_param_joint%alpha(ialpha)! 2017.12.13
     if ( ialphaflag .eq. 2 ) alpha = g_param_joint%alpha_init   ! 2017.12.13
     if ( ialphaflag .eq. 3 ) then  ! 2017.12.13 Minami et al. (2018) cooling
       kmax = 30
       ! call alphaspectralradius_v1(Cd,JJ,BM,alpha,kmax)   ! 2017.12.25 m_spectral.f90
       if (      itype_roughness .eq. 3 ) then            ! MSG 2018.01.18
         call alphaspectralradius_v2(Cd,JJ,RI,alpha,kmax)  ! 2018.01.18 m_spectral.f90
       else if ( itype_roughness .eq. 2 ) then            ! MS  2018.02.04
         call alphaspectralradius_v2(Cd,JJ,BMI,alpha,kmax) ! 2018.01.18 m_spectral.f90
         alpha = alpha * (g_param_joint%beta**2.)          ! 2018.02.04
       else                                               ! 2018.02.04
         call alphaspectralradius_v2(Cd,JJ,BMI,alpha,kmax) ! 2017.12.24 m_spectral.f90
	     end if
       alpha = g_param_joint%gamma*alpha                   ! 2017.12.21
     end if                                               ! 2017.12.11
   end if

   !#[36]## New alpha for cooling strategies
     if ( ite .ge. 2 ) then   ! 2017.12.13
       if ( ialphaflag .eq. 2 .or. ialphaflag .eq. 3 ) then !== cooling strategy
          if ( nrms/nrms0 .gt. 0.9  ) alpha = alpha*(10.**(-1./3.d0))
       end if ! if ialphaflag = 2 or 3
     end if ! 2017.12.13
     if ( ialphaflag .eq. 4 ) then ! Modified version of Grayver et al. (2013)
	     call alphaspectralradius_v2(Cd,JJ,BMI,alpha,30) ! 2018.01.09 m_spectral.f90
       alpha = g_param_joint%gamma*alpha/(1.*ite) ! 2018.01.09
     end if ! 2018.09.01

   !#[37]## obtain new model
    pre_model = h_model ! 2017.06.14 keep previous model
    call getnewmodel_joint(JJ,JJ_mt,g_model_ref,h_model,g_data,h_data,&
    &   g_data_mt,h_data_mt,BMI,CD,CD_mt,alpha,g_param_joint) ! 2022.01.05
    nrms0 = nrms                                 !  2017.12.25

    80 continue                                    !  2017.12.20

  end if ! ip = 0

  CALL MPI_BARRIER(mpi_comm_world, errno)
  CALL MPI_BCAST(iflag,1,MPI_INTEGER4,0,mpi_comm_world,errno)
  if (iflag .eq. 1) goto 100

 !#[38]## Share new model
   CALL SHAREMODEL(h_model,ip) ! see m_shareformpi.f90, 2017.06.05

   if ( ip .eq. 0) call watchstop(t_watch0)
   if ( ip .eq. 0) write(*,'(a,i3,a,f8.4,a)') " ======= iteration step =",ite," END!! ======= Time=",t_watch0%time," [min]"!2020.09.17

end do ! iteration loop============================================  iteration end!

 100 continue
 if ( ip .eq. 0 ) close(21) ! 2017.05.18
 if ( ip .eq. 0 ) then
  call watchstop(t_watch)   ! 2017.09.11
  write(*,'(a,g15.7,a,f15.7,a)') "### END alpha=",alpha," Time=",t_watch%time," [min]" !2022.1.2
 end if

end do ! alpha loop end! 2017.09.08

999 continue ! 2022.10.14

 CALL MPI_FINALIZE(errno) ! 2017.06.05


 10 format(a,i3,a,i3,a,f8.3,a) !2022.01.04

end program inversion_joint ! 2021.12.25

!############################################# LINKglobalmodel2surface
 subroutine linkglobalmodel2surface(g_model,g_surface)
  use surface_type
  use modelroughness
  use modelpart
 implicit none
 type(model),  intent(in)    :: g_model
 type(surface),intent(inout) :: g_surface(5) ! 4sides and bottom
 integer(4) :: i,j,ii,iele,nmodel_global
 integer(4),allocatable,dimension(:) :: stack,item

 !# model count at surface
 do i=1,5 ! surface loop
  do j=1,nmodel_global
    do ii=stack(j-1)+1,stack(j)   
     iele = item(ii)! global element id
      

    end do
  end do
 end do

 return
 end 
!#############################################
! coded on 2021.09.15
subroutine calrespmt(resp5,resp_mt,omega,ip) ! 2022.12.05
  use outresp
  use constants ! dmu,pi
  implicit none
  real(8),       intent(in)    :: omega
  integer(4),    intent(in)    :: ip     ! 2022.12.05
  type(respdata),intent(in)    :: resp5(5,2) ! 1 for ex, 2 for ey polarization
  type(respmt),  intent(inout) :: resp_mt
  complex(8),    allocatable   :: be5_ex(:,:),be5_ey(:,:) ! be5 = bx,by,bz,ex,ey
  complex(8)                   :: a,b,c,d,iunit=(0.d0,1.d0)
  complex(8)                   :: det,z(2,2),bi(2,2),e(2,2)
  integer(4)                   :: i,j,nobs
  real(8)                      :: coef,amp,phase
  
  nobs = resp_mt%nobs
  allocate(be5_ex(5,nobs),be5_ey(5,nobs))
  ! calculate impedance Z = E/B [mV/km]/[nT]
  ! (Ex_ex Ex_ey) = (Zxx Zxy)(Bx_ex Bx_ey)
  ! (Ey_ex Ey_ey) = (Zyx Zyy)(By_ex By_ey)
  !  [E]=Z[B]
  !  Z = [E][B]^-1
  
  !# set bxyzexy_ex and bxyzexy_ey
  write(*,*) "nobs",nobs,"ip",ip
  do i=1,5
!    write(*,*) "allocated(resp5(i,1)%ftobs)",allocated(resp5(i,1)%ftobs)
!    write(*,*) "allocated(resp5(i,2)%ftobs)",allocated(resp5(i,2)%ftobs)
!    write(*,'(a,i3,i8,a,i3)') "i,size(resp5(i,1)%ftobs)",i,size(resp5(i,1)%ftobs),"ip",ip
!    write(*,'(a,i3,i8,a,i3)') "i,size(resp5(i,2)%ftobs)",i,size(resp5(i,2)%ftobs),"ip",ip
   do j=1,nobs
    be5_ex(i,j)=resp5(i,1)%ftobs(j) ! ex polarization
    be5_ey(i,j)=resp5(i,2)%ftobs(j) ! ey plarization
   end do
  end do
  
  write(*,*)
  !# calculate impedance
  do j=1,nobs
   a = be5_ex(1,j) ! Bx_ex
   c = be5_ex(2,j) ! By_ex ! 2022.01.12
   b = be5_ey(1,j) ! Bx_ey ! 2022.01.12
   d = be5_ey(2,j) ! By_ey
   det = a*d - b*c
   bi(1,1:2)=(/ d, -b/)
   bi(2,1:2)=(/ -c, a/)
   bi = bi/det
   e(1,1:2)=(/be5_ex(4,j), be5_ey(4,j)/) ! (ex1, ex2)
   e(2,1:2)=(/be5_ex(5,j), be5_ey(5,j)/) ! (ey1, ey2)
   z = matmul(e,bi)
   resp_mt%zxx(j) = z(1,1) ! [mV/km]/[nT]
   resp_mt%zxy(j) = z(1,2)
   resp_mt%zyx(j) = z(2,1)
   resp_mt%zyy(j) = z(2,2)
  ! rho and pha, rhoa = mu/omega*|Z|**2.
  coef = dmu/omega*1.d+6
   resp_mt%rhoxx(j) = coef*amp(z(1,1))**2. ! [Ohm.m]
   resp_mt%rhoxy(j) = coef*amp(z(1,2))**2. ! [Ohm.m]
   resp_mt%rhoyx(j) = coef*amp(z(2,1))**2. ! [Ohm.m]
   resp_mt%rhoyy(j) = coef*amp(z(2,2))**2. ! [Ohm.m]
  !
   resp_mt%phaxx(j) = phase(z(1,1))
   resp_mt%phaxy(j) = phase(z(1,2))
   resp_mt%phayx(j) = phase(z(2,1))
   resp_mt%phayy(j) = phase(z(2,2))
  end do
  
  write(*,'(a,i2)') " ### CALRESPMT        END !! ###  ip =",ip  ! 2022.12.05
  return
  end

!#############################################
  subroutine declareinversiontype(ijoint,ierr) ! 2022.10.14
  implicit none
  integer(4),intent(in)  :: ijoint
  integer(4),intent(inout) :: ierr
  if (ijoint .eq. 1 ) then
    write(*,*) "START ACTIVE INVERSION!!"
  else if (ijoint .eq. 2 ) then
    write(*,*) "START MT INVERSION!!"
  else if (ijoint .eq. 3 ) then
    write(*,*) "START ACTIVE + MT JOINT INVERSION!!"
  else 
   write(*,*) "GEGEGE ijoint =",ijoint,"should be 1: ACTIVE, 2:MT, or 3: Joint"
   ierr = 1
  end if

  return
  end
!#############################################
!# copied from ../src_3DMT/n_ebfem_3DMT.f90
subroutine CALOBSEBCOMP_MT(fs,nline,nsr,omega,coeffobs,resp5,ip) ! 2022.12.05
  use matrix
  use outresp
  implicit none
  real(8),              intent(in)    :: omega
  integer(4),           intent(in)    :: nline, nsr
  integer(4),           intent(in)    :: ip ! 2022.12.05
  complex(8),           intent(inout) :: fs(nline,nsr) ! 2017.07.11
  type(real_crs_matrix),intent(in)    :: coeffobs(2,3)
  type(respdata),       intent(inout) :: resp5(5,nsr)                 !2017.07.11
  integer(4)                          :: i  ! 2017.07.11
  
  do i=1, nsr ! 2017.07.11
   CALL CALOBSRESP_3DMT(fs(:,i),nline,coeffobs(2,1),resp5(1,i)  ) !bx,fp deleted 2021.09.15
   CALL CALOBSRESP_3DMT(fs(:,i),nline,coeffobs(2,2),resp5(2,i)  ) !by,fp deleted 2021.09.15
   CALL CALOBSRESP_3DMT(fs(:,i),nline,coeffobs(2,3),resp5(3,i)  ) !bz,fp deleted 2021.09.15
   fs(:,i)= - (0.d0,1.d0)*omega*fs(:,i) !E= -i*omega*A
   CALL CALOBSRESP_3DMT(fs(:,i),nline,coeffobs(1,1),resp5(4,i)  ) !ex,fp deleted 2021.09.15
   CALL CALOBSRESP_3DMT(fs(:,i),nline,coeffobs(1,2),resp5(5,i)  ) !ey,fp deleted 2021.09.15
  end do      ! 2017.07.11
  
  write(*,'(a,i2)') " ### CALOBSEBCOMP_MT  END !! ###  ip =",ip ! 2022.12.05
  
  return
  end

!#################################################### PREPAREPT
!# PT [nobs,nline]*ncomp is generated 2021.12.30
subroutine PREPAREPT_JOINT(coeffobs,coeffobs_mt,g_param_joint,PT,PT_mt) ! 2021.12.30
 use matrix
 use param_jointinv
 implicit none
 type(param_joint),       intent(in)  :: g_param_joint
 type(real_crs_matrix),   intent(in)  :: coeffobs(2,3)
 type(real_crs_matrix),   intent(in)  :: coeffobs_mt(2,3)
 type(reaL_crs_matrix),   intent(out) :: PT(5)        ! 2018.10.04
 type(reaL_crs_matrix),   intent(out) :: PT_mt(4)     ! 2021.12.30
 integer(4)                           :: iflag_comp(5)
 integer(4)                           :: i,j,icomp

 !#[1]## set
 iflag_comp = g_param_joint%iflag_comp
 do i=1,5
  PT(i)%ntot = 0
 end do

 !#[2]##Prepare PT: calculate nrow
  icomp = 0
  do i=2,1,-1 !#  B, E loop
  do j=1,3 !# x,y,z comp loop
   if ( icomp .eq. 5 ) exit    ! 2018.10.05 coeffobs(i,j) -> i=1:E,i=2:B, j=1-3: x,y,z component 
   icomp = icomp + 1
   if ( iflag_comp(icomp) .eq. 0 ) cycle
   CALL DUPLICATE_CRSMAT(coeffobs(i,j),PT(icomp)) ! PT(1:5) [nobs,nline] , m_matrix.f90
  end do
  end do

 !#[3]##Prepare PT_mt: calculate nrow
  CALL DUPLICATE_CRSMAT(coeffobs_mt(2,1),PT_mt(1)) !bx PT(1:5) [nobs,nline] , m_matrix.f90
  CALL DUPLICATE_CRSMAT(coeffobs_mt(2,2),PT_mt(2)) !by PT(1:5) [nobs,nline] , m_matrix.f90
  CALL DUPLICATE_CRSMAT(coeffobs_mt(1,1),PT_mt(3)) !ex PT(1:5) [nobs,nline] , m_matrix.f90
  CALL DUPLICATE_CRSMAT(coeffobs_mt(1,2),PT_mt(4)) !ey PT(1:5) [nobs,nline] , m_matrix.f90

 return
 end
!#################################################### READREFINITCOND
!# Coded on 2018.06.21
subroutine READREFINITCOND(r_cond,h_cond,g_param_joint,g_mesh) ! 2018.10.04
 use param_jointinv ! 2018.06.21
 use modelpart
 use param
 use mesh_type   ! 2018.10.04
 implicit none
 type(mesh),           intent(in)    :: g_mesh          ! 2018.10.04
 type(param_cond),     intent(inout) :: r_cond
 type(param_cond),     intent(inout) :: h_cond
 type(param_joint),intent(in)        :: g_param_joint   ! 2018.06.21
 character(50)                       :: modelfile,connectfile
 integer(4)                          :: icondflag_ini, icondflag_ref
 real(8)                             :: sigmahomo       ! 2018.10.04

 !#[0]##
 icondflag_ref = g_param_joint%icondflag_ref         ! 2018.06.21
 icondflag_ini = g_param_joint%icondflag_ini         ! 2018.06.21

 !#[ref cond]
 if     ( icondflag_ref .eq. 0 ) then                ! 2018.10.04
     sigmahomo = g_param_joint%sigmahomo_ref            ! 2018.10.04
     CALL SETCOND(r_cond,g_mesh,sigmahomo)              ! 2018.10.04
   elseif ( icondflag_ref .eq. 1 ) then ! condfile       2018.03.18
     r_cond%condfile = g_param_joint%g_refcondfile      ! 2017.07.19
     CALL READCOND(r_cond)                              ! 2017.07.19
     elseif ( icondflag_ref .eq. 2 ) then ! modelfile
     modelfile   = g_param_joint%g_refmodelfile
     connectfile = g_param_joint%g_refmodelconn
     CALL READMODEL2COND(r_cond,connectfile,modelfile)
   else                                                ! 2018.10.04
     write(*,*) "GEGEGE! icondflag_ref",icondflag_ref   ! 2018.10.04
     stop                                               ! 2018.10.04
 end if

 !#[ini cond]
 if     ( icondflag_ini .eq. 0 ) then                ! 2018.10.04
   sigmahomo = g_param_joint%sigmahomo_ini            ! 2018.10.04
   CALL SETCOND(h_cond,g_mesh,sigmahomo)              ! 2018.10.04
 elseif ( icondflag_ini .eq. 1 ) then ! condfile       2018.03.18
   h_cond%condfile = g_param_joint%g_initcondfile     ! see m_param_inv.f90
   CALL READCOND(h_cond)                              ! see m_param.f90
 elseif ( icondflag_ini .eq. 2 ) then ! modelfile
   modelfile   = g_param_joint%g_initmodelfile
   connectfile = g_param_joint%g_initmodelconn
   CALL READMODEL2COND(h_cond,connectfile,modelfile)
 else                                                ! 2018.10.04
   write(*,*) "GEGEGE! icondflag_ini",icondflag_ini   ! 2018.10.04
   stop                                               ! 2018.10.04
 end if

 !# set nphys1
 CALL SETNPHYS1INDEX2COND(g_mesh,r_cond)            ! ref     cond; see below
 CALL SETNPHYS1INDEX2COND(g_mesh,h_cond)            ! initial cond; see below

 return
 end
!##################################################### SETCOND
!# coded on 2018.10.04
!# case where sigmahomo is given
subroutine setcond(h_cond,g_mesh,sigmahomo)
 use mesh_type
 use param     ! 2018.10.05
 implicit none
 type(mesh),       intent(in)    :: g_mesh
 real(8),          intent(in)    :: sigmahomo
 type(param_cond), intent(inout) :: h_cond
 integer(4)                      :: nphys2,ntet,i
 integer(4),allocatable,dimension(:,:) :: n4flag
 real(8),   allocatable,dimension(:)   :: rho, sigma

 !# set
 ntet   = g_mesh%ntet
 allocate( n4flag(ntet,2) )
 n4flag = g_mesh%n4flag

 !# cal culate nphys2
 nphys2 = 0
 do i=1,ntet
  if ( n4flag(i,2) .ge. 2 ) nphys2 = nphys2 + 1
 end do
 write(*,*) "nphys2=",nphys2

 !# set output
 allocate(h_cond%sigma(nphys2))
 allocate(h_cond%rho(  nphys2))
 allocate(h_cond%index(nphys2))
 h_cond%nphys2   = nphys2
 h_cond%rho(:)   = 1.d0/sigmahomo
 h_cond%sigma(:) = sigmahomo

 return
 end
!##################################################### READMODEL2COND
!# coded on 2018.06.21
subroutine readmodel2cond(r_cond,connectfile,modelfile)
 use modelpart
 use param
 implicit none
 type(param_cond), intent(inout)     :: r_cond
 character(50),    intent(in)        :: connectfile
 character(50),    intent(in)        :: modelfile
 integer(4)                          :: nphys2,nmodel,nmodel2
 integer(4)                          :: i
 integer(4),allocatable,dimension(:) :: id,ele2model
 real(8),   allocatable,dimension(:) :: rho

 !#[1]## read connectfile
  open(1,file=connectfile,status='old',err=90)
  read(1,*) nphys2,nmodel
  allocate(id(nphys2),ele2model(nphys2))
  do i=1,nphys2
   read(1,'(2i10)') id(i),ele2model(i)
  end do
  close(1)

 !#[2]## read model
 allocate(rho(nmodel))
 open(1,file=modelfile,status='old',err=80) ! 20127.12.21
  read(1,*,err=81) nmodel2
  if ( nmodel .ne. nmodel2 ) then
   write(*,*) "GEGEGE nmodel",nmodel,"nmodel2",nmodel2
   stop
  end if
  do i=1,nmodel
   read(1,*) rho(i)
  end do
 close(1)

 !#[3]## gen r_cond
 allocate(r_cond%rho(  nphys2) )
 allocate(r_cond%sigma(nphys2) )
 allocate(r_cond%index(nphys2)) ! 2018.03.20
 r_cond%nphys2 = nphys2
 r_cond%sigma  = -9999
 do i=1,nphys2
  r_cond%rho(i)=rho(ele2model(i))
  if (abs(r_cond%rho(i)) .gt. 1.d-10 ) r_cond%sigma(i) = 1.d0/r_cond%rho(i)
 end do

 return
  90 continue
  write(*,*) "File is not exist",connectfile
  stop
  80 continue
  write(*,*) "File is not exist",modelfile
  stop
  81 continue
  write(*,*) "File is not exist",modelfile,"line",i
  stop
 end

!######################################### OUTRMS
!# rms -> misfit on 2017.12.22
!# coded on 2017.09.08
subroutine OUTRMS(idev,ite,nrms,misfit,alpha,rough1,rough2,icflag)
 implicit none
 integer(4),intent(in) :: idev,ite
 integer(4),intent(in) :: icflag ! 0 : not converged, 1: converged
 real(8),   intent(in) :: misfit,nrms,alpha
 real(8),   intent(in) :: rough1 ! roughness using usual model
 real(8),   intent(in) :: rough2 ! roughness using transformed model (only iboundtype = 2)

 if (ite .eq. 1) write(idev,'(a)') "iteration#     nRMS            RMS          alpha         Roughness   tRoughness"!2020.09.29
 if      ( icflag .eq. 0 ) then
     write(idev, '(i10,5g15.7)') ite,nrms, misfit, alpha,rough1,rough2 ! 2017.09.04
     write(*,*) "Not converged yet..." ! 2022.10.31
 else if ( icflag .eq. 1 ) then
     write(*,*) "Converged!!with nrms =",nrms ! normalized rms 2017.09.08
     write(idev,'(i10,5g15.7,a)') ite, nrms, misfit, alpha,rough1,rough2,"Converged!"! 2017.09.04
 end if
 return
 end
!######################################### OUTRMS
!# rms -> misfit on 2017.12.22
!# coded on 2017.09.08
subroutine OUTRMS_MT(idev,ite,nrms_mt,misfit_mt,alpha,rough1,rough2,icflag)
  implicit none
  integer(4),intent(in) :: idev,ite
  integer(4),intent(in) :: icflag ! 0 : not converged, 1: converged
  real(8),   intent(in) :: misfit_mt,nrms_mt,alpha
  real(8),   intent(in) :: rough1 ! roughness using usual model
  real(8),   intent(in) :: rough2 ! roughness using transformed model (only iboundtype = 2)
  
  if (ite .eq. 1) write(idev,'(a)') "iteration#     nRMSmt          RMSmt        alpha         Roughness   tRoughness"!2020.09.29
  if      ( icflag .eq. 0 ) then
       write(idev, '(i10,5g15.7)') ite,nrms_mt, misfit_mt, alpha,rough1,rough2 ! 2017.09.04
  else if ( icflag .eq. 1 ) then
       write(*,*) "Converged!!with nrms_mt =",nrms_mt ! normalized rms 2017.09.08
       write(idev,'(i10,5g15.7,a)') ite, nrms_mt, misfit_mt, alpha,rough1,rough2,"Converged!"! 2017.09.04
  end if
  return
  end
  

!######################################### OUTOBSFILESINV 2018.10.05
!# modified on 2018.10.05
!# copied from src_inv_mpi/n_inv_mpi.f90 2017.09.03
subroutine OUTOBSFILESINV(g_param,g_param_joint,nsr_inv,sparam,tresp,nfreq,ite,ialpha)
use param_jointinv ! 2017.09.03
use param
use outresp
implicit none
integer(4),                intent(in)    :: ite,nfreq,ialpha ! 2017.09.11
integer(4),                intent(in)    :: nsr_inv
type(respdata),            intent(in)    :: tresp(5,nsr_inv,nfreq)!2017.07.14
type(param_forward),       intent(in)    :: g_param
type(param_source),        intent(in)    :: sparam          ! 2017.07.14
type(param_joint),   intent(in)    :: g_param_joint   ! 2017.07.14
integer(4), allocatable, dimension(:)    :: srcindex        ! 2017.07.14
real(8),    allocatable, dimension(:)    :: freq            ! 2017.07.14
character(2)                             :: num             ! 2017.07.14
integer(4)                               :: nh, nsi,nso     ! 2018.10.05
integer(4)                               :: i,j,k,l,nobs    ! 2017.07.14
character(100)                           :: ampname,phaname ! 2017.09.03
character(50)                            :: head,site, sour ! 2017.07.14
logical,allocatable,dimension(:,:,:,:,:) :: data_avail      ! 2018.10.05
integer(4)                               :: idat,ialphaflag ! 2017.09.11
character(1)                             :: num2            ! 2017.09.11
integer(4)                               :: icomp
integer(4)                               :: iflag_comp(5)

!#[0]## set
 nobs       = g_param%nobs
 head       = g_param_joint%outputfolder     ! 2017.07.25
 nh         = len_trim(head)
 allocate(srcindex(nsr_inv), freq(nfreq) )   ! 2017.07.14
 allocate(data_avail(2,5,nfreq,nobs,nsr_inv))! 2018.10.05
 data_avail = g_param_joint%data_avail       ! 2017.09.03 see m_param_inv.f90
 srcindex   = g_param_joint%srcindex         ! 2017.07.14
 freq       = g_param%freq                   ! 2017.07.14
 ialphaflag = g_param_joint%ialphaflag       ! 2017.09.11
 iflag_comp = g_param_joint%iflag_comp       ! 2018.10.05
 if ( ialphaflag .eq. 1 ) then      ! L-curve  2017.09.11
    write(num2,'(i1)') ialpha                ! 2017.09.11
    head = head(1:nh)//"a"//num2(1:1)//"/"   ! 2017.09.11
    nh = len_trim(head)                      ! 2017.09.11
 end if

 write(num,'(i2.2)')  ite                    ! 2017.07.14
 write(*,*) "nsr_inv",nsr_inv
 write(*,*) "nfreq",  nfreq
 write(*,*) "nobs",   nobs

!#[2]##
 do l=1,nobs
  site  = g_param%obsname(l)
  nsi   = len_trim(site)

  do k=1,nsr_inv
   sour     = sparam%sourcename(srcindex(k))
   nso      = len_trim(sour)

   do icomp = 1,5                        ! 2018.10.05
    if ( iflag_comp(icomp) .eq. 0) cycle ! 2018.10.05

    ampname = head(1:nh)//site(1:nsi)//"_"//sour(1:nso)//"_"//comp(icomp)//"amp"//num//".dat"
    phaname = head(1:nh)//site(1:nsi)//"_"//sour(1:nso)//"_"//comp(icomp)//"pha"//num//".dat"
    open(31,file=ampname) ! 2017.09.04
    open(32,file=phaname) ! 2017.09.04

    do i=1,nfreq
      idat = 0 ! 2017.09.04
      if ( data_avail(1,icomp,i,l,k)) idat = 1
       write(31,110) freq(i),idat,tresp(icomp,k,i)%ftobsamp(l) ! amp of bz 2018.10.05
      idat = 0  ! 2017.09.04
      if ( data_avail(2,icomp,i,l,k)) idat = 1 ! pha 2017.09.04
       write(32,110) freq(i),idat,tresp(icomp,k,i)%ftobsphase(l)! pha of bz 2018.10.05
    end do

    close(31) ! 2017.09.04
    close(32) ! 2017.09.04
   end do     ! 2018.10.05
  end do
 end do

 write(*,*) "### OUTOBSFILESINV END!! ###"

110 format(g15.7,i5,g15.7)
return
end

!######################################### OUTOBSFILESINV 2018.10.05
!# coded on 2022.01.02
subroutine OUTOBSFILESINV_MT(g_param,g_param_joint,resp_mt,nfreq,ite,ialpha)
  use param_jointinv ! 2017.09.03
  use param_mt
  use outresp
  implicit none
  integer(4),                intent(in)    :: ite,nfreq,ialpha! 2017.09.11
  type(respmt),              intent(in)    :: resp_mt(nfreq)  ! 2022.01.02
  type(param_forward_mt),    intent(in)    :: g_param         ! 2022.01.02
  type(param_joint),         intent(in)    :: g_param_joint   ! 2017.07.14
  real(8)                                  :: freq            ! 2017.07.14
  character(2)                             :: num             ! 2017.07.14
  integer(4)                               :: nhead, nsite    ! 2018.10.05
  integer(4)                               :: i,j,k,l,nobs    ! 2017.07.14
  character(100)                           :: filename1       ! 2022.01.22
  character(100)                           :: filename2       ! 2022.01.22
  character(50)                            :: head,site       ! 2017.07.14
  logical,allocatable,dimension(:,:,:,:)   :: data_avail      ! 2022.01.02
  integer(4)                               :: idat,ialphaflag ! 2017.09.11
  character(1)                             :: num2            ! 2017.09.11
  integer(4)                               :: icomp
  
  !#[0]## set
   nobs       = g_param%nobs
   head       = g_param_joint%outputfolder     ! 2017.07.25
   nhead      = len_trim(head)
   allocate(data_avail(2,4,nfreq,nobs) )        ! 2022.01.02
   data_avail = g_param_joint%data_avail_mt     ! 2017.09.03 see m_param_inv.f90
   ialphaflag = g_param_joint%ialphaflag        ! 2017.09.11
   if ( ialphaflag .eq. 1 ) then      ! L-curve  2017.09.11
      write(num2,'(i1)') ialpha                 ! 2017.09.11
      head  = head(1:nhead)//"a"//num2(1:1)//"/"! 2017.09.11
      nhead = len_trim(head)                    ! 2017.09.11
   end if
  
   write(num,'(i2.2)')  ite                    ! 2017.07.14
   write(*,*) "nfreq",  nfreq
   write(*,*) "nobs",   nobs
  
  !#[2]##
 do l=1,nobs
  site  = g_param%obsname(l)
  nsite = len_trim(site)
  filename1 = head(1:nhead)//site(1:nsite)//"_MT"//num//".dat"     ! 2022.01.02
  filename2 = head(1:nhead)//site(1:nsite)//"_MT_imp"//num//".dat" ! 2022.01.02
   open(31,file=filename1)
   open(32,file=filename2) ! 2021.12.15

   do i=1,nfreq
    freq = g_param%freq(i)
    write(31,'(9g15.7)') freq,resp_mt(i)%rhoxx(l),resp_mt(i)%phaxx(l),&
    &                         resp_mt(i)%rhoxy(l),resp_mt(i)%phaxy(l),&
    &                         resp_mt(i)%rhoyx(l),resp_mt(i)%phayx(l),&
    &                         resp_mt(i)%rhoyy(l),resp_mt(i)%phayy(l)
    write(32,'(9g15.7)') freq,real(resp_mt(i)%zxx(l)),imag(resp_mt(i)%zxx(l)),& ! 2021.12.15
    &                         real(resp_mt(i)%zxy(l)),imag(resp_mt(i)%zxy(l)),& ! 2021.12.15
    &                         real(resp_mt(i)%zyx(l)),imag(resp_mt(i)%zyx(l)),& ! 2021.12.15
    &                         real(resp_mt(i)%zyy(l)),imag(resp_mt(i)%zyy(l))    ! 2021.12.15
   end do

   close(31)
   close(32) ! 2021.12.15
  end do

  
   write(*,*) "### OUTOBSFILESINV_MT END!! ###"
  
  110 format(g15.7,i5,g15.7)
  return
  end
!############################################
!# coded on 2017.06.05
subroutine SET_ICCG_VAR(ntet,nline,g_line,A,ip)
use iccg_var_takuto
use line_type
implicit none
integer(4),         intent(in)    :: ntet,nline,ip
type(line_info),    intent(in)    :: g_line
type(global_matrix),intent(inout) :: A
integer(4),allocatable,dimension(:,:) :: table_dof
integer(4) :: dofn = 1

!#[1]## set table_dof
  allocate( table_dof(nline,dofn))
  CALL SET_TABLE_DOF(dofn,nline,table_dof,nline)! table_dof is generated see forward_2DTMinv.f90

!#[2]## set allocate A
  CALL set_iccg_var7_dofn(dofn,6,nline,nline,nline,ntet,g_line%n6line,table_dof,A,ip)

return
end

!############################################
!# modified on 2017.09.03 for multiple sources
!# coded on 2017.06.05
subroutine SENDBZAPRESULT_AP(g_apdm,gt_apdm,nobs,nfreq_act,nfreq_act_ip,g_freq_joint,nsr,ip,np,g_param_joint)!2018.10.05
use matrix
use shareformpi_joint ! 2021.12.25
use jacobian_joint    ! 2017.06.08
use param_jointinv    ! 2018.10.08
use freq_mpi_joint    ! 2022.10.20
implicit none
integer(4),             intent(in)    :: nobs,nfreq_act,nfreq_act_ip,ip,np
integer(4),             intent(in)    :: nsr                   ! 2017.09.03
type(freq_info_joint),  intent(in)    :: g_freq_joint          ! 2022.10.20
type(amp_phase_dm),     intent(in)    :: g_apdm(nfreq_act_ip)  ! 2017.09.03
type(param_joint),      intent(in)    :: g_param_joint         ! 2018.10.08
type(amp_phase_dm),     intent(inout) :: gt_apdm(nfreq_act)    ! 2017.09.03
integer(4)                            :: i,ifreq,ip_from,errno
integer(4)                            :: icomp                 ! 2018.10.05
integer(4)                            :: k                     ! 2017.09.03
integer(4),          dimension(5)     :: iflag_comp            ! 2018.10.08

!#[0]## set
  iflag_comp = g_param_joint%iflag_comp ! 2018.10.08

!#[1]## share ampbz
 do i=1,nfreq_act
!  if (mod(i,np) .eq. 1 ) ip_from = -1
!  ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1
  ip_from  = g_freq_joint%ip_from_act(i) ! 2022.10.20
  ifreq    = g_freq_joint%if_g2l_act(i)  ! 2022.10.20
  

  do k=1,nsr      ! 2017.09.03

  do icomp = 1,5  ! 2018.10.05

  if ( iflag_comp(icomp) .eq. 0)  cycle ! 2018.10.08

  !# dampdm(icomp,k) is a matrix of [nobs,nmodel] for icomp component and k-th source
  if ( ip .eq. ip_from )  gt_apdm(i)%dampdm(icomp,k) = g_apdm(ifreq)%dampdm(icomp,k)!2018.10.05
  if ( ip .eq. ip_from )  gt_apdm(i)%dphadm(icomp,k) = g_apdm(ifreq)%dphadm(icomp,k)!2018.10.05

  !# share     2017.09.03
  call sharerealcrsmatrix(gt_apdm(i)%dampdm(icomp,k),ip_from,ip) ! m_shareformpi_joint.f90
  call sharerealcrsmatrix(gt_apdm(i)%dphadm(icomp,k),ip_from,ip) ! m_shareformpi_joint.f90

  end do     ! 2018.10.05 component loop

  end do     ! 2017.09.03

 end do

!#[2]##
  if (ip .eq. 0) write(*,*) "### SENDBZAPRESULT END!! ###"

return
end
!############################################
!# modified on 2017.09.03 for multiple sources
!# coded on 2017.06.05
subroutine SENDRESULTINV_MT(g_mtdm,gt_mtdm,nobs_mt,nfreq_mt,nfreq_mt_ip,g_freq_joint,ip,np,g_param_joint)!2018.10.05
  use matrix
  use shareformpi_joint ! 2021.12.25
  use jacobian_joint    ! 2017.06.08
  use param_jointinv    ! 2018.10.08
  use freq_mpi_joint    ! 2022.10.20
  implicit none
  integer(4),       intent(in)    :: nobs_mt,nfreq_mt,nfreq_mt_ip,ip,np ! 2022.10.20
  type(freq_info_joint),intent(in) :: g_freq_joint
  type(param_joint),intent(in)    :: g_param_joint         ! 2018.10.08
  type(mt_dm),      intent(in)    :: g_mtdm(nfreq_mt_ip)      ! 2017.09.03
  type(mt_dm),      intent(inout) :: gt_mtdm(nfreq_mt)        ! 2017.09.03
  integer(4)                      :: i,ifreq,ip_from,errno
  integer(4)                            :: k                     ! 2017.09.03
    
  !#[1]## share ampbz
   do i=1,nfreq_mt
!    if (mod(i,np) .eq. 1 ) ip_from = -1
!    ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1
   
   ip_from  = g_freq_joint%ip_from_mt(i) ! 2022.10.20
   ifreq    = g_freq_joint%if_g2l_mt(i)  ! 2022.10.20
      
    if ( ip .eq. ip_from )  gt_mtdm(i)%dzxxdm = g_mtdm(ifreq)%dzxxdm!2018.10.05
    if ( ip .eq. ip_from )  gt_mtdm(i)%dzxydm = g_mtdm(ifreq)%dzxydm!2018.10.05
    if ( ip .eq. ip_from )  gt_mtdm(i)%dzyxdm = g_mtdm(ifreq)%dzyxdm!2018.10.05
    if ( ip .eq. ip_from )  gt_mtdm(i)%dzyydm = g_mtdm(ifreq)%dzyydm!2018.10.05
  
    !# share     2017.09.03
    call sharecomplexcrsmatrix(gt_mtdm(i)%dzxxdm,ip_from,ip) ! m_shareformpi_joint.f90
    call sharecomplexcrsmatrix(gt_mtdm(i)%dzxydm,ip_from,ip) ! m_shareformpi_joint.f90
    call sharecomplexcrsmatrix(gt_mtdm(i)%dzyxdm,ip_from,ip) ! m_shareformpi_joint.f90
    call sharecomplexcrsmatrix(gt_mtdm(i)%dzyydm,ip_from,ip) ! m_shareformpi_joint.f90
  
   end do
  
  !#[2]##
    if (ip .eq. 0) write(*,*) "### SENDRESULTINV_MT END!! ###"
  
  return
  end
  
!#############################################
!# modified on 2022.10.20
!# coded on 2017.06.05
subroutine SENDRECVRESULT(resp5,tresp,ip,np,nfreq_act,nfreq_act_ip,g_freq_joint,nsr) !
use outresp
use shareformpi_joint ! 2021.12.25
use freq_mpi_joint    ! 2022.10.20
implicit none
!include 'mpif.h'
integer(4),    intent(in)    :: ip,np,nfreq_act,nfreq_act_ip
integer(4),    intent(in)    :: nsr                   ! 2017.09.03
type(freq_info_joint),intent(in) :: g_freq_joint      ! 2022.10.20
type(respdata),intent(in)    :: resp5(5,nsr,nfreq_act_ip) ! 2017.09.03
type(respdata),intent(inout) :: tresp(5,nsr,nfreq_act)    ! 2017.09.03
integer(4)                   :: errno,i,j,k,l,ip_from,ifreq_ip ! 2022.10.20

!#[1]##
do i=1,nfreq_act

  ip_from  = g_freq_joint%ip_from_act(i) ! 2022.10.20
  ifreq_ip = g_freq_joint%if_g2l_act(i)  ! 2022.10.20

  if ( ip .eq. ip_from ) then
   do k=1,nsr; do j=1,5 ! 20200807
    tresp(j,k,i) = resp5(j,k,ifreq_ip) ! 20200807
   end do ; end do     ! 20200807
  end if


  do j=1,5

   do k=1,nsr ! 2017.09.03
     call sharerespdata(tresp(j,k,i),ip_from) ! see m_shareformpi_ap 2017.09.03
   end do     ! 2017.09.03

  end do

end do

if ( ip .eq. 0 ) write(*,*) "### SENDRECVRESULT END!! ###"

return
end
!#############################################
!2022.10.20
subroutine SENDRECVIMP(imp_mt,timp_mt,ip,np,nfreq_mt,nfreq_mt_ip,g_freq_joint)!2017.09.03
  use outresp
  use shareformpi_joint ! 2021.12.25
  use freq_mpi_joint    ! 2022.10.20
  implicit none
  !include 'mpif.h'
  integer(4),    intent(in)    :: ip,np,nfreq_mt,nfreq_mt_ip
  type(freq_info_joint),intent(in) :: g_freq_joint  ! 2022.10.20
  type(respmt),  intent(in)    ::  imp_mt(nfreq_mt_ip) ! 2022.10.20
  type(respmt),  intent(inout) :: timp_mt(nfreq_mt)    ! 2022.10.20
  integer(4)                   :: errno,i,j,k,l,ip_from,ifreq             ! 2020.08.06
  
  !#[1]##
   do i=1,nfreq_mt
!    if (mod(i,np) .eq. 1 ) ip_from = -1
!    ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1

  ip_from  = g_freq_joint%ip_from_mt(i) ! 2022.10.20
  ifreq    = g_freq_joint%if_g2l_mt(i)  ! 2022.10.20

  
  if ( ip .eq. ip_from ) then
      timp_mt(i) = imp_mt(ifreq) ! 20200807
    end if
  
  
!  do j=1,5
!     do k=1,nsr ! 2017.09.03
  
  !if ( ip .eq. ip_from .and. j .eq. 3 ) then ! 20200806
  !  write(*,*) "-- before sharing --SENDRECVRESULT --"
  !  write(*,*)"ip_from",ip_from,"tresp i",i,"k",k,"j",j
  !  write(*,'(a,4g15.7)') "resp5 amp",(resp5(j,k,ifreq)%ftobsamp(l),l=1,resp5(j,k,ifreq)%nobs) ! 20200807
  !  write(*,'(a,4g15.7)') "resp5 phase",(resp5(j,k,ifreq)%ftobsphase(l),l=1,resp5(j,k,ifreq)%nobs) ! 20200807
  !  write(*,'(a,4g15.7)') "tresp amp",(tresp(j,k,i)%ftobsamp(l),l=1,tresp(j,k,i)%nobs)
  !  write(*,'(a,4g15.7)') "tresp phase",(tresp(j,k,i)%ftobsphase(l),l=1,tresp(j,k,i)%nobs)
  !  write(*,*) "--"
  ! end if
  
     call shareimpdata(timp_mt(i),ip_from) ! see m_shareformpi_ap 2022.01.02
  
  !if ( ip .eq. 0 .and. j .eq. 3 ) then ! 20200806
  ! write(*,*) "-- after sharing --SENDRECVRESULT --"
  ! write(*,*)"ip=0 tresp i",i,"k",k,"j",j
  ! write(*,'(a,4g15.7)') "tresp amp",(tresp(j,k,i)%ftobsamp(l),l=1,tresp(j,k,i)%nobs)
  ! write(*,'(a,4g15.7)') "tresp phase",(tresp(j,k,i)%ftobsphase(l),l=1,tresp(j,k,i)%nobs)
  ! write(*,*) "--"
  !end if
  
 !    end do     ! 2017.09.03
 !   end do
   end do
  
  if ( ip .eq. 0 ) write(*,*) "### SENDRECVIMP END!! ###" ! 2022.01.02
  
  return
  end

!#############################################
!# modified on 2017.07.13 to include nsr
!# coded on 2017.05.31
subroutine ALLOCATERESP(g_param,nsr,resp5,tresp,ip,nfreq,nfreq_ip)!2017.07.13
use outresp
use param
implicit none
type(param_forward),intent(in)    :: g_param
integer(4),         intent(in)    :: nsr ! 2017.07.13
integer(4),         intent(in)    :: nfreq,ip,nfreq_ip
type(respdata),     intent(inout) :: resp5(5,nsr,nfreq_ip) ! 2017.07.13
type(respdata),     intent(inout) :: tresp(5,nsr,nfreq)    ! 2017.07.13
integer(4) :: i,j,k,nobs

nobs = g_param%nobs

do j=1,nfreq_ip
 do i=1,5
  do k=1,nsr ! 2017.07.13
   CALL ALLOCATERESPDATA(nobs,resp5(i,k,j)) ! 2017.07.13
  end do     ! 2017.07.13
 end do
end do

do j=1,nfreq
 do i=1,5
  do k=1,nsr ! 2017.07.13
   CALL ALLOCATERESPDATA(nobs,tresp(i,k,j)) ! 2017.07.13
  end do     ! 2017.07.13
 end do
end do

if( ip .eq. 0) write(*,*) "### ALLOCATERESP END!! ###"
return
end
!#############################################
!# copied from ../src_3DMT/ALLOCATERESP on 2021.12.30
subroutine ALLOCATERESP_MT(nobs,nsr,resp,resp_mt,ip,nfreq)
  use outresp
  use param
  implicit none
  integer(4),         intent(in)    :: nobs
  integer(4),         intent(in)    :: nsr ! 2017.07.11
  integer(4),         intent(in)    :: nfreq,ip
  type(respdata),     intent(inout) :: resp(5,nsr,nfreq) !2017.07.11
  type(respmt),       intent(inout) :: resp_mt(nfreq)     !2021.09.14
  integer(4)                        :: i,j,k
  
  do j=1,nfreq
   CALL ALLOCATERESPMT(  nobs,resp_mt(j)    ) ! 2021.09.14 m_outresp.f90
   do i=1,5
    do k=1,nsr ! 2017.07.11
     CALL ALLOCATERESPDATA(nobs,resp(  i,k,j)) ! 2017.07.11
    end do     ! 2017.07.11
   end do
  end do
  
  if( ip .eq. 0) write(*,*) "### ALLOCATERESP END!! ###"
  return
  end
!#############################################
!# modified on 2018.10.05 for multiple components
!# modified on 2017.09.03 to include multiple sources
!# coded on 2017.05.31
subroutine CALOBSEBCOMP(fp,fs,nline,nsr,omega,coeffobs,resp5,g_param_joint)!2017.09.03
use matrix
use outresp
use param_jointinv ! 2018.10.05
implicit none
type(param_joint),intent(in)    :: g_param_joint   ! 2018.10.05
real(8),                intent(in)    :: omega
integer(4),             intent(in)    :: nline
integer(4),             intent(in)    :: nsr             ! 2017.09.03
complex(8),             intent(in)    :: fp(nline,nsr),fs(nline,nsr)! 2018.10.05
type(real_crs_matrix),  intent(in)    :: coeffobs(2,3)
type(respdata),         intent(inout) :: resp5(5,nsr)    ! 2017.09.03
integer(4)                            :: i,j,isr,icomp         ! 2017.09.03
integer(4),             dimension(5)  :: iflag_comp      ! 2018.10.05
complex(8), allocatable,dimension(:,:):: fpin,fsin       ! 20200729

allocate(fpin(nline,nsr),fsin(nline,nsr)) ! 20200729

!# set
iflag_comp = g_param_joint%iflag_comp ! 2018.10.08

!# calculate
icomp = 0 ! 1:Bx, 2: By, 3: Bz, 4: Ex, 5: Ey
do i=2,1,-1
 do j=1,3
 if ( icomp .eq. 5 ) exit
 icomp = icomp + 1

 if ( iflag_comp(icomp) .eq. 0 ) cycle

 if ( icomp .le. 3 ) then
  fpin(:,:) = fp(:,:)
  fsin(:,:) = fs(:,:)
 else if ( icomp .ge. 4) then
  fpin = - (0.d0,1.d0)*omega*fp ! E= -i*omega*A
  fsin = - (0.d0,1.d0)*omega*fs !
 end if

 do isr= 1, nsr ! 2017.09.03
  CALL CALOBSRESP(fpin(:,isr),fsin(:,isr),nline,coeffobs(i,j),resp5(icomp,isr))!bz 2018.10.08
 end do       ! 2017.09.03

 end do
end do

return
end subroutine

!############################################## subroutine CALRMS_AP
!# Coded 2017.06.08
subroutine CALRMS_AP(g_data,h_data,Cd,misfit,nrms)
use param_jointinv
use matrix
use caltime ! 2017.12.22
type(real_crs_matrix),intent(in)     :: Cd
type(data_vec_ap),    intent(in)     :: g_data ! obs 2018.10.08 phase are possibly modified
type(data_vec_ap),    intent(in)     :: h_data ! cal 2018.10.08
real(8),              intent(out)    :: misfit,nrms
real(8),allocatable,dimension(:)     :: dvec1,dvec2,dvec
integer(4)                           :: ndat,ndat1,ndat2,i
real(8)                              :: d ! 2017.09.06
type(watch) :: t_watch ! 2017.12.22

call watchstart(t_watch) ! 2017.12.22
!#[0]## set
ndat  = g_data%ndat
ndat1 = h_data%ndat
ndat2 = Cd%ntot
if (ndat .ne. ndat1 .or. ndat .ne. ndat2 ) goto 99
allocate(dvec1(ndat),dvec2(ndat),dvec(ndat))
dvec1 = g_data%dvec
dvec2 = h_data%dvec
dvec  = dvec2 - dvec1

!# check
if (.false.) then
 write(*,'(5x,a)') " obs          |  cal         | d^2           | d^2/e^2"
 write(*,*) "ndat1"
 do i=1,ndat1
  d = dvec2(i)-dvec1(i)
  write(*,'(i5,4g15.7)') i,dvec1(i),dvec2(i),d**2.d0,d**2.d0/Cd%val(i)
 end do
end if

!#[1]## cal rms
misfit = 0.d0
do i=1,ndat
 misfit = misfit + dvec(i)**2.d0 / Cd%val(i) ! dvec(i)**2/err(i)**2
end do
nrms = sqrt(misfit/dble(ndat)) ! 2017.12.22
misfit = 0.5d0 * misfit ! 2017.12.22

write(*,'(a,g15.7)') "            RMS =",misfit ! 2020.09.18
write(*,'(a,g15.7)') " Normarized RMS =",nrms   ! 2020.09.18
call watchstop(t_watch) ! 2017.12.22
!write(*,'(a,f9.4,a)') " ### CALRMS_AP END!! ### Time=",t_watch%time," [min]"! 2020.09.18
write(*,'(a)') " ### CALRMS_AP END!! ###" ! Time=",t_watch%time," [min]"! 2020.09.29

return
99 continue
write(*,*) "GEGEGE! ndat",ndat,".ne. ndat1",ndat1,"or .ne. ndat2",ndat2
stop

end subroutine
!############################################## subroutine CALRMS_MT
!# Coded 2022.01.04
subroutine CALRMS_MT(g_data_mt,h_data_mt,Cd_mt,misfit_mt,nrms_mt)
  use param_jointinv
  use matrix
  use caltime ! 2017.12.22
  type(real_crs_matrix),intent(in)     :: Cd_mt
  type(data_vec_mt),    intent(in)     :: g_data_mt ! obs
  type(data_vec_mt),    intent(in)     :: h_data_mt ! cal 
  real(8),              intent(out)    :: misfit_mt,nrms_mt
  real(8),allocatable,dimension(:)     :: dvec_mt1,dvec_mt2,dvec_mt
  integer(4)                           :: ndat_mt,ndat_mt1,ndat_mt2,i
  real(8)                              :: d ! 2017.09.06
  type(watch) :: t_watch ! 2017.12.22
  
  call watchstart(t_watch) ! 2017.12.22
  !#[0]## set
  ndat_mt  = g_data_mt%ndat_mt
  ndat_mt1 = h_data_mt%ndat_mt
  ndat_mt2 = Cd_mt%ntot
  if (ndat_mt .ne. ndat_mt1 .or. ndat_mt .ne. ndat_mt2 ) goto 99
  allocate(dvec_mt1(ndat_mt),dvec_mt2(ndat_mt),dvec_mt(ndat_mt))
  dvec_mt1 = g_data_mt%dvec_mt
  dvec_mt2 = h_data_mt%dvec_mt
  dvec_mt  = dvec_mt2 - dvec_mt1
  
  !# check
  if (.false.) then
   write(*,'(5x,a)') " obs          |  cal         | d^2           | d^2/e^2"
   write(*,*) "ndat1"
   do i=1,ndat_mt1
    d = dvec_mt2(i)-dvec_mt1(i)
    write(*,'(i5,4g15.7)') i,dvec_mt1(i),dvec_mt2(i),d**2.d0,d**2.d0/Cd_mt%val(i)
   end do
  end if
  
  !#[1]## cal rms
  misfit_mt = 0.d0
  do i=1,ndat_mt
   misfit_mt = misfit_mt + dvec_mt(i)**2.d0 / Cd_mt%val(i) ! dvec(i)**2/err(i)**2
  end do
  nrms_mt = sqrt(misfit_mt/dble(ndat_mt)) ! 2017.12.22
  misfit_mt = 0.5d0 * misfit_mt ! 2017.12.22
  
  write(*,'(a,g15.7)') "            RMS_MT =",misfit_mt ! 2020.09.18
  write(*,'(a,g15.7)') " Normarized RMS_MT =",nrms_mt   ! 2020.09.18
  call watchstop(t_watch) ! 2017.12.22

  write(*,'(a)') " ### CALRMS_mT END!! ###" ! Time=",t_watch%time," [min]"! 2020.09.29
  
  return
  99 continue
  write(*,*) "GEGEGE! ndat_mt",ndat_mt,".ne. ndat_mt1",ndat_mt1,"or .ne. ndat_mt2",ndat_mt2
  stop
  
  end subroutine
!############################################## subroutine GENDVEC_MT
!# coded on 2021.01.04
subroutine GENDVEC_MT(g_param_joint,timp_mt,nfreq_mt,h_data_mt) !2018.10.08
  use param_jointinv ! 2016.06.08
  use outresp     ! 2016.06.08
  implicit none
  integer(4),                    intent(in)     :: nfreq_mt
  type(param_joint),             intent(in)     :: g_param_joint ! 2017.09.04
  type(respmt),                  intent(in)     :: timp_mt(nfreq_mt)! 2022.01.04
  type(data_vec_mt),             intent(inout)  :: h_data_mt     ! cal data
  logical,allocatable, dimension(:,:,:,:)       :: data_avail_mt    ! 2018.10.05
  integer(4) :: iobs,i1,i2,ii,ifreq,i,isr,icomp,l,nobs_mt         ! 2017.09.03
  complex(8) :: c
  !#[0]## set amp and phase of bz
    nobs_mt       = g_param_joint%nobs_mt        ! 2022.01.04
    allocate(data_avail_mt(2,4,nobs_mt,nfreq_mt))   ! 2022.01.04
    data_avail_mt = g_param_joint%data_avail_mt  ! 2022.01.04
    
  !#[1]## cal g_data
  ii=0
  do ifreq=1,nfreq_mt
    do iobs=1,nobs_mt
     do icomp =1,4  ! zxx,zxy,zyx,zyy
      if (icomp .eq. 1) c=timp_mt(ifreq)%zxx(iobs)
      if (icomp .eq. 2) c=timp_mt(ifreq)%zxy(iobs)
      if (icomp .eq. 3) c=timp_mt(ifreq)%zyx(iobs)
      if (icomp .eq. 4) c=timp_mt(ifreq)%zyy(iobs)
      do l=1,2   ! real,phase
       if ( data_avail_mt(l,icomp,iobs,ifreq) ) then ! 2017.09.04
        ii = ii + 1
        if ( l .eq. 1) h_data_mt%dvec_mt(ii)  = real(c) ! 2018.10.05
        if ( l .eq. 2) h_data_mt%dvec_mt(ii)  = imag(c) ! 2018.10.05
       end if
      end do ! l loop
     end do ! icomp     2018.10.05
    end do  ! nobs
  end do    ! nfreq     2017.09.03
  
  !#[2]## out
  if (.false.) then
   do i=1,h_data_mt%ndat_mt
    write(*,*) i,"h_data_mt%dvec=",h_data_mt%dvec_mt(i)
   end do
  end if
  
  write(*,*) "### GENDVEC_MT END!! ###"
  return
  end subroutine
!############################################## subroutine GENDVEC_AP
!# modified for multiple component on 2018.10.05
!# modified for multiple sources   on 2017.09.04
!# modified on 2017.06.08
subroutine GENDVEC_AP(g_param_joint,nsr_inv,tresp,nfreq,g_data_ap,g_data) !2018.10.08
use param_jointinv ! 2016.06.08
use outresp     ! 2016.06.08
implicit none
integer(4),                          intent(in)     :: nfreq
integer(4),                          intent(in)     :: nsr_inv       ! 2017.09.04
type(param_joint),             intent(in)     :: g_param_joint ! 2017.09.04
type(data_vec_ap),                   intent(in)     :: g_data        ! 2018.10.08 obs data
type(respdata),                      intent(in)     :: tresp(5,nsr_inv,nfreq) ! 2017.09.03
type(data_vec_ap),                   intent(inout)  :: g_data_ap     ! cal data
real(8),         allocatable, dimension(:,:,:,:)    :: ampbe,phabe   ! 2018.10.05
logical,         allocatable, dimension(:,:,:,:,:)  :: data_avail    ! 2018.10.05
integer(4) :: iobs,i1,i2,ii,ifreq,i,isr,l,nobs         ! 2017.09.03
real(8)    :: amp,phase
integer(4) :: icomp         ! 2018.10.05
integer(4) :: iflag_comp(5) ! 2018.10.05
real(8)    :: pha0,pha1     ! 2018.10.08

!#[0]## set amp and phase of bz
  nobs       = g_param_joint%nobs              ! 2017.09.04
  allocate(data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2017.07.14
  data_avail = g_data_ap%data_avail            ! 2017.09.03
  allocate( ampbe(5,nobs,nfreq,nsr_inv) )      ! 2018.10.05
  allocate( phabe(5,nobs,nfreq,nsr_inv) )      ! 2018.10.05
  iflag_comp = g_param_joint%iflag_comp        ! 2018.10.05

  do i=1,nfreq
   do isr = 1,nsr_inv ! 2017.09.03
    do icomp = 1,5    ! 2018.10.05
     if (iflag_comp(icomp) .eq. 0) cycle                     ! 2018.10.05
     ampbe(icomp,:,i,isr) = tresp(icomp,isr,i)%ftobsamp(:)   ! 2017.09.03
     phabe(icomp,:,i,isr) = tresp(icomp,isr,i)%ftobsphase(:) ! 2017.09.03
    end do            ! 2018.10.05
   end do             ! 2017.09.03
  end do

!write(*,'(2i5,g15.7)') ((iobs,ifreq,ampbz(iobs,ifreq),ifreq=1,nfreq),iobs=1,nobs)

!#[1]## cal g_data
ii=0
do ifreq=1,nfreq
 do isr=1,nsr_inv ! 2017.09.03
  do iobs=1,nobs
   do icomp =1,5  ! 2018.10.05
    do l=1,2   ! amp,phase 2017.09.03
    if ( data_avail(l,icomp,ifreq,iobs,isr) ) then ! 2017.09.04
    ii = ii + 1
    if (l .eq. 1) then       ! amp   2017.09.03
!     write(*,'(a,4i5,1x,a,g15.7)') "icomp,iobs,ifreq,isr",icomp,iobs,ifreq,isr,"amp",ampbe(icomp,iobs,ifreq,isr)!2020.09.18
     g_data_ap%dvec(ii)  = dlog10(ampbe(icomp,iobs,ifreq,isr)) ! 2018.10.05
    else if ( l .eq. 2) then ! phase 2017.09.03
     pha1                = phabe(icomp,iobs,ifreq,isr)  ! 2018.10.08
     pha0                = g_data%dvec(ii)              ! 2018.10.08
     if ( abs( pha0 - pha1) .gt. abs(pha0 - (pha1+360.d0)) ) then      ! 2018.10.08
      pha1 = pha1 + 360.d0                                             ! 2018.10.08
     else if ( abs( pha0 - pha1) .gt. abs(pha0 - (pha1-360.d0)) ) then ! 2018.10.08
      pha1 = pha1 - 360.d0       ! 2018.10.08
     end if                      ! 2018.10.08
     g_data_ap%dvec(ii)  = pha1  ! 2018.10.08
    end if
   end if
   end do ! amp phase 2017.09.03
   end do ! icomp     2018.10.05
  end do  ! nobs
 end do   ! nsr_inv   2017.09.03
end do    ! nfreq     2017.09.03

!#[2]## out
if (.false.) then
 do i=1,g_data_ap%ndat
  write(*,*) i,"g_data_ap%dvec=",g_data_ap%dvec(i)
 end do
end if

write(*,*) "### GENDVEC_AP END!! ###"
return
end subroutine

!############################################## subroutine GENCD_AP
! Modified on 2017.09.04
! Coded on May 13, 2017
subroutine GENCD(g_data_ap,g_data_mt,CD,CD_mt) ! 2021.12.29
use matrix
use param_jointinv ! 2017.09.04
implicit none
type(data_vec_ap),      intent(in)  :: g_data_ap
type(data_vec_mt),      intent(in)  :: g_data_mt ! 2021.12.29
type(real_crs_matrix),  intent(out) :: CD, CD_mt ! 2021.12.29
real(8),allocatable,dimension(:)    :: error,error_mt ! 2021.12.29
integer(4) :: ndat,i
integer(4) :: ndat_mt

!#[1]# set
ndat    = g_data_ap%ndat
ndat_mt = g_data_mt%ndat_mt
write(*,*) "ndat_active=",ndat,"ndat_mt=",ndat_mt ! 2021.12.29
allocate(error(ndat),error_mt(ndat_mt)) ! 2021.12.29
error    = g_data_ap%error
error_mt = g_data_mt%error_mt ! 2021.12.29

!#[2]## gen Cd : diagonal matrix
CD%nrow  = ndat   ; CD_mt%nrow  = ndat_mt  
CD%ncolm = ndat   ; CD_mt%ncolm = ndat_mt   
CD%ntot  = ndat   ; CD_mt%ntot  = ndat_mt   
allocate(CD%stack(0:ndat),CD%item(ndat),CD%val(ndat))
allocate(CD_mt%stack(0:ndat_mt),CD_mt%item(ndat_mt),CD_mt%val(ndat_mt)) ! 2021.12.29
CD%stack(0)=0
do i=1,ndat
 CD%val(i)   = (error(i))**2.d0
 CD%stack(i) = i
 CD%item(i)  = i
end do
!
CD_mt%stack(0)=0
do i=1,ndat_mt
  CD_mt%val(i)   = (error_mt(i))**2.d0
  CD_mt%stack(i) = i
  CD_mt%item(i)  = i
 end do
 
write(*,*) "### GENCD END!! ###" ! 2021.12.29
return
end


!############################################## subroutine setNPHYS1INDEX2COND
!Coded on 2017.05.12
subroutine SETNPHYS1INDEX2COND(g_mesh,g_cond)
use param
use mesh_type
implicit none
type(mesh),      intent(in)    :: g_mesh
type(param_cond),intent(inout) :: g_cond
integer(4) :: nphys1,nphys2,i,ntet

  nphys2        = g_cond%nphys2 ! # of elements in 2nd physical volume (land)
  nphys1        = g_mesh%ntet - nphys2
  ntet          = g_mesh%ntet
  g_cond%nphys1 = nphys1
  g_cond%ntet   = ntet
  write(*,*) "nphys1=",nphys1,"nphys2=",nphys2,"ntet=",g_mesh%ntet
  do i=1,nphys2
   g_cond%index(i) = nphys1 +i ! element id for whole element space
  end do

return
end

!###################################################################
! copied from n_ebfem_bxyz.f90 on 2017.05.10
subroutine GENXYZMINMAX(em_mesh,g_param)
use param ! 2016.11.20
use mesh_type
implicit none
type(mesh),            intent(inout) :: em_mesh
type(param_forward),   intent(inout) :: g_param
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
real(8) :: xyz(3,em_mesh%node),xyzminmax(6)
integer(4) :: i
xyz = em_mesh%xyz ! normal
xmin=xyz(1,1) ; xmax=xyz(1,1)
ymin=xyz(2,1) ; ymax=xyz(2,1)
zmin=xyz(3,1) ; zmax=xyz(3,1)

do i=1,em_mesh%node
 xmin=min(xmin,xyz(1,i))
 xmax=max(xmax,xyz(1,i))
 ymin=min(ymin,xyz(2,i))
 ymax=max(ymax,xyz(2,i))
 zmin=min(zmin,xyz(3,i))
 zmax=max(zmax,xyz(3,i))
end do

xyzminmax(1:6)=(/xmin,xmax,ymin,ymax,zmin,zmax/)

!# set output
g_param%xyzminmax = xyzminmax
em_mesh%xyzminmax = xyzminmax    ! 2021.12.27

write(*,*) "### GENXYZMINMAX END!! ###"
return
end
!########################################### OUTMODEL
! modified on 2017.12.20
subroutine OUTMODEL(g_param_joint,g_model,g_mesh,ite,ialpha)! 2017.09.11
use mesh_type
use modelpart
use param_jointinv ! 2017.08.31
use caltime     ! 2017.12.25
implicit none
integer(4),              intent(in) :: ite,ialpha       ! 2017.09.11
type(param_joint),       intent(in) :: g_param_joint    ! 2017.08.31
type(model),             intent(in) :: g_model
type(mesh),              intent(in) :: g_mesh
character(50)                       :: modelfile, head   ! 2017.12.20
character(50)                       :: connectfile      ! 2017.12.20
integer(4)                          :: i,j,nphys2,nmodel,npoi,nlin,ntri,ishift
real(8),   allocatable,dimension(:) :: rho_model, logrho_model
integer(4),allocatable,dimension(:) :: ele2model,index
integer(4)   :: ialphaflag,nhead  ! 2017.09.11
character(2) :: num
character(1) :: num2              ! 2017.09.11
type(watch)  :: t_watch           ! 2017.12.25

call watchstart(t_watch) ! 2017.12.25
!#[0]## set
head         = g_param_joint%outputfolder ! 2017.09.03
nhead        = len_trim(head)             ! 2017.09.11
nphys2       = g_model%nphys2
nmodel       = g_model%nmodel
allocate(index(nphys2),rho_model(nmodel),ele2model(nphys2))
allocate(logrho_model(nmodel))
index        = g_model%index
rho_model    = g_model%rho_model
logrho_model = g_model%logrho_model
ele2model    = g_model%ele2model
npoi         = g_mesh%npoi
nlin         = g_mesh%nlin
ntri         = g_mesh%ntri
ialphaflag   = g_param_joint%ialphaflag            ! 2017.09.11
 if ( ialphaflag .eq. 1 .and. ialpha .ne. 0 ) then ! L-curve 2017.09.11
  write(num2,'(i1)')  ialpha ! 2017.09.11
  head      = head(1:nhead)//"a"//num2(1:1)//"/"
  nhead     = len_trim(head)
 end if

!#[1]## when ite = 0, output connection from nmodel to nphys2
 if ( ite .eq. 0 ) then
  connectfile=head(1:nhead)//"model_connect.dat"
  open(1,file=connectfile)
  write(1,'(2i10)') nphys2,nmodel
  ishift = npoi + nlin + ntri
  do i=1,nphys2
   write(1,'(2i10)')  ishift+index(i),ele2model(i)
  end do
  close(1)
 end if

!#[3]## output rho
 write(num,'(i2.2)') ite
 modelfile =head(1:nhead)//"model"//num(1:2)//".dat" ! 2017.12.20
 open(1,file=modelfile)       ! 2017.12.20
  write(1,'(i10)') nmodel
  do j=1,nmodel
!   if (ite .eq. 2 ) write(*,*) "j",j,"logrho_model(j)",logrho_model(j) ! 2017.12.22 check
   write(1,'(e15.7)') 10**logrho_model(j) !2017.12.20
  end do
 close(1)

call watchstop(t_watch) ! 2017.12.25
write(*,'(a,f9.4,a)') " ### OUTMODEL END!! ### Time =",t_watch%time," [min]"!2020.09.18
return
end
!########################################### OUTCOND
! Output folder is changed on 2017.09.04
! Coded on 2017.05.18
subroutine OUTCOND(g_param_joint,g_cond,g_mesh,ite,ialpha) ! 2017.09.11
use mesh_type
use param
use param_jointinv ! 2017.09.04
implicit none
integer(4),             intent(in) :: ite,ialpha    ! 2017.09.11
type(param_joint),intent(in) :: g_param_joint ! 2017.09.04
type(param_cond),       intent(in) :: g_cond
type(mesh),             intent(in) :: g_mesh
character(50) :: condfile,head ! 2017.07.25
integer(4) :: j, nphys2,nmodel,npoi,nlin,ntri,ishift,nphys1
real(8),allocatable,dimension(:)    :: rho,sigma
integer(4),allocatable,dimension(:) :: index
integer(4)   :: nhead ! 2017.09.11
character(2) :: num
character(1) :: num2       ! 2017.09.11
integer(4)   :: ialphaflag ! 2017.09.11

!#[0]## set
head         = g_param_joint%outputfolder ! 2017.09.04
nhead        = len_trim(head) ! 2017.09.11s
nphys1       = g_cond%nphys1
nphys2       = g_cond%nphys2
allocate(index(nphys2),rho(nphys2),sigma(nphys2))
index        = g_cond%index
rho          = g_cond%rho
sigma        = g_cond%sigma
npoi         = g_mesh%npoi
nlin         = g_mesh%nlin
ntri         = g_mesh%ntri
ialphaflag   = g_param_joint%ialphaflag ! 2017.09.11
write(*,*) "nphys1=",nphys1
write(*,*) "nphys2=",nphys2
write(num,'(i2.2)') ite
 if ( ialphaflag .eq. 1 .and. ialpha .ne. 0 ) then ! L-curve 2017.09.11
  write(num2,'(i1)')  ialpha ! 2017.09.11
  head      = head(1:nhead)//"a"//num2(1:1)//"/"
  nhead     = len_trim(head)
 end if

!#[1]## output rho
 condfile = head(1:nhead)//"cond"//num(1:2)//".msh"
 open(1,file=condfile)

 !# standard info
! CALL MESHOUT(1,g_mesh)

 write(1,'(a)') "$MeshFormat"     ! 2017.09.13
 write(1,'(a)') "2.2 0 8"        ! 2017.09.13
 write(1,'(a)') "$EndMeshFormat"  ! 2017.09.13
 write(1,'(a)') "$ElementData"
 write(1,'(a)') "1"
 write(1,'(a)') '"A rho model view"'
 write(1,'(a)') "1"
 write(1,'(a)') "0.0"
 write(1,'(a)') "3"
 write(1,'(a)') "0"
 write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
 write(1,'(i10)') nphys2
 ishift = npoi + nlin + ntri
 do j=1,nphys2
!  write(*,*) "j=",j,"nphys2=",nphys2,"ele2model(j)=",ele2model(j),"nmodel=",nmodel
  write(1,*) ishift+index(j),1.d0/sigma(j)
 end do
 write(1,'(a)') "$EndElementData"
close(1)

write(*,*) "### OUTCOND END!! ###"
return
end
!####################################################### PREPOBSCOEFF
! copied from ../solver/n_ebfem_bxyz.f90 on 2017.05.14
! copied from ../FEM_node/n_bzfem.f90
! adjusted to edge-FEM code
subroutine PREPOBSCOEFF(g_param,h_mesh,l_line,coeffobs)
use mesh_type
use line_type
use param
use matrix
use fem_edge_util
implicit none
type(mesh),           intent(in)  :: h_mesh
type(line_info),      intent(in)  :: l_line
type(param_forward),  intent(in)  :: g_param
type(real_crs_matrix),intent(out) :: coeffobs(2,3)     ! (1,(1,2,3)) for edge (x,y,z)
real(8) :: x3(3),a4(4),r6(6),len(6),w(3,6),elm_xyz(3,4),v
real(8) :: coeff(3,6),coeff_rot(3,6)
real(8) :: wfe(3,6) ! face basis function * rotation matrix
integer(4) :: iele,i,ii,j,k,l,jj,n6(6),ierr

!#[1]# allocate coeffobs
do i=1,2 ; do j=1,3
 coeffobs(i,j)%nrow=g_param%nobs
 coeffobs(i,j)%ncolm=l_line%nline ! added on 2017.05.15
 coeffobs(i,j)%ntot=g_param%nobs * 6 ! for lines of tetrahedral mesh
 allocate(coeffobs(i,j)%stack(0:coeffobs(i,j)%nrow))
 allocate( coeffobs(i,j)%item(  coeffobs(i,j)%ntot))
 allocate(  coeffobs(i,j)%val(  coeffobs(i,j)%ntot))
 coeffobs(i,j)%stack(0)=0
end do ; end do

!#[2]# find element and set values to coeffobs
ii=0
do i=1,g_param%nobs
 if (g_param%lonlatflag .eq. 2 ) then ! xyobs is already set
  x3(1:3)=(/g_param%xyzobs(1,i),g_param%xyzobs(2,i),g_param%xyzobs(3,i)/) ! [km]
  call FINDELEMENT0(x3,h_mesh,iele,a4) ! see m_mesh_type.f90
!  write(*,*) "x3(1:3,i)=",x3(1:3)
!  write(*,*) "ieleobs(i)=",iele
!  write(*,*) "coeff(i,1:4)=",coeff(1:4)
  do j=1,2
   do k=1,3
   coeffobs(j,k)%stack(i)=coeffobs(j,k)%stack(i-1) + 6
   end do
  end do
  do j=1,4
    elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j))
  end do
  CALL EDGEBASISFUN(elm_xyz,x3,w  ,len,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
  CALL FACEBASISFUN(elm_xyz,x3,wfe,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
  do j=1,6
   ! coeff is values for line-integrated value
   coeff(1:3,j)     =   w(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
   coeff_rot(1:3,j) = wfe(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
   r6(j)=abs(l_line%n6line(iele,j))*1.d0
  end do
  n6(1:6)=(/1,2,3,4,5,6/)
  call SORT_INDEX(6,n6,r6) !sort n4 index by r4 : see sort_index.f90
  do j=1,6
   jj=abs(l_line%n6line(iele,n6(j))) ! jj is global node id
   ii=ii+1                  ! ii is entry id for coeffobs matrix
   do l=1,3 ! l for component
    coeffobs(1,l)%item(ii)=jj
    coeffobs(2,l)%item(ii)=jj
    coeffobs(1,l)%val(ii)=    coeff(l,n6(j)) ! for x component
    coeffobs(2,l)%val(ii)=coeff_rot(l,n6(j)) ! for x component
   end do
  end do
 else
  write(*,*) "GEGEGE! in PREPOBSCOEFF"
  write(*,*) "g_param%lonlatflag",g_param%lonlatflag,"should be 2 here."
  stop
 end if
end do

!write(*,*) "### PREPOBSCOEFF END!! ###"   ! commented out on 2017.09.04
return

end
!####################################################### PREPOBSCOEFF
! copied from ../solver/n_ebfem_bxyz.f90 on 2021.09.15
subroutine PREPOBSCOEFF_MT(g_param,h_mesh,l_line,coeffobs,ip) ! 2022.01.05
  use mesh_type
  use line_type
  use param_mt ! 2021.12.15
  use matrix
  use fem_edge_util
  implicit none
  type(mesh),              intent(in)  :: h_mesh
  type(line_info),         intent(in)  :: l_line
  type(param_forward_mt),  intent(in)  :: g_param ! 2021.12.15
  integer(4),              intent(in)  :: ip      ! 2022.01.05
  type(real_crs_matrix),   intent(out) :: coeffobs(2,3)!(1,(1,2,3)) for edge (x,y,z)
  real(8) :: x3(3),a4(4),r6(6),len(6),w(3,6),elm_xyz(3,4),v
  real(8) :: coeff(3,6),coeff_rot(3,6)
  real(8) :: wfe(3,6) ! face basis function * rotation matrix
  integer(4) :: iele,i,ii,j,k,l,jj,n6(6),ierr
  
  !#[1]# allocate coeffobs
  do i=1,2 ; do j=1,3
   coeffobs(i,j)%nrow=g_param%nobs
   coeffobs(i,j)%ntot=g_param%nobs * 6 ! for lines of tetrahedral mesh
   allocate(coeffobs(i,j)%stack(0:coeffobs(i,j)%nrow))
   allocate( coeffobs(i,j)%item(  coeffobs(i,j)%ntot))
   allocate(  coeffobs(i,j)%val(  coeffobs(i,j)%ntot))
   coeffobs(i,j)%stack(0)=0
  end do ; end do
  
  !#[2]# find element and set values to coeffobs
  ii=0
  do i=1,g_param%nobs
   if (g_param%lonlatflag .eq. 2 ) then ! xyobs is already set
    x3(1:3)=(/g_param%xyzobs(1,i),g_param%xyzobs(2,i),g_param%xyzobs(3,i)/) ! [km]
    call FINDELEMENT0(x3,h_mesh,iele,a4) ! see m_mesh_type.f90
  !  write(*,*) "x3(1:3,i)=",x3(1:3)
  !  write(*,*) "ieleobs(i)=",iele
  !  write(*,*) "coeff(i,1:4)=",coeff(1:4)
    do j=1,2
     do k=1,3
     coeffobs(j,k)%stack(i)=coeffobs(j,k)%stack(i-1) + 6
     end do
    end do
    do j=1,4
      elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j))
    end do
    CALL EDGEBASISFUN(elm_xyz,x3,w  ,len,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
    CALL FACEBASISFUN(elm_xyz,x3,wfe,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
    do j=1,6
     ! coeff is values for line-integrated value
     coeff(1:3,j)     =   w(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
     coeff_rot(1:3,j) = wfe(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
     r6(j)=abs(l_line%n6line(iele,j))*1.d0
    end do
    n6(1:6)=(/1,2,3,4,5,6/)
    call SORT_INDEX(6,n6,r6) !sort n4 index by r4 : see sort_index.f90
    do j=1,6
     jj=abs(l_line%n6line(iele,n6(j))) ! jj is global node id
     ii=ii+1                  ! ii is entry id for coeffobs matrix
     do l=1,3 ! l for component
      coeffobs(1,l)%item(ii)=jj
      coeffobs(2,l)%item(ii)=jj
      coeffobs(1,l)%val(ii)=    coeff(l,n6(j)) ! for x component
      coeffobs(2,l)%val(ii)=coeff_rot(l,n6(j)) ! for x component
     end do
    end do
   else
    write(*,*) "GEGEGE! in PREPOBSCOEFF_MT"
    write(*,*) "g_param%lonlatflag",g_param%lonlatflag,"should be 2 here."
    stop
   end if
  end do
  
  if (ip .eq. 0) write(*,*) "### PREPOBSCOEFF_MT END!! ###" ! 2021.01.05
  return
  
  end
  
!######################################## ALLOCATERESP
! modified on 2017.05.18
subroutine PREPRESPFILES(g_param,files,resp5,nfreq)
use outresp
use param
implicit none
integer(4),         intent(in) :: nfreq ! 2017.05.18
type(param_forward),intent(in) :: g_param
type(obsfiles),intent(inout)   :: files
type(respdata),intent(inout)   :: resp5(5,nfreq)
integer(4) :: nobs,i,j

!#[1]## set
nobs = g_param%nobs

!#[2]## make files and allocate
CALL MAKEOBSFILES(g_param,files)
!CALL OPENOBSFILES(files)  !commented out on 2017.05.18

return
end
!######################################## CALOBSRESP
!# Coded on Nov. 21, 2015
!# This calculates the output b fields and output results
subroutine CALOBSRESP_3DMT(ft,nline,coeffobs,resp) ! fp, fs -> ft 2021.09.15
  use matrix
  use outresp
  implicit none
  integer(4),           intent(in)    :: nline
  complex(8),           intent(in)    :: ft(nline) ! 2021.09.15
  type(real_crs_matrix),intent(in)    :: coeffobs ! see m_matrix.f90
  type(respdata),       intent(inout) :: resp   ! see m_outresp.f90
  complex(8),allocatable,dimension(:) :: ftobs
  real(8) :: amp,phase
  integer(4) :: i
  allocate(ftobs(resp%nobs)) ! 2021.09.15
  
  !#[1]## generate btotal
  
  !#[2]## calculate bp,bs,bt at observation points
  CALL mul_matcrs_cv(coeffobs,ft(1:nline),nline,ftobs) ! see m_matrix.f90
  
  !#[3]## cal b3 comp and output
  do i=1,resp%nobs
   resp%ftobsamp(i)  =amp  (ftobs(i)) ! amp of bz
   resp%ftobsphase(i)=phase(ftobs(i)) ! phase of bz
   resp%ftobs(i)     =ftobs(i)        ! 2021.09.15
  end do
  
!  write(*,*) "### CALOBSRESP_3DMT END!! ###"
  return
  end
!######################################## CALOBSRESP
!# Coded on Nov. 21, 2015
!# This calculates the output b fields and output results
subroutine CALOBSRESP(fp,fs,doftot,coeffobs,resp)
use matrix
use outresp
implicit none
integer(4),           intent(in)    :: doftot
complex(8),           intent(in)    :: fp(doftot)
complex(8),           intent(in)    :: fs(doftot)
type(real_crs_matrix),intent(in)    :: coeffobs ! see m_matrix.f90
type(respdata),       intent(inout) :: resp   ! see m_outresp.f90
complex(8),allocatable,dimension(:) :: ft
complex(8),allocatable,dimension(:) :: fpobs,fsobs,ftobs
real(8) :: amp,phase
integer(4) :: i
allocate(ft(doftot))
allocate(fpobs(resp%nobs))
allocate(fsobs(resp%nobs),ftobs(resp%nobs))

!#[1]## generate btotal
ft(1:doftot)=fs(1:doftot)+fp(1:doftot) ! [nT*km] or [mV/km * km]
!write(*,*) "bt is created!"

!#[2]## calculate bp,bs,bt at observation points
CALL mul_matcrs_cv(coeffobs,fp(1:doftot),doftot,fpobs) ! see m_matrix.f90
CALL mul_matcrs_cv(coeffobs,fs(1:doftot),doftot,fsobs) ! see m_matrix.f90
CALL mul_matcrs_cv(coeffobs,ft(1:doftot),doftot,ftobs) ! see m_matrix.f90
!write(*,*) "bpobs,bsobs,btobs is created!"

!#[3]## cal b3 comp and output
do i=1,resp%nobs
 resp%fpobsamp(i)  =amp  (fpobs(i)) ! amplitude of bz primary
 resp%fpobsphase(i)=phase(fpobs(i)) ! phase of primary fields
 resp%ftobsamp(i)  =amp  (ftobs(i)) ! amp of bz
 resp%ftobsphase(i)=phase(ftobs(i)) ! phase of bz
 resp%fsobsamp(i)  =amp  (fsobs(i)) ! amp of bz
 resp%fsobsphase(i)=phase(fsobs(i)) ! phase of bz
end do

!write(*,*) "### CALOBSRESP END!! ###"
return
end
!######################################## function phase
function phase(c) ! [deg]
implicit none
complex(8),intent(in) :: c
real(8) :: phase
real(8),parameter :: pi=4.d0*datan(1.d0), r2d=180.d0/pi
 phase=datan2(dimag(c),dreal(c))*r2d
 return
end
!######################################## function amp
function amp(c)
implicit none
complex(8),intent(in) :: c
real(8) :: amp
 amp=dsqrt(dreal(c)**2.d0 + dimag(c)**2.d0)
 return
end

!################################################################# PREPZSRCOBS
! modified on 2017.07.11 to include multi source
!# Coded on 2017.02.21
subroutine PREPZSRCOBS(h_mesh,g_param,s_param)
use param
use mesh_type
use triangle
implicit none
type(mesh),         intent(inout)     :: h_mesh  ! deallocated at the end 2017.05.15
type(param_forward),intent(inout)     :: g_param
type(param_source), intent(inout)     :: s_param
type(grid_list_type)                  :: glist
integer(4)                            :: nobs,nx,ny
real(8),   allocatable,dimension(:,:) :: xyzobs,xyz
integer(4),allocatable,dimension(:,:) :: n3k
real(8),   allocatable,dimension(:)   :: znew
real(8)    :: a3(3)
integer(4) :: iele,n1,n2,n3,j,k,ntri
integer(4) :: nsr                                ! 2017.07.14
real(8)    :: xyzminmax(6),zorigin
real(8),   allocatable,dimension(:,:) :: xs1,xs2 ! 2017.07.14

!#[0]## cal xyzminmax of h_mesh
  CALL GENXYZMINMAX(h_mesh,g_param)

!#[1]## set
nsr       = s_param%nsource     ! 2017.07.14
allocate(xs1(3,nsr),xs2(3,nsr)) ! 2017.07.14
allocate(xyz(3,h_mesh%node),n3k(h_mesh%ntri,3))
allocate(xyzobs(3,g_param%nobs))
allocate(znew(g_param%nobs))
nobs      = g_param%nobs
xyz       = h_mesh%xyz    ! triangle mesh
n3k       = h_mesh%n3
ntri      = h_mesh%ntri
xyzobs    = g_param%xyzobs
xs1       = s_param%xs1
xs2       = s_param%xs2
xyzminmax = g_param%xyzminmax


!#[2]## cal z for nobsr
nx=300;ny=300
CALL allocate_2Dgrid_list(nx,ny,ntri,glist)   ! see m_mesh_type.f90
CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
CALL classifytri2grd(h_mesh,glist)   ! classify ele to glist,see


!#[3] search for the triangle including (x1,y1)
do j=1,nobs
    call findtriwithgrid(h_mesh,glist,xyzobs(1:2,j),iele,a3)
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    znew(j) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xyzobs(3,j)
!
    write(*,10) " Obs # ",j, g_param%obsname(j)    ! 2020.09.29
    write(*,11) " x,y =",xyzobs(1,j)," , ",xyzobs(2,j)," [km]" ! 2022.10.14
    write(*,11) "   z =",xyzobs(3,j)," ->",znew(j),    " [km]" ! 2020.09.17
end do

!#[3-2]## source z
 do k=1,nsr                                                               ! 2017.07.14
!#  start point
    call findtriwithgrid(h_mesh,glist,xs1(1:2,k),iele,a3)                 ! 2017.07.14
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs1(3,k) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs1(3,k) ! 2017.07.14
    !
    write(*,9 ) " Src # = ",k    ! 2020.09.17
    write(*,*)   "Start point:"  ! 2020.09.17
    write(*,11) " x,y =",xs1(1,k),        " , ",xs1(2,k)," [km]"   ! 2020.09.17
    write(*,11) "   z =",s_param%xs1(3,k)," ->",xs1(3,k)," [km]"   ! 2020.09.17

!#  end point
    call findtriwithgrid(h_mesh,glist,xs2(1:2,k),iele,a3)                 ! 2017.07.14
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs2(3,k) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs2(3,k) ! 2017.07.14
    !
    write(*,*)   "End point:"    ! 2020.09.17
    write(*,11) " x,y =",xs2(1,k),        " , ",xs2(2,k)," [km]" ! 2020.09.17
    write(*,11) "   z =",s_param%xs2(3,k)," ->",xs2(3,k)," [km]" ! 2020.09.17
end do
9  format(a,i3)
10 format(a,i3,a)
11 format(a,f8.3,a,f8.3,a) ! 2020.09.17

!#[4]## set znew to xyz_r
    g_param%xyzobs(3,1:nobs) = znew(1:nobs)
    s_param%xs1(3,:) = xs1(3,:)                 ! 2017.07.14
    s_param%xs2(3,:) = xs2(3,:)                 ! 2017.07.14

!#[5]## kill mesh for memory 2017.05.15
    call killmesh(h_mesh) ! see m_mesh_type.f90

write(*,*) "### PREPZSRCOBS  END!! ###"
return
end
