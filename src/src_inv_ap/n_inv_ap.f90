! Modified on 2017.08.31 for multiple sources
! Coded on 2017.06.08
! to use amp and phase of bz
program inversion_ap
!
use param
use param_apinv      ! 2017.08.31
use mesh_type
use line_type
use face_type
use matrix
use modelpart
use constants        ! added on 2017.05.14
use iccg_var_takuto  ! added on 2017.05.14
use outresp
use jacobian_ap      ! added on 2017.06.08
use shareformpi_ap   ! added on 2017.06.03
use freq_mpi         ! added on 2017.06.07FF
use caltime          ! added on 2017.09.11
use spectral         ! added on 2017.12.12
use modelroughness   ! added on 2017.12.12
!
implicit none
type(param_forward)    :: g_param      ! see m_param.f90
type(param_source)     :: sparam       ! see m_param.f90
type(param_cond)       :: g_cond       ! see m_param.f90 : goal conductivity
type(param_cond)       :: h_cond       ! see m_param.f90 : initial structure
type(param_cond)       :: r_cond       ! see m_param.f90 : ref 2017.07.19
type(param_apinversion):: g_param_apinv! see m_param_apinv.f90
type(mesh)             :: g_mesh       ! see m_mesh_type.f90
type(mesh)             :: h_mesh       ! z file; see m_mesh_type.f90
type(line_info)        :: g_line       ! see m_line_type.f90
type(face_info)        :: g_face       ! see m_line_type.f90
type(modelpara)        :: g_modelpara  ! see m_modelpart.f90
type(model)            :: g_model_ref  ! see m_modelpart.f90
type(model)            :: g_model_ini  ! see m_modelpart.f90 2017.08.31
type(model)            :: h_model      ! see m_modelpart.f90 2017.05.17
type(model)            :: pre_model    ! see m_modelpart.f90, 2017.08.31
type(data_vec_ap)      :: g_data       ! observed data  ; see m_param_apinv.f90
type(data_vec_ap)      :: h_data       ! calculated data; see m_param_apinv.f90
type(freq_info)        :: g_freq       ! see m_freq_mpi.f90
type(real_crs_matrix)  :: RTR,CD       ! see m_matrix.f90 2017.12.25
type(real_crs_matrix)  :: BMI          ! see m_matrix.f90 2017.12.13
type(real_crs_matrix)  :: BM           ! see m_matrix.f90 2017.12.25
type(real_crs_matrix)  :: PT(5)        ! (1:5)[nobs,nlin] ! 2018.10.04, Bx,By,Bz,Ex,Ey
type(real_crs_matrix)  :: R,RI         ! see m_matrix.f90 2017.12.18
type(global_matrix)    :: A            ! see m_iccg_var_takuto.f90
type(real_crs_matrix)  :: coeffobs(2,3)! see m_matrix.f90 ; 1 for edge, 2 for face
integer(4)                                  :: ite,i,j,errno
integer(4)                                  :: nline, ntet, nobs, ndat ! 2017.08.31
integer(4)                                  :: nsr_inv, nmodel         ! 2017.09.04
integer(4)                                  :: nmodelactive            ! 2018.06.25
real(8)                                     :: omega, freq_ip
integer(4)                                  :: nfreq, nfreq_ip         ! for MPI
real(8)                                     :: nrms, nrms0, alpha      ! 2017.09.08
real(8)                                     :: misfit    ! 2017.12.22
complex(8),allocatable,dimension(:,:)       :: fp,fs     ! (nline,nsr_inv) 2017.08.31
type(obsfiles)                              :: files     ! see m_outresp.f90
type(respdata),allocatable,dimension(:,:,:) :: resp5     ! 2017.08.31resp5(5,nsr,nfreq_ip)
type(respdata),allocatable,dimension(:,:,:) :: tresp     ! 2017.08.31 tresp(5,nsr,nfreq)
integer                                     :: access    ! 2017.05.15
type(complex_crs_matrix)                    :: ut(5)     ! 2018.10.04, Bx,By,Bz,Ex,Ey
type(amp_phase_dm),allocatable,dimension(:) :: g_apdm    ! 2017.06.07
type(amp_phase_dm),allocatable,dimension(:) :: gt_apdm   ! 2017.06.07
type(real_crs_matrix)                       :: JJ ! Jacobian matrix, 2017.05.17
real(8)                                     :: rough1,rough2 ! 2017.12.13
integer(4)                                  :: ialpha    ! 2017.07.19
real(8)                                     :: nrms_init ! 2017.07.19
character(50)                               :: head      ! 2017.07.25
character(1)                                :: num2      ! 2017.09.11
type(watch)                                 :: t_watch   ! 2017.09.11
type(watch)                                 :: t_watch0  ! 2018.03.02
real(8)                                     :: nrms_ini  ! 2018.06.25
real(8)                                     :: frms      ! 2018.06.25
integer(4),    allocatable,dimension(:,:)   :: n4        ! 2022.01.27
integer(4)                                  :: node      ! 2022.01.27
!
integer(4) :: ip,np, itemax = 20, iflag
integer(4) :: kmax ! maximum lanczos procedure 2017.12.13
integer(4) :: nalpha, ialphaflag  ! 2017.09.08
integer(4) :: itype_roughness     ! 2017.12.13
!# explanation ##
 !#---------------------------------------------------------
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

!#[0]## read parameters
  CALL READPARAM(g_param,sparam,g_cond) ! include READCOND for initial model
  CALL READPARAINVAP(g_param_apinv,g_modelpara,g_param,sparam,g_data)!m_param_apinv.f90 2018.06.26

!#[1]## read mesh
  CALL READMESH_TOTAL(g_mesh,g_param%g_meshfile)
  if ( access( g_param%z_meshfile, " ") .eq. 0 ) then! if exist, 2017.05.15
    CALL READMESH_TOTAL(h_mesh,g_param%z_meshfile)
    CALL PREPZSRCOBS(h_mesh,g_param,sparam)          ! see below, include kill h_mesh
  end if
  CALL GENXYZMINMAX(g_mesh,g_param)                  ! see below
  CALL READLINE(g_param%g_lineinfofile,g_line)       ! see m_line_type.f90
!  CALL READFACE(g_param_apinv%g_faceinfofile,g_face) ! see m_face_type.f90

!#[2]## Make face of 3D mesh
  nline   = g_line%nline      ! 2022.01.27
  ntet    = g_mesh%ntet       ! 2022.01.27
  node    = g_mesh%node       ! 2022.01.27
  nfreq   = g_param%nfreq     ! 2022.01.27
  n4      = g_mesh%n4         ! allocate n4(ntet,4) here 2020.10.04
  CALL MKFACE(  g_face,node,ntet,4, n4) ! make g_line        2022.01.27
  CALL MKN4FACE(g_face,node,ntet,   n4) ! make g_line%n6line 2022.01.27
  !  CALL READFACE(g_param_joint%g_faceinfofile,g_face) ! see m_face_type.f90
  CALL FACE2ELEMENT(g_face)                          ! cal g_face%face2ele(1:2,1:nface)

!#[2]## prepare initial and ref cond, note sigma_air of ref and init is from g_cond
  CALL DUPLICATESCALARCOND(g_cond,r_cond)            ! see m_param.f90,  2017.08.31
  CALL DUPLICATESCALARCOND(g_cond,h_cond)            ! see m_param.f90,  2017.08.31
  CALL deallocatecond(g_cond)                        ! see m_param.f90,  2017.08.31
  CALL READREFINITCOND(r_cond,h_cond,g_param_apinv,g_mesh)  ! 2018.10.04 homo init is added

!#[3]## generate model space
  CALL genmodelspace(g_mesh,g_modelpara,g_model_ref,g_param,h_cond)! m_modelpart.f90
  CALL modelparam(g_model_ref)                       ! 2017.08.31 see m_modelpart.f90
  CALL OUTMODEL(g_param_apinv,g_model_ref,g_mesh,0,0)! 2017.09.11 for model parameter
  CALL assignmodelrho(r_cond,g_model_ref)            ! 2017.08.31 renew g_model_ref
  g_model_ini = g_model_ref                          ! 2017.08.31
  CALL assignmodelrho(h_cond,g_model_ini)            ! 2017.08.31

!#[4]## cal BMI for SM
  itype_roughness = g_param_apinv%itype_roughness  ! 2017.12.12 see m_param_apinv.f90
  if     ( itype_roughness == 1 )  then  ! SM: smoothest model
   CALL GENBMI_SM(BM,BMI,g_face,g_mesh,g_model_ref)  ! 2017.06.14 m_modelroughenss.f90
  elseif ( itype_roughness == 3 )  then  ! MSG: minimum support gradient, 2017.12.18
   CALL GENRI_MSG(R,RI,g_face,g_mesh,g_model_ref) ! 2017.12.25 m_modelroughenss.f90
  end if

!#[5]## cal Cd
 CALL GENCD_AP(g_data,CD)                          ! see below

end if !#################################################################   ip = 0 end

!#[3]## share params, mesh, and line (see m_shareformpi.f90)
  CALL shareapinv(g_param,sparam,h_cond,g_mesh,g_line,g_param_apinv,g_model_ini,ip) ! 2018.10.04

!#[6]## cal Pt: matrix for get data from simulation results
  CALL PREPOBSCOEFF(g_param,g_mesh,g_line,coeffobs) ! for x,y,z component, see below
  CALL PREPAREPT(coeffobs,g_param_apinv,PT)         ! 2018.10.04 see below

!#[7]## Preparation for global stiff matrix
  nline        = g_line%nline
  nobs         = g_param%nobs
  ntet         = g_mesh%ntet
  nmodel       = g_model_ini%nmodel             ! 2018.06.25
  nmodelactive = g_model_ini%nmodelactive       ! 2018.06.25
  nsr_inv      = g_param_apinv%nsr_inv          ! 2017.08.31
  CALL SET_ICCG_VAR(ntet,nline,g_line,A,ip)     ! see below, 2017.06.05
!  write(*,'(a,2i6,a,i3)') "nmodel,nmodelactive",nmodel,nmodelactive,"ip",ip ! 2018.06.26

!#[8]## set frequency
  call SETFREQIP(g_param,ip,np,g_freq) ! see m_freq_mpi.f90, 2017.06.07
  nfreq    = g_freq%nfreq
  nfreq_ip = g_freq%nfreq_ip

!#[4]## allocate respdata and open output files for each observatory
  allocate( resp5(5,nsr_inv,nfreq_ip)  )         ! 2017.08.31
  allocate( tresp(5,nsr_inv,nfreq   )  )         ! 2017.08.31
  CALL ALLOCATERESP(g_param,nsr_inv,resp5,tresp,ip,nfreq,nfreq_ip) ! 2017.07.13

!  CALL PREPRESPFILES(g_param,files,resp5,nfreq) ! 2017.05.18

  CALL initializedatavec(g_param_apinv,h_data) ! m_param_apinv.f90, 2017.08.31
  ndat     = g_data%ndat                              ! 2017.08.31

  allocate(g_apdm(nfreq_ip))
     call allocateapdm(nobs,nfreq_ip,nmodelactive,nsr_inv,g_apdm ) ! 2018.06.25
  allocate(gt_apdm(nfreq)  )
     call allocateapdm(nobs,nfreq,   nmodelactive,nsr_inv,gt_apdm) ! 2018.06.25
  allocate( fp(nline,nsr_inv),fs(nline,nsr_inv) )      ! 2017.08.31

!============================================================ alpha loop start
 ialphaflag = g_param_apinv%ialphaflag ! 2017.09.11
 ialpha     = 1
 nalpha     = 1
 if ( ialphaflag .eq. 1 ) nalpha = g_param_apinv%nalpha ! L-curve 2017.09.08
! write(*,*) "ip",ip,"nalpha=",nalpha

 do ialpha = 1,nalpha ! nalpha > 1 when L-curve 2017.09.08

 h_model     = g_model_ini                        ! initial model, 2017.08.31
 iflag       = 0                                  ! 2017.09.11,  1: end of iteration
 alpha       = 0.d0                               ! 2017.12.13

 if ( ip .eq. 0 ) then

  call watchstart(t_watch)                                    ! 2017.09.11
!  write(*,*) " watch start alpha = ",alpha                   ! 2017.09.11

  head = g_param_apinv%outputfolder                           ! 2017.08.31
  if ( ialphaflag .eq. 1 ) then                      !  L-curve 2017.09.11
   write(num2,'(i1)') ialpha                                  ! 2017.09.11
   head = head(1:len_trim(head))//"a"//num2(1:1)//"/"         ! 2017.08.31
  end if
  open(21,file=head(1:len_trim(head))//"rms.dat") ! 2017.08.31

 end if ! ip = 0 end

!============================================================ iteration loop start
 do ite = 1,itemax

  if ( ip .eq. 0 ) call watchstart(t_watch0)                   ! 2017.03.02
! if ( ip .eq. 0 ) CALL OUTCOND(g_param_apinv,h_cond,g_mesh,ite,ialpha)!commented out 2017.21.20
if (  ip .eq. 0  ) write(*,*) ""
 if ( ip .eq. 0  ) write(*,'(a,i3,a)') " ======= iteration step =",ite," START!! ======" ! 2020.09.29
 if ( ip .eq. 0  ) CALL OUTMODEL(g_param_apinv,h_model,g_mesh,ite,ialpha) ! 2017.12.20
 
 CALL MPI_BARRIER(mpi_comm_world, errno) ! 2020.09.29
 call model2cond(h_model,h_cond)  ! see m_modelpart.f90 2017.12.22
 !call showcond(h_cond,1)

!============================================================= freq loop start
  do i=1,nfreq_ip
   freq_ip = g_freq%freq_ip(i)
   if ( freq_ip .lt. 0. ) cycle
   write(*,'(a,i3,a,i3,a,f9.4,a)') " ip =",ip,"/",np," freq =",freq_ip," [Hz] start!!" ! 2020.09.29
   omega=2.d0*pi*freq_ip

!#[9]## forward calculateion to get d(m)
!  if ( ip .le. 7 ) then
!  write(*,'(a,i5,a,f15.7)') "IP",ip,"freq_ip",freq_ip
!  write(*,'(a,i10,a,i5)') "nline=",nline,"nsr_inv",nsr_inv
!  write(*,'(a,3f15.7)') "sparam%xs1",sparam%xs1
!  write(*,'(a,3f15.7)') "sparam%xs2",sparam%xs2
!  write(*,'(a,f15.7)')  "omega",omega
!  write(*,'(a,2f15.7)')  "h_cond%rho(1,50)=",h_cond%rho(1),h_cond%rho(50)
!  write(*,'(a,i8,a,i8,a,f20.7)') "PT%nrow",PT%nrow,"PT%ncolm",PT%ncolm,"PT%val(10)",PT%val(10)
!  end if
  CALL forward_apinv(A,g_mesh,g_line,nline,nsr_inv,fs,&
      & omega,sparam,g_param,g_param_apinv,h_cond,PT,ut,ip,np)!2020.09.18

!#[11]## output resp to obs file
  fp=0.d0 ! 2018.06.25
  CALL CALOBSEBCOMP(fp,fs,nline,nsr_inv,omega,coeffobs,resp5(:,:,i),g_param_apinv)!2018.10.05 below

!  if ( ip .le. 7 ) then
! write(*,'(a)')"-------------------------------------------------------"
!  write(*,'(a,i2,a,i2,a,f9.4,a)') " ip =",ip,"/",np," freq_ip =",freq_ip," [Hz]" ! 2020.09.18
!  write(*,'(a)') "coeffobs(2,3)%val(1:ntot)"
!  write(*,'(6f15.7)') coeffobs(2,3)%val(1:coeffobs(2,3)%ntot)
!  write(*,'(a)') "fs(coeffobs(2,3)%item(ntot),1)"
!  write(*,'(12f8.4)') (fs(coeffobs(2,3)%item(j),1),j=1,coeffobs(2,3)%ntot)
!  do j=1,5                         ! 2018.10.08
!   if ( g_param_apinv%iflag_comp(j) .eq. 0 ) cycle ! 2018.10.08
!   write(*,*) "comp(j): ",comp(j)   ! 2018.10.08
!   write(*,'(a,i5,a,4f15.7)') "i",i,"resp5(3,1,i)%ftobsamp(1:4)",resp5(j,1,i)%ftobsamp(1:4)
!   write(*,'(a,i5,a,4f15.7)') "i",i,"resp5(3,1,i)%ftobspha(1:4)",resp5(j,1,i)%ftobsphase(1:4)
!  end do                           ! 2018.10.08
!  end if

!#[10]## generate d|amp|dm and d(pha)dm for jacobian
  CALL genjacobian1(nobs,nline,nsr_inv,ut,fs,PT,h_model,g_mesh,g_line,&
                  &  omega,g_apdm(i),g_param_apinv,ip,np) ! 2020.09.18

 end do

do i=1,nfreq_ip ! 20200806
 freq_ip = g_freq%freq_ip(i)
  if ( ip .le. 7 ) then
 ! write(*,'(a)')"-------------------------------------------------------"
!  write(*,'(a,i5,1x,a,f15.7)') "ip",ip,"freq_ip=",freq_ip,"after genjacobian1"
!  write(*,'(a)') "coeffobs(2,3)%val(1:ntot)"
!  write(*,'(6f15.7)') coeffobs(2,3)%val(1:coeffobs(2,3)%ntot)
!  write(*,'(a)') "fs(coeffobs(2,3)%item(ntot),1)"
!  write(*,'(12f8.4)') (fs(coeffobs(2,3)%item(j),1),j=1,coeffobs(2,3)%ntot)
!  do j=1,5                         ! 2018.10.08
!   if ( g_param_apinv%iflag_comp(j) .eq. 0 ) cycle ! 2018.10.08
!   write(*,*) "comp(j): ",comp(j)   ! 2018.10.08
!   write(*,'(a,i5,a,4f15.7)') "i",i,"resp5(3,1,i)%ftobsamp(1:4)",resp5(j,1,i)%ftobsamp(1:4)
!   write(*,'(a,i5,a,4f15.7)') "i",i,"resp5(3,1,i)%ftobspha(1:4)",resp5(j,1,i)%ftobsphase(1:4)
!  end do                           ! 2018.10.08
  end if
end do

!============================================================= freq loop end

!#[12]## send result to ip = 0
  CALL SENDRECVRESULT(resp5,tresp,ip,np,nfreq,nfreq_ip,nsr_inv) ! see below 2017.09.03
  CALL SENDBZAPRESULT_AP(g_apdm,gt_apdm,nobs,nfreq,nfreq_ip,nsr_inv,ip,np,g_param_apinv) !2018.10.05

!   write(*,*) "[in 12] ip",ip,"freq_ip=",freq_ip
!   write(*,'(a,i5,a,4f15.7)') "i",i,"tresp(3,1,i)%ftobsamp(1:4)",tresp(3,1,i)%ftobsamp(1:4)
!   write(*,'(a,i5,a,4f15.7)') "i",i,"tresp(3,1,i)%ftobsphase(1:4)",tresp(3,1,i)%ftobsphase(1:4)

  CALL MPI_BARRIER(mpi_comm_world, errno) ! 2017.09.03

!#[13]## OUTPUT RESULTS for each iteration step
  if ( ip .eq. 0 ) then
   CALL OUTOBSFILESINV(g_param,g_param_apinv,nsr_inv,sparam,tresp,nfreq,&
                     & ite,ialpha)!2017.09.11

!# gen dvec and cal misfit term ! 2017.12.22
   CALL GENDVEC_AP(g_param_apinv,nsr_inv,tresp,nfreq,h_data,g_data) ! 2018.10.08
   CALL CALRMS_AP(g_data,h_data,Cd,misfit,nrms)   ! cal rms, 2017.09.08
   if ( ite .eq. 1 ) nrms_ini = nrms              ! 2018.01.25

!# cal roughness 2018.01.22
   call CALROUGHNESS(g_param_apinv,h_model,g_model_ref,rough1,rough2,BM,R) !m_modelroughenss

!# check RMS and set alpha
    !#[1]## output misfit and nrms
    frms = g_param_apinv%finalrms  ! 2018.06.25
    if ( frms .lt. 0.1  ) frms=1.0 ! 2018.06.25
    if ( nrms .lt. frms ) then     ! 2018.06.25
     call OUTRMS(21,ite,nrms,misfit,alpha,rough1,rough2,1) !     converged  2017.09.08
     iflag = 1 ;  goto 80                           ! iflag = 1 means "end" 2017.12.20
    else
     call OUTRMS(21,ite,nrms,misfit,alpha,rough1,rough2,0) ! Not converged  2017.09.08
    end if
    if ( ite .ge. 2 .and. nrms_ini .lt. nrms ) then ! 2018.01.25 stop when nrms is larger
     iflag = 1 ; goto 80                            ! 2018.01.25
    end if                                          ! 2018.01.25

!    #[Option]## if nrms increased, replace the reference model
     if ( g_param_apinv%iflag_replace .eq. 1 ) then ! 2018.06.26
     if ( ite .eq. 1 ) nrms0 = nrms/0.9             ! initial     2018.06.26
     if ( nrms/nrms0 .gt. 1.0  ) then               !             2018.06.26
      if (g_param_apinv%ialphaflag    .eq. 2 ) then ! cooling strategy 2017.07.19
      g_model_ref = pre_model
      h_model     = pre_model
      goto 80
      end if
     end if; end if

!-------------------------------------------------------- followings are for next model
    ! ialphaflag = 1 : L-curve method
    ! ialphaflag = 2 : cooling with given alpha0
    ! ialphaflag = 3 : cooling with automatic alpha0 (Minami et al. 2018)

    !#[1]## Generate BMI for next model : see m_modelroughenss.f90
    if (      itype_roughness .eq. 2 ) then ! MS:  Minimum support
     call GENBMI_MS( h_model,g_model_ref,g_param_apinv,BM,BMI,ite,ialphaflag)! 2017.12.25
    else if ( itype_roughness .eq. 3 ) then ! MSG: Minimum support gradient
     call GENBMI_MSG(h_model,g_model_ref,g_param_apinv,R,RI,BM,BMI,ite,ialphaflag) ! 2017.12.25
    end if

    !#[2]## generate JJ (jacobian)
    call genjacobian2(g_param_apinv,nfreq,gt_apdm,JJ) !m_jacobian_ap.f90 2017.12.11
    if ( g_param_apinv%ioutlevel .eq. 1 ) then  ! 2018.06.25
     CALL OUTJACOB(g_param_apinv,JJ,ite,h_model,g_param,sparam)!m_jacobian_ap.f90 2018.06.25
    end if                                      ! 2018.06.25

    if ( g_param_apinv%iboundflag .eq. 2) then ! 2018.01.22 J -> J'
     call TRANSJACOB(g_param_apinv,JJ,h_model) ! 2018.01.22 m_jacobian_ap.f90
    end if

    !#[2]# Initial alpha
    if ( ite == 1 ) then   ! 2017.12.13
     if ( ialphaflag .eq. 1 ) alpha = g_param_apinv%alpha(ialpha)! 2017.12.13
     if ( ialphaflag .eq. 2 ) alpha = g_param_apinv%alpha_init   ! 2017.12.13
     if ( ialphaflag .eq. 3 ) then  ! 2017.12.13 Minami et al. (2018) cooling
      kmax = 30

!      call alphaspectralradius_v1(Cd,JJ,BM,alpha,kmax)   ! 2017.12.25 m_spectral.f90
       if (      itype_roughness .eq. 3 ) then            ! MSG 2018.01.18
        call alphaspectralradius_v2(Cd,JJ,RI,alpha,kmax)  ! 2018.01.18 m_spectral.f90
       else if ( itype_roughness .eq. 2 ) then            ! MS  2018.02.04
        call alphaspectralradius_v2(Cd,JJ,BMI,alpha,kmax) ! 2018.01.18 m_spectral.f90
        alpha = alpha * (g_param_apinv%beta**2.)          ! 2018.02.04
       else                                               ! 2018.02.04
        call alphaspectralradius_v2(Cd,JJ,BMI,alpha,kmax) ! 2017.12.24 m_spectral.f90
	 end if
      alpha = g_param_apinv%gamma*alpha                   ! 2017.12.21
     end if                                               ! 2017.12.11
    end if

    !#[3]## New alpha for cooling strategies
    if ( ite .ge. 2 ) then   ! 2017.12.13
     if ( ialphaflag .eq. 2 .or. ialphaflag .eq. 3 ) then !== cooling strategy
      if ( nrms/nrms0 .gt. 0.9  ) alpha = alpha*(10.**(-1./3.d0))
     end if ! if ialphaflag = 2 or 3
    end if ! 2017.12.13
    if ( ialphaflag .eq. 4 ) then ! Modified version of Grayver et al. (2013)
	call alphaspectralradius_v2(Cd,JJ,BMI,alpha,30) ! 2018.01.09 m_spectral.f90
      alpha = g_param_apinv%gamma*alpha/(1.*ite) ! 2018.01.09
    end if ! 2018.09.01

    !#[5]## obtaine new model
    pre_model = h_model ! 2017.06.14 keep previous model
    call getnewmodel(JJ,g_model_ref,h_model,g_data,h_data,BMI,CD,alpha,g_param_apinv)!m_jacobian_ap.f90

    nrms0 = nrms                                 !  2017.12.25

  80 continue                                    !  2017.12.20

  end if ! ip = 0

   CALL MPI_BARRIER(mpi_comm_world, errno)
   CALL MPI_BCAST(iflag,1,MPI_INTEGER4,0,mpi_comm_world,errno)
   if (iflag .eq. 1) goto 100

   CALL SHAREMODEL(h_model,ip) ! see m_shareformpi.f90, 2017.06.05

   if ( ip .eq. 0) call watchstop(t_watch0)
   if ( ip .eq. 0) write(*,'(a,i3,a,f8.4,a)') " ======= iteration step =",ite," END!! ======= Time=",t_watch0%time," [min]"!2020.09.17

end do ! iteration loop============================================  iteration end!

 100 continue
 if ( ip .eq. 0 ) close(21) ! 2017.05.18
 if ( ip .eq. 0 ) then
  call watchstop(t_watch)   ! 2017.09.11
  write(*,*) "### END alpha=",alpha,"Time=",t_watch%time,"[min]" !2017.09.11
 end if

end do ! alpha loop end! 2017.09.08

 CALL MPI_FINALIZE(errno) ! 2017.06.05

end program inversion_ap

!#################################################### PREPAREPT
!# PT [nobs*ncomp,nline] is generated 2018.10.05
subroutine PREPAREPT(coeffobs,g_param_apinv,PT)
use matrix
use param_apinv
implicit none
type(param_apinversion), intent(in)  :: g_param_apinv
type(real_crs_matrix),   intent(in)  :: coeffobs(2,3)
type(reaL_crs_matrix),   intent(out) :: PT(5)        ! 2018.10.04
integer(4),dimension(5)              :: iflag_comp
integer(4) :: i,j,icomp

!#[1]## set
 iflag_comp = g_param_apinv%iflag_comp
 do i=1,5
  PT(i)%ntot = 0
 end do

!#[2]## calculate nrow
  icomp = 0
  do i=2,1,-1 !#  B, E loop
  do j=1,3 !# x,y,z comp loop
   if ( icomp .eq. 5 ) exit    ! 2018.10.05
   icomp = icomp + 1
   if ( iflag_comp(icomp) .eq. 0 ) cycle
   CALL DUPLICATE_CRSMAT(coeffobs(i,j),PT(icomp)) ! PT(1:5) [nobs,nline] , m_matrix.f90
  end do
  end do

return
end
!#################################################### READREFINITCOND
!# Coded on 2018.06.21
subroutine READREFINITCOND(r_cond,h_cond,g_param_apinv,g_mesh) ! 2018.10.04
use param_apinv ! 2018.06.21
use modelpart
use param
use mesh_type   ! 2018.10.04
implicit none
type(mesh),           intent(in)    :: g_mesh          ! 2018.10.04
type(param_cond),     intent(inout) :: r_cond
type(param_cond),     intent(inout) :: h_cond
type(param_apinversion),intent(in)  :: g_param_apinv   ! 2018.06.21
character(50)                       :: modelfile,connectfile
integer(4)                          :: icondflag_ini, icondflag_ref
real(8)                             :: sigmahomo       ! 2018.10.04

!#[0]##
icondflag_ref = g_param_apinv%icondflag_ref         ! 2018.06.21
icondflag_ini = g_param_apinv%icondflag_ini         ! 2018.06.21

!#[ref cond]
if     ( icondflag_ref .eq. 0 ) then                ! 2018.10.04
 sigmahomo = g_param_apinv%sigmahomo_ref            ! 2018.10.04
 CALL SETCOND(r_cond,g_mesh,sigmahomo)              ! 2018.10.04
elseif ( icondflag_ref .eq. 1 ) then ! condfile       2018.03.18
 r_cond%condfile = g_param_apinv%g_refcondfile      ! 2017.07.19
 CALL READCOND(r_cond)                              ! 2017.07.19
elseif ( icondflag_ref .eq. 2 ) then ! modelfile
 modelfile   = g_param_apinv%g_refmodelfile
 connectfile = g_param_apinv%g_refmodelconn
 CALL READMODEL2COND(r_cond,connectfile,modelfile)
else                                                ! 2018.10.04
 write(*,*) "GEGEGE! icondflag_ref",icondflag_ref   ! 2018.10.04
 stop                                               ! 2018.10.04
end if

!#[ini cond]
if     ( icondflag_ini .eq. 0 ) then                ! 2018.10.04
 sigmahomo = g_param_apinv%sigmahomo_ini            ! 2018.10.04
 CALL SETCOND(h_cond,g_mesh,sigmahomo)              ! 2018.10.04
elseif ( icondflag_ini .eq. 1 ) then ! condfile       2018.03.18
 h_cond%condfile = g_param_apinv%g_initcondfile     ! see m_param_inv.f90
 CALL READCOND(h_cond)                              ! see m_param.f90
elseif ( icondflag_ini .eq. 2 ) then ! modelfile
 modelfile   = g_param_apinv%g_initmodelfile
 connectfile = g_param_apinv%g_initmodelconn
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
else if ( icflag .eq. 1 ) then
     write(*,*) "Converged!!with nrms =",nrms ! normalized rms 2017.09.08
     write(idev,'(i10,5g15.7,a)') ite, nrms, misfit, alpha,rough1,rough2,"Converged!"! 2017.09.04
end if
return
end


!######################################### OUTOBSFILESINV 2018.10.05
!# modified on 2018.10.05
!# copied from src_inv_mpi/n_inv_mpi.f90 2017.09.03
subroutine OUTOBSFILESINV(g_param,g_param_apinv,nsr_inv,sparam,tresp,nfreq,ite,ialpha)
use param_apinv ! 2017.09.03
use param
use outresp
implicit none
integer(4),                intent(in)    :: ite,nfreq,ialpha ! 2017.09.11
integer(4),                intent(in)    :: nsr_inv
type(respdata),            intent(in)    :: tresp(5,nsr_inv,nfreq)!2017.07.14
type(param_forward),       intent(in)    :: g_param
type(param_source),        intent(in)    :: sparam          ! 2017.07.14
type(param_apinversion),   intent(in)    :: g_param_apinv   ! 2017.07.14
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
 head       = g_param_apinv%outputfolder     ! 2017.07.25
 nh         = len_trim(head)
 allocate(srcindex(nsr_inv), freq(nfreq) )   ! 2017.07.14
 allocate(data_avail(2,5,nfreq,nobs,nsr_inv))! 2018.10.05
 data_avail = g_param_apinv%data_avail       ! 2017.09.03 see m_param_inv.f90
 srcindex   = g_param_apinv%srcindex         ! 2017.07.14
 freq       = g_param%freq                   ! 2017.07.14
 ialphaflag = g_param_apinv%ialphaflag       ! 2017.09.11
 iflag_comp = g_param_apinv%iflag_comp       ! 2018.10.05
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
  CALL SET_TABLE_DOF(dofn,nline,table_dof,nline)! table_dof is generated

!#[2]## set allocate A
  CALL set_iccg_var7_dofn(dofn,6,nline,nline,nline,ntet,g_line%n6line,table_dof,A,ip)

return
end

!############################################
!# modified on 2017.09.03 for multiple sources
!# coded on 2017.06.05
subroutine SENDBZAPRESULT_AP(g_apdm,gt_apdm,nobs,nfreq,nfreq_ip,nsr,ip,np,g_param_apinv)!2018.10.05
use matrix
use shareformpi_ap ! 2017.06.08
use jacobian_ap    ! 2017.06.08
use param_apinv    ! 2018.10.08
implicit none
integer(4),             intent(in)    :: nobs,nfreq,nfreq_ip,ip,np
integer(4),             intent(in)    :: nsr                   ! 2017.09.03
type(amp_phase_dm),     intent(in)    :: g_apdm(nfreq_ip)      ! 2017.09.03
type(param_apinversion),intent(in)    :: g_param_apinv         ! 2018.10.08
type(amp_phase_dm),     intent(inout) :: gt_apdm(nfreq)        ! 2017.09.03
integer(4)                            :: i,ifreq,ip_from,errno
integer(4)                            :: icomp                 ! 2018.10.05
integer(4)                            :: k                     ! 2017.09.03
integer(4),          dimension(5)     :: iflag_comp            ! 2018.10.08

!#[0]## set
  iflag_comp = g_param_apinv%iflag_comp ! 2018.10.08

!#[1]## share ampbz
 do i=1,nfreq
  if (mod(i,np) .eq. 1 ) ip_from = -1
  ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1

  do k=1,nsr      ! 2017.09.03

  do icomp = 1,5  ! 2018.10.05

  if ( iflag_comp(icomp) .eq. 0)  cycle ! 2018.10.08

  !# dampdm(icomp,k) is a matrix of [nobs,nmodel] for icomp component and k-th source
  if ( ip .eq. ip_from )  gt_apdm(i)%dampdm(icomp,k) = g_apdm(ifreq)%dampdm(icomp,k)!2018.10.05
  if ( ip .eq. ip_from )  gt_apdm(i)%dphadm(icomp,k) = g_apdm(ifreq)%dphadm(icomp,k)!2018.10.05

  !# share     2017.09.03
  call sharerealcrsmatrix(gt_apdm(i)%dampdm(icomp,k),ip_from,ip) ! m_shareformpi_ap.f90
  call sharerealcrsmatrix(gt_apdm(i)%dphadm(icomp,k),ip_from,ip) ! m_shareformpi_ap.f90

  end do     ! 2018.10.05 component loop

  end do     ! 2017.09.03

 end do

!#[2]##
  if (ip .eq. 0) write(*,*) "### SENDBZAPRESULT END!! ###"

return
end

!#############################################
!# modified on 2017.09.03 to include nsr
!# coded on 2017.06.05
subroutine SENDRECVRESULT(resp5,tresp,ip,np,nfreq,nfreq_ip,nsr)!2017.09.03
use outresp
use shareformpi_ap
implicit none
!include 'mpif.h'
integer(4),    intent(in)    :: ip,np,nfreq,nfreq_ip
integer(4),    intent(in)    :: nsr                   ! 2017.09.03
type(respdata),intent(in)    :: resp5(5,nsr,nfreq_ip) ! 2017.09.03
type(respdata),intent(inout) :: tresp(5,nsr,nfreq)    ! 2017.09.03
integer(4)                   :: errno,i,j,k,l,ip_from,ifreq             ! 2020.08.06

!#[1]##
 do i=1,nfreq
  if (mod(i,np) .eq. 1 ) ip_from = -1
  ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1

if ( ip .eq. ip_from ) then
   do k=1,nsr; do j=1,5 ! 20200807
    tresp(j,k,i) = resp5(j,k,ifreq) ! 20200807
   end do ; end do     ! 20200807
  end if


do j=1,5
   do k=1,nsr ! 2017.09.03

!if ( ip .eq. ip_from .and. j .eq. 3 ) then ! 20200806
!  write(*,*) "-- before sharing --SENDRECVRESULT --"
!  write(*,*)"ip_from",ip_from,"tresp i",i,"k",k,"j",j
!  write(*,'(a,4g15.7)') "resp5 amp",(resp5(j,k,ifreq)%ftobsamp(l),l=1,resp5(j,k,ifreq)%nobs) ! 20200807
!  write(*,'(a,4g15.7)') "resp5 phase",(resp5(j,k,ifreq)%ftobsphase(l),l=1,resp5(j,k,ifreq)%nobs) ! 20200807
!  write(*,'(a,4g15.7)') "tresp amp",(tresp(j,k,i)%ftobsamp(l),l=1,tresp(j,k,i)%nobs)
!  write(*,'(a,4g15.7)') "tresp phase",(tresp(j,k,i)%ftobsphase(l),l=1,tresp(j,k,i)%nobs)
!  write(*,*) "--"
! end if

   call sharerespdata(tresp(j,k,i),ip_from) ! see m_shareformpi_ap 2017.09.03

!if ( ip .eq. 0 .and. j .eq. 3 ) then ! 20200806
! write(*,*) "-- after sharing --SENDRECVRESULT --"
! write(*,*)"ip=0 tresp i",i,"k",k,"j",j
! write(*,'(a,4g15.7)') "tresp amp",(tresp(j,k,i)%ftobsamp(l),l=1,tresp(j,k,i)%nobs)
! write(*,'(a,4g15.7)') "tresp phase",(tresp(j,k,i)%ftobsphase(l),l=1,tresp(j,k,i)%nobs)
! write(*,*) "--"
!end if

   end do     ! 2017.09.03
  end do
 end do

if ( ip .eq. 0 ) write(*,*) "### SENDRECVRESULT END!! ###"

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
!# modified on 2018.10.05 for multiple components
!# modified on 2017.09.03 to include multiple sources
!# coded on 2017.05.31
subroutine CALOBSEBCOMP(fp,fs,nline,nsr,omega,coeffobs,resp5,g_param_apinv)!2017.09.03
use matrix
use outresp
use param_apinv ! 2018.10.05
implicit none
type(param_apinversion),intent(in)    :: g_param_apinv   ! 2018.10.05
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
iflag_comp = g_param_apinv%iflag_comp ! 2018.10.08

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
use param_apinv
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
!############################################## subroutine GENDVEC_AP
!# modified for multiple component on 2018.10.05
!# modified for multiple sources   on 2017.09.04
!# modified on 2017.06.08
subroutine GENDVEC_AP(g_param_apinv,nsr_inv,tresp,nfreq,g_data_ap,g_data) !2018.10.08
use param_apinv ! 2016.06.08
use outresp     ! 2016.06.08
implicit none
integer(4),                          intent(in)     :: nfreq
integer(4),                          intent(in)     :: nsr_inv       ! 2017.09.04
type(param_apinversion),             intent(in)     :: g_param_apinv ! 2017.09.04
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
  nobs       = g_param_apinv%nobs              ! 2017.09.04
  allocate(data_avail(2,5,nfreq,nobs,nsr_inv)) ! 2017.07.14
  data_avail = g_data_ap%data_avail            ! 2017.09.03
  allocate( ampbe(5,nobs,nfreq,nsr_inv) )      ! 2018.10.05
  allocate( phabe(5,nobs,nfreq,nsr_inv) )      ! 2018.10.05
  iflag_comp = g_param_apinv%iflag_comp        ! 2018.10.05

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
subroutine GENCD_AP(g_data_ap,CD) ! 2017.08.31
use matrix
use param_apinv ! 2017.09.04
implicit none
type(data_vec_ap),      intent(in)  :: g_data_ap
type(real_crs_matrix),  intent(out) :: CD
real(8),allocatable,dimension(:)    :: error
integer(4) :: ndat,i

!#[1]# set
ndat   = g_data_ap%ndat
write(*,*) "ndat=",ndat
allocate(error(ndat))
error  = g_data_ap%error

!#[2]## gen Cd : diagonal matrix
CD%nrow  = ndat
CD%ncolm = ndat
CD%ntot  = ndat
allocate(CD%stack(0:ndat),CD%item(ndat),CD%val(ndat))
CD%stack(0)=0
do i=1,ndat
 CD%val(i)   = (error(i))**2.d0
 CD%stack(i) = i
 CD%item(i)  = i
end do

write(*,*) "### GENCD_AP END!! ###"
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
type(mesh),intent(in) :: em_mesh
type(param_forward),intent(inout) :: g_param
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

write(*,*) "### GENXYZMINMAX END!! ###"
return
end
!########################################### OUTMODEL
! modified on 2017.12.20
subroutine OUTMODEL(g_param_apinv,g_model,g_mesh,ite,ialpha)! 2017.09.11
use mesh_type
use modelpart
use param_apinv ! 2017.08.31
use caltime     ! 2017.12.25
implicit none
integer(4),              intent(in) :: ite,ialpha       ! 2017.09.11
type(param_apinversion), intent(in) :: g_param_apinv    ! 2017.08.31
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
head         = g_param_apinv%outputfolder ! 2017.09.03
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
ialphaflag   = g_param_apinv%ialphaflag            ! 2017.09.11
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
subroutine OUTCOND(g_param_apinv,g_cond,g_mesh,ite,ialpha) ! 2017.09.11
use mesh_type
use param
use param_apinv ! 2017.09.04
implicit none
integer(4),             intent(in) :: ite,ialpha    ! 2017.09.11
type(param_apinversion),intent(in) :: g_param_apinv ! 2017.09.04
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
head         = g_param_apinv%outputfolder ! 2017.09.04
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
ialphaflag   = g_param_apinv%ialphaflag ! 2017.09.11
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

!###################################################################
! Coded on Feb 15, 2016
subroutine set_table_dof(dofn,nodtot,table_dof,doftot_ip)
implicit none
	integer(4),intent(in) :: dofn, nodtot, doftot_ip
	integer(4),intent(out) :: table_dof(nodtot,dofn)
	integer(4) :: i,j,dof_id
	dof_id=0
	do i=1,nodtot
       do j=1,dofn
	  dof_id = dof_id + 1
	  table_dof(i,j)= dof_id
       end do
	end do
	if (dof_id .ne. doftot_ip ) then
       write(*,*) "GEGEGE dof_id=",dof_id,"doftot_ip=",doftot_ip
	 stop
	end if
!	write(*,*) "### SET TABLE_DOF END!! ###" ! commented out on 2017.09.04
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
    write(*,11) " x,y =",xyzobs(1,j)," , ",xyzobs(1,j)," [km]" ! 2020.09.17
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
