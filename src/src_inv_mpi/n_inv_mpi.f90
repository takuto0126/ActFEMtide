! modified on 2017.07.13 for multiple sources
! Coded on 2017.06.02
!
program inversion_csem
!
use param
use param_inv
use mesh_type
use line_type
use face_type
use matrix
use modelpart
use constants        ! added on 2017.05.14
use iccg_var_takuto  ! added on 2017.05.14
use outresp
use jacobian         ! added on 2017.05.16
use shareformpi      ! added on 2017.06.03
use freq_mpi         ! added on 2017.06.07
use caltime          ! added on 2018.03.16
!
implicit none
type(param_forward)   :: g_param       ! see m_param.f90
type(param_source)    :: sparam        ! see m_param.f90
type(param_cond)      :: g_cond        ! see m_param.f90 : goal conductivity
type(param_cond)      :: h_cond        ! see m_param.f90 : initial structure
type(param_cond)      :: r_cond        ! see m_param.f90 : ref 2017.07.19
type(param_inversion) :: g_param_inv   ! see m_param_inv.f90
type(mesh)            :: g_mesh        ! see m_mesh_type.f90
type(mesh)            :: h_mesh        ! z file; see m_mesh_type.f90
type(line_info)       :: g_line        ! see m_line_type.f90
type(face_info)       :: g_face        ! see m_line_type.f90
type(modelpara)       :: g_modelpara   ! see m_modelpart.f90
type(model)           :: g_model_ref   ! see m_modelpart.f90
type(model)           :: g_model_ini   ! see m_modelpart.f90 2017.07.19
type(model)           :: h_model       ! see m_modelpart.f90 (2017.05.17)
type(model)           :: pre_model     ! see m_modelpart.f90, 2017.07.07
type(data_vec)        :: g_data        ! observed data  ; see m_param_inv.f90
type(data_vec)        :: h_data        ! calculated data; see m_param_inv.f90
type(freq_info)       :: g_freq        ! see m_freq_mpi.f90
type(real_crs_matrix) :: DQTDQ,CM,CD   ! see m_matrix.f90
type(real_crs_matrix) :: PT            ! [nobs,nlin]
type(global_matrix)   :: A             ! see m_iccg_var_takuto.f90
type(real_crs_matrix) :: coeffobs(2,3) ! see m_matrix.f90 ; 1 for edge, 2 for face
integer(4)                                  :: ite,i,j,errno
integer(4)                                  :: nline, ntet, nobs, ndat ! 2017.07.14
integer(4)                                  :: nsr_inv                 ! 2017.07.13
real(8)                                     :: omega, freq_ip
integer(4)                                  :: nfreq,nfreq_ip ! for MPI
real(8)                                     :: rms, rms0, alpha
complex(8),    allocatable,dimension(:,:)   :: fp,fs     ! (nline,nsr_inv)2017.07.13
type(obsfiles)                              :: files     ! see m_outresp.f90
type(respdata),allocatable,dimension(:,:,:) :: resp5     !resp5(5,nsr,nfrq_ip)2017.07.13
type(respdata),allocatable,dimension(:,:,:) :: tresp     !tresp(5,nsr,nfreq)  2017.07.13
integer                                     :: access    ! 2017.05.15
type(complex_crs_matrix)                    :: ut        ! 2017.05.16
type(amp_phase_dm),allocatable,dimension(:) :: g_apdm    ! 2017.07.14
type(amp_phase_dm),allocatable,dimension(:) :: gt_apdm   ! 2017.07.14
type(real_crs_matrix)                       :: JJ        ! Jacobian matrix, 2017.05.17
real(8)                                     :: roughness ! 2017.07.19
integer(4)                                  :: ialpha    ! 2017.07.19
real(8)                                     :: rms_init  ! 2017.07.19
character(50)                               :: head      ! 2017.07.25
type(watch)                                 :: t_watch   ! 2018.03.16
type(watch)                                 :: t_watch0  ! 2018.03.16
real(8)                                     :: frms      ! 2018.06.14
!
integer(4) :: ip,np, itemax = 15, iflag = 0
!#---------------------------------------------------------
!# Algorithm
!# Phi(m) = F(m) + alpha*R(m)
!# F(m)     = 1/2*(d(m)-d_obs)^T Cd^-1 (d(m)-d_obs)
!# dF/dm    = J(m)^T Cd^-1 (d(m)-d_obs)
!# d2F/dm2  = J(m)^T Cd^-1 J(m)
!# R(m)     = 1/2*(m - m_ref)^T(DQ)^T(DQ)(m-m_ref)
!# dR/dm    = (DQ)^T(DQ)(m-m_ref)
!# d2R/dm2  = (DQ)^T(DQ)
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
!# [2] gen Cm from DQTDQ (by BiCG ?)
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
 write(*,*) "ip=",ip,"np=",np

 if ( ip .eq. 0) then !################################################# ip = 0

!#[0]## read parameters
  CALL READPARAM(g_param,sparam,g_cond) ! include READCOND for initial model
  CALL READPARAINV(g_param_inv,g_modelpara,g_param,sparam,g_data) ! 2018.03.20

!#[1]## read mesh
  CALL READMESH_TOTAL(g_mesh,g_param%g_meshfile)
  if ( access( g_param%z_meshfile, " ") .eq. 0 ) then ! if exist, 2017.05.15
    CALL READMESH_TOTAL(h_mesh,g_param%z_meshfile)
    CALL PREPZSRCOBS(h_mesh,g_param,sparam)        ! see below, include kill h_mesh
  end if
  CALL GENXYZMINMAX(g_mesh,g_param)                ! see below
  CALL READLINE(g_param%g_lineinfofile,g_line)     ! see m_line_type.f90
  CALL READFACE(g_param_inv%g_faceinfofile,g_face) ! see m_face_type.f90
  CALL FACE2ELEMENT(g_face)                     ! cal g_face%face2ele(1:2,1:nface)

!#[2]## read initial(same as reference) model
  CALL DUPLICATESCALARCOND(g_cond,r_cond)          ! see m_param.f90,  2017.07.19
  CALL DUPLICATESCALARCOND(g_cond,h_cond)          ! see m_param.f90,  2017.05.14
  CALL deallocatecond(g_cond)                      ! see m_param.f90,  2017.07.19
  CALL READREFINITCOND(r_cond,h_cond,g_param_inv)  ! 2018.03.18
  CALL SETNPHYS1INDEX2COND(g_mesh,r_cond)          ! ref     cond; see below
  CALL SETNPHYS1INDEX2COND(g_mesh,h_cond)          ! initial cond; see below

!#[3]## generate model space
  CALL genmodelspace(g_mesh,g_modelpara,g_model_ref,g_param,h_cond)! m_modelpart.f90
  CALL modelparam(g_model_ref)                     ! 2017.08.28 see m_modelpart.f90
  CALL OUTMODEL(g_param_inv,g_model_ref,g_mesh,0)  ! 2018.03.17

!#[4]## cal Cm
  CALL GENDQ_AP(DQTDQ,g_face,g_mesh,g_model_ref)   ! see below
  CALL GENCM(DQTDQ,CM)                             ! see below

!#[5]##
  CALL assignmodelrho(r_cond,g_model_ref)          ! 2017.07.19
  g_model_ini = g_model_ref                        ! 2017.07.19
  CALL assignmodelrho(h_cond,g_model_ini)          ! 2017.07.19

!#[6]## cal Cd
  CALL GENCD(g_data,CD)                            ! see below

end if !#################################################################   ip = 0

!#[3]## share params, mesh, and line (see m_shareformpi.f90)
  CALL shareinv(g_param,sparam,h_cond,g_mesh,g_line,g_face,g_param_inv,&
		 &  g_model_ini,ip) !2017.07.19 g_model_ref -> h_model

!#[6]## cal Pt: matrix for get data from simulation results
  CALL PREPOBSCOEFF(g_param,g_mesh,g_line,coeffobs) ! for x,y,z component, see below
  CALL DUPLICATE_CRSMAT(coeffobs(2,3),PT) ! PT [nobs, nline] , see m_matrix.f90

!#[7]## Preparation for global stiff matrix
  nline   = g_line%nline
  nobs    = g_param%nobs
  ntet    = g_mesh%ntet
  nsr_inv = g_param_inv%nsr_inv ! 2017.07.13
  CALL SET_ICCG_VAR(ntet,nline,g_line,A,ip) ! see below, 2017.06.05

!#[8]## set frequency
  call SETFREQIP(g_param,ip,np,g_freq) ! see m_freq_mpi.f90, 2017.06.07
  nfreq    = g_freq%nfreq
  nfreq_ip = g_freq%nfreq_ip

!#[4]## allocate respdata and open output files for each observatory
  allocate( resp5(5,nsr_inv,nfreq_ip)  )              ! 2017.07.13
  allocate( tresp(5,nsr_inv,nfreq   )  )              ! 2017.07.13
  CALL ALLOCATERESP(g_param,nsr_inv,resp5,tresp,ip,nfreq,nfreq_ip) ! 2017.07.13

!  CALL PREPRESPFILES(g_param,files,resp5,nfreq)      ! 2017.05.18

 CALL initializedatavec(g_param_inv,h_data) ! m_param_inv.f90, 2017.07.14
  ndat     = g_data%ndat                              ! 2017.07.14

 allocate(g_apdm(nfreq_ip)  )                         ! 2017.07.14
     call allocateapdm(nobs,nfreq_ip,nsr_inv,g_apdm ) ! 2017.07.14
 allocate(gt_apdm(nfreq  )  )                         ! 2017.07.14
     call allocateapdm(nobs,nfreq, nsr_inv,  gt_apdm) ! 2017.07.14

 allocate( fp(nline,nsr_inv),fs(nline,nsr_inv) )      ! 2017.07.13

!#[8]## iteration loop start
 h_model = g_model_ini           ! initial model,       2017.07.19
! if ( ip .eq. 0) CALL OUTMODEL(g_param_inv,h_model,g_mesh,1)    ! 2017.07.25
 head = g_param_inv%outputfolder                                 ! 2017.07.25
 if ( ip .eq. 0) open(21,file=head(1:len_trim(head))//"rms.dat") ! 2017.07.25

!============================================================ iteration loop start
 ialpha = 1
 if ( ip .eq. 0 ) then
  if (g_param_inv%ialphaflag .eq. 1 ) alpha = g_param_inv%alpha(ialpha) ! 2017.07.19
  if (g_param_inv%ialphaflag .eq. 2 ) alpha = g_param_inv%alpha_init    ! 2017.07.19
 end if
! write(*,*) "alpha",alpha

 do ite = 1,itemax

  if ( ip .eq. 0 ) CALL watchstart(t_watch0)                     ! 2017.03.16
  if ( ip .eq. 0 ) CALL OUTMODEL(g_param_inv,h_model,g_mesh,ite) ! 2017.12.20
  call model2cond(h_model,h_cond)  ! see m_modelpart.f90
!  if ( ip .eq. 0 ) CALL OUTCOND(g_param_inv,h_cond,g_mesh,ite)!commented out 2018.03.16
  !call showcond(h_cond,1)

!============================================================= freq loop start
  do i=1,nfreq_ip
   freq_ip = g_freq%freq_ip(i)
   if ( freq_ip .lt. 0. ) cycle
   write(*,*) "ip=",ip,"frequency =",freq_ip,"[Hz]"
   omega=2.d0*pi*freq_ip

!#[9]## forward calculateion to get d(m)
  CALL forward_inv(A,g_mesh,g_line,nline,nsr_inv,fs,&
		    & omega,sparam,g_param,g_param_inv,h_cond,PT,ut)!2017.07.14

!#[10]## generate d|bz|/dm and bzamp for jacobian
  CALL genjacobian1(nobs,nline,nsr_inv,ut,fs,PT,h_model,g_mesh,g_line,&
                  & omega,g_apdm(i))            !2017.07.14

!#[11]## output resp to obs file
  CALL CALOBSEBCOMP(fp,fs,nline,nsr_inv,omega,coeffobs,resp5(:,:,i)) ! 2017.07.13 see below

 end do
!============================================================= freq loop end

!#[12]## send result to ip = 0
  CALL SENDRECVRESULT(resp5,tresp,ip,np,nfreq,nfreq_ip,nsr_inv)  ! 2017.07.13
  CALL SENDBZAPRESULT(g_apdm,gt_apdm,nobs,nfreq,nfreq_ip,nsr_inv,ip,np)!2017.07.14

  CALL MPI_BARRIER(mpi_comm_world, errno)

!#[13]## OUTPUT RESULTS for each iteration step
  if ( ip .eq. 0 ) then
   CALL OUTOBSFILESINV(g_param,g_param_inv,nsr_inv,sparam,tresp,nfreq,ite)!2017.07.14

!# gen dvec and cal rms
   CALL GENDVEC(gt_apdm,nfreq,h_data) ! gen data vector
   CALL CALRMS(g_data,h_data,Cd,rms)  ! cal rms, 2017.05.17

 !# conversion check
   call CALROUGHNESS(DQTDQ,h_model,g_model_ref,roughness)
   frms = g_param_inv%finalrms    ! 2018.06.14
   if ( frms .lt. 0.1  ) frms=1.0 ! 2018.06.14
   if ( rms .lt. frms  ) then     ! 2018.06.14
    write(*,*) "Converged!!with rms =",rms
    write(21,'(i10,3g15.7,a)') ite, rms, alpha, roughness,"Converged!"! 2017.06.14
    iflag = 1
    goto 70  ! end
   end if
   write(21,'(i10,3g15.7)') ite, rms, alpha, roughness               ! 2017.07.19

 !# h_model is revised for the next iteration
   if ( ite .eq. 1 ) rms0     = rms/0.9 ! initial alpha will be used
!   if ( ite .eq. 1 ) rms_init = rms/0.9
   write(*,*) "rms =",rms
   write(*,*) "rms0=",rms0
   if ( rms/rms0 .gt. 0.9 .and. g_param_inv%ialphaflag == 2) alpha = alpha*10**(-0.3d0)

   !# replace model_ref if rms > rms0 ! 2017.06.14
   if ( rms/rms0 .gt. 1.0 ) then ! if there are no movement 2017.07.19
     if (g_param_inv%ialphaflag .eq. 2 ) then ! cooling strategy 2017.07.19
!	alpha = alpha*10**(0.3d0)  ! keep alpha still 2017.07.20
      g_model_ref = pre_model
      h_model     = pre_model
      goto 80
!     else if (g_param_inv%ialphaflag .eq. 1 ) then ! L-curve strategy 2017.07.19
!      ialpha = ialpha + 1
!	if ( ialpha .gt. g_param_inv%nalpha) goto 70 ! end
!      alpha    = g_param_inv%alpha(ialpha)
!      h_model  = g_model_ini    ! 2017.07.19
!	rms0     = rms_init       ! 2017.07.19
!	goto 90
!     else
!      write(*,*) "GEGEGE"
!	iflag = 1
!	goto 70 ! end
     end if
   end if

 !# h_model is revised for the next iteration
   call genjacobian2(g_param_inv,nfreq,gt_apdm,JJ,g_model_ref) ! 2017.07.14
   pre_model = h_model ! 2017.06.14 keep previous model
   if ( g_param_inv%ioutlevel .eq. 1 ) then  ! 2018.03.29
    CALL OUTJACOB(g_param_inv,JJ,ite,h_model,g_param,sparam)!m_jacobian.f90 2018.03.29
   end if                                    ! 2018.03.29
   call getnewmodel(JJ,g_model_ref,h_model,g_data,h_data,CM,CD,alpha)   ! m_jacobian.f90, 2017.05.17

!90 continue  ! 2017.07.19

 !# output new model
80 continue
!   CALL OUTMODEL(g_param_inv,h_model,g_mesh,ite+1)   ! 2017.09.19

70 continue
  end if ! ip = 0

  CALL MPI_BARRIER(mpi_comm_world, errno)
  CALL MPI_BCAST(iflag,1,MPI_INTEGER4,0,mpi_comm_world,errno)
  if (iflag .eq. 1) goto 100

  CALL SHAREMODEL(h_model,ip) ! see m_shareformpi.f90, 2017.06.05

   if ( ip .eq. 0 ) rms0 = rms
   if ( ip .eq. 0)  call watchstop(t_watch0)
   if ( ip .eq. 0)  write(*,*) "### END ite=",ite,"Time=",t_watch0%time,"[min]"!2018.03.16

end do ! iteration loop

 100 continue
 if ( ip .eq. 0 ) close(21) ! 2017.05.18

 CALL MPI_FINALIZE(errno) ! 2017.06.05

end program inversion_csem

!#################################################### READREFINITCOND
!# Coded on 2018.03.18
subroutine READREFINITCOND(r_cond,h_cond,g_param_inv) ! 2018.03.18
use param_inv
use modelpart
use param
implicit none
type(param_cond),     intent(inout) :: r_cond
type(param_cond),     intent(inout) :: h_cond
type(param_inversion),intent(in)    :: g_param_inv ! 2018.03.19
character(50)                       :: modelfile,connectfile
integer(4)                          :: icondflag_ini, icondflag_ref

!#[0]##
icondflag_ref = g_param_inv%icondflag_ref
icondflag_ini = g_param_inv%icondflag_ini

!#[ref cond]
if ( icondflag_ref .eq. 1 ) then ! condfile        2018.03.18
 r_cond%condfile = g_param_inv%g_refcondfile      ! 2017.07.19
 CALL READCOND(r_cond)                            ! 2017.07.19
elseif ( icondflag_ref .eq. 2 ) then ! modelfile
 modelfile   = g_param_inv%g_refmodelfile
 connectfile = g_param_inv%g_refmodelconn
 CALL READMODEL2COND(r_cond,connectfile,modelfile)
end if

!#[ini cond]
if ( icondflag_ini .eq. 1 ) then ! condfile        2018.03.18
 h_cond%condfile = g_param_inv%g_initcondfile     ! see m_param_inv.f90
 CALL READCOND(h_cond)                            ! see m_param.f90
elseif ( icondflag_ini .eq. 2 ) then ! modelfile
 modelfile   = g_param_inv%g_initmodelfile
 connectfile = g_param_inv%g_initmodelconn
 CALL READMODEL2COND(h_cond,connectfile,modelfile)
end if

return
end
!##################################################### READMODEL2COND
!# coded on 2018.03.18
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
 r_cond%sigma=-9999
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
!######################################### CALROUGHNESS 2017.07.19
!# nmodelactive is introduced on 2018.03.20
subroutine CALROUGHNESS(DQTDQ,h_model,g_model_ref,roughness) ! 2017.07.19
use matrix
use modelpart
implicit none
type(real_crs_matrix),intent(in)  :: DQTDQ
type(model),          intent(in)  :: h_model
type(model),          intent(in)  :: g_model_ref
real(8),              intent(out) :: roughness
integer(4)                        :: nmodel,i,ii, nmodelactive ! 2018.03.20
real(8),allocatable,dimension(:)  :: dmodel,dmodel2
integer(4),allocatable,dimension(:) :: iactive      ! 2018.03.20
!#[1]## set
nmodel       = h_model%nmodel
nmodelactive = h_model%nmodelactive ! 2018.03.20
allocate(iactive(nmodel)) ! 2018.03.20
iactive       = h_model%iactive
allocate(dmodel(nmodelactive),dmodel2(nmodelactive))

ii=0
do i=1,nmodel
if (iactive(i) .eq. 1 ) then
 ii = ii + 1
 dmodel(ii) = h_model%logrho_model(i) - g_model_ref%logrho_model(i)
end if
end do

!#[2]## cal roughness
call mul_matcrs_rv(DQTDQ,dmodel,nmodelactive,dmodel2) ! 2018.03.20
roughness = 0.d0
do i=1,nmodelactive ! 2018.03.20
 roughness = roughness + dmodel(i)*dmodel2(i)
end do

write(*,*) "### CALROUGHNESS END!! ###"

return
end

!######################################### OUTOBSFILESINV 2017.05.18
!# output folder is changed on 2017.07.25
!# modified on 2017.07.14 for multiple source
!# coded on 2017.05.18
subroutine OUTOBSFILESINV(g_param,g_param_inv,nsr_inv,sparam,tresp,nfreq,ite)!2017.07.14
use param_inv
use param
use outresp
implicit none
integer(4),              intent(in) :: ite,nfreq
integer(4),              intent(in) :: nsr_inv
type(respdata),          intent(in) :: tresp(5,nsr_inv,nfreq)!2017.07.14
type(param_forward),     intent(in) :: g_param
type(param_source),      intent(in) :: sparam         ! 2017.07.14
type(param_inversion),   intent(in) :: g_param_inv    ! 2017.07.14
integer(4),allocatable,dimension(:) :: srcindex       ! 2017.07.14
real(8),   allocatable,dimension(:) :: freq           ! 2017.07.14
character(2)   :: num                       ! 2017.07.14
integer(4)     :: nhead                     ! 2017.07.14
integer(4)     :: i,j,k,l,nobs,nsour,nsite  ! 2017.07.14
character(50)  :: filename,head,site, sour  ! 2017.07.14
logical,allocatable,dimension(:,:,:) :: data_avail ! 2017.07.14

!#[0]## set
 nobs     = g_param%nobs
 head     = g_param_inv%outputfolder       ! 2017.07.25
 nhead    = len_trim(head)
 allocate(srcindex(nsr_inv), freq(nfreq) ) ! 2017.07.14
 srcindex = g_param_inv%srcindex           ! 2017.07.14
 write(num,'(i2.2)') ite                   ! 2017.07.14
 freq     = g_param%freq                   ! 2017.07.14
 allocate(data_avail(nfreq,nobs,nsr_inv))  ! 2017.07.14
 data_avail = g_param_inv%data_avail ! 2017.07.14 see m_param_inv.f90

 write(*,*) "nsr_inv",nsr_inv
 write(*,*) "nfreq",nfreq
 write(*,*) "nobs",nobs
 write(*,*) "data_avail",data_avail(:,:,:)

!#[2]##
 do l=1,nobs
  site  = g_param%obsname(l)
  nsite = len_trim(site)

  do k=1,nsr_inv
   sour     = sparam%sourcename(srcindex(k))
   nsour    = len_trim(sour)
   filename = head(1:nhead)//site(1:nsite)//"_"//sour(1:nsour)//"_"//num(1:2)//".dat"
   open(31,file=filename)

   do i=1,nfreq
    if ( data_avail(i,l,k)) then ! 2017.07.14 if data is used for inversion
     write(31,'(11g15.7)') freq(i),(tresp(j,k,i)%ftobsamp(l),&
     &                           tresp(j,k,i)%ftobsphase(l),j=1,5)
    else
     write(31,'(g15.7,a)') freq(i),"0. 0. 0. 0. 0. 0. 0. 0. 0. 0."!2017.07.14
    end if
   end do

   close(31)
  end do
 end do

 write(*,*) "### OUTOBSFILESINV END!! ###"

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
  CALL SET_TABLE_DOF(dofn,nline,table_dof,nline)

!#[2]## set allocate A
  CALL set_iccg_var7_dofn(dofn,6,nline,nline,nline,ntet,g_line%n6line,table_dof,A,ip)

return
end

!############################################
!# modified on 2017.07.14 for multiple sources
!# coded on 2017.06.05
subroutine SENDBZAPRESULT(g_apdm,gt_apdm,nobs,nfreq,nfreq_ip,nsr,ip,np)
use matrix
use shareformpi
use jacobian   ! 2017.06.07
implicit none
integer(4),           intent(in)    :: nobs,nfreq,nfreq_ip,ip,np
integer(4),           intent(in)    :: nsr ! 2017.07.14
type(amp_phase_dm),   intent(in)    :: g_apdm(nfreq_ip)      ! 2017.06.07
type(amp_phase_dm),   intent(inout) :: gt_apdm(nfreq)        ! 2017.06.07
integer(4) :: i,ifreq,ip_from,errno
integer(4) :: k ! 2017.07.14

!#[1]## share ampbz
 do i=1,nfreq
  if (mod(i,np) .eq. 1 ) ip_from = -1
  ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1
!  write(*,*) "i=",i,"ip=",ip,"ampbz",g_apdm(ifreq)%ampbz

  if ( ip .eq. ip_from)gt_apdm(i)%ampbz(:,:)  = g_apdm(ifreq)%ampbz(:,:) !2017.07.14
  call MPI_BCAST(gt_apdm(i)%ampbz,nobs*nsr,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)! 2017.07.14

  do k=1,nsr ! 2017.07.17
   if ( ip .eq. ip_from) gt_apdm(i)%dampbzdm(k)=g_apdm(ifreq)%dampbzdm(k)!2017.07.17
   call sharerealcrsmatrix(gt_apdm(i)%dampbzdm(k),ip_from,ip) !2017.07.14
  end do      ! 2017.07.14

 end do

!#[check]## 
 if (.false.) then
 if (ip .eq. 0) then
 do k=1,nsr
 do i=1,nfreq
  write(*,*) "i",i,"ntot=",gt_apdm(i)%dampbzdm(k)%ntot
 end do
 end do
 write(*,*) "ampbz="
 do i=1,nfreq
  do k=1,nobs
  write(*,*) "ifreq iobs",i,k,"ampbz=",gt_apdm(i)%ampbz(k,1)
  end do
 end do
 end if
 end if

!#[2]##
  if (ip .eq. 0) write(*,*) "### SENDBZAPRESULT END!! ###"

return
end

!#############################################
!# modified on 2017.07.13 to include nsr
!# coded on 2017.06.05
subroutine SENDRECVRESULT(resp5,tresp,ip,np,nfreq,nfreq_ip,nsr)
use outresp
use shareformpi
implicit none
!include 'mpif.h'
integer(4),    intent(in)    :: ip,np,nfreq,nfreq_ip
integer(4),    intent(in)    :: nsr                   ! 2017.07.13
type(respdata),intent(in)    :: resp5(5,nsr,nfreq_ip) ! 2017.07.13
type(respdata),intent(inout) :: tresp(5,nsr,nfreq)    ! 2017.07.13
integer(4) :: errno,i,j,k,ip_from,ifreq               ! 2017.07.13

!write(*,*) "nfreq,nfreq_ip,nsr",nfreq,nfreq_ip,nsr,"ip",ip

!#[1]##
 do i=1,nfreq
  if (mod(i,np) .eq. 1 ) ip_from = -1
  ip_from = ip_from + 1 ; ifreq = (i-1)/np + 1
  if ( ip .eq. ip_from ) tresp(1:5,1:nsr,i) = resp5(1:5,1:nsr,ifreq)!2017.07.13
  do j=1,5
   do k=1,nsr ! 2017.07.13
   call sharerespdata(tresp(j,k,i),ip_from) ! see m_shareformpi 2017.07.13
   end do     ! 2017.07.13
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
!# modified on 2017.07.14 to include multiple sources
!# coded on 2017.05.31
subroutine CALOBSEBCOMP(fp,fs,nline,nsr,omega,coeffobs,resp5)!2017.07.13
use matrix
use outresp
implicit none
real(8),              intent(in)    :: omega
integer(4),           intent(in)    :: nline
integer(4),           intent(in)    :: nsr             ! 2017.07.14
complex(8),           intent(inout) :: fp(nline,nsr),fs(nline,nsr)! 2017.07.14
type(real_crs_matrix),intent(in)    :: coeffobs(2,3)
type(respdata),       intent(inout) :: resp5(5,nsr)    ! 2017.07.14
integer(4)                          :: i               ! 2017.07.14

!CALL CALOBSRESP(fp,fs,nline,coeffobs(2,1),resp5(1)  ) !bx
!CALL CALOBSRESP(fp,fs,nline,coeffobs(2,2),resp5(2)  ) !by
do i= 1, nsr ! 2017.07.14
 CALL CALOBSRESP(fp(:,i),fs(:,i),nline,coeffobs(2,3),resp5(3,i)  ) !bz 2017.07.14
end do       ! 20917.07.14
!fp= - (0.d0,1.d0)*omega*fp ! E= -i*omega*A
!fs= - (0.d0,1.d0)*omega*fs !
!CALL CALOBSRESP(fp,fs,nline,coeffobs(1,1),resp5(4)  ) ! ex
!CALL CALOBSRESP(fp,fs,nline,coeffobs(1,2),resp5(5)  ) ! ey

return
end

!############################################## subroutine CALRMS
!# Coded 2017.05.18
subroutine CALRMS(g_data,h_data,Cd,rms)
use param_inv
use matrix
type(real_crs_matrix),intent(in)  :: Cd
type(data_vec),       intent(in)  :: g_data,h_data
real(8),              intent(out) :: rms
real(8),allocatable,dimension(:)  :: dvec1,dvec2,dvec
integer(4)                        :: ndat,ndat1,ndat2,i
real(8)                           :: d2

!#[0]## set
ndat  = g_data%ndat
ndat1 = h_data%ndat
ndat2 = Cd%ntot
if (ndat .ne. ndat1 .or. ndat .ne. ndat2 ) goto 99
allocate(dvec1(ndat),dvec2(ndat),dvec(ndat))
dvec1 = g_data%dvec
dvec2 = h_data%dvec
dvec = dvec2 - dvec1

!#[1]## cal rms
rms = 0.d0
do i=1,ndat
 d2 = dvec(i)**2.d0 ! 2017.07.25
 write(*,10) "i",i,"rms",d2," d2/err2",d2/Cd%val(i) ! 2017.07.25
 rms = rms + d2/ Cd%val(i) ! dvec(i)**2/err(i)**2
end do
 rms = sqrt(rms/dble(ndat))

10 format(a,i10,a,g15.7,a,g15.7)

write(*,*) "rms=",rms
write(*,*) "### CALRMS END!! ###"

return
99 continue
write(*,*) "GEGEGE! ndat",ndat,".ne. ndat1",ndat1,"or .ne. ndat2",ndat2
stop

end
!############################################## subroutine GENDVEC
!# modified on 2017.07.14 for multisource inversion
subroutine GENDVEC(g_apdm,nfreq,g_data) ! 2017.07.14
use param_inv
use jacobian ! 2017.06.07
implicit none
integer(4),           intent(in)     :: nfreq
type(amp_phase_dm),   intent(in)     :: g_apdm(nfreq)
type(data_vec),       intent(inout)  :: g_data      ! 2017.07.14
real(8),allocatable,dimension(:,:,:) :: ampbz       ! 2017.07.14
logical,allocatable,dimension(:,:,:) :: data_avail  ! 2017.07.14
integer(4)                           :: iobs,ii,ifreq,i,isr ! 2017.07.14
integer(4)                           :: nobs,nsr_inv ! 2017.07.14
integer(4)                           :: icount     ! 2017.07.14

!#[0]## set
  nobs       = g_data%nobs    ! 2017.07.14
  nsr_inv    = g_data%nsr_inv ! 2017.07.14
  allocate(data_avail(nfreq,nobs,nsr_inv)) ! 2017.07.14
  allocate(     ampbz(nfreq,nobs,nsr_inv)) ! 2017.07.14
  data_avail = g_data%data_avail           ! 2017.07.14
  do i=1,nfreq
   ampbz(i,1:nobs,1:nsr_inv) = g_apdm(i)%ampbz(1:nobs,1:nsr_inv)
  end do

!#[1]## cal g_data
ii = 0
do ifreq=2,nfreq
 do isr=1,nsr_inv ! 2017.07.14
  do iobs=1,nobs
   if ( data_avail(ifreq,iobs,isr)) then ! 2017.07.14
   ii = ii + 1
   g_data%dvec(ii)  = ampbz(ifreq,iobs,isr)/ampbz(1,iobs,isr) ! 2017.07.17
   end if
 end do
end do
end do ! 2017.07.14

if (ii .ne. g_data%ndat) then !2017.07.14 for check
 write(*,*) "GEGEGE icount",ii,"is not equal to ndat",g_data%ndat
end if

!#[2]## out
do ii=1,g_data%ndat
 write(*,*) ii,"g_data%dvec=",g_data%dvec(ii)
end do

write(*,*) "### GENDVEC END!! ###"
return
end

!############################################## subroutine GENCD
! Coded on May 13, 2017
subroutine GENCD(g_data,CD)
use matrix
use param_inv
implicit none
type(data_vec),       intent(in)  :: g_data
type(real_crs_matrix),intent(out) :: CD
real(8),allocatable,dimension(:)  :: error
integer(4) :: ndat,i

!#[1]# set
ndat   = g_data%ndat
write(*,*) "ndat=",ndat
allocate(error(ndat))
error  = g_data%error

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

write(*,*) "### GENCD END!! ###"
return
end

!############################################## subroutine GENCM
subroutine GENCM(DQTDQ,CM)
use matrix
implicit none
type(real_crs_matrix),intent(inout) :: DQTDQ
type(real_crs_matrix),intent(out)   :: CM
type(real_crs_matrix)               :: crsin,crsout
real(8),allocatable,dimension(:,:)  :: x
integer(4) :: nmodel,i,j,j1,j2,ntot
real(8)   :: threshold = 1.d-10, epsilon
integer(8),allocatable,dimension(:) :: stack,item

nmodel = DQTDQ%nrow
ntot   = DQTDQ%ntot
allocate(stack(0:nmodel),item(ntot))
stack  = DQTDQ%stack
item   = DQTDQ%item

!#[1]## add small values to diagonals (Usui et al., 2017)
epsilon =1.d-3 ! small value
do i=1,nmodel
 do j=stack(i-1)+1,stack(i)
 if ( item(j) .ge. i ) then
  DQTDQ%val(j) = DQTDQ%val(j) + epsilon
end if
 end do
end do
!write(*,*) "DQTDQ=" ! 2017.07.18
!call realcrsout(CM) ! 2017.07.18

!#[2]## CM^-1 -> CM

write(*,*) "nmodel check =",nmodel
allocate(x(nmodel,nmodel))
call solveCM(nmodel,DQTDQ,x)

!#[3]##
call conv_full2crs(x,nmodel,nmodel,CM,threshold)


!# check
if ( .false. ) then
call mulreal_crs_crs_crs(DQTDQ,CM,crsout)
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
do i=1,crsout%nrow
! if (crsout%stack(i)-crsout%stack(i-1) .ne. 0) then
 write(*,*) i,"# of content",crsout%stack(i)-crsout%stack(i-1)!,"ele2model=",ele2model(i)
 j1=crsout%stack(i-1)+1;j2=crsout%stack(i)
 write(*,'(5g15.7)') (crsout%item(j),j=j1,j2)
 write(*,'(5g15.7)') (crsout%val(j),j=j1,j2)
 write(*,*) ""
! end if
end do
end if

 write(*,*) "### GENCM END!! ###"
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

  nphys2 = g_cond%nphys2 ! # of elements in 2nd physical volume (land)
  nphys1 = g_mesh%ntet - nphys2
  ntet   = g_mesh%ntet
  g_cond%nphys1 = nphys1
  g_cond%ntet   = ntet
  write(*,*) "nphys1=",nphys1,"nphys2=",nphys2,"ntet=",g_mesh%ntet
  do i=1,nphys2
   g_cond%index(i) = nphys1 +i
  end do

return
end
!############################################## subroutine genDQ_AP
! Copied from GENDQ in ../src_inv_ap/n_inv_ap.f90
!#  DQ : [nphys2,nmodel]
!#  D :  [nphys2,nphys2]
!#  Q :  [nphys2,nmodel]
subroutine GENDQ_AP(DQTDQ,g_face,g_mesh,g_model) ! 2018.03.19
use face_type
use matrix
use mesh_type
use modelpart
use fem_util ! 2017.05.18 for calculations of face areas
implicit none
type(face_info),      intent(inout)   :: g_face
type(mesh),           intent(in)      :: g_mesh
type(model),          intent(in)      :: g_model     ! 2018.03.16
type(real_crs_matrix),intent(out)     :: DQTDQ
type(real_crs_matrix)                 :: D, Q,DQT,DQ,crsout,QT,QTDQ,crs,DQ2
type(real_ccs_matrix)                 :: DQTCCS,QTCCS
integer(4),allocatable,dimension(:,:) :: n4face,n4flag,n4
integer(4),allocatable,dimension(:)   :: index, ele2model
integer(4),allocatable,dimension(:)   :: iactive   ! 2018.03.18
integer(4),allocatable,dimension(:)   :: icount
real(8),allocatable,dimension(:,:)    :: band
integer(4),allocatable,dimension(:,:) :: band_ind
real(8),   allocatable,dimension(:,:) :: xyz
integer(4) :: ncolm,nrow,iface,iele,nphys1,nphys2,ntot,nc,n5(5)
integer(4) :: nface,ntet,i,j,nmodel,icele,j1,j2,node,k
real(8)    :: r5(5)
real(8)    :: elm_xyz(3,4),a4(4)    ! 2017.05.18
integer(4) :: icombine,nmodelactive ! 2018.03.16
integer(4) :: ii,jj                 ! 2018.03.18

!#[0]## set
  nface        = g_face%nface
  ntet         = g_face%ntet
  node         = g_mesh%node
  nphys1       = g_model%nphys1
  nphys2       = g_model%nphys2  ! # of elements in land
  allocate(n4face(ntet,4),n4flag(ntet,2),n4(ntet,4),xyz(3,node))
  allocate(index(nphys2),ele2model(nphys2))
  n4face       = g_face%n4face
  index        = g_model%index ! element id for nphys2
  nmodel       = g_model%nmodel
  ele2model    = g_model%ele2model
  n4           = g_mesh%n4            ! 2017.05.18
  xyz          = g_mesh%xyz           ! 2017.05.18
  n4flag       = g_mesh%n4flag        ! 2017.06.18
  icombine     = g_model%icombine     ! 2018.03.19
  write(*,*) "GENDQ_AP nmodel",nmodel
  allocate( iactive(nmodel) )      ! 2018.03.18
  iactive      = g_model%iactive      ! 0 not active, 1 active 2018.03.18
  nmodelactive = g_model%nmodelactive ! 2018.03.18

!#[2]## gen band matrix for D = Nele*Nele

  !#[2-1]## count only the land element connection
  ncolm = 5
  nrow  = nphys2
  allocate(icount(nrow),band(ncolm,nrow),band_ind(ncolm,nrow))
  ! count (iflag = 1) and assign (iflag = 2) loop
  do i=1,nphys2 ! element loop
   icele=index(i) ! element id for whole element group
   icount(i) = 1
   band_ind(1,i) = i ! element id within nphys2
   do j=1,4               ! 2017.05.18
    elm_xyz(:,j) = xyz(:,n4(icele,j))
   end do
   call area4(elm_xyz,a4) ! 2017.05.18
   do j=1,4     ! face loop
    iface=n4face(icele,j)
    if (iface .gt. 0 ) iele=g_face%face2ele(2, iface) ! face is outward
    if (iface .lt. 0 ) iele=g_face%face2ele(1,-iface) ! face is inward
!    if ( icele .eq. 671195 .or. icele .eq. 671196) then
!     write(*,*) "elm_xyz(1:3,1:4)="
!     write(*,'(3g15.7)') (elm_xyz(1:3,k),k=1,4)
!     write(*,*) "icele",icele,"iface",iface,&
!&    "face2ele(1:2,iface)",g_face%face2ele(1:2,abs(iface)),"nphys1",nphys1
!    end if
    if ( iele .eq. 0 ) goto 100 ! there are no neighbor element in land
    if ( iele .gt. nphys1 ) then ! neighboring element is in land
      icount(i) = icount(i) + 1
      band(1,i) = band(1,i) + 1. ! add 1 to the center element
      band(icount(i),i) = -1.    ! add -1 for the neighboring element
!      band(1,i) = band(1,i) + a4(j)      ! face area; 2017.05.18
!      band(icount(i),i)    = -a4(j)      ! face area; 2017.05.18
	band_ind(icount(i),i) = iele - nphys1 ! element id within nphys2
    end if
    100 continue
   end do
  end do

  !#[2-2]## sort and store as crs matrix, D
   D%nrow  = nphys2
   D%ncolm = nphys2
   allocate(D%stack(0:nrow))
   D%stack(0)=0
   do i=1,nrow
    D%stack(i) = D%stack(i-1) + icount(i)
   end do
   ntot = D%stack(nrow)
   D%ntot = ntot
   D%nrow = nrow
   allocate(D%item(ntot),D%val(ntot))
   do i=1,nrow
    n5(1:5) = (/ 1,2,3,4,5 /)
    nc=icount(i) ; r5(1:nc)= band_ind(1:nc,i)*1.d0
    CALL sort_index(nc,n5(1:nc),r5(1:nc))
    do j=1,nc
     D%item(D%stack(i-1)+j) = band_ind(n5(j),i)
     D%val( D%stack(i-1)+j) = band(n5(j),i)
    end do
   end do

!#[3]## generate crsmatrix Q
   Q%nrow  = nphys2
   Q%ncolm = nmodel
   Q%ntot  = nphys2
   allocate(Q%stack(0:nphys2),Q%item(nphys2))
   allocate(Q%val(nphys2))
   Q%stack(0)=0
   do i=1,nphys2
    Q%stack(i)= i
    Q%item(i) = ele2model(i)
    Q%val(i)  = 1.d0
   end do
   write(*,*) "## Q is generated ##"

!#[4]## calculate DQ [nphys2,nmodel]

  call mulreal_crs_crs_crs(D,Q,DQ)     ! see m_matrix.f90
  call trans_crs2ccs(Q,QTCCS)          ! 2017.06.16
  call conv_ccs2crs(QTCCS,QT)          ! 2017.06.16
  call mulreal_crs_crs_crs(QT,DQ,QTDQ) ! 2017.06.16 for connectivity of models
  DQ = QTDQ ! 2017.06.16
  write(*,*) "## DQ is generated! ##"

  !# 2017.05.19
  !# make weights of the interface jump between models the same
  !# regardress of # of interface triangles, like Usui et al. (2017)
  if (.true.) then ! 2017.05.19
   do i=1,DQ%nrow
    j1 = DQ%stack(i) - DQ%stack(i-1)
    do j2 = DQ%stack(i-1)+1,DQ%stack(i)
     DQ%val(j2) = -1.d0
     if ( DQ%item(j2) .eq. i ) DQ%val(j2) = 1.d0*(j1-1)
     if ( j1 .eq. 1 ) then
      write(*,*) "GEGEGE j1=1!! i=",i
	stop
     end if
    end do
   end do
  end if

!#[5]## When icombine = 2 ; modify DQ and g_model 2018.03.16
  !# update g_model and DQ
  if ( icombine .eq. 2 ) then ! 2018.03.16
   DQ2%nrow  = nmodelactive
   DQ2%ncolm = nmodelactive
   allocate(DQ2%stack(0:DQ2%nrow))
   DQ2%stack(0)=0
   ii=0
   do i=1,nmodel
    if ( iactive(i) .eq. 1 ) then ! only when active
     ii = ii + 1
     DQ2%stack(ii) = DQ2%stack(ii-1)
     do j=DQ%stack(i-1)+1,DQ%stack(i)
      if ( iactive(DQ%item(j)) .eq. 1 ) DQ2%stack(ii) = DQ2%stack(ii) + 1
     end do
    end if
   end do
   if ( ii .ne. nmodelactive ) goto 99
   DQ2%ntot  = DQ2%stack(DQ2%nrow)
   allocate(DQ2%item(DQ2%ntot) )
   allocate(DQ2%val( DQ2%ntot) )
   ii = 0 ; jj = 0
   do i=1,nmodel
    if ( iactive(i) .eq. 1 ) then
     ii = ii + 1
     do j=DQ%stack(i-1)+1,DQ%stack(i)
      if ( iactive(DQ%item(j)) .eq. 1 ) then
	 jj = jj + 1
       DQ2%item(jj) = DQ%item(j)
       DQ2%val( jj) = DQ%val(j)
      end if
     end do
    end if
   end do
   DQ = DQ2 !# DQ is updated from [nmodel * nmodel] -> [nmodel-1 * nmodel-1]
  end if   ! 2018.03.16

!#[6]## calculate DQTDQ
  call trans_crs2ccs(DQ,DQTCCS)    ! DQ [ntet,nmodel]
  write(*,*) "## DQTCCS is generated! ##"
  call conv_ccs2crs(DQTCCS,DQT)
  write(*,*) "## DQT is generated! ##"
  call mulreal_crs_crs_crs(DQT,DQ,DQTDQ)
  write(*,*) "## DQTDQ is generated! ##"

!# output for check
if (.false.) then
crsout = DQ
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
do i=1,crsout%nrow
! if (crsout%stack(i)-crsout%stack(i-1) .ne. 0) then
 write(*,*) i,"# of content",crsout%stack(i)-crsout%stack(i-1)!,"ele2model=",ele2model(i)
 j1=crsout%stack(i-1)+1;j2=crsout%stack(i)
 write(*,'(5g15.7)') (crsout%item(j),j=j1,j2)
 write(*,'(5g15.7)') (crsout%val(j),j=j1,j2)
 write(*,*) ""
! end if
end do
end if

if (.false.) then
crsout = D
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
do i=1,crsout%nrow
! if (crsout%stack(i)-crsout%stack(i-1) .ne. 0) then
 write(*,*) i,"# of content",crsout%stack(i)-crsout%stack(i-1)!,"ele2model=",ele2model(i)
 j1=crsout%stack(i-1)+1;j2=crsout%stack(i)
 write(*,'(5g15.7)') (crsout%item(j),j=j1,j2)
 write(*,'(5g15.7)') (crsout%val(j),j=j1,j2)
 write(*,*) ""
! end if
end do
end if


 write(*,*) "### GENDQ END! ###"
return
99 continue
write(*,*) "GEGEGE! nmodelactive",nmodelactive,"ii",ii
stop
end

!############################################## subroutine genDQ
!#  DQ : [nphys2,nmodel]
!#  D :  [nphys2,nphys2]
!#  Q :  [nphys2,nmodel]
subroutine GENDQ(DQTDQ,g_face,g_mesh,g_model)
use face_type
use matrix
use mesh_type
use modelpart
use fem_util ! 2017.05.18 for calculations of face areas
implicit none
type(face_info), intent(inout) :: g_face
type(mesh),      intent(in)    :: g_mesh
type(model),     intent(in)    :: g_model
type(real_crs_matrix),intent(out) :: DQTDQ
type(real_crs_matrix)  :: D, Q,DQT,DQ,crsout
type(real_ccs_matrix)  :: DQTCCS
integer(4),allocatable,dimension(:,:) :: n4face,n4flag,n4
integer(4),allocatable,dimension(:)   :: index, ele2model
integer(4),allocatable,dimension(:)   :: icount
real(8),allocatable,dimension(:,:)    :: band
integer(4),allocatable,dimension(:,:) :: band_ind
real(8),   allocatable,dimension(:,:) :: xyz
integer(4) :: ncolm,nrow,iface,iele,nphys1,nphys2,ntot,nc,n5(5)
integer(4) :: nface,ntet,i,j,nmodel,icele,j1,j2,node
real(8)    :: r5(5)
real(8)    :: elm_xyz(3,4),a4(4) ! 2017.05.18

!#[0]## set
  nface     = g_face%nface
  ntet      = g_face%ntet
  node      = g_mesh%node
  nphys1    = g_model%nphys1
  nphys2    = g_model%nphys2  ! # of elements in land
  allocate(n4face(ntet,4),n4flag(ntet,2),n4(ntet,4),xyz(3,node))
  allocate(index(nphys2),ele2model(nphys2))
  n4face    = g_face%n4face
  index     = g_model%index ! element id for nphys2
  nmodel    = g_model%nmodel
  ele2model = g_model%ele2model
  n4        = g_mesh%n4     ! 2017.05.18
  xyz       = g_mesh%xyz    ! 2017.05.18

!#[2]## gen band matrix for D = Nele*Nele

  !#[2-1]## count only the land element connection
  ncolm = 5
  nrow  = nphys2
  allocate(icount(nrow),band(ncolm,nrow),band_ind(ncolm,nrow))
  ! count (iflag = 1) and assign (iflag = 2) loop
  do i=1,nphys2 ! element loop
   icele=index(i) ! element id for whole element group
   icount(i) = 1
   band_ind(1,i) = i ! element id within nphys2
   do j=1,4               ! 2017.05.18
    elm_xyz(:,j) = xyz(:,n4(icele,j))
   end do
   call area4(elm_xyz,a4) ! 2017.05.18
   do j=1,4     ! face loop
    iface=n4face(icele,j)
    if (iface .gt. 0 ) iele=g_face%face2ele(1, iface) ! face is outward
    if (iface .lt. 0 ) iele=g_face%face2ele(2,-iface) ! face is inward
    if ( iele .eq. 0 ) goto 100 ! there are no neighbor element in land
    if ( iele .gt. nphys1 ) then ! neighboring element is in land
      icount(i) = icount(i) + 1
!      band(1,i) = band(1,i) + 1 ! add 1 to the center element
!      band(icount(i),i) = -1    ! add -1 for the neighboring element
      band(1,i) = band(1,i) + a4(j)      ! face area; 2017.05.18
      band(icount(i),i)    = -a4(j)      ! face area; 2017.05.18
	band_ind(icount(i),i) = iele - nphys1 ! element id within nphys2
    end if
    100 continue
   end do
  end do

  !#[2-2]## sort and store as crs matrix, D
   n5(1:5) = (/ 1,2,3,4,5 /)
   D%nrow  = nphys2
   D%ncolm = nphys2
   allocate(D%stack(0:nrow))
   D%stack(0)=0
   do i=1,nrow
    D%stack(i) = D%stack(i-1) + icount(i)
   end do
   ntot = D%stack(nrow)
   D%ntot = ntot
   D%nrow = nrow
   allocate(D%item(ntot),D%val(ntot))
   do i=1,nrow
    nc=icount(i) ; r5(1:nc)= band_ind(1:nc,i)*1.d0
    CALL sort_index(nc,n5(1:nc),r5(1:nc))
    do j=1,nc
     D%item(D%stack(i-1)+j) = band_ind(n5(j),i)
     D%val( D%stack(i-1)+j) = band(n5(j),i)
    end do
   end do

!#[3]## generate crsmatrix Q
   Q%nrow  = nphys2
   Q%ncolm = nmodel
   Q%ntot=nphys2
   allocate(Q%stack(0:nphys2),Q%item(nphys2))
   allocate(Q%val(nphys2))
   Q%stack(0)=0
   do i=1,nphys2
    Q%stack(i)= i
    Q%item(i) = ele2model(i)
    Q%val(i)  = 1
   end do

!#[4]## calculate DQ [nphys2,nmodel]
  call mulreal_crs_crs_crs(D,Q,DQ) ! see m_matrix.f90
  write(*,*) "## DQ is generated! ##"

  !# 2017.05.19
  !# make weights of the interface jump between models the same throught whole the region
  if (.true.) then ! 2017.05.19
   do i=1,DQ%nrow
    j1 = DQ%stack(i) - DQ%stack(i-1)
    do j2 = DQ%stack(i-1)+1,DQ%stack(i)
     DQ%val(j2) = -1.d0
     if (DQ%item(j2) .eq. i) DQ%val(j2) = 1.d0*(j1-1)
    end do
   end do
  end if

!#[5]## calculate DQTDQ
  call trans_crs2ccs(DQ,DQTCCS)    ! DQ [ntet,nmodel]
  write(*,*) "## DQTCCS is generated! ##"
  call conv_ccs2crs(DQTCCS,DQT)
  write(*,*) "## DQT is generated! ##"
  call mulreal_crs_crs_crs(DQT,DQ,DQTDQ)
  write(*,*) "## DQTDQ is generated! ##"

!# output for check
if (.false.) then
crsout = DQTDQ
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
do i=1,crsout%nrow
! if (crsout%stack(i)-crsout%stack(i-1) .ne. 0) then
 write(*,*) i,"# of content",crsout%stack(i)-crsout%stack(i-1)!,"ele2model=",ele2model(i)
 j1=crsout%stack(i-1)+1;j2=crsout%stack(i)
 write(*,'(5g15.7)') (crsout%item(j),j=j1,j2)
 write(*,'(5g15.7)') (crsout%val(j),j=j1,j2)
 write(*,*) ""
! end if
end do
end if

 write(*,*) "### GENDQ END! ###"
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
!# modified on 2018.03.16 to generate model_connect.dat
! OUTPUT folder is changed on 2017.07.25
subroutine OUTMODEL(g_param_inv,g_model,g_mesh,ite) ! 2017.07.25
use mesh_type
use modelpart
use param_inv
use caltime   ! 2018.03.16 ../common/m_caltime.f90
implicit none
integer(4),            intent(in)   :: ite
type(param_inversion), intent(in)   :: g_param_inv ! 2017.07.25
type(model),           intent(in)   :: g_model
type(mesh),            intent(in)   :: g_mesh
character(50)                       :: modelfile, head, connectfile
integer(4)                          :: i,j, nphys2,nmodel
integer(4)                          :: npoi,nlin,ntri,ishift
integer(4)                          :: nhead, icombine ! 2018.03.17
real(8),   allocatable,dimension(:) :: logrho_model
integer(4),allocatable,dimension(:) :: ele2model,index
character(2)                        :: num
type(watch)                         :: t_watch    ! 2018.03.16

call watchstart(t_watch)   ! 2018.03.16

!#[0]## set
head         = g_param_inv%outputfolder ! 2017.07.25
nphys2       = g_model%nphys2
nmodel       = g_model%nmodel
allocate(index(nphys2),ele2model(nphys2))
! allocate(rho_model(nmodel))     ! 2017.08.28
allocate(logrho_model(nmodel))
index        = g_model%index
!rho_model   = g_model%rho_model  ! 2017.08.28
logrho_model = g_model%logrho_model
ele2model    = g_model%ele2model
npoi         = g_mesh%npoi
nlin         = g_mesh%nlin
ntri         = g_mesh%ntri
nhead        = len_trim(head)           ! 2018.03.16

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
 open(1,file=modelfile)              ! 2017.12.20
 write(1,'(i10)') nmodel
 do j=1,nmodel
  write(1,'(e15.7)') 10**logrho_model(j) !2017.12.20
 end do
 close(1)

call watchstop(t_watch) ! 2018.03.16
write(*,'(a,g15.7,a)') "### OUTMODEL END!! Time=",t_watch%time,"[min] ###"
return
end
!########################################### OUTCOND
! Output folder is changed on 2017.07.25
! Coded on 2017.05.18
subroutine OUTCOND(g_param_inv,g_cond,g_mesh,ite) ! 2017.07.25
use mesh_type
use param
use param_inv
implicit none
integer(4),           intent(in) :: ite
type(param_inversion),intent(in) :: g_param_inv ! 2017.07.25
type(param_cond),     intent(in) :: g_cond
type(mesh),           intent(in) :: g_mesh
character(50) :: condfile,head ! 2017.07.25
integer(4) :: j, nphys2,nmodel,npoi,nlin,ntri,ishift,nphys1
real(8),allocatable,dimension(:)    :: rho,sigma
integer(4),allocatable,dimension(:) :: index
character(2) :: num

!#[0]## set
head         = g_param_inv%outputfolder ! 2017.07.25
nphys1       = g_cond%nphys1
nphys2       = g_cond%nphys2
allocate(index(nphys2),rho(nphys2),sigma(nphys2))
index        = g_cond%index
rho          = g_cond%rho
sigma        = g_cond%sigma
npoi         = g_mesh%npoi
nlin         = g_mesh%nlin
ntri         = g_mesh%ntri
write(*,*) "nphys1=",nphys1
write(*,*) "nphys2=",nphys2
write(num,'(i2.2)') ite

!#[1]## output rho
 condfile = head(1:len_trim(head))//"cond"//num(1:2)//".msh"
 open(1,file=condfile)

 !# standard info
! CALL MESHOUT(1,g_mesh)

 write(1,'(a)') "$MeshFormat"   ! 2017.09.19
 write(1,'(a)') "2.2 0 8"       ! 2017.09.19
 write(1,'(a)') "$EndMeshFormat"! 2017.09.19
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

!write(*,*) "### PREPOBSCOEFF END!! ###"
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
!	write(*,*) "### SET TABLE_DOF END!! ###"
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
integer(4),           intent(in) :: doftot
complex(8),           intent(in) :: fp(doftot)
complex(8),           intent(in) :: fs(doftot)
type(real_crs_matrix),intent(in) :: coeffobs ! see m_matrix.f90
type(respdata),       intent(inout):: resp   ! see m_outresp.f90
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

write(*,*) "### CALOBSRESP END!! ###"
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
    write(*,*) "xyobs(1:2,j)=",xyzobs(1:2,j)
    write(*,*) j,"/nobs",xyzobs(3,j),"->",znew(j),"[km]"
end do

!#[3-2]## source z
 do k=1,nsr                                                               ! 2017.07.14
!#  start point
    call findtriwithgrid(h_mesh,glist,xs1(1:2,k),iele,a3)                 ! 2017.07.14
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs1(3,k) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs1(3,k) ! 2017.07.14
    !
    write(*,*) "xs1(1:2,k)=",xs1(1:2,k)                                   ! 2017.07.14
    write(*,*) "z",s_param%xs1(3,k),"->",xs1(3,k),"[km]"                  ! 2017.07.14

!#  end point
    call findtriwithgrid(h_mesh,glist,xs2(1:2,k),iele,a3)                 ! 2017.07.14
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs2(3,k) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs2(3,k) ! 2017.07.14
    !
    write(*,*) "xs2(1:2,k)=",xs2(1:2,k)                                   ! 2017.07.14
    write(*,*) "z",s_param%xs2(3,k),"->",xs2(3,k),"[km]"                  ! 2017.07.14
end do

!#[4]## set znew to xyz_r
    g_param%xyzobs(3,1:nobs) = znew(1:nobs)
    s_param%xs1(3,:) = xs1(3,:)                 ! 2017.07.14
    s_param%xs2(3,:) = xs2(3,:)                 ! 2017.07.14

!#[5]## kill mesh for memory 2017.05.15
    call killmesh(h_mesh) ! see m_mesh_type.f90

write(*,*) "### PREPZSRCOBS END ###"
return
end
