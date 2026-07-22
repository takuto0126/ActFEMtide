!# 2017.12.12
module modelroughness
implicit none
!# itype_roughness = 1 (Smoothest model)
! subroutine GENDQ
! subroutine GENCM(RTR,CM,g_mesh,g_model)
!
!# itype_roughness = 2 (Minimum Support)
! subroutine GENCM0_MS


contains
!############################################## calroughness
!# 2018.01.22
subroutine CALROUGHNESS(g_param_joint,h_model,g_model_ref,rough1,rough2,BM,R)
use param_jointinv
use modelpart
use jacobian_joint ! 2022.01.05
use caltime
implicit none
type(param_joint),intent(in)  :: g_param_joint
type(model),            intent(in)  :: h_model
type(model),            intent(in)  :: g_model_ref
type(real_crs_matrix),  intent(in)  :: R        ! used only when MSG
type(real_crs_matrix),  intent(in)  :: BM       ! used only when MS
real(8),                intent(out) :: rough1   ! usual roughness
real(8),                intent(out) :: rough2   ! roughness using transformed model
type(model)                         :: ki,ki_ref
integer(4)                          :: itype_roughness,iboundflag
type(watch)                         :: t_watch

!#[1]## set
 rough1=0.d0  ! 2018.06.26
 rough2=0.d0  ! 2018.06.26
 call watchstart(t_watch)
 itype_roughness = g_param_joint%itype_roughness
 iboundflag      = g_param_joint%iboundflag      !0:off,1:cut,2:trnsformed model

!#[2]## cal roughness
   if (     itype_roughness .eq. 1 ) then ! SM: Smoothest model
    CALL CALROUGHNESS_SM(BM,h_model,g_model_ref,rough1 ) ! 2017.12.25 m_modelroughness.
   else if (itype_roughness .eq. 2 ) then ! MS: Minimum support
    CALL CALROUGHNESS_MS(g_param_joint,h_model,g_model_ref,rough1 ) ! 2017.12.13
   else if (itype_roughness .eq. 3 ) then ! MSG: Minimum support gradient
    CALL CALROUGHNESS_MSG(R,g_param_joint,h_model,g_model_ref,rough1 )! 2017.12.13
   end if

!#[3]## cal roughness of
  if (iboundflag .eq. 2 ) then
   ki     = h_model
   ki_ref = g_model_ref
   call TRANSMODELX(ki,    g_param_joint,1) ! model -> ki
   call TRANSMODELX(ki_ref,g_param_joint,1) ! model -> ki
   !#
   if (     itype_roughness .eq. 1 ) then ! SM: Smoothest model
    CALL CALROUGHNESS_SM(BM,ki,ki_ref,rough2 )              ! 2018.01.22 m_modelroughness.
   else if (itype_roughness .eq. 2 ) then ! MS: Minimum support
    CALL CALROUGHNESS_MS(g_param_joint,ki,ki_ref,rough2 )   ! 2018.01.22
   else if (itype_roughness .eq. 3 ) then ! MSG: Minimum support gradient
    CALL CALROUGHNESS_MSG(R,g_param_joint,ki,ki_ref,rough2 )! 2018.01.22
   end if
   !#
  end if
  call watchstop(t_watch)
!  write(*,'(a,f8.4,a)') " ### CALROUGHNESS END!! ###  Time=",t_watch%time," [min]" ! 2020.09.18
  write(*,'(a)') " ### CALROUGHNESS END!! ###" ! Time=",t_watch%time," [min]" ! 2020.09.18
return
end
!############################################## subroutine GENCM0_MS
!# iboundflag is added on 2018.02.04
!# Coded on 2017.12.12
subroutine GENBMI_MS(g_model,g_model_ref,g_param_joint,BM,BMI,it,ialphaflag)
use modelpart    ! m_modelpart.f90
use matrix       ! m_matrix.f90
use param_jointinv  ! m_param_jointinv.f90
use jacobian_joint  ! m_jacobian_joint.f90 for TRANSMODELX 2018.02.04
use caltime      ! 2017.12.22
implicit none
type(param_joint),intent(in)    :: g_param_joint
integer(4),             intent(in)    :: it
integer(4),             intent(in)    :: ialphaflag ! 2017.12.25
type(model),            intent(in)    :: g_model
type(model),            intent(in)    :: g_model_ref
type(real_crs_matrix),  intent(inout) :: BM,BMI
real(8)                               :: beta       ! beta for Minimum support
integer(4)                            :: i,nmodel
real(8),  allocatable,dimension(:)    :: logrhomodel1,logrhomodel2
integer(4)                            :: iboundflag ! 2018.02.04
type(model)                           :: ki,ki_ref  ! 2018.02.04
integer(4)                            :: BMIdenomflag ! 2018.02.05
type(watch)                           :: t_watch

 call watchstart(t_watch) ! 2017.12.22
 !#[1]## set
 beta         = g_param_joint%beta
 nmodel       = g_model%nmodel
 allocate(logrhomodel1(nmodel),logrhomodel2(nmodel))
 iboundflag   = g_param_joint%iboundflag !0:off,1:cut,2:trnsformed model, 2018.02.04
 BMIdenomflag = g_param_joint%BMIdenomflag ! 1:model, 2:chi are used in BMI dnominator


 !#[2]## generate CM = [BM]^-1
 if ( it .eq. 0 ) then
  BMI%nrow  = nmodel
  BMI%ncolm = nmodel
  BMI%ntot  = nmodel
  allocate(BMI%stack(0:nmodel),BMI%item(nmodel),BMI%val(nmodel))
  BMI%stack(0)=0
 end if

 !#[3]## when ibound = 2, transfer model to ki, 2018.02.04
  if (iboundflag .eq. 2 .and. BMIdenomflag .eq. 2 ) then  ! 2018.02.04
   ki     = g_model                         ! 2018.02.04
   ki_ref = g_model_ref                     ! 2018.02.04
   call TRANSMODELX(ki,    g_param_joint,1) ! model -> ki m_jacobian_joint.f90
   call TRANSMODELX(ki_ref,g_param_joint,1) ! model -> ki m_jacobian_joint.f90
   logrhomodel1 = ki%logrho_model           ! 2018.02.04
   logrhomodel2 = ki_ref%logrho_model       ! 2018.02.04
  else
   logrhomodel1 = g_model%logrho_model      ! 2018.02.04
   logrhomodel2 = g_model_ref%logrho_model  ! 2018.02.04
  end if

 !#[4]## calculate BM and BMI
 do i=1,nmodel
  BMI%item(i)  = i
  BMI%stack(i) = i
  BMI%val(i)   = (((logrhomodel1(i)-logrhomodel2(i))**2. + beta**2.)**2.)/(beta**2.)
 end do

 if (ialphaflag .eq. 3 ) then ! Minami method 2017.12.25
 if ( it .eq. 1 ) then ! 2017.12.25 for initial alpha by spectral method
  BM=BMI
  do i=1,nmodel
   BM%val(i) = 1.d0/BMI%val(i)
  end do
 end if
 if ( it .eq. 2 ) then           ! BM is required only for the first iteration 2017.12.25
  call deallocate_real_crsmat(BM) ! 2017.12.25
 end if ! 2017.12.25
 end if ! 2017.12.25

 call watchstop(t_watch) ! 2017.12.22
 write(*,'(a,g15.7,a)') "### GENCM0_MS END!! Time=",t_watch%time,"[min] ###" ! 2017.12.22
return
end
!######################################### GENBMI_MSG 2017.12.13
!# nmodelactive is introduced 2018.06.25
!# iboudnflag is added on 2018.02.04
!# modified on 2017.12.25
!# coded on 2017.12.13
subroutine GENBMI_MSG(g_model,g_model_ref,g_param_joint,R,RI,BM,BMI,it,ialphaflag)
use modelpart    ! m_modelpart.f90
use matrix       ! m_matrix.f90
use param_jointinv  ! m_param_jointinv.f90
use caltime      ! m_caltime.f90     2017.12.18
use jacobian_joint  ! m_jacobian_joint.f90 2018.02.04
implicit none
type(param_joint),intent(in)    :: g_param_joint
integer(4),             intent(in)    :: it
integer(4),             intent(in)    :: ialphaflag ! 2017.12.25
type(model),            intent(in)    :: g_model
type(model),            intent(in)    :: g_model_ref
type(real_crs_matrix),  intent(in)    :: R,RI       ! 2017.12.25
type(real_crs_matrix),  intent(inout) :: BM,BMI
real(8)                               :: beta   ! beta for Minimum support
integer(4)                            :: i,nmodel
real(8),  allocatable,dimension(:)    :: logrhomodel1,logrhomodel2
type(model)                           :: ki, ki_ref   ! 2018.02.04
integer(4)                            :: iboundflag   ! 2018.02.04
integer(4)                            :: BMIdenomflag ! 2018.02.05
real(8)                               :: m2
type(real_crs_matrix)                 :: W1I,W2I,RIW1I
type(real_crs_matrix)                 :: W1,W2,RW2
integer(4)                            :: nmodelactive  ! 2018.06.25
integer(4)                            :: icombine      ! 2018.06.25
integer(4),allocatable,dimension(:)   :: iactive       ! 2018.06.25
integer(4)                            :: kk            ! 2018.06.25
type(watch) :: t_watch

!#[1]## set
 call watchstart(t_watch) ! 2017.12.18
 beta         = g_param_joint%beta
 nmodel       = g_model%nmodel
 iboundflag   = g_param_joint%iboundflag   ! 0:off,1:cut,2:trnsformed model, 2018.02.04
 BMIdenomflag = g_param_joint%BMIdenomflag ! 1:model, 2:chi are used in BMI 2018.02.05
 nmodelactive = g_model%nmodelactive ! 2018.06.25
 allocate( iactive(nmodel))          ! 2018.06.25
 iactive      = g_model%iactive      ! 2018.06.25
 icombine     = g_model%icombine

 allocate(logrhomodel1(nmodel),logrhomodel2(nmodel))

!#[2]## when ibound = 2, transfer model to ki, 2018.02.04
  if (iboundflag .eq. 2 .and. BMIdenomflag .eq. 2 ) then ! 2018.02.04
   ki     = g_model                         ! 2018.02.04
   ki_ref = g_model_ref                     ! 2018.02.04
   call TRANSMODELX(ki,    g_param_joint,1) ! model -> ki m_jacobian_joint.f90
   call TRANSMODELX(ki_ref,g_param_joint,1) ! model -> ki m_jacobian_joint.f90
   logrhomodel1 = ki%logrho_model           ! 2018.02.04
   logrhomodel2 = ki_ref%logrho_model       ! 2018.02.04
  else
   logrhomodel1 = g_model%logrho_model      ! 2018.02.04
   logrhomodel2 = g_model_ref%logrho_model  ! 2018.02.04
  end if


!#[3]## calculate W1I and W2I
 allocate(W1I%stack(0:nmodel),W1I%item(nmodelactive),W1I%val(nmodelactive)) ! 2018.06.25
 W1I%nrow=nmodelactive ; W1I%ncolm=nmodelactive ; W1I%ntot=nmodelactive     ! 2017.06.25
 W1I%stack(0)=0  ! 2017.12.25
 W2I=W1I
 kk=0            ! 2018.06.25
 do i=1,nmodel
  if ( iactive(i) .eq. 0 ) cycle      ! 2018.06.25
  kk=kk+1                             ! 2018.06.25
  W1I%stack(kk)=kk ; W2I%stack(kk)=kk ! 2018.06.25
  W1I%item(kk) =kk ; W2I%item(kk) =kk ! 2018.06.25
  m2 = (logrhomodel2(i) - logrhomodel1(i) )**2.d0
  W1I%val(kk) = ((m2 + beta**2.)**(3./2.))/(beta**2.) ! 2018.06.25
  W2I%val(kk) = sqrt(m2 + beta**2.)                   ! 2018.06.25
 end do

 !#[4]## reflect W1I and W2I to BMI
 call mulreal_crs_crs_crs(RI,W1I,RIW1I)
 call mulreal_crs_crs_crs(W2I,RIW1I,BMI)

 !#[5]## BM for initial alpha by spectral method 2017.12.25
 if ( ialphaflag .eq. 3 ) then ! for Minami 2018 method
  if ( it .eq. 1 ) then
   W1=W1I       ! 2017.12.25
   W2=W2I       ! 2017.12.25
   do i=1,nmodel
    W1%val(i)=1.d0/W1I%val(i)
    W2%val(i)=1.d0/W2I%val(i)
   end do
   call mulreal_crs_crs_crs(R,W2,RW2)  ! 2017.12.25
   call mulreal_crs_crs_crs(W1,RW2,BM) ! 2017.12.25
  end if
  if ( it .eq. 2 ) then            ! BM is required only for first iteration 2017.12.25
   call deallocate_real_crsmat(BM) ! 2017.12.25
  end if                           ! 2017.12.25
 end if

 call watchstop(t_watch)
 write(*,'(a,g15.7,a)') "### GENBMI_MSG END!! Time=",t_watch%time,"[min] ###"
return
end
!######################################### CALROUGHNESS_MSG 2017.12.13
!# nmodelactive is introduced on 2018.06.25
!# confirmed on 2017.12.18
!# coded on 2017.12.13
!# roughness for minimum support
subroutine CALROUGHNESS_MSG(R,g_param_joint,g_model,g_model_ref,rough)
use param_jointinv
use modelpart
use matrix
use outerinnerproduct
use caltime           ! 2017.12.19
implicit none
real(8),                intent(out)  :: rough
type(real_crs_matrix),  intent(in)   :: R
type(model),            intent(in)   :: g_model
type(model),            intent(in)   :: g_model_ref
type(param_joint),intent(in)   :: g_param_joint
real(8),    allocatable,dimension(:) :: logrhomodel1
real(8),    allocatable,dimension(:) :: logrhomodel2
real(8),    allocatable,dimension(:) :: q,qq
integer(4)   :: nmodel,i
integer(4)   :: nmodelactive                   ! 2018.06.25
integer(4)   :: kk                             ! 2018.06.25
integer(4),allocatable,dimension(:) :: iactive ! 2018.06.25
real(8)      :: beta
real(8)      :: m
type(watch)  :: t_watch ! 2017.12.19

!#[1] set
 call watchstart(t_watch) ! 2017.12.22
 nmodel       = g_model%nmodel
 nmodelactive = g_model%nmodelactive ! 2018.06.25
 allocate(iactive(nmodel))           ! 2018.06.25
 iactive      = g_model%iactive      ! 2018.06.25
 allocate( logrhomodel1(nmodel),logrhomodel2(nmodel))
 allocate( q(nmodelactive), qq(nmodelactive) ) ! 2018.06.25
 logrhomodel1 = g_model%logrho_model
 logrhomodel2 = g_model_ref%logrho_model
 beta         = g_param_joint%beta

 kk=0 ! 2018.06.25
 do i=1,nmodel
  if ( iactive(i) .eq. 0) cycle ! 2018.06.25
  kk=kk+1                       ! 2018.06.25
  m = logrhomodel1(i) - logrhomodel2(i)
  q(kk)   = m/(sqrt(m**2. + beta**2.)) ! BM is the inverse of CM 2018.06.25
 end do

!#[2]## cal rough1 (true roughness)
 call mul_matcrs_rv(R,q,nmodelactive,qq) ! 2018.06.25
 rough = inner_n(q,qq,nmodelactive)      ! 2018.06.25
 rough = rough*0.5d0

 call watchstop(t_watch)
 write(*,'(a,g15.7,a)') "### CALROUGHNESS_MSG END!! Time=",t_watch%time,"[min] ###" ! 2017.12.22

return
end

!######################################### CALROUGHNESS_MS 2017.12.13
!# nmodelactive is introduced 2018.06.25
!# coded on 2017.12.13
!# roughness for minimum support
subroutine CALROUGHNESS_MS(g_param_joint,g_model,g_model_ref,rough1)
use param_jointinv
use modelpart
use caltime  ! 2017.12.22
implicit none
real(8),                intent(out)  :: rough1 ! true
type(model),            intent(in)   :: g_model
type(model),            intent(in)   :: g_model_ref
type(param_joint),intent(in)   :: g_param_joint
real(8),    allocatable,dimension(:) :: logrhomodel1
real(8),    allocatable,dimension(:) :: logrhomodel2
type(real_crs_matrix)                :: BM
integer(4)                           :: nmodel,i
integer(4),allocatable,dimension(:)  :: iactive      ! 2018.06.25
real(8)      :: beta
real(8)      :: m2
type(watch) :: t_watch  ! see m_caltime.f90 2017.12.22

 call watchstart(t_watch) ! see m_caltime.f90 2017.12.22

!#[1] set
 nmodel       = g_model%nmodel
 allocate( iactive(nmodel))          ! 2018.06.25
 iactive      = g_model%iactive      ! 2018.06.25
 allocate( logrhomodel1(nmodel),logrhomodel2(nmodel))
 logrhomodel1 = g_model%logrho_model
 logrhomodel2 = g_model_ref%logrho_model
 beta         = g_param_joint%beta

!#[2]## cal rough1 (true roughness)
 rough1=0.d0
 do i=1,nmodel
  if ( iactive(i) .eq. 0 ) cycle ! 2018.06.25
  m2=(logrhomodel1(i)-logrhomodel2(i))**2.d0
  rough1 = rough1 + m2/(m2+beta**2.d0)
 end do
  rough1 = 0.5d0*rough1

 call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12
 write(*,'(a,g15.7,a)') "### CALROUGHNESS_MS END!! Time=",t_watch%time,"[min] ###"

return
end

!######################################### CALROUGHNESS_SM 2017.07.19
!# nmodelactive is introduced on 2018.06.25
!# moved from n_inv_ap.f90 on 2017.12.12
subroutine CALROUGHNESS_SM(BM,h_model,g_model_ref,roughness) ! 2017.12.12
use matrix
use modelpart
use caltime   ! 2017.12.22
implicit none
type(real_crs_matrix),intent(in)   :: BM
type(model),          intent(in)   :: h_model
type(model),          intent(in)   :: g_model_ref
real(8),              intent(out)  :: roughness
integer(4)                         :: nmodel,i
real(8),allocatable,  dimension(:) :: dmodel,dmodel2
integer(4)                         :: nmodelactive,kk ! 2018.06.25
integer(4),allocatable,dimension(:):: iactive         ! 2018.06.25
type(watch) :: t_watch             ! see m_caltime.f90 2017.12.22

 call watchstart(t_watch) ! see m_caltime.f90 2017.12.22

!#[1]## set
nmodel       = h_model%nmodel
nmodelactive = h_model%nmodelactive                   ! 2018.06.25
write(*,'(a,2i5)') " [CALROUGHNESS_SM] nmodel,nmodelactive",nmodel,nmodelactive
allocate( dmodel(nmodelactive),dmodel2(nmodelactive)) ! 2018.06.25
allocate( iactive(nmodel) )                           ! 2018.06.25
iactive      = h_model%iactive                        ! 2018.06.25

!#[2]# generate dmodel
kk=0       ! 2018.06.25
do i=1,nmodel
if ( iactive(i) .ne. 1 ) cycle ! 2018.06.25
 kk=kk+1   ! 2018.06.25
 dmodel(kk) = h_model%logrho_model(i) - g_model_ref%logrho_model(i) ! 2018.06.25
end do     ! 2018.06.25

!#[2]## cal roughness
call mul_matcrs_rv(BM,dmodel,nmodelactive,dmodel2) ! 2018.06.25
roughness = 0.d0
do i=1,nmodelactive                               ! 2018.06.25
 roughness = roughness + dmodel(i)*dmodel2(i)    ! 2018.06.25
end do

call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12
write(*,'(a,f8.4,a)') " ### CALROUGHNESS_SM END!!  Time=",t_watch%time," [min]" ! 2020.09.18

return
end


!############################################## subroutine GENRI_MSG
!# nmode is converted to nmodelactive on 2018.06.25
!# GENDQ is included on 2017.12.25
!# coded on 2017.12.18
subroutine GENRI_MSG(R,RI,g_face,g_mesh,g_model_ref) ! 2017.12.25
use matrix
use solveCM_ap! 2017.12.18
use caltime   ! 2017.12.22
use face_type ! 2017.12.25
use mesh_type ! 2017.12.25
use modelpart ! 2017.12.25
implicit none
type(real_crs_matrix),  intent(out)   :: R
type(real_crs_matrix),  intent(out)   :: RI
type(face_info),        intent(inout) :: g_face     ! 2017.12.25
type(mesh),             intent(in)    :: g_mesh     ! 2017.12.25
type(model),            intent(in)    :: g_model_ref! 2017.12.25
type(real_crs_matrix)                 :: RE         ! R with epsilon 2017.12.25
type(real_crs_matrix)                 :: RTR        ! This is not used 2017.12.25
integer(4),parameter                  :: itype_roughness = 3 ! MSG 2017.12.25
!# internal variables
type(PARDISO_PARAM)                   :: B        ! 2017.12.18
real(8)                               :: epsilon
integer(4)                            :: i,j
integer(4)                            :: nmodelactive ! 2018.06.25
type(watch) :: t_watch  ! 2017.12.22

call watchstart(t_watch) ! 2017.12.22
write(*,*) "### GENRI_MSG START!! ###" !2017.12.28

!#[0]##
 CALL GENDQ(R,RTR,g_face,g_mesh,g_model_ref,itype_roughness) ! see m_modelroughness.f90
 ! when itye_roughenss = 3, RTR is not allocated 2017.12.25

!#[1]## add small values to diagonals of R (Usui et al., 2017)
 RE = R  ! 2017.12.25
 epsilon =1.d-2 ! small value
 do i=1,RE%nrow
  do j=RE%stack(i-1)+1,RE%stack(i)
   if ( RE%item(j) .eq. i ) then ! only diagonal
   RE%val(j) = RE%val(j) + epsilon                 ! 2017.08.31
   end if
  end do
 end do

!#[3]## calculate R inverse
 nmodelactive = R%nrow                    ! 2018.06.25
 write(*,*) "nmodel check =",nmodelactive ! 2018.06.25
 call PARDISOphase1(RE,B)                 ! see m_solveCM_ap.f90 2017.12.25
 call PARDISOphase2(B)                    ! see m_solveCM_ap.f90
 call PARDISOphase3(B,nmodelactive,RI)    ! see m_solveCM_ap.f90 2018.06.25
 call PARDISOphase4(B)                    ! see m_solveCM_ap.f90

 call watchstop(t_watch)
 write(*,'(a,g15.7,a)') "### GENRI_MSG END!! Time=",t_watch%time,"[min] ###"

return
end
!############################################## subroutine GENBMI_SM
!# GENDQ is included on 2017.12.25
subroutine GENBMI_SM(BM,BMI,g_face,g_mesh,g_model) ! 2017.12.25
use matrix
use modelpart ! 2017.06.14
use mesh_type ! 2017.06.14
use fem_util  ! 2017.06.14
use solveCM_ap! 2017.06.16
use caltime   ! 2017.12.22
use face_type ! 2017.12.25
implicit none
type(real_crs_matrix), intent(out)   :: BM      ! 2017.12.25
type(real_crs_matrix), intent(out)   :: BMI     ! 2017.12.25
type(mesh),            intent(in)    :: g_mesh  ! 2017.06.14
type(model),           intent(in)    :: g_model ! 2017.06.14
type(face_info),       intent(inout) :: g_face     ! 2017.12.25
type(PARDISO_PARAM)                  :: B       ! 2017.06.16
type(real_crs_matrix)                :: R       ! 2017.12.25
type(real_crs_matrix)                :: RTR     ! 2017.12.25
type(real_crs_matrix)                :: crsin,crsout ! 2017.12.25
integer(4)                           :: nmodel,i,j,j1,j2,ntot
real(8)                              :: epsilon
real(8),allocatable,dimension(:)     :: volmodel    ! 2017.06.14
type(real_crs_matrix)                :: model2ele ! 2017.06.14
real(8)                              :: v,xx(3,6), elm_xyz(3,4) ! 2017.06.14
integer(8),allocatable,dimension(:)  :: stack,item
integer(4),parameter                 :: itype_roughness = 1 ! Smoothest Model 2017.12.25
type(watch) :: t_watch ! 2017.12.22

call watchstart(t_watch) ! 2017.12.22

!#[0]## GEN RTR
  CALL GENDQ(R,RTR,g_face,g_mesh,g_model,itype_roughness) ! see m_modelroughness.f90
  !# only RTR come out, and R is not allocated when itype_roughness = 1
  BM = RTR ! 2017.12 25

!#[1]## set
 nmodel = RTR%nrow
 ntot   = RTR%ntot
 allocate(stack(0:nmodel),item(ntot))
 stack  = RTR%stack
 item   = RTR%item
 allocate(volmodel(nmodel))      ! 2017.06.14
 model2ele = g_model%model2ele   ! 2017.06.14

!#[2]## calculate volume of each model 2017.06.14
 if ( .false. ) then ! only when the volume of models will be used 2017.12.18
 do i=1,nmodel
  volmodel(i)=0.d0
  do j=model2ele%stack(i-1)+1,model2ele%stack(i) ! element id
   do j1=1,4
    elm_xyz(1:3,j1)=g_mesh%xyz(1:3,g_mesh%n4(model2ele%item(j),j1))
   end do
   call calxmn(elm_xyz,xx)
   call volume(xx(:,1), xx(:,5), -xx(:,4), v) ! volume of this element [km^3]
   volmodel(i) = volmodel(i) + v
  end do
 end do
 end if

!#[3]## add small values to diagonals (Usui et al., 2017)
 epsilon =1.d-2 ! small value
 do i=1,nmodel
  do j=stack(i-1)+1,stack(i)
  if ( item(j) .eq. i ) then ! only diagonal
   RTR%val(j) = RTR%val(j) + epsilon                 ! 2017.08.31
   !  write(*,*) "i",i,"RTR%val(j)=",RTR%val(j),"volmodel(i)",volmodel(i)
   !  RTR%val(j) = RTR%val(j) + 10.*(volmodel(i)**3.d0) ! commented out 2017.08.31
  end if
  end do
 end do

!#[4]## calculate BMI

!write(*,*) "nmodel check =",nmodel
 call PARDISOphase1(RTR,B)   ! see m_solveCM_ap.f90
 call PARDISOphase2(B)       ! see m_solveCM_ap.f90
 call PARDISOphase3(B,nmodel,BMI) ! see m_solveCM_ap.f90 2017.12.25
 call PARDISOphase4(B)       ! see m_solveCM_ap.f90


!# check---------------------------------------------------------------------------
if ( .false. ) then
call mulreal_crs_crs_crs(RTR,BMI,crsout)
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
stop
end if

 call watchstop(t_watch) ! 2017.12.22
 write(*,'(a,f9.4,a)') " ### GENBMI_SM END!! ### Time=",t_watch%time," [min]"!2020.09.18
return
end

!############################################## subroutine genDQ ""
!#  DQ : [nphys2,nmodel]
!#  D :  [nphys2,nphys2]
!#  Q :  [nphys2,nmodel]
!#  output is RTR when ityperoughenss = 1 (SM) 2017.12.25
!#  output is R   when ityperoughness = 3 (MSG) 2017.12.25
!#  When R is output, RTR is empty (not allocated), vise versa. 2017.12.25
subroutine GENDQ(R,RTR,g_face,g_mesh,g_model,itype_roughness) ! 2017.12.14
use face_type
use matrix
use mesh_type
use modelpart
use fem_util ! 2017.05.18 for calculations of face areas
use caltime  ! 2017.09.06 see m_caltime.f90
implicit none
type(face_info),        intent(inout) :: g_face
integer(4),             intent(in)    :: itype_roughness ! 2017.12.14
type(mesh),             intent(in)    :: g_mesh
type(model),            intent(in)    :: g_model
type(real_crs_matrix),  intent(out)   :: R                ! 2017.12.14
type(real_crs_matrix),  intent(out)   :: RTR              ! 2017.12.25
type(real_crs_matrix)                 :: D,Q,DQT,DQ,crsout,QT,QTDQ,crs, DQ2 ! 2018.06.22
type(real_ccs_matrix)                 :: DQTCCS,QTCCS
integer(4),allocatable,dimension(:,:) :: n4face,n4flag,n4
integer(4),allocatable,dimension(:)   :: index, ele2model
integer(4),allocatable,dimension(:)   :: iactive   ! 2018.06.22
integer(4),allocatable,dimension(:)   :: icount
real(8),   allocatable,dimension(:,:) :: band
integer(4),allocatable,dimension(:,:) :: band_ind
real(8),   allocatable,dimension(:,:) :: xyz
integer(4)  :: ncolm,nrow,iface,iele,nphys1,nphys2,ntot,nc,n5(5)
integer(4)  :: nface,ntet,i,j,nmodel,icele,j1,j2,node,k
real(8)     :: r5(5)
real(8)     :: elm_xyz(3,4),a4(4)    ! 2017.05.18
type(watch) :: t_watch,t_watch_total ! 2017.09.06
!# element difference matrix H 2017.12.18
type(real_crs_matrix)     :: H            ! 2017.12.18
integer(4)                :: nface_nphys2 ! 2017.12.18
integer(4)                :: ifacecount   ! 2017.12.18
integer(4)                :: nface_model  ! 2017.12.18
integer(4)                :: icombine
integer(4)                :: nmodelactive ! 2018.03.16
integer(4)                :: ii,jj        ! 2018.03.18

call watchstart(t_watch_total) ! 2017.12.22

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
  n4flag    = g_mesh%n4flag ! 2017.06.18
  icombine     = g_model%icombine     ! 2018.06.22
  write(*,'(a,i7)') " GENDQ_AP nmodel",nmodel ! 2020.09.18
  allocate( iactive(nmodel) )         ! 2018.06.22
  iactive      = g_model%iactive      ! 0 not active, 1 active 2018.06.22
  nmodelactive = g_model%nmodelactive ! 2018.06.22

!#[2]## gen band matrix for D = Nele*Nele

  !#[2-1]## count only the land element connection
  ncolm = 5
  nrow  = nphys2
  allocate(icount(nrow),band(ncolm,nrow),band_ind(ncolm,nrow))
  ! count (iflag = 1) and assign (iflag = 2) loop
  nface_nphys2 = 0 ! 2017.12.18
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
    nface_nphys2 = nface_nphys2 + ( icount(i) -1 ) ! 2017.12.18
  end do

  !#[2-2]## generate element difference matrix H=[nface*nele] 2017.12.18
  !# [2-2] is not tested yet. This routine is found not required for MSG 2017.12.18
   if ( .false. ) then ! only for MSG
     !if ( itype_roughness .eq. 3) then
     !#[2-2-1]# allocate
     nface_nphys2 = nface_nphys2 / 2 ! eliminating double countintg
     write(*,*) "nface_nphys2=",nface_nphys2
     allocate(H%stack(0:nface_nphys2),H%item(nface_nphys2*2),H%val(nface_nphys2*2))
     H%stack(0)=0
     H%nrow   = nface_nphys2
     H%ncolm  = nphys2
     H%ntot   = nface_nphys2 * 2
     !#[2-2-2]# generate H
     ifacecount = 0
     do i=1,nphys2
       icele = band_ind(1,i)
       do j=2,icount(i)
         iele = band_ind(j,i)
         if ( icele .lt. iele ) then ! to avoid double conting
           ifacecount = ifacecount + 1
           H%stack(ifacecount:nface_nphys2) = H%stack(ifacecount:nface_nphys2) + 1
           H%item(2*(ifacecount-1)+1 ) = icele
           H%item(2*(ifacecount-1)+2 ) = iele
           H%val( 2*(ifacecount-1)+1 ) = -1.d0
           H%val( 2*(ifacecount-1)+2 ) =  1.d0
	       end if
        end do
      end do
      if ( ifacecount .ne. nface_nphys2 ) then
        write(*,*) "GEGEGE ifacecount",ifacecount,"nface_nphys2",nface_nphys2
        stop
      end if!--------------------------------------------  [2-2] for MSG end
   end if


  !#[2-4]## sort and store as crs matrix, D=[nphys2*nmodel]
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
   write(*,'(a)') " ## D      is generated!! ##" ! 2020.09.18

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
   write(*,'(a)') " ## Q      is generated!! ##" ! 2020.09.18

!#[4]## calculate DQ [nphys2,nmodel]
  call watchstart(t_watch)        ! 2017.09.06
  call mulreal_crs_crs_crs(D,Q,DQ,iflag_mkl=1) ! see m_matrix.f90 iflag_mkl=1 added 2022.01.27
  call watchstop(t_watch)         ! 2017.09.06
  write(*,'(a,f9.4,a)') " ## DQ     is generated!! ## Time =",t_watch%time,"[min]" ! 2020.09.18

  call trans_crs2ccs(Q,QTCCS)     ! 2017.06.16
  call conv_ccs2crs(QTCCS,QT)     ! 2017.06.16

  call watchstart(t_watch)        ! 2017.09.06
  call mulreal_crs_crs_crs(QT,DQ,QTDQ) ! 2017.06.16 for connectivity of models
  call watchstop(t_watch)        ! 2017.09.06
  write(*,'(a,f9.4,a)') " ## QTDQ   is generated!! ## Time =",t_watch%time,"[min]" ! 2020.09.18

  DQ = QTDQ ! 2017.06.16
  write(*,'(a)') " ## DQ     is generated!! ##"  ! 2020.09.18

  !# 2017.05.19
  !# make powerless the weights from # of interfaces between models
  !# like R in Usui et al. (2017)
  if (.true.) then ! 2017.12.14
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
  write(*,'(a)') " ## R      is generated!! ##"  ! 2020.09.18

!#[5]## When icombine = 2 ; modify DQ and g_model 2018.06.22
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
  end if   ! 2018.06.22


  !# for MSG, only R is required       ! 2017.12.14
  if ( itype_roughness .eq. 3 ) then   ! 2017.12.18
   R = DQ                              ! 2017.12.14
   goto 88                             ! 2017.12.22
  end if                               ! 2017.12.18

!#[5]## calculate RTR
  call trans_crs2ccs(DQ,DQTCCS)    ! DQ [ntet,nmodel]
  write(*,'(a)') " ## DQTCCS is generated!! ##" ! 2020.09.18
  call conv_ccs2crs(DQTCCS,DQT)
  write(*,'(a)') " ## DQT    is generated!! ##" ! 2020.09.18
  call mulreal_crs_crs_crs(DQT,DQ,RTR)
  write(*,'(a)') " ## RTR    is generated!! ##" ! 2020.09.18

!# output for check ========================================================== for check
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
!# output for check ========================================================== for check

88 continue
 call watchstop(t_watch_total) ! 2017.12.22
 write(*,'(a,f9.4,a)') " ### GENDQ END!! ### Time=",t_watch_total%time," [min]" ! 2020.09.18

return

99 continue
write(*,*) "GEGEGE! nmodelactive",nmodelactive,"ii",ii
stop

end

!### GENDQ END

end module modelroughness
