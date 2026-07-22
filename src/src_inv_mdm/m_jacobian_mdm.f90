!# Jacobian matrix for ACTIVE responses
!# coded on May 13, 2017
!#     | J2 |
!# J = | J3 | [Ndat,Nmodel]
!#     | vv |
!#     | JN | Ji is for i-th frequency
!#
!# Ji [Nobs,Nmodel]; j = obs index, k = model index
!# Ji[j,k]  = {d/dmk(|bz_i,j|)|bz_1,j| - |bz_i,j|d/dmk(|bz_1,j|)}/|bz_1,j|^2
!#
!# d/dmk(|bz_i,j|) = |bz_i,j| Re{(1/PT*Al)*PT*d/dmk(Al) }
!#
!# Au = P
!# P [nedge, nobs]
!# u [nedge, nobs]


module jacobian_mdm
use matrix
implicit none

!# 2017.06.07
 type amp_phase_dm
  integer(4) :: nobs     ! for each freqency
  integer(4) :: nsr_inv  ! # of sources,     2017.07.14
  real(8),              allocatable,dimension(:,:) :: ampbz    ! (nobs,nsr) 2017.07.14
  type(real_crs_matrix),allocatable,dimension(:)   :: dampbzdm ! (nsr) 2017.07.14
 end type

contains

!###########################
!# modified on 2017.07.14 for multiple sources
!# coded on 2017.06.07
subroutine allocateapdm(nobs,nfreq,nsr_inv,g_apdm) ! 2017.07.14
implicit none
integer(4),        intent(in)    :: nobs,nfreq
integer(4),        intent(in)    :: nsr_inv         ! 2017.07.14
type(amp_phase_dm),intent(inout) :: g_apdm(nfreq)
integer(4) :: i

do i=1,nfreq
 g_apdm(i)%nobs     = nobs
 g_apdm(i)%nsr_inv  = nsr_inv   ! 2017.07.14
 allocate( g_apdm(i)%ampbz(nobs,nsr_inv) ) ! 2017.07.14
 allocate( g_apdm(i)%dampbzdm(  nsr_inv) )   ! 2017.07.14
end do

return
end
!########################### gen jacobian
!# Modified on 2017.07.14 for multiple sources
!# Coded on 2017.05.16
subroutine genjacobian1(nobs,nline,nsr_inv,ut,bs,PT,g_model,h_mesh,l_line,omega,g_apdm)
use matrix            ! see m_matrix.f90
use modelpart         ! see m_modelpart.f90
use line_type         ! see m_line_type.f90
use constants         ! see m_constants.f90
use outerinnerproduct ! see m_outerinnerproduct.f90
use fem_util          ! for volume, intv, (see m_fem_utiil.f90 )
use fem_edge_util     ! see fem_edge_util.f90
use mesh_type         ! see m_mesh_type.f90
implicit none
integer(4),              intent(in)    :: nobs,nline
integer(4),              intent(in)    :: nsr_inv          ! # of sources 2017.07.14
type(complex_crs_matrix),intent(in)    :: ut               ! [nobs,nedge]
complex(8),              intent(in)    :: bs(nline,nsr_inv)! 2017.07.14
type(real_crs_matrix),   intent(in)    :: PT   ! [nobs,nline]
type(model),             intent(in)    :: g_model
type(line_info),         intent(in)    :: l_line
type(mesh),              intent(in)    :: h_mesh
real(8),                 intent(in)    :: omega
type(amp_phase_dm),      intent(inout) :: g_apdm             ! 2017.07.14
type(real_crs_matrix)                  :: dbzdm(nsr_inv)     ! 2017.07.14 (g_apdm%dampdm)
real(8)                                :: ampbz(nobs,nsr_inv)! 2017.09.03 (g_apdm&ampbz)
complex(8)                             :: bz(nobs,nsr_inv)   ! 2017.07.14
integer(4)                             :: nmodel,nphys2
integer(4)                             :: imodel,i,j,k,l,m,n,ii,jj
integer(4)                             :: iele,idirection(6),n6line(6)
integer(4),allocatable,dimension(:)    :: stack,item
real(8),   allocatable,dimension(:)    :: logrho_model
complex(8),allocatable,dimension(:,:)  :: AAL,dbzobs       ! 2017.07.14
real(8),   allocatable,dimension(:,:,:):: dbzdmfull        ! 2017.07.14
real(8)    :: threshold =1.d-10
!# for dA/dm
complex(8) :: S1(6,6)
real(8)    :: elm_xyz(3,4),yy,gn(3,4),v
real(8)    :: RM,sigma ! model = log10(rho)
complex(8) :: dBBdm, iunit=(0.d0,1.d0)
real(8), parameter  :: L0=1.d+3  ! [m]  scale length
integer(4) :: tstart,tend,t_rate,t_max  ! 2017.07.21
real(8)    :: t1,t2,t3,time             ! 2017.07.21
complex(8),allocatable,dimension(:,:) :: utfull  ! 2017.07.21
complex(8) :: z                         ! 2017.07.21

!#[0]## set
 nmodel       = g_model%nmodel
 nphys2       = g_model%nphys2
 allocate(stack(0:nmodel),item(nphys2))
 allocate(logrho_model(nmodel))
 stack        = g_model%model2ele%stack
 item         = g_model%model2ele%item  ! element id for whole element space
 logrho_model = g_model%logrho_model
 !
 allocate(utfull(ut%ncolm,ut%nrow)) ! 2017.07.21
 utfull = 0.d0   ! 2017.07.21
 do i=1,ut%nrow  ! 2017.07.21
  do j=ut%stack(i-1)+1,ut%stack(i)  ! 2017.07.21
   utfull(ut%item(j),i) = ut%val(j) ! 2017.07.21
  end do         ! 2017.07.21
 end do          ! 2017.07.21

!#[1]## cal bz and ampbz : complex
 do i=1,nsr_inv  ! 2017.07.14
  call mul_matcrs_cv(PT,bs(:,i),nline,bz(:,i)) ! 2017.07.14
 end do          ! 2017.07.14
 do j=1,nsr_inv  ! 2017.07.14
  do i=1,nobs
   ampbz(i,j) = sqrt(real(bz(i,j))**2.d0 + imag(bz(i,j))**2.d0) ! 2017.07.14
  end do
 end do      ! 2017.07.14

!#[2]## cal ut*dA/dm*Al
! allocate(AAL(nline,nsr_inv))
 allocate(AAL(6,nsr_inv))       ! 2017.07.21
 allocate(dbzobs(nobs,nsr_inv)) ! 2017.07.14
 allocate(dbzdmfull(nobs,nmodel,nsr_inv))          ! 2017.07.14

 t1 = 0.0 ![min] 2017.07.21
 t2 = 0.0 ![min] 2017.07.21
 t3 = 0.0 ![min] 2017.07.21
 do imodel=1,nmodel
  AAL = 0.d0    ! 2017.07.21
  dbzobs = 0.d0 ! 2017.07.21
  call system_clock(tstart) ! 2017.07.21
  do jj=stack(imodel-1)+1,stack(imodel) ! element loop for i-th model
    iele=item(jj)

   !#[2-1]## ! check the direction of edge, compared to the defined lines
   idirection(1:6)=1
   n6line(1:6) = l_line%n6line(iele,1:6)
   do j=1,6
     if ( n6line(j) .lt. 0 ) idirection(j)=-1
     n6line(j) = n6line(j)*idirection(j)
   end do
   do j=1,4
    elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j)) ! [km]
   end do
   call gradnodebasisfun(elm_xyz,gn,v) ! see fem_util.f90

   !#[2-2]## Second term from i * omega * mu * sigma * int{ sigma w cdot w }dv {Bsl}
   !# [4-1] ## assemble coefficient for i * omega*
   if ( h_mesh%n4flag(iele,1) .ne. 2 ) goto 99 ! in the case of not in land
   RM    = g_model%logrho_model(imodel) ! log10(rho)
   sigma = 10**(-RM)                     ! sigma = 1/10**(log10(rho)) = 10**(-M)
!   write(*,*) "imodel=",imodel,"sigma=",sigma,"RM=",RM

   !#[2-3]##
   !# BB     = iunit*omega*dmu*10^(-M)  (where M=log10(rho))
   !# dBB/dm = iunit*omega*dmu*{10^-M)}*(-log10)
   dBBdm =iunit*omega*dmu*sigma*L0**2.d0*(-log(10.d0))

   !#[2-4] ## assemble scheme No.2 ( analytical assembly)
   S1(:,:)=(0.d0,0.d0)
   do i=1,6
     do j=1,6
      k=kl(i,1);l=kl(i,2) ; m=kl(j,1) ; n=kl(j,2) ! gn*gn [km^-2], intv [km^3], yy[km]
      yy =     intv(k,m,v)*inner(gn(:,l), gn(:,n))   & ! first term
     &	-  intv(k,n,v)*inner(gn(:,l), gn(:,m))   & ! second term
     &      -  intv(l,m,v)*inner(gn(:,k), gn(:,n))   & ! third term
     &      +  intv(l,n,v)*inner(gn(:,k), gn(:,m))     ! forth term
      S1(i,j)= yy*idirection(i)*idirection(j)*dBBdm  ! S1 [km*rad/s*S/m]
     end do
   end do
   AAL = 0.d0
   do i=1,6
    ii = n6line(i)
    do j=1,6
!     write(*,101) "iele",iele,"i,j,",i,j,"bs(n6line)=",bs(n6line(j)),"S1(i,j)",S1(i,j)
     do k=1,nsr_inv   ! 2017.07.14
!      AAL(ii,k)=AAL(ii,k) + S1(i,j)*bs(n6line(j),k) ! (dA/dm)*Al 2017.07.14
       AAL(i,k)=AAL(i,k) + S1(i,j)*bs(n6line(j),k) ! (dA/dm)*Al 2017.07.14
     end do           ! 2017.07.14
    end do
   end do
   !# uT * (-AAL)
     do i=1,ut%nrow   ! 2017.07.21
      do k=1,nsr_inv  ! 2017.07.21
	do j=1,6        ! 2017.07.21
	 ii = n6line(j) ! 2017.07.21
	 dbzobs(i,k) = dbzobs(i,k) + utfull(ii,i)*(-AAL(j,k)) ! 2017.07.21
	end do          ! 2017.07.21
	end do          ! 2017.07.21
     end do           ! 2017.07.21
  end do ! element loop (iele) for i-th model
  call system_clock(tend,t_rate,t_max)     ! 2017.07.21
  call caltime(tstart,tend,t_rate,t_max,t1)! 2017.07.21

  call system_clock(tstart)                ! 2017.07.21
  !#[2-4]## cal ut*dA/dm*Al
  !# d/dm(|bz(iobs)|) = |bz(iobs)| Re{1/bz(iobs)*uT*(-dA/dm)*Al}
  !AAL(:,:) = - AAL(:,:) ! commented out on 2017.07.21

  do k=1,nsr_inv  ! 2017.07.14
!   call mul_matcrscomp_cv(ut,AAL(:,k),nline,dbzobs(:,k)) ! ut*AAL (ut) 2017.07.21
   do i=1,nobs
    dbzdmfull(i,imodel,k)=ampbz(i,k)*dreal(dbzobs(i,k)/bz(i,k))   ! 2017.07.14
!   write(*,100) "imodel=",imodel,"iobs=",i,"dbzobs=",dbzobs(i,k),&
!    & "bz=",bz(i,k),"dbzdmfull", dbzdmfull(i,imodel,k)
   end do
  end do       ! nsr loop, 2017.07.14

  call system_clock(tend,t_rate,t_max)       ! 2017.07.21
  call caltime(tstart,tend,t_rate,t_max,t2)  ! 2017.07.21

 end do  ! model loop   (imodel)
! write(*,*) "t1=",t1,"[min]"
! write(*,*) "t2=",t2,"[min]"
! write(*,*) "### GENJACOBIAN1 model loop end! ###"


!#[3]## full to real crs matrix
 call system_clock(tstart)            ! 2017.07.21
 do k=1,nsr_inv   ! 2017.07.14
  call conv_full2crs(dbzdmfull(:,:,k),nobs,nmodel,dbzdm(k),threshold) ! 2017.07.14
 end do       ! 2017.07.14
 call system_clock(tend,t_rate,t_max)     ! 2017.07.21
 call caltime(tstart,tend,t_rate,t_max,t3)! 2017.07.21
! write(*,*) "t3=",t3,"[min]"

! write(*,*) "dbzdm:" ! realcrs [nobs,nmodel]
! call realcrsout(dbzdm)

!#[4]## set output  2017.06.07
 g_apdm%nobs         = nobs
 g_apdm%ampbz        = ampbz  ! (nobs,nsr) 2017.07.14
 do k=1,nsr_inv                 ! 2017.07.14
  g_apdm%dampbzdm(k) = dbzdm(k) ! 2017.07.14
 end do                         ! 2017.07.14

 write(*,*) "### GENJACOBIAN1 END!! ###"
return

!# error #
99 continue
   write(*,*) "GEGEGE air region is included in nmodel. imodel=",imodel
   write(*,*) "iele=",iele,"h_mesh%n4flag(iele,1)=",h_mesh%n4flag(iele,1)
   do j=1,4
    write(*,*) "node=",h_mesh%n4(iele,j),"elm_xyz=",elm_xyz(1:3,j)
   end do
stop

!# format #
100 format(a,i3,a,i3,a,2g15.7,a,2g15.7,a,g15.7)
101 format(a,i3,a,2i3,a,2g15.7,a,2g15.7)

end
!############################################################# caltime
! 2017.07.21
subroutine caltime(t1,t2,t_rate,t_max,time)
implicit none
integer(4),intent(in)    :: t1,t2,t_rate,t_max
real(8),   intent(inout) :: time ! [min]
real(8)                  :: diff
  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  time = time + diff/dble(t_rate)/60.d0

return
end
!############################################################# genjacobian2
!# Modified on 2017.07.14 for multiple sources
!# Coded on 2017.05.17
subroutine genjacobian2(g_param_inv,nfreq,g_apdm,JJ,h_model,h_model_dm) ! 2017.10.31
use param_mdminv ! 2017.10.31
use matrix
implicit none
type(param_inversion),  intent(in)    :: g_param_inv     ! 2017.07.14
integer(4),             intent(in)    :: nfreq
type(amp_phase_dm),     intent(in)    :: g_apdm(nfreq)   ! 2017.06.07
type(model)      ,      intent(in)    :: h_model,h_model_dm
type(real_crs_matrix),  intent(out)   :: JJ(2,2)          ! Jacobian matrix [ndat,nmodel]
!=====================================================================
real(8),   allocatable,dimension(:,:,:)  :: ampbz(:,:,:)    ! (nobs,nsr,nfreq) 2017.07.14
type(real_crs_matrix),allocatable,dimension(:,:)   :: dampbzdm   ![nobs,nmodel](nsr,nfreq)
real(8),   allocatable, dimension(:)  :: dmbz,dmbzr      ! d|bz|/dm[1:nmodel](iobs,ifreq)
integer(4),allocatable, dimension(:)  :: stackr,itemr    ! for reference frequency
real(8),   allocatable, dimension(:)  :: valr
integer(4),allocatable, dimension(:)  :: stack,item      ! 2017.07.14
real(8),   allocatable, dimension(:)  :: val             ! 2017.07.14
real(8),   allocatable, dimension(:,:):: JP              ! full matrix of jacobian
integer(4)                            :: ndat            ! 2017.07.14
integer(4)                            :: nsr_inv         ! 2017.07.14
integer(4)                            :: nobs            ! 2017.07.14
integer(4)                            :: i,j,nmodel      ! 2017.07.14
integer(4)                            :: ntotr,ntot      ! 2017.07.17
integer(4)                            :: iobs,ifreq,isr  ! 2017.07.14
integer(4)                            :: icount          ! 2017.07.14
real(8)                               :: threshold = 1.d-10
integer(4)                            :: ndat1,ndat2     ! 2017.10.31
integer(4)                            :: ndatii          ! 2017.11.01
logical,allocatable,dimension(:,:,:)  :: data_avail      ! 2017.10.31
integer(4)                            :: ii,jm,imodel    ! 2017.10.31
integer(4)                            :: ishift,jshift   ! 2017.10.31
integer(4)                            :: nmodel1,nmodel2 ! 2017.10.31
integer(4),allocatable,dimension(:)   :: dm2modelptr     ! 2017.10.31

!#[0]## set
 ndat    = g_param_inv%ndat               ! 2017.07.14
 ndat1   = g_param_inv%ndat1              ! 2017.10.31
 ndat2   = g_param_inv%ndat2              ! 2017.10.31
 nobs    = g_param_inv%nobs               ! 2017.07.14
 nsr_inv = g_param_inv%nsr_inv            ! 2017.07.14
 nmodel1 = h_model%nmodel                 ! 2017.10.31
 nmodel2 = h_model_dm%nmodel              ! 2017.10.31
 allocate( dm2modelptr(nmodel2) )         ! 2017.10.31
 dm2modelptr = h_model_dm%dm2modelptr     ! 2017.10.31
 allocate( ampbz(nobs,nsr_inv,nfreq) )    ! 2017.07.14
 allocate( dampbzdm(nsr_inv,nfreq)   )    ! 2017.07.14
 do i=1,nfreq    ! 2017.06.07
  do j=1,nsr_inv ! 2017.07.14
   ampbz(1:nobs,j,i) = g_apdm(i)%ampbz(1:nobs,j) ! 2017.07.14
   dampbzdm(j,i)     = g_apdm(i)%dampbzdm(j)     ! 2017.07.14
  end do                                         ! 2017.07.14
 end do
 allocate(data_avail(nfreq,nobs,nsr_inv)) ! 2017.10.31


!#[2]## assemble J
!#  d(iobs)   = |bz(ifreq)|/|bz(i_reffreq)|
!#  JJ(ifreq)[iobs,imodel] = {(d/dm|bz_i|)(ifreq,iobs,imodel)*|bz(ifreq_ref,iobs)|
!#                          - (d/dm|bz(i_reffreq,iobs,imodel))|bz(ifreq,    iobs)|}
!#                           /|bz(ifreq_ref,iobs)|^2

 write(*,*) "nfreq",nfreq,"nsr_inv",nsr_inv,"nobs",nobs

!#
! d  = | y1 |  m = | m0 |
!      | y2 |      | dm |
! JP = | JP(1,1) JP(1,2) |
!      | JP(2,1) JP(2,2) |
! JP(1,1) = dy1/dm0 ; JP(1,2) = dy1/ddm
! JP(2,1) = dy2/dm0 ; JP(2,2) = dy2/ddm
!
 allocate( dmbz(nmodel1),dmbzr(nmodel1) )
 do ii=1,2   ! 2017.10.31 row loop
  if ( ii .eq. 1) data_avail     = g_param_inv%data_avail1 ! 2017.10.31
  if ( ii .eq. 2) data_avail     = g_param_inv%data_avail2 ! 2017.10.31
  if ( ii .eq. 1) ndatii         = ndat1      ! 2017.11.01
  if ( ii .eq. 2) ndatii         = ndat2      ! 2017.11.01
 do jm=1,2   ! 2017.10.31 colmn loop
  if ( jm .eq. 1 ) nmodel        = nmodel1    ! 2017.11.01
  if ( jm .eq. 2 ) nmodel        = nmodel2    ! 2017.11.01
  allocate( JP(ndatii,nmodel)           )     ! 2017.11.01

 icount = 0           ! 2017.07.14
 do ifreq = 2, nfreq  ! 2017.07.14
  do isr =1,nsr_inv    ! 2017.07.14
   !
   ntotr   = g_apdm(1)%dampbzdm(isr)%ntot     ! 2017.07.14
   allocate( stackr(0:nobs)             )
   allocate( itemr(ntotr), valr(ntotr)  ) ! for reference freq
   stackr = dampbzdm(isr,    1)%stack            ! for reference freq
   itemr  = dampbzdm(isr,    1)%item             ! for reference freq
   valr   = dampbzdm(isr,    1)%val              ! for reference freq
!   write(*,*) "ntotr",g_apdm(1)%dampbzdm(isr)%ntot
!   write(*,*) "ifreq,isr",ifreq,isr
!   write(*,*) "size(item)",size(itemr)
!   write(*,*) "size(val)",size(valr)
   !
   ntot = dampbzdm(isr,ifreq)%ntot         ! 2017.07.17
   allocate( stack(0:nobs)              )  ! 2017.07.14
   allocate( item( ntot), val( ntot)  ) ! for current freq 2017.07.17
   stack  = dampbzdm(isr,ifreq)%stack      ! 2017.07.14
   item   = dampbzdm(isr,ifreq)%item       ! 2017.07.14
   val    = dampbzdm(isr,ifreq)%val        ! 2017.07.14
 !  write(*,*) "ntot",g_apdm(ifreq)%dampbzdm(isr)%ntot
 !  write(*,*) "ifreq,isr",ifreq,isr
 !  write(*,*) "size(item)",size(item)
 !  write(*,*) "size(val)",size(val)
   do iobs=1,nobs    ! 2017.07.14
    if (data_avail(ifreq,iobs,isr)) then ! 2017.07.14
    icount = icount + 1                ! 2017.07.14
    !# set d|bz|/dm for reference and iobs
    dmbzr = 0.d0 ;    dmbz  = 0.d0
    do j = stackr(iobs-1)+1 , stackr(iobs)
     dmbzr( itemr(j) ) = valr(j)   ! d|bz|/dm for iobs [1: nmodel]
    end do
    do j = stack( iobs-1)+1 , stack( iobs)
     dmbz( item(j)   ) = val( j)
    end do
    do i=1,nmodel
     imodel = i                               ! 2017.10.31
     if ( jm .eq. 2 ) imodel = dm2modelptr(i) ! 2017.10.31
     JP(icount,i)=&
     &(dmbz(imodel)*ampbz(iobs,isr,1)- ampbz(iobs,isr,ifreq)*dmbzr(imodel))&
                &                /(ampbz(iobs,isr,1)**2.d0) ! 2017.10.31
!    write(*,*) "ifreq,imodel",ifreq,i,"JP(icount,i)=",JP(icount,i),"ampbz(iobs,1)",ampbz(iobs,isr,1),ampbz(iobs,isr,ifreq),dmbz(i),dmbzr(i)
    end do ! nmodel loop
    end if ! only for inv data     2017.07.14
   end do  ! obs loop              2017.07.14
   deallocate(stack, item, val )!  2017.07.17
   deallocate(stackr,itemr,valr)!  2017.07.17
  end do   ! source loop           2017.07.14
  end do   ! frequency loop        2017.07.14

  if ( ii .eq. 1 .and. jm .eq. 2) JP = -1./2.*JP    ! 2017.11.01
  if ( ii .eq. 1 .and. jm .eq. 2) JP =  1./2.*JP    ! 2017.11.01


 if (icount .ne. ndatii) then                       ! 2017.11.01
  write(*,*) "GEGEGE ndat=",ndatii,"icount=",icount ! 2017.11.01
  stop
 end if

!#[2]## full 2 crs matrix
 CALL conv_full2crs(JP,ndatii,nmodel,JJ(ii,jm),threshold)  ! 2017.11.01
 !write(*,*) "JJ="
 !call realcrsout(JJ)

   deallocate(JP)       ! 2017.11.01
  end do !2017.10.31 ii row loop
  end do !2017.10.31 jj colmn loop

 write(*,*) "### GENJACOBIAN2 END!! ###"
return
end

!########################################################## getnewmodel
!# coded on 2017.11.01
!#  [C]beta = X
!#  where
!#  [C] = [J*Cm*JT + alpha*Cd]
!#   X  = d_obs - d(M) - J(M - M_ref)
!# ---
!#  M_ite+1 = M_ref + Cm*JT*beta
!#  |
!#  | 2017.11.01
!#  v
!#  [C]beta' = X     (beta'=alpha*beta, CM'=CM/alpha)
!#  where
!#  [C] = [J*Cm'*JT + Cd]
!#   X  = d_obs - d(M) - J(M - M_ref)
!# ---
!#  M_ite+1 = M_ref + Cm'*JT*beta'
subroutine getnewmodel(JJ_P,g_model_ref,h_model,h_model_dm,&
                     & g_data,h_data,CM,CD,alpha1,alpha2) ! 2017.11.01
use modelpart
use matrix
use param_mdminv ! 2017.10.31
implicit none
type(real_crs_matrix),intent(in)    :: JJ_P(2,2)   ! [ndat,nmodel] 2017.11.01
type(model),          intent(in)    :: g_model_ref ! reference model
type(model),          intent(inout) :: h_model     ! old -> new
type(model),          intent(inout) :: h_model_dm  ! old -> new    2017.11.01
type(data_vec),       intent(in)    :: g_data ! observed
type(data_vec),       intent(in)    :: h_data ! calculated
type(real_crs_matrix),intent(in)    :: CM     ! model covariance matrix
type(real_crs_matrix),intent(in)    :: CD     ! data covariance matrix
real(8),              intent(in)    :: alpha1  ! 2017.11.01
real(8),              intent(in)    :: alpha2  ! 2017.11.01
type(real_crs_matrix)               :: C,JCMJT,CMJT2,crsout
type(real_ccs_matrix)               :: CMJT,JT
real(8),allocatable,dimension(:)    :: X, dobs, dcal,JM,beta   ! [ndat]
real(8),allocatable,dimension(:)    :: model1,model_ref,dmodel ! [nmodel]
integer(4)                          :: ndat,ndat1,ndat2
integer(4)                          :: nmodel,nmodel1,nmodel2,i
type(real_crs_matrix)               :: JJ,JJ1,JJ2,ACM    ! 2017.11.01
integer(4)                          :: ntot1,ntot2       ! 2017.11.01

!#[0]## set
  ndat      = g_data%ndat       ! 2017.11.01
  ndat1     = g_data%ndat1      ! 2017.11.01
  ndat2     = g_data%ndat2      ! 2017.11.01
  nmodel1   = h_model%nmodel    ! 2017.11.01 former
  nmodel2   = h_model_dm%nmodel ! 2017.11.01 latter
  nmodel    = nmodel1 + nmodel2 ! 2017.11.01
  allocate(X(ndat),dobs(ndat),dcal(ndat),JM(ndat))
  dobs      = g_data%dvec
  dcal      = h_data%dvec
  allocate( model1(nmodel),model_ref(nmodel) )        ! 2017.11.01
  allocate( dmodel(nmodel)                   )        ! 2017.11.01
  model1(1:nmodel1)        = h_model%logrho_model     ! 2017.11.01
  model1(nmodel1+1:nmodel) = h_model_dm%logrho_model  ! 2017.11.01
  model_ref(:)             = 0.d0                     ! 2017.11.01
  model_ref(1:nmodel1)     = g_model_ref%logrho_model ! 2017.11.01

!#[1]## assemble JJ and CM, and assemble  JCMJT
 call lateralcomb_real_crs_mat(JJ_P(1,1),JJ_P(1,2),JJ1) ! 2017.11.01
 call lateralcomb_real_crs_mat(JJ_P(2,1),JJ_P(2,2),JJ2) ! 2017.11.01
 JJ%nrow  = ndat
 JJ%ncolm = nmodel
 JJ%ntot  = JJ1%ntot + JJ2%ntot
 allocate(JJ%stack(0:ndat),JJ%item(JJ%ntot),JJ%val(JJ%ntot)) ! 2017.11.01
 JJ%stack(0)=0
 call combine_real_crs_mat(JJ,JJ1,      1) ! 2017.11.01
 call combine_real_crs_mat(JJ,JJ2,ndat1+1) ! 2017.11.01 JJ is assembled!
 call deallocate_real_crsmat(JJ1)          ! 2017.11.01
 call deallocate_real_crsmat(JJ2)          ! 2017.11.01
 call trans_crs2ccs(JJ,JT)

 !#[2]## reflect alpha1 and alpha2 in CM to generate ACM
 ACM   = CM                ! 2017.11.01
 ntot1 = CM%stack(nmodel1) ! 2017.11.01
 ntot2 = CM%stack(nmodel2) ! 2017.11.01
 ACM%val(      1:ntot1) = ACM%val(      1:ntot1)/alpha1 ! 2017.11.01
 ACM%val(ntot1+1:ntot2) = ACM%val(ntot1+1:ntot2)/alpha2 ! 2017.11.01

 !write(*,*) "CM="
 !call realcrsout(CM)
 !write(*,*) "JJ="
 !call realcrsout(JJ)
 call mulreal_crs_ccs_ccs(ACM,JT,CMJT)    ! 2017.11.01
 call mulreal_crs_ccs_crs(JJ,CMJT,JCMJT)
 call conv_ccs2crs(CMJT,CMJT2)            ! [nmodel,ndat]

 !write(*,*) "JCMJT:"
 !call realcrsout(JCMJT)
 !write(*,*) "ACD:"
 !call realcrsout(ACD)
 call add_crs_crs_crs(JCMJT,CD,C) ! C = JCM'JT + Cd
! write(*,*) "C:"
! call realcrsout(C)


!#[2]## gen X
!#   X  = d_obs - d(M) + J(M - M_ref) ! 2017.07.20
 dmodel = model1 - model_ref
! write(*,*) "model1, model_ref, dmodel"
! do i=1,nmodel
!  write(*,'(3g15.7)') model1(i),model_ref(i),dmodel(i)
! end do
 call mul_matcrs_rv(JJ,dmodel,nmodel,JM)
! write(*,*) "dobs=",dobs
! write(*,*) "dcal=",dcal
! write(*,*) "JM=",JM
 X = dobs - dcal + JM                 ! 2017.07.20
! write(*,*) "X=",X

!#[3]## cal beta
!# [C]beta = X
!# M_ite+1 = M_ref + Cm*JT*beta
 allocate(beta(ndat))                     ! [ndat,1]
 call solvebeta(C,X,ndat,beta)            ! [ndat,ndat]*[ndat,1] = [ndat,1]
 call mul_matcrs_rv(CMJT2,beta,ndat,dmodel) ! [nmodel,ndat]*[ndat,1]=[nmodel,1]

 !#[4]## output
  h_model%logrho_model(1:nmodel1)    = dmodel(1:nmodel1) + model_ref(1:nmodel1)
  h_model_dm%logrho_model(1:nmodel2) = dmodel(nmodel1+1:nmodel)


 write(*,*) "### GETNEWMODEL END!! ###"
return
end


end module
