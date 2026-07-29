!# Modified on 2017.09.03 for multisource inversion
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


module jacobian_apimp ! 2021.12.24
use matrix
implicit none

!# 2017.06.07
 type amp_phase_dm
  integer(4) :: nobs          ! for each freqency
  integer(4) :: nmodelactive  ! 2018.06.25
  integer(4) :: nsr_inv  ! # of sources,     2017.08.31
!  real(8),  allocatable,dimension(:) :: ampbz
  type(real_crs_matrix),allocatable,dimension(:,:) :: dampdm ! damp/dm (5,nsr) 2018.10.05
  type(real_crs_matrix),allocatable,dimension(:,:) :: dphadm ! dpha/dm (5,nsr) 2018.10.05
 end type



contains

!###########################
!# nmode -> nmodelactive on 2018.06.25
!# generated on 2017.08.31 for multiple sources
subroutine allocateapdm(nobs,nfreq,nmodelactive,nsr_inv,g_apdm) ! 2018.06.25
implicit none
integer(4),        intent(in)    :: nobs,nfreq
integer(4),        intent(in)    :: nmodelactive    ! 2017.09.04
integer(4),        intent(in)    :: nsr_inv         ! 2017.07.14
type(amp_phase_dm),intent(inout) :: g_apdm(nfreq)
integer(4) :: i

do i=1,nfreq
 g_apdm(i)%nobs         = nobs
 g_apdm(i)%nsr_inv      = nsr_inv        ! 2017.07.14
 g_apdm(i)%nmodelactive = nmodelactive   ! 2018.06.25
 allocate( g_apdm(i)%dampdm(5,nsr_inv) ) ! 2018.10.05
 allocate( g_apdm(i)%dphadm(5,nsr_inv) ) ! 2018.10.05
end do

return
end
!########################### Transxtomodel
!# coded on 2018.01.22
!# iflag = 1 : model to X
!# iflag = 2 : X to model
subroutine TRANSMODELX(h_model,g_param_apinv,iflag)
use param_apinv       ! see m_param_apinv.f90
use caltime
implicit none
integer(4)             , intent(in)     :: iflag
type(param_apinversion), intent(in)     :: g_param_apinv! see m_param_apinv.f90
type(model),             intent(inout)  :: h_model
integer(4) :: nmodel,i
real(8),allocatable,dimension(:) :: logrho_model
real(8) :: p_param
real(8) :: logrho_upper,logrho_lower
type(watch) :: t_watch

!#[1]## set
 call watchstart(t_watch)
 nmodel = h_model%nmodel
 allocate(logrho_model(nmodel))
 logrho_model = h_model%logrho_model
 p_param      = g_param_apinv%p_param
 logrho_upper = g_param_apinv%logrho_upper
 logrho_lower = g_param_apinv%logrho_lower

!#[2]## transform
 if ( iflag .eq. 1 ) then ! model -> X
  do i=1,nmodel
  logrho_model(i) = 1./p_param*&
  & log( (logrho_model(i)-logrho_lower)/(logrho_upper-logrho_model(i)))
  end do
 end if

 if ( iflag .eq. 2 ) then ! X to model
  do i=1,nmodel
  logrho_model(i) = (logrho_lower + logrho_upper*exp(p_param*logrho_model(i)))&
                                      & /(1.d0 + exp(p_param*logrho_model(i)))
  end do
 end if

!#[3]## set output
 h_model%logrho_model = logrho_model

 call watchstop(t_watch) ! 2018.01.23
 write(*,'(a,g15.7,a)') "### TRANSMODELX       END!! Time=",t_watch%time,"[min]"
return
end

!########################### Transferjacob 2018.01.22
!# 2018.01.22
subroutine transjacob(g_param_apinv,JJ,h_model)
use caltime           ! see m_caltime.f90 2017.09.06
use param_apinv       ! see m_param_apinv.f90
use modelpart            !
implicit none
type(real_crs_matrix),   intent(inout)  :: JJ          ! 2017.06.07
type(param_apinversion), intent(in)     :: g_param_apinv! see m_param_apinv.f90
type(model),             intent(in)     :: h_model
real(8),    allocatable,dimension(:)    :: logrho_model,dmdx
integer(4)                              :: i,ii,ntot,nrow,ncolm,nmodel ! 2018.06.25
integer(4)                              :: nmodelactive ! 2018.06.25
integer(4), allocatable,dimension(:)    :: iactive   ! 2018.06.25
integer(4), allocatable,dimension(:)    :: ptrnmodel ! 2018.06.25
real(8)     :: p_param, mi
real(8)     :: logrho_upper
real(8)     :: logrho_lower
type(watch) :: t_watch

!#[0]## set
 call watchstart(t_watch) ! 2018.01.23
 nmodel       = h_model%nmodel
 nmodelactive = h_model%nmodelactive ! 2018.06.25
 allocate( iactive(nmodel))          ! 2018.06.25
 iactive      = h_model%iactive      ! 2018.06.25
 allocate( ptrnmodel(nmodelactive))  ! 2018.06.25
 allocate(logrho_model(nmodel),dmdx(nmodel))
 logrho_model = h_model%logrho_model
 p_param      = g_param_apinv%p_param
 logrho_upper = g_param_apinv%logrho_upper
 logrho_lower = g_param_apinv%logrho_lower
 ntot         = JJ%ntot
 nrow         = JJ%nrow
 ncolm        = JJ%ncolm
 if ( ncolm .ne. nmodelactive ) then ! 2018.06.25
 write(*,*) "GEGEGE! nrow",nrow,"nmodelactive",nmodelactive,"nmodel",nmodel ! 2018.06.26
 stop
 end if

!#[1]## set ptrnmodel 2018.06.25
 ii = 0
 do i=1,nmodel
  if ( iactive(i) .eq. 1 ) then
   ii=ii+1
   ptrnmodel(ii) = i
  end if
 end do

!#[2]## gen dm/dx
 do i=1,nmodel
  mi   = logrho_model(i)
  dmdx(i) = p_param*(logrho_upper - mi)*(mi-logrho_lower)/(logrho_upper - logrho_lower)
 end do

!#[3]## J-> J'
 do i=1,ntot
  JJ%val(i)=JJ%val(i)*dmdx(ptrnmodel(JJ%item(i))) ! 2018.06.25
 end do

 call watchstop(t_watch) ! 2018.01.23
 write(*,'(a,g15.7,a)') "### TRANSJACOB END!! ### Time=",t_watch%time,"[min]"

return
end
!########################### gen jacobian
!# nmodelactive is introduced on 2018.06.22
!# Modified on 2017.09.03
!# Coded on 2017.05.16
subroutine genjacobian1(nobs,nline,nsr_inv,ut,bs,PT,g_model,h_mesh,l_line,omega,&
                      & g_apdm,g_param_apinv,ip,np) !2020.09.18
use matrix            ! see m_matrix.f90
use modelpart         ! see m_modelpart.f90
use line_type         ! see m_line_type.f90
use constants         ! see m_constants.f90
use outerinnerproduct ! see m_outerinnerproduct.f90
use fem_util          ! for volume, intv, (see m_fem_utiil.f90 )
use fem_edge_util     ! see fem_edge_util.f90
use mesh_type         ! see m_mesh_type.f90
use caltime           ! see m_caltime.f90 2017.09.06
use param_apinv       ! see m_param_apinv.f90  2018.10.05
implicit none
integer(4),              intent(in)     :: ip,np            ! 2020.09.18
type(param_apinversion), intent(in)     :: g_param_apinv    ! 2018.10.05
integer(4),              intent(in)     :: nobs,nline
integer(4),              intent(in)     :: nsr_inv          ! # of sources 2017.09.03
type(complex_crs_matrix),intent(in)     :: ut(5)            ! [nobs,nline]  2018.10.05
type(complex_crs_matrix)                :: utt              ! 2018.10.05
complex(8),              intent(in)     :: bs(nline,nsr_inv)!        2017.09.03
type(real_crs_matrix),   intent(in)     :: PT(5)            ! [nobs,nline] 2018.10.05
type(model),             intent(in)     :: g_model
type(line_info),         intent(in)     :: l_line
type(mesh),              intent(in)     :: h_mesh
real(8),                 intent(in)     :: omega
type(amp_phase_dm),      intent(inout)  :: g_apdm           ! 2017.06.07
type(real_crs_matrix)                   :: dampdm(nsr_inv)  ! 2017.09.03 (->g_apdm)
type(real_crs_matrix)                   :: dphadm(nsr_inv)  ! 2017.09.03 (->g_apdm)
complex(8),     dimension(nobs,nsr_inv) :: be               ! 2018.10.05
integer(4)                              :: nmodel,nphys2
integer(4)                              :: imodel,i,j,k,l,m,n,ii,jj,kk!2018.06.21
integer(4)                              :: iele,idirection(6),n6line(6)
integer(4),allocatable,dimension(:)     :: stack,item
real(8),   allocatable,dimension(:)     :: logrho_model
complex(8),allocatable,dimension(:,:)   :: AAL,dbeobs       ! 2017.09.03
real(8),   allocatable,dimension(:,:,:) :: dampdmfull       ! 2017.09.03
real(8),   allocatable,dimension(:,:,:) :: dphadmfull       ! 2017.09.03
real(8)                                 :: threshold =1.d-10
integer(4)                              :: nmodelactive     ! 2018.06.21
integer(4),allocatable,dimension(:)     :: iactive          ! 2018.06.21
integer(4),            dimension(5)     :: iflag_comp       ! 2018.10.05
integer(4)                              :: icomp            ! 2018.10.05
!# for dA/dm
complex(8)                              :: S1(6,6)
real(8)                                 :: elm_xyz(3,4),yy,gn(3,4),v
real(8)                                 :: RM,sigma         ! model = log10(rho)
complex(8)                              :: dBBdm, iunit=(0.d0,1.d0)
real(8),   parameter                    :: L0=1.d+3         ! [m]  scale length
complex(8),allocatable,dimension(:,:)   :: utfull           ! 2017.09.03
complex(8)                              :: z                ! 2017.09.03
type(watch)                             :: t_watch          ! see m_caltime.f90 2017.09.06

 call watchstart(t_watch) ! see m_caltime.f90 2017.12.22

!#[0]## set
 nmodel       = g_model%nmodel
 nmodelactive = g_model%nmodelactive ! 2018.06.21
 allocate( iactive(nmodel) )         ! 2018.06.21
 iactive      = g_model%iactive      ! 2018.06.21
 nphys2       = g_model%nphys2
 allocate(stack(0:nmodel),item(nphys2))
 allocate(logrho_model(nmodel))
 stack        = g_model%model2ele%stack
 item         = g_model%model2ele%item  ! element id for whole element space
 logrho_model = g_model%logrho_model
 iflag_comp   = g_param_apinv%iflag_comp ! 2018.10.05
 allocate(utfull(nline,nobs))        ! 2018.10.05
 allocate( AAL(6,nsr_inv) )
 allocate( dbeobs(nobs,nsr_inv)            )       ! 2018.10.05
 allocate( dampdmfull(nobs,nmodelactive,nsr_inv) ) ! 2018.06.22
 allocate( dphadmfull(nobs,nmodelactive,nsr_inv) ) ! 2018.06.22

!#[1]## component loop start          2018.10.05
 do icomp = 1,5                     ! 2018.10.05
  if (iflag_comp(icomp) .eq. 0)cycle! 2018.10.05

!#[1]## cal bx,by,bz,ex,ey dependent on iflag_comp
  utt    = ut(icomp)                 ! 2018.10.05
  utfull = 0.d0                      ! 2017.09.03
  do   i = 1,utt%nrow                ! 2017.09.03
   do j=utt%stack(i-1)+1,utt%stack(i)! 2017.09.03
    utfull(utt%item(j),i) = utt%val(j) ! 2017.09.03
   end do                            ! 2017.09.03
  end do                             ! 2017.09.03
  if ( icomp .ge. 4 ) utfull = -iunit*omega*utfull ! only for Ex,Ey, 2018.10.05

!#[1]## cal either of bx,by,bz,ex,ey dependent on icomp
  do i=1,nsr_inv                     ! 2017.09.03
   call mul_matcrs_cv(PT(icomp),bs(:,i),nline,be(:,i)) ! 2018.10.05
  end do
  if ( icomp .ge. 4 ) be     = -iunit*omega*be        ! only for Ex,Ey, 2018.10.05

!#[2]## cal ut*dA/dm*Al

 kk=0 ! 2018.06.21
 do imodel=1,nmodel
  if ( iactive(imodel) .ne. 1 ) cycle ! 2018.06.22
  kk = kk + 1                         ! 2018.06.22
  AAL    = 0.d0                       ! 2017.09.03
  dbeobs = 0.d0                       ! 2018.10.05 either of bx,by,bz,ex,ey
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
   if ( h_mesh%n4flag(iele,1) .lt. 2 ) goto 99 ! in the case of not in land
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
   AAL = 0.d0  ! 2017.09.03
   do i=1,6
    ii = n6line(i)
    do j=1,6
!     write(*,101) "iele",iele,"i,j,",i,j,"bs(n6line)=",bs(n6line(j)),"S1(i,j)",S1(i,j)
     do k=1,nsr_inv    ! 2017.09.03
     AAL(i,k)=AAL(i,k) + S1(i,j)*bs(n6line(j),k) ! (dA/dm)*Al 2017.09.03
     end do            ! 2017.09.03
    end do
   end do
   !# uT * (-AAL)
     do i=1,utt%nrow  ! 2018.10.05
      do k=1,nsr_inv  ! 2017.09.03
	do j=1,6        ! 2017.09.03
	 ii = n6line(j) ! 2017.09.03
	 dbeobs(i,k) = dbeobs(i,k) + utfull(ii,i)*(-AAL(j,k)) ! 2017.07.21
	end do          ! 2017.09.03
	end do          ! 2017.09.03
     end do           ! 2017.09.03
  end do ! element loop (iele) for i-thmodel

  call watchstart(t_watch) ! 2017.09.06
  !#[2-4]## cal ut*dA/dm*Al
  !# d/dm(|bz(iobs)|) = |bz(iobs)| Re{1/bz(iobs)*uT*(-dA/dm)*Al}

  do k=1,nsr_inv  ! 2017.09.03
   do i=1,nobs
    dampdmfull(i,kk,k)=1./log(10.d0)*dreal(dbeobs(i,k)/be(i,k)) ! 2018.10.05
    dphadmfull(i,kk,k)=180./pi      *dimag(dbeobs(i,k)/be(i,k)) ! 2018.10.05
!   write(*,100) "imodel=",imodel,"iobs=",i,"dbeobs=",dbeobs(i),&
!    & "bz=",bz(i),"dbzdmfull", dbzdmfull(i,imodel)
   end do
  end do       ! nsr loop, 2017.09.03

 end do  ! model loop   (imodel)

!#[3]## full to real_crs_matrix
 do k=1,nsr_inv   ! 2017.09.03
  call conv_full2crs(dampdmfull(:,:,k),nobs,nmodelactive,dampdm(k),threshold)!2017.06.22
  call conv_full2crs(dphadmfull(:,:,k),nobs,nmodelactive,dphadm(k),threshold)!2018.06.22
 end do           ! 2017.09.03

! write(*,*) "dbzdm:" ! realcrs [nobs,nmodel]
! call realcrsout(dbzdm)

!#[4]## set output  2017.06.07
 g_apdm%nobs     = nobs
 do k=1,nsr_inv ! 2017.09.03
  g_apdm%dampdm(icomp,k) = dampdm(k) ! 2018.10.05
  g_apdm%dphadm(icomp,k) = dphadm(k) ! 2018.10.05
 end do

 end do ! comp loop end 2018.10.05

 call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12

!write(*,'(a,i2,a,i2,a,f9.4,a)') " ### GENJACOBIAN1   END !! ###  ip =",ip," /",np," Time =",t_watch%time," [min]" !2020.09.18
write(*,'(a,i2,a,i2)') " ### GENJACOBIAN1   END !! ###  ip =",ip !" /",np," Time =",t_watch%time," [min]" !2020.09.18
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

!############################################################# genjacobian2
!# modified on 2018.10.05 for multiple components
!# nmodelactive is introduced on 2018.06.25
!# modified on 2017.09.04 for multiple source inversion
!# Coded on 2017.05.17
subroutine genjacobian2(g_param_apinv,nfreq,g_apdm,JJ) ! 2018.06.25
use param_apinv ! 2017.09.04
use matrix
use caltime     ! 2017.12.22
implicit none
type(param_apinversion),   intent(in)    :: g_param_apinv   ! 2017.09.04
integer(4),                intent(in)    :: nfreq           ! 2017.09.04
type(amp_phase_dm),        intent(in)    :: g_apdm(nfreq)   ! 2017.06.07
type(real_crs_matrix),     intent(out)   :: JJ              ! Jacobian matrix [ndat,nmodel]
type(real_crs_matrix),allocatable,dimension(:,:,:,:)   :: dapdm       !(2,5,nsr,nfreq) [nobs,nmodel]
logical,              allocatable,dimension(:,:,:,:,:) :: data_avail  ! 2018.10.05
integer(4)                               :: ii,i,nmodel,ntotr,nsr_inv, nobs    ! 2017.09.04
integer(4)                               :: j,k,l,nsft,lsft,ncount,ndat        ! 2017.09.04
integer(4)                               :: nmodelactive,icomp    ! 2018.10.05
integer(4),                dimension(5)  :: iflag_comp            ! 2018.10.05
type(watch) :: t_watch  ! see m_caltime.f90 2017.12.22

 call watchstart(t_watch) ! see m_caltime.f90 2017.12.22

!#[0]## set
 ndat           = g_param_apinv%ndat          ! 2017.09.04
 nobs           = g_param_apinv%nobs          ! 2017.09.04
 nsr_inv        = g_param_apinv%nsr_inv       ! 2017.09.04
 nmodelactive   = g_apdm(1)%nmodelactive      ! 2018.06.25
 iflag_comp     = g_param_apinv%iflag_comp    ! 2018.10.05
 allocate(dapdm(2,5,nsr_inv,nfreq))           ! 2018.10.05
 allocate(data_avail(2,5,nfreq,nobs,nsr_inv)) !2018.10.05
 data_avail = g_param_apinv%data_avail      ! 2017.09.04
 do k=1,nsr_inv ! 2017.09.04
  do i=1,nfreq  ! 2017.06.07
   do icomp = 1,5 ! 2018.10.05
    if ( iflag_comp(icomp) .eq. 0 ) cycle          ! 2018.10.05
    dapdm(1,icomp,k,i) = g_apdm(i)%dampdm(icomp,k) ! 2017.10.05
    dapdm(2,icomp,k,i) = g_apdm(i)%dphadm(icomp,k) ! 2018.10.05
   end do         ! 2018.10.05
  end do
 end do

!#[1]## count ntotr; total number of elements in JJ
 ntotr = 0
 do i = 1,nfreq
  do k= 1, nsr_inv  ! 2017.09.04
   do j = 1,nobs    ! 2017.09.04
    do icomp = 1,5  ! 2018.10.05
     do l = 1,2     ! 2017.09.04
      if (data_avail(l,icomp,i,j,k)) then    ! 2018.10.05
       ntotr = ntotr + dapdm(l,icomp,k,i)%stack(j) - dapdm(l,icomp,k,i)%stack(j-1) ! 2018.10.05
      end if
     end do
    end do          ! 2018.10.05
   end do
  end do
 end do

!#[2]## assemble J
 JJ%nrow  = ndat         ! 2017.09.04
 JJ%ntot  = ntotr        ! 2017.09.04
 JJ%ncolm = nmodelactive ! 2018.06.25
 allocate(JJ%stack(0:JJ%nrow),JJ%item(ntotr),JJ%val(ntotr))
 JJ%stack(:)=0
 ii = 0            ! 2017.09.04
 do i =1,nfreq
  do k=1,nsr_inv   ! 2017.09.04
   do j = 1,nobs   ! 2017.09.04
    do icomp = 1,5 ! 2018.10.5
     do l = 1,2    ! 2017.09.04
      if ( data_avail(l,icomp,i,j,k) ) then  ! 2018.10.05
       ii = ii + 1 ! current row index
	 ncount = dapdm(l,icomp,k,i)%stack(j) - dapdm(l,icomp,k,i)%stack(j-1)
	 JJ%stack(ii) = JJ%stack(ii-1) + ncount ! 2017.09.04
	 nsft = JJ%stack(ii-1)                  ! 2017.09.04
	 lsft = dapdm(l,icomp,k,i)%stack(j-1)         ! 2017.09.04
	 JJ%item(nsft+1:nsft+ncount) = dapdm(l,icomp,k,i)%item(lsft+1:lsft+ncount)
	 JJ%val( nsft+1:nsft+ncount) = dapdm(l,icomp,k,i)%val( lsft+1:lsft+ncount)
     end if
    end do ! l     loop 2017.09.04
    end do ! icomp loop 2018.10.05
   end do  ! j     loop 2017.09.04
  end do   ! k     loop 2017.09.04
 end do    ! i     loop 2017.09.04
 if ( ii .ne. ndat) then
  write(*,*) "GEGEGE ii",ii,"should be ndat",ndat,"!!!"
  stop
 end if

!#[3]## set output
 call watchstop(t_watch)  ! see m_caltime.f90 2017.12.12
 write(*,'(a,f8.4,a)') " ### GENJACOBIAN2 END!! ### Time =",t_watch%time," [min]"!2020.09.18
return
end

!######################################################### OUTJACOB
!# modified for multiple components on 2018.10.05
!# coded on 2018.06.25
subroutine OUTJACOB(g_param_apinv,JJ,ite,g_model,g_param,sparam)
use param_apinv
use matrix
use modelpart
use outresp  ! for free_unit ./common/m_outresp.f90
use param
implicit none
type(real_crs_matrix),    intent(in)        :: JJ
type(param_apinversion),  intent(in)        :: g_param_apinv
type(param_source),       intent(in)        :: sparam
integer(4),               intent(in)        :: ite
type(model),              intent(in)        :: g_model
type(param_forward),      intent(in)        :: g_param
character(50)                               :: head,outfile,site,sour
character(2)                                :: num,nf
integer(4)                                  :: i,j,k,l,m,n,ii,nh ! 2018.06.25
logical,   allocatable,dimension(:,:,:,:,:) :: data_avail        ! 2018.10.05
integer(4)                                  :: nfreq,nobs,nsr_inv
integer(4)                                  :: nmodel,nmodelactive
integer(4),allocatable,dimension(:)         :: iactive      ! 2018.06.25
integer(4),allocatable,dimension(:)         :: ptrnmodel    ! 2018.06.25
integer(4)                                  :: idev,nsi,nso ! 2018.10.05
integer(4),allocatable,dimension(:,:,:,:,:) :: icount       ! 2018.10.05
real(8),   allocatable,dimension(:,:,:)     :: JP           ! 2018.06.25
real(8),   allocatable,dimension(:)         :: freq
integer(4),allocatable,dimension(:)         :: srcindex
integer(4)                                  :: icomp        ! 2018.10.05
integer(4),            dimension(5)         :: iflag_comp   ! 2018.10.05

!#[0]## set
 nfreq        = g_param_apinv%nfreq           ! 2018.06.25
 nobs         = g_param_apinv%nobs            ! 2018.06.25
 nsr_inv      = g_param_apinv%nsr_inv         ! 2018.06.25
 head         = g_param_apinv%outputfolder    ! 2018.06.25
 nmodel       = g_model%nmodel                ! 2018.06.25
 nmodelactive = g_model%nmodelactive          ! 2018.06.25
 allocate( iactive(  nmodel)      )           ! 2018.06.25
 iactive      = g_model%iactive               ! 2018.06.25
 allocate( ptrnmodel(nmodelactive))           ! 2018.06.25
 nh           = len_trim(head)
 allocate(srcindex(nsr_inv), freq(nfreq) )    ! 2017.07.14
 srcindex = g_param_apinv%srcindex            ! 2018.06.25
 write(num,'(i2.2)') ite                      ! 2017.07.14
 freq         = g_param%freq                  ! 2017.07.14
 allocate( data_avail(2,5,nfreq,nobs,nsr_inv))! 2018.10.05
 data_avail   = g_param_apinv%data_avail      ! 2018.06.25
 if ( nmodelactive .ne. JJ%ncolm ) then
  write(*,*) "GEGEGE nmodelactive",nmodelactive,"JJ%ncolm",JJ%ncolm
 end if

!#[1]## set ptractive
 ii=0              ! 2018.06.25
 do i=1,nmodel     ! 2018.06.25
  if ( iactive(i) .eq. 1 ) then ! 2018.06.25
   ii = ii + 1     ! 2018.06.25
   ptrnmodel(ii)=i ! 2018.06.25
  end if           ! 2018.06.25
 end do            ! 2018.06.25
 iflag_comp = g_param_apinv%iflag_comp ! 2018.10.05

!#[2]## output
 allocate(JP(2,nfreq,nmodel))   ! 2018.06.25

!#[3]## generate icount
 ii = 0 ! 2018.06.25
 allocate(icount(2,5,nfreq,nobs,nsr_inv)) ! 2018.10.05
 do i=1,nfreq
 do k=1,nsr_inv
 do l=1,nobs
  do icomp = 1,5                      ! 2018.10.05
  do m=1,2                            ! 2018.06.25
  if (data_avail(m,icomp,i,l,k)) then ! 2018.10.05
   ii = ii + 1
   icount(m,icomp,i,l,k) = ii         ! 2018.10.05
  end if
 end do                               ! 2018.06.25
 end do                               ! 2018.10.05
 end do
 end do
 end do

!#[4]## generate JP and output
 do icomp = 1,5                       ! 2018.10.05
  if (iflag_comp(icomp) .eq. 0) cycle ! 2018.10.05

  do k=1,nsr_inv
   do l=1,nobs

   site    = g_param%obsname(l)
   nsi   = len_trim(site)
   sour    = sparam%sourcename(srcindex(k))
   nso   = len_trim(sour)
   outfile = head(1:nh)//"Jacob"//site(1:nsi)//"_"//sour(1:nso)//"_"//comp(icomp)//num(1:2)//".dat"!2018.10.05
   idev    = free_unit()   ! see m_outresp.f90

   !# gen JP
   JP = 0.d0 ! 2018.06.25
   do m=1,2  ! 1:amp, 2:phase
   do i=1,nfreq
    if ( data_avail(m,icomp,i,l,k) ) then        ! 2018.10.05
     ii = icount(m,icomp,i,l,k)                  ! 2018.10.05
     do j=JJ%stack(ii-1)+1,JJ%stack(ii)
      JP(m,i,ptrnmodel(JJ%item(j))) = JJ%val(j)  ! 2018.06.25
     end do
    end if
   end do
   end do    ! 2018.06.25

   !# output JP
   open(idev,file=outfile)
    write(idev,'(8g15.7)') freq(1:nfreq)
    write(idev,*) nmodel
    do m=1,2      ! 2018.06.25 1:amp, 2:phase
     do n=1,nmodel
      write(idev,'(8g15.7)') JP(m,1:nfreq,n)
     end do
    end do        ! 2018.06.25
   close(idev)

   end do
  end do
 end do           ! icomp loop   2018.10.05

return
end
!########################################################## getnewmodel
!# nmodelactive is introduced on 2018.06.25
!# coded on 2017.05.17
!#  [C]beta = X
!#  where
!#  [C] = [J*Cm*JT + alpha*Cd]
!#   X  = d_obs - d(M) - J(M - M_ref)
!# ---
!#  M_ite+1 = M_ref + Cm*JT*beta
!# 
!#
subroutine getnewmodel(JJ,g_model_ref,h_model,g_data,h_data,CM,CD,alpha,g_param_apinv)
use modelpart
use matrix
use param_apinv ! 2017.06.09
use caltime ! 2017.12.22
implicit none
type(param_apinversion),intent(in)  :: g_param_apinv
real(8),              intent(in)    :: alpha
type(real_crs_matrix),intent(in)    :: JJ          ! [ndat,nmodel]
type(model),          intent(in)    :: g_model_ref ! reference model
type(model),          intent(inout) :: h_model     ! old -> new
type(data_vec_ap),    intent(in)    :: g_data      ! observed
type(data_vec_ap),    intent(in)    :: h_data      ! calculated
type(real_crs_matrix),intent(in)    :: CM          ! model covariance matrix
type(real_crs_matrix),intent(in)    :: CD          ! data covariance matrix
type(real_crs_matrix)               :: C,JCMJT,ACD,crsout      ! 2018.01.23
type(real_crs_matrix)               :: CMJT, JT    ! 2018.01.23
real(8),allocatable,dimension(:)    :: X, dobs, dcal,JM,beta   ! [ndat]
real(8),allocatable,dimension(:)    :: model1,model_ref,dmodel ! [nmodel]
integer(4)                          :: ndat,nmodel,i,ii ! 2018.06.25
integer(4)                          :: nmodelactive ! 2018.06.25
integer(4),allocatable,dimension(:) :: iactive      ! 2018.06.25
integer(4)                          :: iboundflag
integer(4)                          :: icheck = 0
type(model)                         :: ki,kiref
type(model)                         :: m0,m_ref
type(watch) :: t_watch ! 2017.12.22
type(watch) :: t_watch1 ! 2018.01.23

 call watchstart(t_watch) ! 2017.12.22
 write(*,'(a,f9.4)') " alpha =",alpha," is adopted in GETNEWMODEL"  ! 2020.09.29

!#[0]## set
  call watchstart(t_watch1)
  ndat         = g_data%ndat
  nmodel       = g_model_ref%nmodel
  nmodelactive = g_model_ref%nmodelactive ! 2018.06.25
  allocate( iactive(nmodel) )             ! 2018.06.25
  iactive      = g_model_ref%iactive      ! 2018.06.25
  allocate(X(ndat),dobs(ndat),dcal(ndat),JM(ndat))
  dobs         = g_data%dvec
  dcal         = h_data%dvec
  allocate( model1(   nmodelactive) )     ! 2018.06.25
  allocate( model_ref(nmodelactive) )     ! 2018.06.25
  allocate( dmodel(   nmodelactive) )     ! 2018.06.25
  !#
  iboundflag = g_param_apinv%iboundflag   ! 0:off, 1:on, 2:transfer
  call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [0]   END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[1-1]## when transformation is used 2018.01.22
 call watchstart(t_watch1)
 if ( iboundflag .eq. 2 ) then   ! 2018.06.25
  ki    = h_model                ! 2018.06.25
  kiref = g_model_ref            ! 2018.06.25
  call TRANSMODELX(ki,   g_param_apinv,1) ! model -> X
  call TRANSMODELX(kiref,g_param_apinv,1) ! model -> X
  m0    = ki                     ! 2018.06.25
  m_ref = kiref                  ! 2018.06.25
 else
  m0    = h_model                ! 2018.06.25
  m_ref = g_model_ref            ! 2018.06.25
 end if
 call watchstop(t_watch1)
if (icheck .eq. 1)  write(*,10) " ### GETNEWMODEL [1]   END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[1-2]## set model1 and model_ref
  ii=0
  do i=1,nmodel                           ! 2018.06.25
   if (iactive(i) .eq. 1) then            ! 2018.06.25
    ii=ii+1                               ! 2018.06.25
    model1(ii)    = m0%logrho_model(i)    ! 2018.06.25
    model_ref(ii) = m_ref%logrho_model(i) ! 2018.06.25
   end if                                 ! 2018.06.25
  end do                                  ! 2018.06.25
  !#
  if ( ii .ne. nmodelactive) then
   write(*,*) "GEGEGE ii",ii,"nmodelactive",nmodelactive
   stop
  end if

!#[2]## assemble C
!#[2-1]# gen JT
  call watchstart(t_watch1)
  call trans_crs2crs(JJ,JT) ! 2018.01.23 m_matrix.f90
  call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-1] END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[2-2]# cal CMJT = CM*JT
  call watchstart(t_watch1)
  call mulreal_crs_crs_crs(CM,JT,CMJT)
  call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-2] END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[2-3]# cal JCMJT = JJ * CMJT
  call watchstart(t_watch1)
  call mulreal_crs_crs_crs(JJ,CMJT,JCMJT)
  call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-3] END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[2-4]# cal C = JCMJT + alpha*CM
  call watchstart(t_watch1)
  ACD = CD
  ACD%val = Cd%val * alpha
  call add_crs_crs_crs(JCMJT,ACD,C) ! C = JCMJT + alpha*CM
  call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [2-4] END!! Time =",t_watch1%time," [min]"!2020.09.18

! write(*,*) "C:"
! call realcrsout(C)

!#[3]## gen X
!#   X  = d_obs - d(M) + J(M - M_ref)
!#[3-1]##
  call watchstart(t_watch1)
  dmodel = model1 - model_ref
  call mul_matcrs_rv(JJ,dmodel,nmodel,JM)
  X = dobs - dcal + JM           ! 2017.07.20
  call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [3-1] END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[4]## cal beta
!# [C]beta = X
!# M_ite+1 = M_ref + Cm*JT*beta
!#[4-1]## solve beta
 call watchstart(t_watch1)
 allocate(beta(ndat))                     ! [ndat,1]
 call solvebeta(C,X,ndat,beta)            ! [ndat,ndat]*[ndat,1] = [ndat,1]
 call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [4-1] END!! Time =",t_watch1%time," [min]"!2020.09.18

!#[4-2]## generate dmodel (crs format)
 call watchstart(t_watch1)
 call mul_matcrs_rv(CMJT,beta,ndat,dmodel)! [nmodel,ndat]*[ndat,1]=[nmodel,1]
 call watchstop(t_watch1)
if (icheck .eq. 1)  write(*,10) " ### GETNEWMODEL [4-2] END!! Time =",t_watch1%time," [min]"!2020.09.18
 ii = 0 ! 2018.06.25
 do i=1,nmodel
  if ( iactive(i) .eq. 1) then ! 2018.06.25
   ii = ii + 1
   h_model%logrho_model(i) = dmodel(ii) + model_ref(ii)  ! [nmodel,1] 2018.06.25
  end if
 end do

!#[5-1]## cut the model
  if ( iboundflag .eq. 1 ) then  ! 2018.01.18
   call BOUNDMODEL(h_model,g_param_apinv)      ! 2018.01.18
  end if

!#[5-2]## transform X to model
  if ( iboundflag .eq. 2 ) then                !  2018.01.22
   call watchstart(t_watch1)
   !# iflag=1: model -> x, iflag = 2: X to model  2018.01.22
   call TRANSMODELX(h_model,g_param_apinv,2)   !  2018.01.22 X -> Model
   call watchstop(t_watch1)
if (icheck .eq. 1) write(*,10) " ### GETNEWMODEL [5-2] END!! Time =",t_watch1%time," [min]"!2020.09.18
  end if                                       !  2018.01.22

 call watchstop(t_watch) ! 2017.12.22
 write(*,10) " ### GETNEWMODEL END!! ###   Time =",t_watch%time," [min]"!2020.09.18

return

10 format(a,f8.4,a) ! 2020.09.18
end
!#################################################### BOUNDMODEL
!# coded on 2018.01.18
subroutine BOUNDMODEL(h_model,g_param_apinv)
use modelpart
use param_apinv
implicit none
type(param_apinversion),intent(in)    :: g_param_apinv
type(model),            intent(inout) :: h_model
integer(4) :: i
real(8)    :: logrho,logrho_upper,logrho_lower
integer(4) :: icount_low,icount_up

!#[1]## set
logrho_upper = g_param_apinv%logrho_upper
logrho_lower = g_param_apinv%logrho_lower
icount_low = 0
icount_up  = 0

!#[2]## cut off
do i=1,h_model%nmodel
 logrho=h_model%logrho_model(i)
 if ( logrho .gt. logrho_upper ) then
  logrho = logrho_upper
  icount_up = icount_up + 1
 else if ( logrho .lt. logrho_lower ) then
  logrho = logrho_lower
  icount_low = icount_low + 1
 else
  goto 100
 end if
  h_model%logrho_model(i) = logrho
 100 continue
end do

write(*,*) "icount_low=",icount_low," icount_up=",icount_up
write(*,*) "### COUNDMODEL END!! ###"

return
end

end module
