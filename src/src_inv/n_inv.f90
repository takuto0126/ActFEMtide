! Coded on 2017.03.07
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
!
implicit none
type(param_forward)   :: g_param     ! see m_param.f90
type(param_source)    :: sparam      ! see m_param.f90
type(param_cond)      :: g_cond      ! see m_param.f90 : goal conductivity
type(param_cond)      :: h_cond      ! see m_param.f90 : initial structure
type(param_inversion) :: g_param_inv ! see m_param_inv.f90
type(mesh)            :: g_mesh      ! see m_mesh_type.f90
type(mesh)            :: h_mesh      ! z file; see m_mesh_type.f90
type(line_info)       :: g_line      ! see m_line_type.f90
type(face_info)       :: g_face      ! see m_line_type.f90
type(modelpara)       :: g_modelpara ! see m_modelpart.f90
type(model)           :: g_model_ref ! see m_modelpart.f90
type(model)           :: h_model     ! see m_modelpart.f90 (2017.05.17)
type(data_vec)        :: g_data      ! observed data  ; see m_param_inv.f90
type(data_vec)        :: h_data      ! calculated data; see m_param_inv.f90
type(real_crs_matrix) :: DQTDQ,CM,CD ! see m_matrix.f90
type(real_crs_matrix) :: PT          ! [nobs,nlin]
!# for forward modelling
type(global_matrix)   :: A             ! see m_iccg_var_takuto.f90
type(real_crs_matrix) :: coeffobs(2,3) ! see m_matrix.f90 ; 1 for edge, 2 for face
integer(4)                            :: ite, itemax, nphys1,nphys2,i,j,nfreq,ip,dofn
integer(4)                            :: nline,ntet, nobs
real(8)                               :: freq, omega
integer(4),allocatable,dimension(:,:) :: table_dof
real(8)                               :: rms
real(8)                               :: alpha ! hyper parameter
complex(8),allocatable,dimension(:)   :: fp,fs
type(obsfiles)        :: files       ! see m_outresp.f90
type(respdata),       allocatable,dimension(:,:) :: resp5    ! m_outresp.f90, 2017.05.18
integer                                          :: access   ! 2017.05.15
type(complex_crs_matrix)                         :: ut       ! 2017.05.16
real(8),              allocatable,dimension(:,:) :: ampbz    ! 2017.05.16
type(real_crs_matrix),allocatable,dimension(:)   :: dampbzdm ! 2017.05.16
type(real_crs_matrix)                            :: JJ ! Jacobian matrix, 2017.05.17
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
!# [7] solve beta in eq. 1 (by BiCG ?)
!# [8] m_k+1 = m_ref + Cm*J*beta

!#[0]## read parameters
  CALL READPARAM(g_param,sparam,g_cond) ! include READCOND for initial model
  CALL READPARAINV(g_param_inv,g_modelpara,g_param,g_data)

!#[1]## read mesh
  CALL READMESH_TOTAL(g_mesh,g_param%g_meshfile)
  if ( access( g_param%z_meshfile, " ") .eq. 0 ) then ! if exist, 2017.05.15
    CALL READMESH_TOTAL(h_mesh,g_param%z_meshfile)
    CALL PREPZSRCOBS(h_mesh,g_param,sparam)        ! see below, include kill h_mesh
  end if
  CALL GENXYZMINMAX(g_mesh,g_param)                ! see below
  CALL READLINE(g_param%g_lineinfofile,g_line)     ! see m_line_type.f90
  CALL READFACE(g_param_inv%g_faceinfofile,g_face) ! see m_face_type.f90

!#[2]## read initial model
  CALL DUPLICATESCALARCOND(g_cond,h_cond)          ! see m_param.f90,  2017.05.14
  h_cond%condfile = g_param_inv%g_initcondfile     ! see m_param_inv.f90
  CALL READCOND(h_cond)                            ! see m_param.f90
  CALL SETNPHYS1INDEX2COND(g_mesh,g_cond)          ! goal cond;    see below
  CALL SETNPHYS1INDEX2COND(g_mesh,h_cond)          ! initial cond; see below
  ! CALL OUTCOND(h_cond,g_mesh,0)     ! this is available on 2017.05.15, see below
  ! write(*,*) "g_cond"
  ! call showcond(g_cond,1)
  ! write(*,*) "h_cond 1"
  ! call showcond(h_cond,1)

!#[3]## generate model space
  CALL genmodelspace(g_mesh,g_modelpara,g_model_ref,g_param,h_cond)! see m_modelpart.f90
  CALL assignmodelrho(h_cond,g_model_ref)                          ! see m_modelpart.f90

!#[4]## cal Cm
  CALL GENDQ(DQTDQ,g_face,g_mesh,g_model_ref) ! see below
  CALL GENCM(DQTDQ,CM)                    ! see below

!#[5]## cal Cd
  CALL GENCD(g_data,CD)                   ! see below

!#[6]## cal Pt: matrix for get data from simulation results
  CALL PREPOBSCOEFF(g_param,g_mesh,g_line,coeffobs) ! for x,y,z component, see below
  CALL DUPLICATE_CRSMAT(coeffobs(2,3),PT) ! PT [nobs, nline] , see m_matrix.f90

!#[7]## Preparation for global stiff matrix
  nline = g_line%nline
  nfreq = g_param%nfreq
  nobs  = g_param%nobs
  ntet  = g_mesh%ntet
  allocate( fp(nline),fs(nline) )
  ip=0 ; dofn=1
  allocate( table_dof(nline,dofn))
  CALL SET_TABLE_DOF(dofn,nline,table_dof,nline)
  CALL set_iccg_var7_dofn(dofn,6,nline,nline,nline,ntet,g_line%n6line,table_dof,A,ip)

!#[4]## allocate respdata and open output files for each observatory
  allocate(resp5(5,nfreq)) ! 2017.05.18
  CALL PREPRESPFILES(g_param,files,resp5,nfreq) ! 2017.05.18

 CALL initializedatavec(nfreq,nobs,h_data) ! see m_param_inv.f90

 allocate(ampbz(nobs,nfreq),dampbzdm(nfreq))

!#[8]## iteration loop start
 h_model = g_model_ref           ! initial model, 2017.05.18
 CALL OUTMODEL(h_model,g_mesh,1) ! see below
 open(21,file="rms.dat")         ! 2017.05.18
 itemax = 10

do ite = 1,itemax

 write(*,*) "ite=",ite
 call model2cond(h_model,h_cond)  ! see m_modelpart.f90
 CALL OUTCOND(h_cond,g_mesh,ite)  ! 2017.05.18
 !call showcond(h_cond,1)

 do i=1,nfreq

  freq = g_param%freq(i) ;   omega=2.d0*pi*freq ; fp = 0.d0
  write(*,*) "i=",i,"freq=",freq

!#[9]## forward calculateion to get d(m)
  CALL forward_inv(A,g_mesh,g_line,nline,fs,omega,sparam,g_param,h_cond,PT,ut)

!#[10]## generate d|bz|dm and bzamp for jacobian
  CALL genjacobian1(nobs,nline,ut,fs,PT,h_model,g_mesh,g_line,omega,ampbz(:,i),dampbzdm(i))

!#[11]## output resp to obs file
    ! calculate bx,by,bz
!   CALL CALOBSRESP(fp,fs,nline,coeffobs(2,1),resp5(1)  ) !bx, see below
!   CALL CALOBSRESP(fp,fs,nline,coeffobs(2,2),resp5(2)  ) !by
    CALL CALOBSRESP(fp,fs,nline,coeffobs(2,3),resp5(3,i)  ) !bz
!   fp = - (0.d0,1.d0)*omega*fp ! E= -i*omega*A
!   fs = - (0.d0,1.d0)*omega*fs !
    ! calculate ex,ey
!   CALL CALOBSRESP(fp,fs,nline,coeffobs(1,1),resp5(4)  ) ! ex
!   CALL CALOBSRESP(fp,fs,nline,coeffobs(1,2),resp5(5)  ) ! ey

 end do ! freq loop end

 !# output responses
 CALL OUTOBSFILESINV(resp5,nfreq,ite,g_param) ! see m_outresp.f90, 2017.05.18

 !# gen dvec and cal rms
 CALL GENDVEC(ampbz,nfreq,nobs,h_data) ! gen data vector
 CALL CALRMS(g_data,h_data,Cd,rms)     ! cal rms, 2017.05.17
 write(21,*) ite, rms                  ! 2017.05.18

 !# conversion check
 if ( rms .lt. 0.2d0 ) then
  write(*,*) "Converged!!with rms =",rms
  stop
 end if

!# h_model is revised for the next iteration
 call genjacobian2(nobs,nfreq,ampbz,dampbzdm,JJ) ! see m_jacobian matrix
 call getnewmodel(JJ,g_model_ref,h_model,g_data,h_data,CM,CD) ! see m_jacobian.f90, 2017.05.17

 !# output new model
 CALL OUTMODEL(h_model,g_mesh,ite+1)     ! see below

end do ! iteration loop

close(21) ! 2017.05.18

end program inversion_csem

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
 rms = rms + dvec(i)**2.d0 / Cd%val(i) ! dvec(i)**2/err(i)**2
end do
 rms = sqrt(rms/dble(ndat))

write(*,*) "rms=",rms

return
99 continue
write(*,*) "GEGEGE! ndat",ndat,".ne. ndat1",ndat1,"or .ne. ndat2",ndat2
stop

end
!############################################## subroutine GENDVEC
subroutine GENDVEC(ampbz,nfreq,nobs,g_data)
use param_inv
implicit none
integer(4),    intent(in)    :: nobs,nfreq
real(8),       intent(in)    :: ampbz(nobs,nfreq)
type(data_vec),intent(inout) :: g_data
integer(4) :: iobs,ii,ifreq,i
real(8)    :: amp,phase

write(*,'(2i5,g15.7)') ((iobs,ifreq,ampbz(iobs,ifreq),ifreq=1,nfreq),iobs=1,nobs)

!#[1]## cal g_data
do ifreq=2,nfreq
 do iobs=1,nobs
  ii = (ifreq-2)*nobs + iobs
  g_data%dvec(ii)  = ampbz(iobs,ifreq)/ampbz(iobs,1)
 end do
end do

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
real(8),allocatable,dimension(:) :: error
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
type(real_crs_matrix),intent(out) :: CM
type(real_crs_matrix) :: crsin,crsout
real(8),allocatable,dimension(:,:) :: x
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

!#[2]##

write(*,*) "nmodel check =",nmodel
allocate(x(nmodel,nmodel))
call solveCM(nmodel,DQTDQ,x)
!do i=1,nmodel
! write(*,*) i,"x(i,1:5)",x(i,1:5)
!end do

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

!#[1]## gen FACE2ELE
  CALL FACE2ELEMENT(g_face) ! cal g_face%face2ele(1:2,1:nface)

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
subroutine OUTMODEL(g_model,g_mesh,ite)
use mesh_type
use modelpart
implicit none
integer(4),      intent(in) :: ite
type(model),     intent(in) :: g_model
type(mesh),      intent(in) :: g_mesh
character(50) :: condfile
integer(4) :: j, nphys2,nmodel,npoi,nlin,ntri,ishift
real(8),allocatable,dimension(:)    :: rho_model, logrho_model
integer(4),allocatable,dimension(:) :: ele2model,index
character(2) :: num

!#[0]## set
nphys2    = g_model%nphys2
nmodel    = g_model%nmodel
allocate(index(nphys2),rho_model(nmodel),ele2model(nphys2))
allocate(logrho_model(nmodel))
index        = g_model%index
rho_model    = g_model%rho_model
logrho_model = g_model%logrho_model
ele2model    = g_model%ele2model
npoi         = g_mesh%npoi
nlin         = g_mesh%nlin
ntri         = g_mesh%ntri

!#[1]## set file name
 write(num,'(i2.2)') ite
 condfile ="model"//num(1:2)//".msh"

!#[2]## output rho
 open(1,file=condfile)

 !# standard info
 CALL MESHOUT(1,g_mesh)

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
  write(1,'(i10,e15.7)') ishift+index(j),10**logrho_model(ele2model(j))
 end do
 write(1,'(a)') "$EndElementData"
close(1)

write(*,*) "### OUTTETMODEL END!! ###"
return
end
!########################################### OUTCOND
! Coded on 2017.05.18
subroutine OUTCOND(g_cond,g_mesh,ite)
use mesh_type
use param
implicit none
integer(4),      intent(in) :: ite
type(param_cond),intent(in) :: g_cond
type(mesh),      intent(in) :: g_mesh
character(50) :: condfile
integer(4) :: j, nphys2,nmodel,npoi,nlin,ntri,ishift,nphys1
real(8),allocatable,dimension(:)    :: rho,sigma
integer(4),allocatable,dimension(:) :: index
character(2) :: num

!#[0]## set
nphys1    = g_cond%nphys1
nphys2    = g_cond%nphys2
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
 condfile = "cond"//num(1:2)//".msh"
 open(1,file=condfile)

 !# standard info
! CALL MESHOUT(1,g_mesh)

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

write(*,*) "### PREPOBSCOEFF END!! ###"
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
	write(*,*) "### SET TABLE_DOF END!! ###"
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
do j=1,nfreq ! 2017.05.18
do i=1,5
 CALL ALLOCATERESPDATA(nobs,resp5(i,j))
end do
end do

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
!# Copied from ../solver/n_ebfem_bxyz.f90 on 2017.05.15
!# Coded on 2017.02.21
subroutine PREPZSRCOBS(h_mesh,g_param,s_param)
use param
use mesh_type
use triangle
implicit none
type(mesh),         intent(inout) :: h_mesh  ! deallocated at the end 2017.05.15
type(param_forward),intent(inout) :: g_param
type(param_source), intent(inout) :: s_param
type(grid_list_type) :: glist
integer(4) :: nobs,nx,ny
real(8),allocatable,dimension(:,:) :: xyzobs,xyz
integer(4),allocatable,dimension(:,:) :: n3k
real(8),allocatable,dimension(:) :: znew
real(8) :: a3(3)
integer(4) :: iele,n1,n2,n3,j,ntri
real(8) :: xyzminmax(6),zorigin,xs1(3),xs2(3)

!#[0]## cal xyzminmax of h_mesh
  CALL GENXYZMINMAX(h_mesh,g_param)

!#[1]## set
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
!#  start point
    call findtriwithgrid(h_mesh,glist,xs1(1:2),iele,a3)
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs1(3) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs1(3)
    !
    write(*,*) "xs1(1:2)=",xs1(1:2)
    write(*,*) "z",s_param%xs1(3),"->",xs1(3),"[km]"

!#  end point
    call findtriwithgrid(h_mesh,glist,xs2(1:2),iele,a3)
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs2(3) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs2(3)
    !
    write(*,*) "xs2(1:2)=",xs2(1:2)
    write(*,*) "z",s_param%xs2(3),"->",xs2(3),"[km]"

!#[4]## set znew to xyz_r
    g_param%xyzobs(3,1:nobs) = znew(1:nobs)
    s_param%xs1(3) = xs1(3)
    s_param%xs2(3) = xs2(3)

!#[5]## kill mesh for memory 2017.05.15
    call killmesh(h_mesh) ! see m_mesh_type.f90

write(*,*) "### PREPZSRCOBS END ###"
return
end
