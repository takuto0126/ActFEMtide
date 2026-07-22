program calprofile
use post_param
use mesh_type
use modelpart ! 2018.10.24
use param
implicit none
type(mesh)       :: em_mesh
type(param_cond) :: m0_cond
type(param_cond) :: dm_cond
type(profile)    :: m0_profile
type(profile)    :: dm_profile
type(profile)    :: dmreal_profile       ! 2018.10.24
type(profile)    :: m1_profile
type(profile)    :: m2_profile
type(param_post) :: g_parampost
character(50)    :: mshfile
character(50)    :: m0condfile,dmcondfile
character(50)    :: connectfile          ! 2018.10.24
integer(4)       :: ifileflag            ! 2018.10.24
integer(4)       :: nphys2,nmodelfile    ! 2018.10.24
integer(4)       :: itype_process        ! 2018.10.24

!#[0]## set
 CALL postparamread(g_parampost)
 mshfile    = g_parampost%mshfile
 m0condfile = g_parampost%m0condfile
 dmcondfile = g_parampost%dmcondfile
 ifileflag  = g_parampost%ifileflag ! 2018.01.06
 itype_process = g_parampost%itype_process ! 2018.10.24

!#[1]## read msh and cond
 write(*,*) "mshfile",mshfile
 CALL READMESH_TOTAL(em_mesh,mshfile)

!#[2]## read resistivity structure
 !#[2-1]## cond file case
 if ( ifileflag .eq. 1 ) then ! 2018.10.24
  m0_cond%condfile = m0condfile
  dm_cond%condfile = dmcondfile
  CALL READCOND(m0_cond)                            ! 2017.07.19
  CALL SETNPHYS1INDEX2COND(em_mesh,m0_cond)         ! ref     cond; see below
  CALL READCOND(dm_cond)                            ! 2017.07.19
  CALL SETNPHYS1INDEX2COND(em_mesh,dm_cond)         ! ref     cond; see below

 !#[2-2]## model file case
 else if ( ifileflag .eq. 2 ) then ! 2018.10.24
  connectfile = g_parampost%connectfile ! 2018.10.24
  call readmodel2cond(m0_cond,connectfile,g_parampost%m0modelfile) ! 2018.10.24
  call readmodel2cond(dm_cond,connectfile,g_parampost%dmmodelfile) ! 2018.10.24
  call SETNPHYS1INDEX2COND(em_mesh,m0_cond) ! 2018.10.24
  call SETNPHYS1INDEX2COND(em_mesh,dm_cond) ! 2018.10.24
 end if


!#[2]## calculate profile
 CALL CALVALUE(em_mesh,g_parampost,m0_cond,m0_profile)
 CALL CALVALUE(em_mesh,g_parampost,dm_cond,dm_profile)

!#[3]## m1, m2, dm, m0
  CALL PREPOSTPROFIlE(m0_profile,dm_profile,m1_profile,m2_profile,itype_process)!2018.10.24

!#[3]## output profile
 CALL OUTPUTPROFILE(m0_profile,g_parampost%outputfile(1))
 CALL OUTPUTPROFILE(dm_profile,g_parampost%outputfile(2))
 CALL OUTPUTPROFILE(m1_profile,g_parampost%outputfile(3))
 CALL OUTPUTPROFILE(m2_profile,g_parampost%outputfile(4))
 dmreal_profile=m1_profile
 dmreal_profile%value(:) = m2_profile%value(:) - m1_profile%value(:)
 CALL OUTPUTPROFILE(dmreal_profile,g_parampost%outputfile(5))


end program

!##################################
subroutine  PREPOSTPROFIlE(m0_profile,dm_profile,m1_profile,m2_profile,itype_process)
use post_param
implicit none
integer(4),   intent(in)    :: itype_process ! 1:m0(average)+dm, 2: m1 & dm 2018.10.24
type(profile),intent(inout) :: m0_profile    ! -> average
type(profile),intent(inout) :: dm_profile    ! -> after - before 2018.10.24
type(profile),intent(out)   :: m1_profile    ! -> before
type(profile),intent(out)   :: m2_profile    ! -> after
integer(4) :: nlayer, i

!#[1]## set
 nlayer     = m0_profile%nlayer

!#[2]##
 if ( itype_process .eq. 1) then ! 2018.10.24
  m1_profile = m0_profile
  m2_profile = m0_profile
  do i=1,nlayer
   m1_profile%value(i) = 10.**(log10(m0_profile%value(i))- 1/2.*dm_profile%value(i))
   m2_profile%value(i) = 10.**(log10(m0_profile%value(i))+ 1/2.*dm_profile%value(i))
  end do
 else if (itype_process .eq. 2 ) then ! 2018.10.24
  m1_profile = m0_profile
  m2_profile = m0_profile
  do i=1,nlayer
   m2_profile%value(i) = 10.**(log10(m0_profile%value(i)) +       dm_profile%value(i))
   m0_profile%value(i) = 10.**(log10(m1_profile%value(i)) + 1./2.*dm_profile%value(i))
  end do
 else
  write(*,*) "GEGEGE"
  stop
 end if

return
end
!###################################
subroutine OUTPUTPROFILE(g_profile,outputfile)
use post_param
implicit none
type(profile),   intent(in) :: g_profile
character(50),   intent(in) :: outputfile
integer(4)    :: nlayer
real(8),allocatable,dimension(:,:) :: bottoplayer
real(8),allocatable,dimension(:)   :: value
integer(4) :: i

!#[0]# set
nlayer     = g_profile%nlayer
allocate(bottoplayer(2,nlayer))
bottoplayer = g_profile%bottoplayer
allocate(value(nlayer))
value      = g_profile%value

!#[1]# out
open(1,file=outputfile)
 do i=1,nlayer
  write(1,*) value(i),bottoplayer(1,i)
  write(1,*) value(i),bottoplayer(2,i)
 end do
close(1)

return
end
!###################################
subroutine CALVALUE(em_mesh,g_parampost,g_cond,g_profile)
use post_param
use mesh_type
use param
implicit none
type(param_post),intent(in) :: g_parampost
type(profile),   intent(out) :: g_profile
type(mesh),      intent(in)  :: em_mesh
type(param_cond),      intent(in)  :: g_cond
integer(4)                         :: ntet,nphys2,nphys1,node
real(8),allocatable,dimension(:)   :: zlayer
real(8),allocatable,dimensioN(:,:) :: bottoplayer
real(8),allocatable,dimension(:)   :: value
integer(4)                         ::  i,j,ii,k,nlayer
real(8),allocatable,dimension(:,:) :: xyz
integer(4),allocatable,dimension(:,:) :: n4
integer(4),allocatable,dimension(:)   :: index
real(8) :: xo,yo,zbot,ztop
real(8) :: elm_xyz(3,4)
real(8) :: xyz_c(3)
real(8) :: radius = 0.1d0 ! [km]

!#[0]## set
ntet   = em_mesh%ntet
node   = em_mesh%node
nphys1 = g_cond%nphys1
nphys2 = g_cond%nphys2
nlayer = g_parampost%nlayer
allocate( zlayer(nlayer+1) )
zlayer = g_parampost%zlayer
allocate( xyz(3,node),n4(ntet,4) )
xyz    = em_mesh%xyz
n4     = em_mesh%n4
xo     = g_parampost%xyprof(1)
yo     = g_parampost%xyprof(2)
allocate(index(nphys2))
index  = g_cond%index

!#[1]## set topbotlayer
 allocate( bottoplayer(2,nlayer) )
 do i=1,nlayer
  bottoplayer(1,i)=zlayer(i)   ! bot
  bottoplayer(2,i)=zlayer(i+1) ! top
 end do

!#[2]## cal profile values
  allocate(value(nlayer))
  do i=1,nlayer
   write(*,*) "i-th",i,"layer"
   zbot = bottoplayer(1,i)
   ztop = bottoplayer(2,i)
   write(*,*) "zbot",zbot,"ztop",ztop
   do j = 1,nphys2
    ii = index(j)
    do k=1,4
     elm_xyz(:,k) = xyz(:,n4(ii,k))
    end do
     xyz_c(1) =sum(elm_xyz(1,:))/4.
     xyz_c(2) =sum(elm_xyz(2,:))/4.
     xyz_c(3) =sum(elm_xyz(3,:))/4.
    if ( abs(xyz_c(1) - xo) .lt. radius ) then
    if ( abs(xyz_c(2) - yo) .lt. radius ) then
    if ( zbot .lt. xyz_c(3) .and. xyz_c(3) .lt. ztop) then
     value(i) = g_cond%rho(j)
     goto 50
    end if
    end if
    end if
   end do
   write(*,*) "GEGEGE no corresponding.."
   stop
   50 continue
  end do

!#[3]## out
 g_profile%nlayer = nlayer
 allocate( g_profile%bottoplayer(2,nlayer) )
 g_profile%bottoplayer = bottoplayer
 allocate( g_profile%value(nlayer+1))
 g_profile%value = value

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
!##################################################### READMODEL2COND
!# copied from n_inv_ap.f90 on 2018.10.24
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


