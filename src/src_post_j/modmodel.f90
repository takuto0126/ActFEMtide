! 2018.03.23
program modmodel
use mesh_type
use param
use modelpart
use param_modmodel
use cond_type      ! ../src_mesh/m_cond_type.f90
implicit none
type(mesh)       :: h_mesh
type(param_cond) :: h_cond
type(model)      :: g_model
type(param_cmod) :: m_param
type(param_mod)  :: g_param_mod

!#[0]## read param
 CALL readcmodparam(m_param)
 CALL readparammod(g_param_mod)

!#[1]## read mesh
 CALL READMESH_TOTAL(h_mesh,m_param%mshfile)
 CALL READCONDBASE(h_cond,g_param_mod)
 CALL SETNPHYS1INDEX2COND(h_mesh,h_cond)          ! ref     cond; see below


!#[2]## read model
 CALL ASSIGNCOND(h_mesh,m_param,h_cond) ! modify h_cond

!#[3]##
 CALL OUTCOND(m_param%outcondfile,h_cond,h_mesh)

end program
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
!#################################################### READREFINITCOND
!# Coded on 2018.03.18
subroutine READCONDBASE(h_cond,g_param_mod) ! 2018.03.18
use modelpart
use param
use param_modmodel
implicit none
type(param_cond), intent(inout) :: h_cond
type(param_mod),  intent(in)    :: g_param_mod
character(50)                   :: modelfile,connectfile
integer(4)                      :: icondflag

!#[0]##
icondflag = g_param_mod%icondflag

!#[ref cond]
if ( icondflag .eq. 1 ) then ! condfile        2018.03.18
 h_cond%condfile = g_param_mod%condfile      ! 2017.07.19
 CALL READCOND(h_cond)                            ! 2017.07.19
elseif ( icondflag .eq. 2 ) then ! modelfile
 modelfile   = g_param_mod%modelfile
 connectfile = g_param_mod%modelconn
 CALL READMODEL2COND(h_cond,connectfile,modelfile)
end if

return
end
!########################################### OUTCOND
! Output folder is changed on 2018.07.23
! Coded on 2017.05.18
subroutine OUTCOND(condfile,g_cond,g_mesh) ! 2017.07.25
use mesh_type
use param
implicit none
type(param_cond),      intent(in)   :: g_cond
type(mesh),            intent(in)   :: g_mesh
character(50),         intent(in)   :: condfile ! 2018.03.23
integer(4)                          :: j, nphys2,nmodel,npoi,nlin,ntri
integer(4)                          :: ishift,nphys1
real(8),   allocatable,dimension(:) :: rho,sigma
integer(4),allocatable,dimension(:) :: index
character(2) :: num

!#[0]## set
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

!#[1]## output rho
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
  write(1,*) ishift+index(j),rho(j)
 end do
 write(1,'(a)') "$EndElementData"
close(1)

write(*,*) "### OUTCOND END!! ###"
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
!####################################  ASSIGNCOND
subroutine ASSIGNCOND(h_mesh,m_param,r_cond)
use mesh_type
use cond_type
use param
implicit none
type(mesh),            intent(in)     :: h_mesh
type(param_cmod),      intent(in)     :: m_param
type(param_cond),      intent(inout)  :: r_cond
integer(4)                            :: node,ntet,nphys2
integer(4),allocatable,dimension(:,:) :: n4
real(8),   allocatable,dimension(:,:) :: xyz
integer(4)                            :: i,j,k,l,ll
real(8)                               :: cxyz(3)
real(8),               dimension(2)   :: x2,y2,z2
real(8)                               :: ratio,rho_o,rho_n
real(8),  allocatable,dimension(:)    :: rho
integer(4),allocatable,dimension(:)   :: index


!#[1]# set
 ntet   = h_mesh%ntet
 node   = h_mesh%node
 nphys2 = r_cond%nphys2
 allocate(xyz(3,node),n4(ntet,4))
 xyz    = h_mesh%xyz
 n4     = h_mesh%n4
 allocate(rho(nphys2),index(nphys2))
 rho    = r_cond%rho
 index  = r_cond%index

 ![2]# search cuboid for this tetrahedron added on 2017.02.22
 do l=1,nphys2
  ll=index(l)
  cxyz(:)=(xyz(:,n4(ll,1))+xyz(:,n4(ll,2))+xyz(:,n4(ll,3))+xyz(:,n4(ll,4)))/4.d0
  do i=1,m_param%ncuboid
   x2=m_param%g_cuboid(i)%xminmax
   y2=m_param%g_cuboid(i)%yminmax
   z2=m_param%g_cuboid(i)%zminmax
   if (z2(1) .lt. cxyz(3) .and. cxyz(3) .lt. z2(2) .and. &
   &   y2(1) .lt. cxyz(2) .and. cxyz(2) .lt. y2(2) .and. &
   &   x2(1) .lt. cxyz(1) .and. cxyz(1) .lt. x2(2) ) then
     ratio = m_param%g_cuboid(i)%ratio
     rho_o = log10(rho(l))
     rho_n = log10(m_param%g_cuboid(i)%rho)
     rho(l)= 10.d0**(ratio*rho_n + (1.-ratio)*rho_o)  ![Ohm.m]
!     write(*,*) "rho(l)=",rho(l)
!     write(*,*) "x2=",x2,"y2=",y2,"z2=",z2
     goto 110
   end if
  end do
  110 continue
 end do

!#[3]# set output
 r_cond%nphys2 = nphys2
 r_cond%rho    = rho
 write(*,*) "### ASSIGNCOND END!! ###"

 return
end
