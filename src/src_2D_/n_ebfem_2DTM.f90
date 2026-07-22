! coded on 2020.10.28
program ebfem_2DTM
use mesh_type
use param
use matrix
use line_type  ! see m_line_type.f90
use iccg_var_takuto ! see ../solver/m_iccg_var_takuto.f90
use outresp
use surface_type ! see src_2D/m_surface_type.f90
use obs_type     ! see m_obs_type.f90 added on 2016.10.27
use constants, only:pi,dmu
use face_type
implicit none
type(param_forward)                         :: g_param
type(param_source)                          :: sparam
type(param_cond)                            :: g_cond     ! see m_param.f90
type(mesh)                                  :: g_mesh     ! see m_mesh_type.f90
type(mesh)                                  :: h_mesh     ! topography file
type(line_info)                             :: g_line     ! see m_line_type.f90
integer(4)                                  :: nline,nsr,nfreq ! 2018.02.22
integer(4)                                  :: i,j,iresfile,dofn,i_surface
real(8)                                     :: omega, freq ! 2018.02.22
integer(4)                                  :: ixyflag     ! 2017.10.11
type(surface),             dimension(6)     :: g_surface
type(face_info)                             :: g_face
integer(4)                                  :: nface,node,ntri, doftot,ntet
integer(4),allocatable,dimension(:,:)       :: n4
integer(4)                                  :: ip=0,iele
complex(8)                                  :: e(2),bx
real(8)                                     :: rhoa,pha,y,z

!#[0]## read param
 CALL READPARAM(g_param,sparam,g_cond) ! get parameter info and conductivity info

!#[1]## Mesh READ
 CALL READMESH_TOTAL(g_mesh,g_param%g_meshfile) ! 3D mesh
 CALL READMESH_TOTAL(h_mesh,g_param%z_meshfile) ! topography file
 CALL GENXYZMINMAX(g_mesh,g_param)             ! generate xyzminmax 2017.10.12
 if (g_cond%condflag .eq. 1) then !  "1" means conductivity file is given
   g_cond%ntet   = g_mesh%ntet
   g_cond%nphys1 = g_mesh%ntet - g_cond%nphys2
 end if
 CALL SETNPHYS1INDEX2COND(g_mesh,g_cond) ! added on 2017.05.31; see below

!#[2]## Line information
  CALL READLINE(g_param%g_lineinfofile,g_line) ! see m_line_type.f90

!#[3]## Make 6 surface information
  ntet    = g_mesh%ntet     ! 2107.09.08
  node    = g_mesh%node     ! 2017.09.08
  ixyflag = g_param%ixyflag ! 2018.02.22
  nfreq   = g_param%nfreq   ! 2018.02.22
  n4      = g_mesh%n4       ! allocate n4(ntet,4) here 2020.10.04
  CALL MKFACE(  g_face,node,ntet,4, n4) ! make g_line
  CALL MKN4FACE(g_face,node,ntet,   n4) ! make g_line%n6line
  CALL FACE2ELEMENT(g_face)

!#[4]## Extract Faceinfo of 3D mesh for TM mode calculation
!        Face2        ^ y
!        ----         |
! Face3 |    | Face5  |
!        ----         ----> x
!        Face4
  CALL EXTRACT6SURFACES(g_mesh,g_line,g_face,g_surface) ! 2020.10.29 m_surface_type.f90

  CALL FINDBOUNDARYLINE(g_mesh,g_surface)  ! 2021.05.25 m_surface_type.f90

  CALL COND3DTO2D(g_mesh,g_surface,g_cond) ! see m_surfae_type.f90 2021.10.14

!#[5]##
freq  = g_param%freq(1)
freq = 0.1d0
omega = g_param%freq(1)*2.*pi
ip=0

!#[6]## allocate table_dof, A
CALL PREPAOFSURFACE(g_surface(2:5),4,ip) ! see forward_2DTM.f90


do j=2,5 ! except 1 (top) and 6 (bottom) surface

!#[7]## forward_2DTM to obtain g_surface%bs [mV/km*km]=[mV]
 call forward_2DTM(g_surface(j),freq,g_param,g_cond,ip)

 do i=1,g_surface(j)%nline
!  write(*,*)"i",i,"y1 z1",g_surface(j)%x1x2_face(1:2,g_surface(j)%line(1,i))
!  write(*,*)"i",i,"y2 z2",g_surface(j)%x1x2_face(1:2,g_surface(j)%line(2,i))
!write(*,'(a,i5,a,2g15.7,2i5)')"i ",i," bs", g_surface(j)%bs(i),g_surface(j)%line(1:2,i)
 end do
 y=0.d0; z=0.d0
 call searchtri(g_surface(j),y,z,iele) ! in src_2D/m_surface_type.f90
 call E_ele(g_surface(j),iele,y,z,e)   ! in src_2D/m_surface_type.f90
 write(*,*) "E",E
 call B_ele(g_surface(j),omega,iele,y,z,bx) ! in src_2D/m_surface_type.f90
 write(*,*) "bx",bx
 call rhoap_tri(e,bx,omega,rhoa,pha)
 write(*,*) "j",j,"rhoa",rhoa," [Ohm.m] pha",pha,"[deg]"

 call outsurface_vec(g_surface(j),j) ! see below

end do


end program

!###################################################################
! modified for spherical on 2016.11.20
! iflag = 0 for xyz
! iflag = 1 for xyzspherical
subroutine GENXYZMINMAX(em_mesh,g_param)
use param ! 2016.11.20
use mesh_type
implicit none
type(mesh),         intent(inout) :: em_mesh
type(param_forward),intent(inout) :: g_param
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
real(8) :: xyzminmax(6)
real(8), allocatable,dimension(:,:) :: xyz
integer(4) :: i

allocate(xyz(3,em_mesh%node))
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

write(*,*) "xmin,xmax",xmin,xmax ! 2021.10.13
write(*,*) "ymin,ymax",ymin,ymax ! 2021.10.13
write(*,*) "zmin,zmax",zmin,zmax ! 2021.10.13

xyzminmax(1:6)=(/xmin,xmax,ymin,ymax,zmin,zmax/)

!# set output
g_param%xyzminmax = xyzminmax
em_mesh%xyzminmax = xyzminmax   ! 2020.10.31

write(*,*) "### GENXYZMINMAX END!! ###"
return
end
!############################################## subroutine setNPHYS1INDEX2COND
! Modified on 2017.05.31
! Coded on 2017.05.12
subroutine SETNPHYS1INDEX2COND(g_mesh,g_cond)
use param
use mesh_type
implicit none
type(mesh),      intent(in)    :: g_mesh
type(param_cond),intent(inout) :: g_cond
integer(4) :: nphys1,nphys2,i,ntet,icount
integer(4),allocatable,dimension(:,:) :: n4flag

!#[0]## set
  ntet   = g_mesh%ntet
  allocate(n4flag(ntet,2))
  n4flag = g_mesh%n4flag

!#[1]## calculate nphys1 and nphys2
  nphys2 = 0
  do i=1,ntet
   if ( n4flag(i,1) .ge. 2 ) nphys2 = nphys2 + 1 ! count land elements 2017.09.29
  end do

  if (g_cond%condflag .eq. 1 ) then ! check when cond file is given
    if ( g_cond%nphys2 .ne. nphys2) then
     write(*,*) "GEGEGE nphys2=",nphys2,"g_cond%nphys2=",g_cond%nphys2
     stop
    end if
  end if

!#[2]## set nphys1 and nphys2
  nphys1 = ntet - nphys2
  g_cond%nphys1 = nphys1
  g_cond%nphys2 = nphys2 ! # of elements in 2nd physical volume (land)
  g_cond%ntet   = ntet
  write(*,*) "nphys1=",nphys1,"nphys2=",nphys2,"ntet=",g_mesh%ntet

!#[3]## prepare the rho, sigma in the case where file is not given
  ! If the file is given the following is allocated in READCOND in m_param.f90
  if (g_cond%condflag .eq. 0 ) then ! file is not given
   allocate( g_cond%sigma(nphys2) )
   allocate( g_cond%rho(  nphys2) )
   allocate( g_cond%index(nphys2) )
   do i=1,nphys2                          ! 2017.09.29
    if ( n4flag(nphys1+i,1) .le. 1 ) then ! 2017.09.29
     write(*,*) "GEGEGE! i=",i,"n4flag(nphys1+i,1)=",n4flag(nphys1+i,1),"nphys1=",nphys1
     stop                                 ! 2017.09.29
    end if                                ! 2017.09.29
    g_cond%sigma(i) = g_cond%sigma_land(n4flag(nphys1+i,1)-1) ! 2017.09.29
    g_cond%rho(i)   = 1.d0/g_cond%sigma(i)
   end do
  end if

!#[3]## set index
  do i=1,nphys2
   g_cond%index(i) = nphys1 +i
  end do

write(*,*) "### SETNPYS1INDEX2COND ###"
return
end

!======================================================  outsurface
! on 2020.10.31
subroutine outsurface_vec(g_surface,jin)
use mesh_type
use surface_type
implicit none
type(surface),intent(in) :: g_surface
integer(4),intent(in) :: jin
type(mesh)   :: h_mesh
integer(4)   :: iout,i,j
character(1) :: ci
real(8)      :: x3_c(3),e_real(3)
character(50) :: mshfile
complex(8)   :: e(2)
character(1) :: num

!# MESHOUT
iout=1

!# point
h_mesh%node = g_surface%node
allocate(h_mesh%xyz(3,h_mesh%node))
h_mesh%xyz=0.d0
if (g_surface%facetype .eq. "xy") h_mesh%xyz(1:2,:)=g_surface%x1x2_face
if (g_surface%facetype .eq. "yz") h_mesh%xyz(2:3,:)=g_surface%x1x2_face
if (g_surface%facetype .eq. "xz") then
h_mesh%xyz(1,:)=g_surface%x1x2_face(1,:)
h_mesh%xyz(3,:)=g_surface%x1x2_face(2,:)
end if

!# point, line
h_mesh%npoi = 0
h_mesh%nlin = 0

!# n3
h_mesh%ntri = g_surface%ntri
allocate(h_mesh%n3(h_mesh%ntri,3))
h_mesh%n3 = g_surface%n3
allocate(h_mesh%n3flag(h_mesh%ntri,2))
h_mesh%n3flag(:,:)=0

write(num,'(i1)') jin
mshfile="surface_vec"//num//".msh"
open(iout,file=mshfile)
call MESHOUT(iout,h_mesh)

!#[4]## create vector fields
write(1,'(a)') "$ElementData"
write(1,'(a)') "1"
write(1,'(a)') '"A vector view"'
write(1,'(a)') "1"
write(1,'(a)') "0.0"
write(1,'(a)') "3"
write(1,'(a)') "0"
write(1,'(a)') "3"
write(1,*) g_surface%ntri

!#[5]
do i=1,g_surface%ntri
 x3_c=0.d0
 do j=1,3
  x3_c(1:2) = x3_c(1:2) + g_surface%x1x2_face(1:2,g_surface%n3(i,j))/3.
 end do
! write(*,*) "x3_c(1:2)",x3_c(1:2)
 call E_ele(g_surface,i,x3_c(1),x3_c(2),e)
  e_real=0.d0

if ( g_surface%facetype .eq. "xy" ) then
  e_real(1:2)=(/real(e(1)),real(e(2)) /)
  end if

if ( g_surface%facetype .eq. "yz" ) then
 e_real(2:3)=(/real(e(1)),real(e(2)) /)
 end if

 if ( g_surface%facetype .eq. "xz" ) then
 e_real(1:3)=(/real(e(1)),0.d0,real(e(2)) /)
 end if

write(1,*) i, e_real(1:3)

end do

write(1,'(a)') "$EndElementData"

close(iout)

return
end

