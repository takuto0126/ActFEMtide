!# coded on Dec. 14, 2015
!# This program create kandamodel.msh
!# from the kand inversion model
!# Input files
!#  aso3d_wt_model.01  : grid and conductivity file
!#
program gridmdl2msh
use param
use mesh_type
use cond_type ! see m_cond_type.f90
implicit none
type(param_cmod)     :: m_param ! see m_cond_type.f90
character(50)        :: modelfile="../from_Kanda/aso3d_wt4_model.01"
character(50)        :: outmshfile="./kandamodel.msh"
character(50)        :: gridtopofile="kandatopo.msh"
integer(4)           :: nx,ny,nz
real(8),allocatable,dimension(:) :: dx,dy,dz
real(8),allocatable,dimension(:) :: xgrd,ygrd,zgrd
real(8),allocatable,dimension(:) :: modelrho
integer(4)           :: nobs_amt
integer(4),parameter :: nobsmax_amt=100
real(8)              :: xyobs(2,nobsmax_amt)
!### observation points for Kanda model inversion ####
type(mesh)           :: h_mesh
character(50)        :: obsfile_xy="../from_Kanda/site.xy"
character(50)        :: obsfile_loc="../from_Kanda/site.loc"
real(8),parameter    :: ztop_grd=1590.d0 ! [m] actual altitude of top of the Kanda grid
real(8)              :: clatlon(2)  ! lat and lon [deg] of the center point of Kanda grid
!### read 3d tetrahedral mesh for nakadake #############
character(50) :: inmshfile
integer(4),  allocatable, dimension(:,:) :: n4g,n3g,n2g
real (8),    allocatable, dimension(:,:) :: xyzg  ! [km]
integer(4),  allocatable, dimension(:)   :: n1g
integer(4) :: IPL, nodeg,ntetg,ntrig,nling,npoig,nsele
!### output 3d tetrahedral rho model file ##############
real(8),dimension(:),allocatable :: tetrho
!character(50) :: outtetmshfile="nakadake3d_init.msh"

!#[]##
call readcmodparam(m_param) ! see m_cond_type.f90
 inmshfile = m_param%mshfile

!#[0]# read observation points
call readobs(nobsmax_amt,nobs_amt,xyobs,obsfile_xy,obsfile_loc,clatlon)

!#[1]# count and read model file
CALL countmodel(modelfile,nx,ny,nz)
 allocate(dx(nx),dy(ny),dz(nz),modelrho(nx*ny*nz))
 allocate(xgrd(nx+1),ygrd(ny+1),zgrd(nz+1))
CALL readmodel(modelfile,dx,dy,dz,nx,ny,nz,modelrho)

!#[2]# gen xyz grd coordinate from dx, dy, and dz
CALL gengrd(dx,dy,dz,nx,ny,nz,xgrd,ygrd,zgrd)

!#[3]# convert the axis to z: upward, x: eastward, y: northward
CALL adjustcoord(nx,ny,nz,xgrd,ygrd,zgrd,clatlon,ztop_grd) ! change the unit to [km]

!#[4]# output conductivity model to .msh file, keeping hexahedral grids
!CALL outmodelmsh(modelrho,nx,ny,nz,xgrd,ygrd,zgrd,outmshfile,xyobs(:,1:nobs_amt),nobs_amt)

!#[5]# output topography of hexahedral mesh
!CALL outmodeltopo(modelrho,nx,ny,nz,xgrd,ygrd,zgrd,gridtopofile)

!#[6]# read 3d msh tetrahedral mesh
CALL READMESH_TOTAL(h_mesh,inmshfile)

!#[7]# assign conductivity to each tetrahedron from Kanda model
     allocate(tetrho(h_mesh%ntet))
CALL ASSIGNCOND(h_mesh,xgrd,ygrd,zgrd,nx,ny,nz,modelrho,tetrho,m_param)

!#[8]# output the tetrahedral model
CALL OUTTETMODEL(m_param,tetrho,h_mesh)

end program gridmdl2msh
!########################################### outmodeltopo
subroutine outmodeltopo(modelrho,nx,ny,nz,xgrd,ygrd,zgrd,gridtopofile)
implicit none
character(50),intent(in) :: gridtopofile
integer(4),intent(in) :: nx,ny,nz
real(8),intent(in) :: xgrd(nx+1),ygrd(ny+1),zgrd(nz+1),modelrho(nx*ny*nz)
integer(4) :: i,j,k,n1,n2,n3,n4
real(8) :: topo(nx,ny),rho
!#[1]# calculate topography
do j=1,ny
 do i=1,nx
  do k=1,nz
   rho=modelrho((k-1)*(nx*ny)+ny*(i-1)+j)
   if ( rho .lt. 1.d+8 ) then
    topo(i,j)=zgrd(k)
    goto 100
   end if
  end do
  100 continue
 end do
end do

!#[2]# output file for gmsh
open(1,file=gridtopofile)
write(1,'(a)') "$MeshFormat"
write(1,'(a)')  "2.2 0 8"
write(1,'(a)')  "$EndMeshFormat"
write(1,'(a)')  "$Nodes"
write(1,*) (nx+1)*(ny+1)
do i=1,nx+1
 do j=1,ny+1
  write(1,*) (ny+1)*(i-1)+j,xgrd(i),ygrd(j),"0.0"
 end do
end do
write(1,'(a)')  "$EndNodes"
write(1,'(a)')  "$Elements"
write(1,*) nx*ny
 do i=1,nx
  do j=1,ny
   n1=(ny+1)*(i-1)+j
   n2=n1+1
   n3=n1+(ny+1)+1
   n4=n1+(ny+1)
   write(1,*) ny*(i-1)+j," 3 2 0 1 ",n1,n2,n3,n4
  end do
 end do
write(1,'(a)')  "$EndElements"
write(1,'(a)') "$ElementData"
write(1,'(a)') "1"
write(1,'(a)') "Topography"
write(1,'(a)') "1"
write(1,'(a)') "0.0"
write(1,'(a)') "3"
write(1,'(a)') "0"
write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
write(1,*) nx*ny
do i=1,nx
 do j=1,ny
  write(1,*) (i-1)*ny+j,topo(i,j)
 end do
end do
write(1,'(a)') "$EndElementData"
close(1)
write(*,*) "### OUTMODELTOPO END!! ###"
return
end
!########################################### OUTTETMODEL
subroutine OUTTETMODEL(m_param,tetrho,h_mesh)
use mesh_type
use cond_type
implicit none
type(param_cmod),intent(in) :: m_param
type(mesh),intent(in)       :: h_mesh
real(8),intent(in) :: tetrho(h_mesh%ntet)
character(50) :: condfile
integer(4) :: ntetg,nling,npoig,ntrig
integer(4) :: j, ishift,icount
real(8) :: ratio

!#[0]## set
ntetg = h_mesh%ntet
nling = h_mesh%nlin
ntrig = h_mesh%ntri
npoig = h_mesh%npoi
ratio = m_param%ratio
condfile = m_param%outcondfile

!#[1]## output rho
 open(1,file=condfile)
 icount=0
 do j=1,ntetg
  if (h_mesh%n4flag(j,1) .eq. 2) icount=icount+1
 end do
 write(1,'(a)') "$MeshFormat"    ! 2017.09.11
 write(1,'(a)') "2.2 0 8"        ! 2017.09.11
 write(1,'(a)') "$EndMeshFormat" ! 2017.09.11
 write(1,'(a)') "$ElementData"
 write(1,'(a)') "1"
 write(1,'(a)') '"A rho model view"'
 write(1,'(a)') "1"
 write(1,'(a)') "0.0"
 write(1,'(a)') "3"
 write(1,'(a)') "0"
 write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
 write(1,'(i10)') icount
  ishift=nling + ntrig + npoig
 do j=1,ntetg
!    only the physical volume of 2 (land region)
  if ( h_mesh%n4flag(j,1) .eq. 2 ) write(1,*) ishift+j,tetrho(j)*ratio
 end do
write(1,'(a)') "$EndElementData"
close(1)

write(*,*) "### OUTTETMODEL END!! ###"
return
end

!####################################  ASSIGNCOND
subroutine ASSIGNCOND(h_mesh,xgrd,ygrd,zgrd,nx,ny,nz,modelrho,tetrho,m_param)
use mesh_type
use cond_type
implicit none
include "meshpara.f90"
type(mesh),intent(in) :: h_mesh
type(param_cmod),intent(in) :: m_param
integer(4),intent(in) :: nx,ny,nz
real(8),intent(in)  :: xgrd(nx+1),ygrd(ny+1),zgrd(nz+1)
real(8),intent(in)  :: modelrho(nx*ny*nz)
real(8),intent(out) :: tetrho(h_mesh%ntet)
integer(4) :: nodeg,ntetg
integer(4),allocatable,dimension(:,:) :: n4g
real(8),allocatable,dimension(:,:)    :: xyzg
integer(4) :: i,j,k,l, ijk(3), nxyz(3), inum,nn
real(8) :: cxyz(3),x2(2),y2(2),z2(2),ratio,rho_o,rho_n
real(8),allocatable,dimension(:,:) :: xyzgrd
nn=max(nx+1,ny+1,nz+1)
allocate(xyzgrd(nn,3))

!#[1]# Preparation
nxyz(1:3)=(/nx,ny,nz/)
xyzgrd(1:nx+1,1)=xgrd(1:nx+1) !  [km]
xyzgrd(1:ny+1,2)=ygrd(1:ny+1) !  [km]
xyzgrd(1:nz+1,3)=zgrd(1:nz+1) !  [km]
!
ntetg            = h_mesh%ntet
nodeg            = h_mesh%node
allocate(n4g(ntetg,5),xyzg(3,nodeg))
n4g(1:ntetg,2:5) = h_mesh%n4(1:ntetg,1:4)
n4g(1:ntetg,1)   = h_mesh%n4flag(1:ntetg,1)
xyzg             = h_mesh%xyz

!#[2]## assign loop
do l=1,ntetg
 if ( n4g(l,1) .eq. 1 ) then
  tetrho(l)=1.d6  ! in the air
 else
 ![2-1]# center of gravity of this tetrahedron
 cxyz(:)=(xyzg(:,n4g(l,2))+xyzg(:,n4g(l,3))+xyzg(:,n4g(l,4))+xyzg(:,n4g(l,5)))/4.d0


 ![2-2]# search the corresponding grid, ijk
 ijk(:)=0
 do i=1,3
  do j=1,nxyz(i)
   if ((xyzgrd(j,i) .lt. cxyz(i) .and. cxyz(i) .le. xyzgrd(j+1,i)) .or. &
     & (xyzgrd(j+1,i) .lt. cxyz(i) .and. cxyz(i) .le. xyzgrd(j,i)) ) then
   ijk(i)=j
   goto 100
   end if
  end do ! nx, ny, nz loop  [j]
  100 continue ! next component
 end do  ! 3 component loop [i]

 ![2-3]# assign the conductivity to tetrahedron
 if ( ijk(1) .ne. 0 .and. ijk(2) .ne. 0 .and. ijk(3) .ne. 0 ) then
  inum=nx*ny*(ijk(3)-1)+ny*(ijk(1)-1)+ijk(2)
  ! write(*,*) "inum=",inum,"ijk(1:3)=",ijk(1:3)
  tetrho(l)=modelrho(inum)
  if ( tetrho(l) .ge. 1.d+5) tetrho(l)=1.d+5
 else ! the corresponding tet is not found
  if ( n4g(l,1) .eq. 2 ) tetrho(l)=1.d+3 ! in the ground
 end if

 ![2-4]# search cuboid for this tetrahedron added on 2017.02.22
  do i=1,m_param%ncuboid
   x2=m_param%g_cuboid(i)%xminmax
   y2=m_param%g_cuboid(i)%yminmax
   z2=m_param%g_cuboid(i)%zminmax
   if (z2(1) .lt. cxyz(3) .and. cxyz(3) .lt. z2(2) .and. &
   &   y2(1) .lt. cxyz(2) .and. cxyz(2) .lt. y2(2) .and. &
   &   x2(1) .lt. cxyz(1) .and. cxyz(1) .lt. x2(2) ) then
     ratio = m_param%g_cuboid(i)%ratio
     rho_o = log10(tetrho(l))
     rho_n = log10(m_param%g_cuboid(i)%rho)
     tetrho(l)= 10.d0**(ratio*rho_n + (1.-ratio)*rho_o)  ![Ohm.m]
     goto 110
   end if
  end do

 end if

110 continue
end do

write(*,*) "### ASSIGNCOND END!! ###"
return
end
!#####################################  readobs
subroutine readobs(nobsmax_amt,nobs_amt,xyobs,obsfile_xy,obsfile_loc,clatlon)
use param
implicit none
include "meshpara.f90"
character(50),intent(in) :: obsfile_xy, obsfile_loc
integer(4),intent(in) :: nobsmax_amt
integer(4),intent(out) :: nobs_amt
real(8),intent(out) :: xyobs(2,nobsmax_amt) ! [m]
integer(4) :: i,ii
character(30),dimension(:),allocatable :: name
real(8),dimension(:),allocatable :: utm_x,utm_y,utm_x0,utm_y0
real(8),dimension(:),allocatable :: lat,lon,z,clat,clon
real(8) :: xobs(nobsmax_amt),yobs(nobsmax_amt),londev,latdev,utm_x00(1),utm_y00(1)
real(8),intent(out) :: clatlon(2)

!#[1]## obsfile_xy
open(1,file=obsfile_xy)
do i=1,nobsmax_amt
 read(1,*,end=99) xyobs(2,i), xyobs(1,i), ii ! xyobs(2,i):north, xyobs(1,i):east [m]
 !write(*,*) xyobs(2,i), xyobs(1,i), ii ! first is north, second is east coordinate
end do
99 close(1)
nobs_amt=i-1
yobs(1:nobs_amt)=xyobs(2,1:nobs_amt)
xobs(1:nobs_amt)=xyobs(1,1:nobs_amt)

!#[2]## obsfile_loc
allocate( lat(nobs_amt), lon(nobs_amt), z(nobs_amt) )
open(1,file=obsfile_loc)
do i=1,nobs_amt
 read(1,*)lat(i),lon(i),z(i),name
end do
close(1)

!#[3]## calculate the center lon and lat
!# LON LAT <-> UTM are conducted by subroutine UTMGMT
allocate(utm_x(nobs_amt),utm_y(nobs_amt),utm_x0(nobs_amt),utm_y0(nobs_amt) )
allocate(clat(nobs_amt), clon(nobs_amt))
call UTMGMT_N(nobs_amt,lon(1:nobs_amt),lat(1:nobs_amt),&
	    & utm_x(1:nobs_amt),utm_y(1:nobs_amt),"52S",0) ! LONLAT2UTM
utm_x0(1:nobs_amt)= utm_x(1:nobs_amt) - xobs(1:nobs_amt)
utm_y0(1:nobs_amt)= utm_y(1:nobs_amt) - yobs(1:nobs_amt)
call UTMGMT_N(nobs_amt,utm_x0(1:nobs_amt),utm_y0(1:nobs_amt),&
           &clon(1:nobs_amt),clat(1:nobs_amt),"52S",1) ! UTM2LONLAT
do i=1,nobs_amt
 clatlon(1:2)=clatlon(1:2)+(/clat(i),clon(i)/)
 write(*,*) "i=",i,"clat,clon=",clat(i),clon(i)
end do
clatlon(1:2)=clatlon(1:2)/float(nobs_amt)

write(*,*) "clatlon(1:2)=",clatlon(1),clatlon(2)

!#[4]## check standard deviation [m]
CALL UTMGMT_N(1,clatlon(2:2),clatlon(1:1), utm_x00,utm_y00,"52S",0) ! 20200728
londev=0.d0 ; latdev=0.d0
do i=1,nobs_amt
 londev=londev+(utm_x0(i)-utm_x00(1))**2.d0
 latdev=latdev+(utm_y0(i)-utm_y00(1))**2.d0
end do
 londev=dsqrt(londev/float(nobs_amt))
 latdev=dsqrt(latdev/float(nobs_amt))
write(*,*) "standard deviation of clon=",londev,"[m]"
write(*,*) "standard deviation of clat=",latdev,"[m]"

write(*,*) "### readobs END!! ###"
return
end
!#####################################  outmodelmsh
subroutine outmodelmsh(modelrho,nx,ny,nz,xgrd,ygrd,zgrd,outmshfile,xyobs,nobs_amt)
implicit none
integer(4),intent(in) :: nx,ny,nz
real(8),intent(in) :: modelrho(nx*ny*nz)
real(8),intent(in) :: xgrd(nx+1),ygrd(ny+1),zgrd(nz+1)
character(50),intent(in) :: outmshfile
!### info of obs
integer(4) :: nobs_amt
real(8),intent(in) :: xyobs(2,nobs_amt)
!## internal
integer(4) :: i,j,k,n1,n2,n3,n4,n5,n6,n7,n8, iele
open(1,file=outmshfile)
!#[1]## Basical data input to msh file
write(1,'(a)') "$MeshFormat"
write(1,'(a)')  "2.2 0 8"
write(1,'(a)')  "$EndMeshFormat"
write(1,'(a)')  "$Nodes"
write(1,*)  (nx+1)*(ny+1)*(nz+1)+nobs_amt
!After convertaxis, x is east, y: north, z:upward
!Then the order of node and element coord be (((ymin,ymax),xmin,xmax),zmax,zmin)
do k=1,nz+1
 do i=1,nx+1
  do j=1,ny+1
   write(1,'(i10,3f15.7)') (k-1)*(nx+1)*(ny+1)+(i-1)*(ny+1)+j,xgrd(i),ygrd(j),zgrd(k)
  end do
 end do
end do
!# observation sites
do i=1,nobs_amt
 j=(nx+1)*(ny+1)*(nz+1)+i
 write(1,'(i10,2f15.7,1x,a)') j, xyobs(1,i),xyobs(2,i),"0.0"
end do
!#
write(1,'(a)') "$EndNodes"
write(1,'(a)') "$Elements"
write(1,*)  nx*ny*nz + nobs_amt
do k=1,nz
 do i=1,nx
  do j=1,ny
    iele=nx*ny*(k-1)+ny*(i-1)+j
    n1=(nx+1)*(ny+1)*(k-1)+(ny+1)*(i-1)+j
    n2=n1+1
    n3=n1+(ny+1)+1
    n4=n1+(ny+1)
    n5=n1+(nx+1)*(ny+1)
    n6=n2+(nx+1)*(ny+1)
    n7=n3+(nx+1)*(ny+1)
    n8=n4+(nx+1)*(ny+1)
    write(1,'(i8,a6,9i8)') iele, " 5 2 0 ",iele,n1,n2,n3,n4,n5,n6,n7,n8
  end do
 end do
end do
do i=1,nobs_amt
 write(1,'(i10,a10,i10)') nx*ny*nz+i, " 15 2 0 1 ",(nx+1)*(ny+1)*(nz+1)+i
end do
write(1,'(a)') "$EndElements"
!goto 100
!#[2]## create vector fields
write(1,'(a)') "$ElementData"
write(1,'(a)') "1"
write(1,'(a)') '"A conductivity model view"'
write(1,'(a)') "1"
write(1,'(a)') "0.0"
write(1,'(a)') "3"
write(1,'(a)') "0"
write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
write(1,*) nx*ny*nz
do i=1,nx*ny*nz
   write(1,*) i,modelrho(i)
end do
write(1,'(a)') "$EndElementData"
100 continue
close(1)
write(*,*)"### outmodelmsh END!! ###"
return
end
!#####################################  count model
subroutine countmodel(modelfile,nx,ny,nz)
implicit none
character(50),intent(in) :: modelfile
integer(4),intent(out) :: nx,ny,nz
integer(4) :: i
open(1,file=modelfile)
 read(1,*) ! header
 read(1,*) nx,ny,nz,i
close(1)
write(*,*) "### coundmodel END!! ###"
return
end
!#####################################  read model
subroutine readmodel(modelfile,dx,dy,dz,nx,ny,nz,modelrho)
implicit none
character(50),intent(in) :: modelfile
integer(4),intent(in) :: nx,ny,nz
real(8),intent(out) :: dx(nx),dy(ny),dz(nz)
real(8),intent(out) :: modelrho(nx*ny*nz)
integer(4) :: i
open(1,file=modelfile)
 read(1,*) ! header
 read(1,*) ! nx,ny,nz
 read(1,'(7e12.4)') (dx(i),i=1,nx)
 read(1,'(7e12.4)') (dy(i),i=1,ny)
 read(1,'(7e12.4)') (dz(i),i=1,nz)
 read(1,'(e12.4)') (modelrho(i),i=1,nx*ny*nz)
close(1)
! write(*,'(7e12.4)') (dx(i),i=1,nx)
! write(*,'(7e12.4)') (dy(i),i=1,ny)
! write(*,'(7e12.4)') (dz(i),i=1,nz)
! write(*,'(e12.4)') (modelrho(i),i=1,nx*ny*nz)
 write(*,*) "### readmodel END!! ###"
return
end
!##################################### gen grid coordinate
subroutine gengrd(dx,dy,dz,nx,ny,nz,xgrd,ygrd,zgrd)
implicit none
integer(4),intent(in) :: nx,ny,nz
real(8),intent(in) :: dx(nx),dy(ny),dz(nz)
real(8),intent(out) :: xgrd(nx+1),ygrd(ny+1),zgrd(nz+1)
integer(4) :: i,j,k
real(8) :: half
!#[1]## y grd
xgrd(1)=0.d0
do i=2,nx+1
 xgrd(i)=xgrd(i-1)+dx(i-1)
end do
half=xgrd(nx+1)/2.d0
xgrd(1:nx+1)=xgrd(1:nx+1) - half

!#[2]## y grd
ygrd(1)=0.d0
do i=2,ny+1
 ygrd(i)=ygrd(i-1)+dy(i-1)
end do
half=ygrd(ny+1)/2.d0
ygrd(1:ny+1)=ygrd(1:ny+1) - half

!#[3]## z grd
zgrd(1)=0.d0
do i=2,nz+1
 zgrd(i)=zgrd(i-1)+dz(i-1)
end do

!#[4]## check
!write(*,'(7e12.4)') (xgrd(i),i=1,nx+1)
!write(*,'(7e12.4)') (ygrd(i),i=1,ny+1)
!write(*,'(7e12.4)') (zgrd(i),i=1,nz+1)
return
end
!####################################### convertzxis
subroutine adjustcoord(nx,ny,nz,xgrd,ygrd,zgrd,clatlon,ztop_grd)
use param ! 20200728
implicit none
include "meshpara.f90"
integer(4),intent(in) :: nx,ny,nz
real(8),intent(in) :: clatlon(2),ztop_grd ! [deg] and [m]
real(8),intent(inout) :: xgrd(nx+1),ygrd(ny+1),zgrd(nz+1)
real(8),allocatable,dimension(:) :: xgrd_old,ygrd_old,zgrd_old
real(8) :: utm_x0(2),utm_y0(2)
integer(4) :: i
!#[1]## check if nx is equal to ny
if (nx .ne. ny) then
 write(*,*) "GEGEGE! nx .ne. ny. This program cannot deal with not-equal nx and ny!"
 write(*,*) "nx=",nx,"ny=",ny
end if

!#[2]## convert x and y
allocate(xgrd_old(nx+1),ygrd_old(ny+1),zgrd_old(nz+1))
xgrd_old(:)=xgrd(:)/1000.d0 ! [m] -> [km]
ygrd_old(:)=ygrd(:)/1000.d0 ! [m] -> [km]
zgrd_old(:)=zgrd(:)/1000.d0 ! [m] -> [km]
xgrd(:)=ygrd_old(:)
ygrd(:)=-xgrd_old(:) ! on Dec. 15, 2015 minus is added because north south is reversed
zgrd(:)=-zgrd_old(:) ! downward positive to upward positive

!#[3]## shift the origin
!# clatlon are the origin for gridded model
!# latorigin, lonorigin are the origin for msh file
CALL UTMGMT_N(2,(/clatlon(2),lonorigin/),(/clatlon(1),latorigin/),utm_x0,utm_y0,"52S",0)
write(*,*) "utm_x0(1)-utm_x0(2)=",utm_x0(2)-utm_x0(1),"[m]"
write(*,*) "utm_y0(1)-utm_y0(2)=",utm_y0(2)-utm_y0(1),"[m]"
xgrd(:)=xgrd(:)+(utm_x0(1)-utm_x0(2))/1.d3 - 0.04 ! [km] ! Dec. 21, 2015
ygrd(:)=ygrd(:)+(utm_y0(1)-utm_y0(2))/1.d3 + 0.04 ! [km] ! -0.04 is shift from comparison between kandatopo.msh and nakadake2d_topo.msh
zgrd(:)=zgrd(:)+ztop_grd/1000.d0 ! [km]

!write(*,'(7e12.4)') (xgrd(i),i=1,nx+1)
!write(*,'(7e12.4)') (ygrd(i),i=1,ny+1)
!write(*,'(7e12.4)') (zgrd(i),i=1,nz+1)
!write(*,*) "### convertaxis END!! ###"
return
end

