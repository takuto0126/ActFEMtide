! coded on 2020.10.29
module surface_type
use iccg_var_takuto ! see ../solver/m_iccg_var_takuto.f90
use matrix ! 2022.01.27
implicit none

type surface

! line information
! The directions of lines in surface and 3D mesh lines are the same
! ilineface_to_3D(1:nline)
integer(4) :: nline    ! # of lines in the face
integer(4),allocatable,dimension(:,:) :: line ! line(1:2,1:nline_surface)
integer(4),allocatable,dimension(:)   :: ilineface_to_3D ! j=iedgeface_to_3D(i) : i face edge, j 3D edge
integer(4),allocatable,dimension(:)   :: iflag_bound ! iflag_bound(1:nline) 2021.05.25
! iflag_bound : 0 not boundary                     ! 2021.05.25
!               1 top   ; 2 bottom boundary line    ! 2021.05.25
!               3 left  ; 4 right  boundary line    ! 2021.05.25
!

character(2) :: facetype ! "xy", "yz", or "xz"
real(8),allocatable,dimension(:,:) :: x1x2_face ! xy_face(1:2,npoi)

! node information
integer(4) :: node  ! # of points in the face
integer(4) :: node3d ! # of points in the parent 3D mesh 2021.12.30
! inodeface_to_3D(1:npoi), inode3D_to_face(1:npoi)
integer(4),allocatable,dimension(:) :: inodeface_to_3D ! l=iedgeface_to_3D(k) : k face node, l 3D node
integer(4),allocatable,dimension(:) :: inode3D_to_face !
!integer(4),allocatable,dimension(:) :: face_stack      ! l=iedgeface_to_3D(k) : k face node, l 3D node
!integer(4),allocatable,dimension(:) :: face_stack_item ! l=iedgeface_to_3D(k) : k face node, l 3D node


! triangle information
integer(4) :: ntri  ! # of triangles
!# Rule for n3: smallest node number to large and unticlockwise in the following plane
!# Rule for n3line: line i starts with ith node and toward unticlockwise
!
!   <xy surface>           <yz surface>        <xz surface>
!  y       3             z      3            z       3
!  ^    |  |\ r (2)      ^    | |\ r (2)     ^     | |\ r (2)
!  | (3)v  | \ \         | (3)v | \ \        |  (3)v | \ \
!  |     1 ---\ 2        |     1 --\ 2       |      1 -->\ 2
!  |       (1)           |       (1)         |        (1)
!  --------------> x     -------------> y    -----------------> x
!  () is the line direction and number for n3lin3
!
integer(4),allocatable,dimension(:,:) :: n3 ! n3(1:ntri,1:3) node id for the triangle
integer(4),allocatable,dimension(:,:) :: n3line ! n3line(1:ntri,1:3) line id for the triangle in the surface
integer(4),allocatable,dimension(:)   :: ifacetri_to_face ! triangle id in the surface to face id in 3dmesh
integer(4),allocatable,dimension(:)   :: ifacetri_to_tet ! triangle id in the surface to tet id in 3dmesh

! coefficient matrix information
type(global_matrix) :: A   ! see m_iccg_var_takuto.f90
integer(4),allocatable,dimension(:,:) ::table_dof
! conductivity of triangle
real(8),   allocatable,dimension(:)   :: cond ! cond(1:ntri) [S/m]

! boundary condition and solution  2021.05.31
logical,   allocatable,dimension(:)   :: line_bc   ! line_bc(nline)
complex(8),allocatable,dimension(:)   :: Avalue_bc ! Avalue_bc(nline)
complex(8),allocatable,dimension(:)   :: bs        ! bs(nline)[mV/km * km],solution obtained from A*bs=b_bec

!# surface model info
integer(4)            :: nmodel_surface    ! 2022.01.26
integer(4)            :: nmodel_global
type(real_crs_matrix) :: model2ele ! 2022.01.16 only stack and element
integer(4)            :: model_face2global ! [nmodel_surface]
integer(4)            :: model_glabal2face ! [nmodel_global] zero for no corresponding model
end type


contains
!##################################################### added in this file on Oct 13, 2021
!# coded on 2021.06.01
subroutine cond3dto2d(g_mesh,g_surface,g_cond)
use param
use mesh_type
implicit none
type(surface),    intent(inout) :: g_surface(6)
type(mesh),       intent(in)    :: g_mesh
type(param_cond), intent(in)    :: g_cond
integer(4),allocatable,dimension(:,:) :: n4flag
integer(4) :: ntri,ntet,nphys2,iele_tet,i,j

ntet   = g_mesh%ntet
n4flag = g_mesh%n4flag

!# allocate g_surface(j)%cond(:)
do j=2,5
 ntri = g_surface(j)%ntri
 if (.not. allocated(g_surface(j)%cond )) allocate( g_surface(j)%cond(ntri) ) ! 2021.12.31
end do

do j=2,5 ! 5 faces except top surface
 
 ntri = g_surface(j)%ntri
 if ( .not. allocated(g_surface(j)%cond)) allocate(g_surface(j)%cond(ntri))

 do i=1,ntri

 iele_tet = g_surface(j)%ifacetri_to_tet(i)

 if (      n4flag(iele_tet,1) .eq. 1 ) then ! in the air
  g_surface(j)%cond(i) = g_cond%sigma_air
 else if ( n4flag(iele_tet,1) .eq. 2 ) then ! under ground
  if ( g_cond%condflag .eq. 0)g_surface(j)%cond(i)=g_cond%sigma_land(1) ! 2017.09.29
  if ( g_cond%condflag .eq. 1)g_surface(j)%cond(i)=g_cond%sigma(iele_tet - g_cond%nphys1)
 else
  write(*,*) "GEGEGE n4flag=",n4flag(iele_tet,1)
  stop
 end if
 end do ! surface triangle loop

end do ! surface loop

return
end
!############################################################## FINDBOUNDARYLINE
!# 2021.05.25
subroutine FINDBOUNDARYLINE(g_mesh,g_surface)
use mesh_type   ! src_mesh/m_mesh_type.f90
implicit none
type(mesh),   intent(in) :: g_mesh
type(surface),intent(inout) :: g_surface(6)
real(8)    :: xmin,xmax,ymin,ymax,zmin,zmax
real(8)    :: yleft,yright,zbot,ztop
real(8)    :: y1,y2,z1,z2
integer(4) :: i,i1,i2,j

xmin = g_mesh%xyzminmax(1) + 0.1 ! 2021.06.08
xmax = g_mesh%xyzminmax(2) - 0.1 ! 2021.06.08
ymin = g_mesh%xyzminmax(3) + 0.1 ! 2021.06.08
ymax = g_mesh%xyzminmax(4) - 0.1 ! 2021.06.08
zmin = g_mesh%xyzminmax(5) + 0.1 ! 2021.06.08
zmax = g_mesh%xyzminmax(6) - 0.1 ! 2021.06.08

do i=1,6
 allocate(g_surface(i)%iflag_bound(g_surface(i)%nline))
g_surface(i)%iflag_bound(:)=0 ! not boundary (default)
end do

!# search bound
do i=1,6                     ! surface loop
 do j=1,g_surface(i)%nline   ! line loop
  i1 = g_surface(i)%line(1,j)
  i2 = g_surface(i)%line(2,j)
  y1 = g_surface(i)%x1x2_face(1,i1)
  y2 = g_surface(i)%x1x2_face(1,i2)
  z1 = g_surface(i)%x1x2_face(2,i1)
  z2 = g_surface(i)%x1x2_face(2,i2)
  if      (g_surface(i)%facetype .eq. "xy" ) then
   yleft = xmin ; yright = xmax ; zbot = ymin ; ztop = ymax
  else if (g_surface(i)%facetype .eq. "yz" ) then
   yleft = ymin ; yright = ymax ; zbot = zmin ; ztop = zmax
  else if (g_surface(i)%facetype .eq. "xz" ) then
   yleft = xmin ; yright = xmax ; zbot = zmin ; ztop = zmax
  end if

  if (y1 .lt. yleft  .and. y2 .lt. yleft )  g_surface(i)%iflag_bound(j)=3 ! left    bound
  if (y1 .gt. yright .and. y2 .gt. yright)  g_surface(i)%iflag_bound(j)=4 ! right   bound
  if (z1 .lt. zbot   .and. z2 .lt. zbot  )  g_surface(i)%iflag_bound(j)=2 ! bottom  bound
  if (z1 .gt. ztop   .and. z2 .gt. ztop  )  g_surface(i)%iflag_bound(j)=1 ! top     bound
 end do
end do

return
end

!################################################################ EXTRACT6FACES
! coded on 2020.10.29
!
!#[1]## face type
!        Face2        ^ y
!        ----         |
! Face3 |    | Face5  |
!        ----         ----> x
!        Face4
!
! Top     Face 1          z
!         ------          ^
! Face 3  |    | Face 5   |
!         ------          -----> x
! Bottom  Face 6
!

subroutine EXTRACT6SURFACES(g_mesh,l_line,g_face,g_surface)
use mesh_type
use line_type
use face_type
use matrix
use outerinnerproduct
implicit none
type(surface),   intent(out) :: g_surface(6)
type(mesh),      intent(in)  :: g_mesh
type(line_info), intent(in)  :: l_line
type(face_info), intent(in)  :: g_face ! see ../common/m_face_type.f90
logical,   allocatable,dimension(:,:)   :: node_on_surface
integer(4),allocatable,dimension(:,:)   :: inode3D_to_face
integer(4),allocatable,dimension(:,:)   :: inodeface_to_3D
integer(4),allocatable,dimension(:,:)   :: ifacetri_to_face
integer(4),allocatable,dimension(:,:)   :: ifacetri_to_tet
integer(4),allocatable,dimension(:,:)   :: iline3D_to_face
integer(4),allocatable,dimension(:,:)   :: ilineface_to_3D
integer(4),allocatable,dimension(:)     :: face
integer(4),allocatable,dimension(:,:)   :: line_3D
integer(4),allocatable,dimension(:,:)   :: n6line
real(8),   allocatable,dimension(:,:)   :: xyz
real(8),   allocatable,dimension(:,:,:) :: x1x2_face
integer(4),allocatable,dimension(:,:,:) :: n3line,n3,line_surface
integer(4),allocatable,dimension(:)     :: face_stack
integer(4),allocatable,dimension(:,:)   :: face_item
integer(4),allocatable,dimension(:,:)   :: face2ele
character(2) :: facetype
integer(4)   :: i,j,k,l,line(2,3),ntet,node,nline,nface,i_3D,iface_3D,inode1,inode2
integer(4)   :: i1,i2,i3,i1_3D,i2_3D,i3_3D,ifa(2),itet_3D
integer(4)   :: node_surface(6),nline_surface(6),nface_surface(6)
integer(4)   :: icount(6),icount_line(6),iface_count(6)
integer(4)   :: nface_surface_max,node_surface_max,nline_surface_max
real(8)      :: xmin,xmax,ymin,ymax,zmin,zmax,line_f(2)
real(8)      :: x1(3),x2(3),x3(3),x,y,z
integer(4)   :: nsurface=6
!#[0]## set
node       = g_mesh%node
ntet       = l_line%ntet
nline      = l_line%nline
nface      = g_face%nface
xyz        = g_mesh%xyz       ! allocate here 2020.11.01
n6line     = l_line%n6line     ! n6line(ntet,6)     allocate here 2020.11.01
line_3D    = l_line%line       ! line_3D(2,nline),  allocate here 2020.11.01
face_stack = g_face%face_stack ! face_stack(0:node) allocate here
face_item  = g_face%face_item  ! face_item(2,nface) allocate here
face2ele   = g_face%face2ele   ! allocate face2ele here

g_surface(1)%facetype="xy" ! top
g_surface(2)%facetype="xz" ! north
g_surface(3)%facetype="yz" ! west
g_surface(4)%facetype="xz" ! south
g_surface(5)%facetype="yz" ! east
g_surface(6)%facetype="xy" ! bottom

!#[2]## extract line
xmin = g_mesh%xyzminmax(1) +0.1
xmax = g_mesh%xyzminmax(2) -0.1
ymin = g_mesh%xyzminmax(3) +0.1
ymax = g_mesh%xyzminmax(4) -0.1
zmin = g_mesh%xyzminmax(5) +0.1
zmax = g_mesh%xyzminmax(6) -0.1

!#[3]## node group
allocate(node_on_surface(6,node))
node_on_surface(:,:)=.false.
write(*,*) "node=",node
do i=1,node
     x=xyz(1,i) ; y=xyz(2,i) ; z=xyz(3,i)
!     write(*,*) "i",i,"x,y,z",x,y,z
     if   ( z .ge. zmax ) node_on_surface(1,i)=.true. ! Face 1
     if   ( y .ge. ymax ) node_on_surface(2,i)=.true. ! Face 2
     if   ( x .le. xmin ) node_on_surface(3,i)=.true. ! Face 3
     if   ( y .le. ymin ) node_on_surface(4,i)=.true. ! Face 4
     if   ( x .ge. xmax ) node_on_surface(5,i)=.true. ! Face 5
     if   ( z .le. zmin ) node_on_surface(6,i)=.true. ! Face 6
end do    ! line loop

!#[4]## count node on surface
node_surface(:)=0
node_surface_max=0
do i=1,node
do j=1,nsurface ! surface loop
 if ( node_on_surface(j,i) ) then
  node_surface(j) = node_surface(j) + 1
 end if
end do
end do
node_surface_max = maxval(node_surface)
write(*,'(a,6i7)') " # of nodes of each surface",node_surface(1:6) ! 2021.10.13
write(*,*) "node_surface_max=",node_surface_max

!#[5]## allocate node
allocate(inodeface_to_3D(node_surface_max,6))
allocate(inode3D_to_face(node,6))
allocate(x1x2_face(2,node_surface_max,6))

!#[6]## set coordinate of node on surface
icount(:)=0
do i=1,node ! 3D node
 do j=1,nsurface ! surface loop
   facetype = g_surface(j)%facetype
   if (node_on_surface(j,i)) then
    icount(j) = icount(j) + 1 ! face id in jth surface
    inodeface_to_3D(icount(j),j) = i         ! node id in the 3dmesh
    inode3D_to_face(i,        j) = icount(j) ! node id in the surface
   if ( facetype .eq. "xy") x1x2_face(1:2,icount(j),j) = xyz(1:2,i)
   if ( facetype .eq. "yz") x1x2_face(1:2,icount(j),j) = xyz(2:3,i)
   if ( facetype .eq. "xz") x1x2_face(1:2,icount(j),j)=(/xyz(1,i),xyz(3,i)/)
   end if
 end do
end do

if (.false.) then ! 2021.12.22
open(1,file="inode3dtoface.dat")
do i=1,node
if(node_on_surface(1,i)) write(1,*) i,inode3D_to_face(i,1), xyz(1:3,i)
end do
close(1)
open(1,file="inodefaceto3D.dat")
do i=1,node_surface(1)
 write(1,*) i,inodeface_to_3D(i,1)
end do
close(1)
end if ! 2021.21.22

!#[7]## identify face on the surface
nface_surface_max=10000
allocate(ifacetri_to_face(    nface_surface_max,6))
allocate(n3(nface_surface_max,3,6))

iface_count(:)=0
do j=1,nsurface ! surface loop
! open(1,file="n3_3d.dat") ! commented out 2021.12.22
 do i=1,node_surface(j)
  i_3D = inodeface_to_3D(i,j)
  do k=face_stack(i_3D - 1)+1,face_stack(i_3D) ! k=face id in 3Dmesh; face loop for face starting from i_3D (smallest node id in the face)
   iface_3D = k
   inode1 = face_item(1,iface_3D) ! remaining face node id in 3dmesh
   inode2 = face_item(2,iface_3D) ! remaining face node id in 3dmesh
   if ( node_on_surface(j,inode1) .and. node_on_surface(j,inode2) ) then ! face in on surface !!
!    write(1,*) k,i_3D,inode1,inode2 ! commented out 2021.12.22
    iface_count(j)=iface_count(j) + 1
    n3(iface_count(j),1:3,j)=(/i,inode3D_to_face(inode1,j),inode3D_to_face(inode2,j)/)
    ifacetri_to_face(iface_count(j),j) = iface_3D ! face in surface to face in 3d mesh
   end if
 end do
 end do
! close(1) ! commented out 2021.12.22
end do
 nface_surface(:)=iface_count(:)
 write(*,'(a,6i7)') " # of face of each surface",nface_surface(1:6)
 nface_surface_max=maxval(nface_surface)

!#[8]## modify n3 anticlockwise and starint with smallest node id
do j=1,nsurface ! surface loop
 do i=1,nface_surface(j)
  ! make i1 is smallest
    i1 = n3(i,1,j) ; i2 = n3(i,2,j) ; i3 = n3(i,3,j)
   if ( n3(i,1,j) .ge. n3(i,2,j)) then
    i1 = n3(i,2,j) ; i2 = n3(i,1,j)
   end if
   if ( i1        .ge. n3(i,3,j)) then
    i3 = i1 ; i1 = n3(i,3,j)
   end if
   n3(i,1:3,j)=(/i1,i2,i3/)
   i1 = n3(i,1,j) ; i2 = n3(i,2,j) ; i3 = n3(i,3,j)
   x1(3)=0. ; x2(3)=0.
!   write(*,'(a,3i4)') "i1,i2,i3",i1,i2,i3
   x1(1:2)=x1x2_face(1:2,i2,j) - x1x2_face(1:2,i1,j)
   x2(1:2)=x1x2_face(1:2,i3,j) - x1x2_face(1:2,i1,j)
   x3=outer(x1,x2)
   if ( x3(3) .lt. 0. ) n3(i,1:3,j)=(/i1,i3,i2/) ! make the node order unticlockwise
!   write(*,'(a,3i3,a,3i4)')"surface,j,i,k",j,i,k,"n3(k,1:3,j)",n3(k,1:3,j)
  end do
 end do

if (.false.) then ! 2021.12.22
 open(2,file="n3_2.dat")
 j=1
 do i=1,nface_surface(j)
  write(2,*) i,n3(i,1:3,j)
 end do
 close(2)
end if ! 2021.12.22

!#[9]## face to element facing the surface
 allocate(ifacetri_to_tet(nface_surface_max,6))
!open(1,file="face2ele.dat") ! 2021.12.22
 do j=1,nsurface
  do i=1,nface_surface(j)
   iface_3D = ifacetri_to_face(i,j)
   ifa(1:2) = face2ele(1:2,iface_3D)
 !  write(1,*) iface_3D,ifa(1:2)  ! 2021.12.22
   if ( ifa(1) .ne. 0 ) ifacetri_to_tet(i,j) = ifa(1)
   if ( ifa(2) .ne. 0 ) ifacetri_to_tet(i,j) = ifa(2)
   if ( ifa(1) .ne. 0 .and.  ifa(2) .ne. 0) then
    write(*,*) "GEGEGE ifa(1:2)",ifa(1:2)
    stop
   end if
end do
end do
!close(1)  ! 2021.12.22

if (.false.) then ! 2021.12.22
open(1,file="ifacetri_to_tet.dat")
 j=1
 do i=1,nface_surface(j)
write(1,*) i,ifacetri_to_tet(i,j)
 end do
close(1)
end if ! 2021.12.22

!#[10]## generate line and connection
  nline_surface_max = 100000
  allocate(ilineface_to_3D(nline_surface_max,6))
  allocate(line_surface(2,nline_surface_max,6))
  allocate(iline3D_to_face(nline,6))
  icount_line(:) = 0
!  open(1,file="line_3D.dat")
  do i= 1,nline
   do j=1,nsurface
   if (node_on_surface(j,line_3D(1,i)) .and. &
&      node_on_surface(j,line_3D(2,i))) then
  !  write(1,*) i,line_3D(1:2,i)
    icount_line(j)=icount_line(j)+1
    line_surface(1,icount_line(j),j)=inode3D_to_face(line_3D(1,i),j)
    line_surface(2,icount_line(j),j)=inode3D_to_face(line_3D(2,i),j)
    ilineface_to_3D(icount_line(j),j)=i !surface line to line in 3Dmesh
    iline3D_to_face(i,j)=icount_line(j)
    end if
   end do
  end do
 ! close(1)
  nline_surface(:)=icount_line(:)
!  open(1,file="line_surface.dat")
!  write(1,'(3i7)') (i,line_surface(1:2,i,1),i=1,nline_surface(1))
!  close(1)
!  open(1,file="ilineface_to_3D.dat")
!  write(1,'(2i7)') (i,ilineface_to_3D(i,1),i=1,nline_surface(1))
!  close(1)

!#[11]## generate n3line
allocate(n3line(nface_surface_max,3,6))
!open(1,file="n3_3D.dat")  ! 2021.12.22
!  open(2,file="n3line_3D.dat") ! 2021.12.22
  do j=1,nsurface
  do i=1,nface_surface(j)
   i1 = n3(i,1,j)
   i2 = n3(i,2,j)
   i3 = n3(i,3,j)
   i1_3D = inodeface_to_3D(i1,j)
   i2_3D = inodeface_to_3D(i2,j)
   i3_3D = inodeface_to_3D(i3,j)
   itet_3D=ifacetri_to_tet(i,j)
   line(1:2,1) =(/i1_3D,i2_3D/)
   line(1:2,2) =(/i2_3D,i3_3D/)
   line(1:2,3) =(/i3_3D,i1_3D/)
  ! write(1,*) i,i1_3D,i2_3D,i3_3D ! 2021.12.22
   do l=1,3
    do k=1,6
     line_f(1:2)=line_3D(1:2,abs(n6line(itet_3D,k)))
    if ( line(1,l) .eq. line_f(1) .and. line(2,l) .eq. line_f(2))then
!      write(2,*) i,line_3D(1:2,abs(n6line(itet_3D,k))),"+1" ! 2021.12.22
      n3line(i,l,j) = iline3D_to_face(abs(n6line(itet_3D,k)),j)
      goto 10
     end if
    if ( line(1,l) == line_f(2) .and. line(2,l) == line_f(1))then
 !    write(2,*) i,line_3D(1:2,abs(n6line(itet_3D,k))),"-1" ! 2021.12.22
     n3line(i,l,j)= - iline3D_to_face(abs(n6line(itet_3D,k)),j)
     goto 10
    end if
   end do ! k loop
   10 continue
   end do ! l line loop
  end do  ! i face loop
  end do  ! j surafce loop
 ! close(1) ! 2021.12.22
!close(2) ! 2021.12.22

!# set output to g_surface
do j=1,nsurface
 allocate( g_surface(j)%x1x2_face(    2,node_surface(j))   )
 allocate( g_surface(j)%inodeface_to_3D(node_surface(j))   )
 allocate( g_surface(j)%inode3D_to_face(node)              )
 allocate( g_surface(j)%n3line(nface_surface(j),3)         )
 allocate( g_surface(j)%n3(node_surface(j),3)              )
 allocate( g_surface(j)%ilineface_to_3D(nline_surface(j))  )
 allocate( g_surface(j)%line(2,nline_surface(j))           )
 allocate( g_surface(j)%ifacetri_to_face(nface_surface(j)) )
 allocate( g_surface(j)%ifacetri_to_tet( nface_surface(j)) )
 g_surface(j)%node3d    = node             ! 2021.12.30
 g_surface(j)%node      = node_surface(j)
 g_surface(j)%ntri      = nface_surface(j)
 g_surface(j)%nline     = nline_surface(j)
! arrays
 g_surface(j)%x1x2_face = x1x2_face(1:2,1:node_surface(j),j)
 g_surface(j)%n3        = n3(1:nface_surface(j),1:3,j)
 g_surface(j)%n3line    = n3line(1:nface_surface(j),1:3,j)
 g_surface(j)%inodeface_to_3D = inodeface_to_3D(1:node_surface(j),j)
 g_surface(j)%inode3D_to_face = inode3D_to_face(1:node,j)
 g_surface(j)%line      = line_surface(1:2,1:nline_surface(j),j)
 g_surface(j)%ilineface_to_3D = ilineface_to_3D(1:nline_surface(j),j)
 g_surface(j)%ifacetri_to_face = ifacetri_to_face(1:nface_surface(j),j)
 g_surface(j)%ifacetri_to_tet = ifacetri_to_tet(1:nface_surface(j),j)
end do

!# output surface
if (.false.) call outsurface(g_surface) ! see below 2021.12.22 false added

write(*,*) "### EXTRACT6FACES END!! ###" ! 2021.09.16
return
end
!======================================================  outsurface
! on 2020.10.31
subroutine outsurface(g_surface)
use mesh_type
implicit none
type(surface),intent(in) :: g_surface(6)
type(mesh)   :: h_mesh
integer(4)   :: iout,i,j
character(1) :: ci
integer(4)   :: nsurface=6

!# MESHOUT
iout=11
do j=1,nsurface
!write(*,*) "j=",j
write(ci,'(i1)') j

!# point
h_mesh%node = g_surface(j)%node
allocate(h_mesh%xyz(3,h_mesh%node))
h_mesh%xyz=0.d0
if (g_surface(j)%facetype .eq. "xy") h_mesh%xyz(1:2,:)=g_surface(j)%x1x2_face
if (g_surface(j)%facetype .eq. "yz") h_mesh%xyz(2:3,:)=g_surface(j)%x1x2_face
if (g_surface(j)%facetype .eq. "xz") then
 h_mesh%xyz(1,:)=g_surface(j)%x1x2_face(1,:)
 h_mesh%xyz(3,:)=g_surface(j)%x1x2_face(2,:)
end if

!# line
h_mesh%nlin = g_surface(j)%nline
allocate(h_mesh%n2(h_mesh%nlin,2))
allocate(h_mesh%n2flag(h_mesh%nlin,2))
h_mesh%n2(:,1) = g_surface(j)%line(1,:)
h_mesh%n2(:,2) = g_surface(j)%line(2,:)
h_mesh%n2flag(:,1:2)=0

!# n3
h_mesh%ntri = g_surface(j)%ntri
allocate(h_mesh%n3(h_mesh%ntri,3))
h_mesh%n3 = g_surface(j)%n3
allocate(h_mesh%n3flag(h_mesh%ntri,2))
h_mesh%n3flag(:,:)=0

open(iout,file="surface"//ci//".msh")
 write(*,*) "Output Mesh: surface"//ci//".msh" ! 2021.10.13
 call MESHOUT(iout,h_mesh) ! see src_mesh/m_mesh_type.f90
close(iout)

!# geofile
if ( .false. ) then ! 2021.12.22
open(1,file="surface"//ci//".geo")
write(1,*) "lc2=1.0;"
do i=1,h_mesh%node
  write(1,*) "Point(",i,")={",h_mesh%xyz(1,i),",",h_mesh%xyz(2,i),",",h_mesh%xyz(3,i),"};"
end do
!#[2]## Lines
do i=1, h_mesh%nlin
write(1,*) "Line(",i,")={",h_mesh%n2(i,1),",",h_mesh%n2(i,2),"};"
end do
!#[3]## Lineloops and Plane Surface
do i=1,h_mesh%ntri
  write(1,*) "Line Loop(",i,")={",g_surface(j)%n3line(i,1),",",g_surface(j)%n3line(i,2), &
&            ",",g_surface(j)%n3line(i,3),"};"
end do
do i=1,h_mesh%ntri
  write(1,*) "Plane Surface(",i,")={",i,"};"
end do
close(1)
end if             ! 2021.12.22

deallocate(h_mesh%xyz)
deallocate(h_mesh%n2,h_mesh%n3,h_mesh%n2flag,h_mesh%n3flag)

end do ! j loop

return
end

!##################################### 2021.06.01
subroutine searchtri(g_surface,y,z,iele)
use outerinnerproduct
implicit none
type(surface),intent(in)  :: g_surface
real(8),      intent(in)  :: y,z ! obs coordinate
integer(4)  , intent(out) :: iele
integer(4) :: i,j
real(8)    :: yt(3),zt(3),x12(3),x13(3),x23(3),r1(3),r2(3)
real(8)    :: o1(3),o2(3),o3(3)
integer(4),allocatable,dimension(:,:) :: n3 ! n3(ntri,3)

n3 = g_surface%n3

do i=1,g_surface%ntri
 
 do j=1,3
  yt(j) = g_surface%x1x2_face(1,n3(i,j))
  zt(j) = g_surface%x1x2_face(2,n3(i,j))
 end do
 x12(1:3) = (/ 0.d0, yt(2) - yt(1) , zt(2) - zt(1) /)
 x13(1:3) = (/ 0.d0, yt(3) - yt(1) , zt(3) - zt(1) /)
 x23(1:3) = (/ 0.d0, yt(3) - yt(2) , zt(3) - zt(2) /)
 r1( 1:3) = (/ 0.d0, y     - yt(1) , z     - zt(1) /)
 r2( 1:3) = (/ 0.d0, y     - yt(2) , z     - zt(2) /)
 o1 = outer(r1,x13)
 o2 = outer(x12,r1)
 o3 = outer(x23,r2)
! write(*,*) "yt",yt
! write(*,*) "zt",zt
! write(*,'(i4,3(a,f15.7))')i, " o1",o1(1)," o2(1)",o2(1)," o3(1)",o3(1)
 if ( o1(1) .ge. 0. .and. o2(1) .ge. 0. .and. o3(1) .ge. 0. ) then
  iele = i
  goto 100
 end if
end do

write(*,*) "GEGEGE iele is not found for y",y,"z",z
stop

100 continue
!do i=1,3
!write(*,*)"y,z", yt(i),zt(i)
!end do
!write(*,*) "o1(:)",o1
write(*,*) "### EARCHTRI END!! ###"
return
end
!#########################################
! 2021.05.31
subroutine rhoap_tri(e,bx,omega,rhoa,pha)
use constants, only:pi,dmu
implicit none
complex(8),intent(in)  :: e(2) ! [mV/km]
complex(8),intent(in)  :: bx   ! [nT]
real(8),   intent(in)  :: omega
real(8),   intent(out) :: rhoa ! apparent resistivity [Ohm.m]
real(8),   intent(out) :: pha  ![deg]
real(8)    :: r2d
complex(8) :: z

r2d=180./pi
z=e(1)/bx * 1.d+3 ! Z [V/m]/[T] = [mV/km]/[nT] *1.d3
pha  = atan2(imag(z),real(z))*r2d     ! [deg]
rhoa = dmu/omega*cdabs(z)**2. ! [Ohm.m]

return
end
!#########################################
!# E = sum of w [1/km] * El [mV/km * km] = [mV/km]
subroutine E_ele(g_surface,iele,y,z,e)
implicit none
real(8),       intent(in)  :: y,z ! obs coordinate
type(surface), intent(in)  :: g_surface
integer(4),    intent(in)  :: iele
complex(8),    intent(out) :: e(2)  ! electric field [mV/km]
real(8)    :: yt(3),zt(3),w(2)
integer(4) :: i,j,ipoi,id(2,3),iline,idirection,l,m
real(8)    :: lambda(3),gn(2,3),elm_yz(2,3)

do i=1,3 ! node loop
 ipoi =  g_surface%n3(iele,i)
 elm_yz(1:2,i) = g_surface%x1x2_face(1:2,ipoi)
end do
call nodebasisfun_tri(elm_yz,y,z,lambda)
!write(*,*) "lambda",lambda(1:3)
call gradnodebasisfun_tri(elm_yz,gn)
!write(*,*)"gn1",gn(1:2,1)
!write(*,*)"gn2",gn(1:2,2)
!write(*,*)"gn3",gn(1:2,3)

id(1,1:3)=(/1,2,3/) ! l
id(2,1:3)=(/2,3,1/) ! m

e=0.d0
do j=1,3 ! line loop
 iline      = abs(g_surface%n3line(iele,j))
 idirection = 1
 if (g_surface%n3line(iele,j) .lt. 0.) idirection = -1
 l = id(1,j)
 m = id(2,j)
 w = lambda(l)*gn(1:2,m) - lambda(m)*gn(1:2,l)
! write(*,*)"l,m,",l,m,"w",w
! write(*,*)"iline",iline
! write(*,*) "bs",g_surface%bs(iline)
 e(1:2) = e(1:2) + w(1:2)*g_surface%bs(iline)*idirection
end do

return
end


!#########################################
!# B = (i/omega)*rot E : [V/  m]/m  -> [T]
!#                       [mV/km]/km -> [nT]
!# E = sum of w [1/km] * El [mV/km * km] = [mV/km]
!# B = (i/omega)*sum of grcgr [1/km^2] * El [mV/km * km] = [mV/km]/km -> [nT]
subroutine B_ele(g_surface,omega,iele,y,z,bx) ! bx [nT]
implicit none
type(surface),intent(in)  :: g_surface
real(8),      intent(in)  :: y,z
real(8),      intent(in)  :: omega
complex(8),   intent(out) :: bx
complex(8)                :: iunit=(0.d0,1.d0),bxn(3)
real(8)                   :: grcgr_i(3) ! 1/(km)^2
real(8)                   :: elm_yz(2,3)
integer(4)                :: i,ipoi,iline,idirection,i_edge,iele

do i=1,3 ! node loop
 ipoi =  g_surface%n3(iele,i)
 elm_yz(1:2,i) = g_surface%x1x2_face(1:2,ipoi)
end do

bxn=0.d0
do i_edge=1,3 ! line loop
 iline      = abs(g_surface%n3line(iele,i_edge))
 idirection = g_surface%n3line(iele,i_edge)/iline
 call grcgr(i_edge,elm_yz(1,:),elm_yz(2,:),grcgr_i)
 bxn = bxn + 2.*iunit/omega*grcgr_i(:)*g_surface%bs(iline)*idirection ! El [mV/km * km] * [1/km2]-> nT
end do

bx = bxn(1)

return
end

!#########################################
subroutine gradnodebasisfun_tri(elm_yz,gn) ! gn [1/km]
use outerinnerproduct
implicit none
real(8),intent(in)  :: elm_yz(2,3)
real(8),intent(out) :: gn(2,3) ! (y,z) vector for node 1 to 3
real(8) :: ex(3),s,elm_xyz(3,3),gn_n(3,3),xmn(3)
integer(4) :: id(3,3),l,m,n,i_edge
ex=(/ 1.d0,0.d0,1.d0/)
elm_xyz=0.d0
elm_xyz(2:3,:)=elm_yz(1:2,:)

id(1,1:3)=(/1,2,3/)
id(2,1:3)=(/2,3,1/)! m
id(3,1:3)=(/3,1,2/)! n

call area_tri(elm_yz(1,1:3),elm_yz(2,1:3),s)

do l=1,3
 m=id(l,2)
 n=id(l,3)
 xmn = elm_xyz(:,n) - elm_xyz(:,m)
 gn_n(1:3,l)=1./2./s * outer(ex,xmn)
end do

gn(1:2,1:3) = gn_n(2:3,1:3)

return
end
!############################################
subroutine nodebasisfun_tri(elm_yz,y,z,lambda)
implicit none
real(8),intent(in) :: elm_yz(2,3) ! triangle coordinate
real(8),intent(in) :: y,z ! obs coordinate
real(8),intent(out) :: lambda(3)
real(8) :: s,s1,s2,s3,y1(3),z1(3)

call area_tri(elm_yz(1,:),elm_yz(2,:),s)
! s3
y1=(/ elm_yz(1,1), elm_yz(1,2),y/) ; z1=(/elm_yz(2,1),elm_yz(2,2),z/)
call area_tri(y1,z1,s3)
! s1
y1=(/ elm_yz(1,2), elm_yz(1,3),y/) ; z1=(/elm_yz(2,2),elm_yz(2,3),z/)
call area_tri(y1,z1,s1)
! s2
y1=(/ elm_yz(1,1), elm_yz(1,3),y/) ; z1=(/elm_yz(2,1),elm_yz(2,3),z/)
call area_tri(y1,z1,s2)

lambda(1)=s1/s
lambda(2)=s2/s
lambda(3)=s3/s

return
end
!############################################
! coded on 2021.05.25
! gradient cross gradient for triangle
subroutine grcgr(i_edge,y,z,grcgr_i)
use outerinnerproduct
implicit none
integer(4),intent(in)  :: i_edge
real(8),   intent(in)  :: y(3),z(3)
real(8),   intent(out) :: grcgr_i(3)
real(8)    :: r(3,3),s,rnl(3),rmn(3)
integer(4) :: id(3,3),l,m,n,i

id(1,1:3)=(/1,2,3/) ! l,m,n for 1 st edge
id(2,1:3)=(/2,3,1/) ! l,m,n for 2 nd edge
id(3,1:3)=(/3,1,2/) !

l=id(i_edge,1)
m=id(i_edge,2)
n=id(i_edge,3)

r(1,1:3)=(/ 0.d0, y(1), z(1) /)
r(2,1:3)=(/ 0.d0, y(2), z(2) /)
r(3,1:3)=(/ 0.d0, y(3), z(3) /)
rnl =r(l,:) - r(n,:)
rmn =r(n,:) - r(m,:)

! grcgr is in x direction
call area_tri(y,z,s)
grcgr_i(1:3) = -1./4./(s**2.)*outer(rnl,rmn)

return
end
!#############################################
! 2021.05.25
! assume closs product of x1 and x2 point on the side of direction of x3
! volume can be minus when (x1 times x2) cdot x3 < 0
subroutine area_tri(y,z,s)
use outerinnerproduct
implicit none
real(8),intent(in) ::y(3),z(3)
real(8),intent(out) :: s

s=1./2.*abs((y(3)-y(1))*(z(2)-z(1)) - (y(2)-y(1))*(z(3)-z(1)))

return

end subroutine area_tri


!############################################# 2020.10.28
! integral of shape function over triangle
function ints(k, m, s)
implicit none
integer(4),intent(in) :: k, m
real(8),intent(in) :: s ![km^2]
real(8) :: ints
! here assume k = l
! intv =6v*(k!l!m!n!)/(k+l+m+n+3)!
if ( k .ne. m ) ints=s/6.d0 ![km^3]
if ( k .eq. m ) ints=s/12.d0 ![km^3]
return
end function ints

!###########################################
! coded on 2021.05.25
function wiwjdS(i_edge,j_edge,y,z)
implicit none
real(8),   intent(in) :: y(3),z(3)
integer(4),intent(in) :: i_edge,j_edge
integer(4)            :: id(3,3),l,m,n,lp,mp
real(8)               :: wiwjdS

id(1,1:3)=(/1,2,3/)
id(2,1:3)=(/2,3,1/)
id(3,1:3)=(/3,1,2/)

l =id(i_edge,1)
m =id(i_edge,2)
lp=id(j_edge,1)
mp=id(j_edge,2)

wiwjdS= NlNmdS(l,lp,y,z)*grdgr(m,mp,y,z) &! see m_surface_type.f90
&     - NlNmdS(m,lp,y,z)*grdgr(l,mp,y,z) &
&     - NlNmdS(l,mp,y,z)*grdgr(m,lp,y,z) &
&     + NlNmdS(m,mp,y,z)*grdgr(l,lp,y,z)
return
end
!############################################# 2021.05.31

!#########################################
! 2021.05.31
function NlNmdS(l,m,y,z)
implicit none
real(8),   intent(in) :: y(3),z(3)
integer(4),intent(in) :: l,m
real(8)    :: s,NlNmdS

call area_tri(y,z,s)
if ( l .eq. m ) NlNmdS = 1./6. *s
if ( l .ne. m ) NlNmdS = 1./12. *s

return
end

!##########################################
! coded on 2021.05.25
! gradient lambda dot gradient dot lambda for triangle
function grdgr(i,j,y,z)
use outerinnerproduct
implicit none
real(8),   intent(in)  :: y(3),z(3)
integer(4),intent(in)  :: i,j
real(8)    :: grdgr
real(8)    :: rnmi(3),rnmj(3)
real(8)    :: elm_yz(2,3),gn(2,3)

elm_yz(1,1:3)=y(1:3)
elm_yz(2,1:3)=z(1:3)

call gradnodebasisfun_tri(elm_yz,gn)

rnmi(1:3)=(/ 0.d0, gn(1,i), gn(2,i) /)
rnmj(1:3)=(/ 0.d0, gn(1,j), gn(2,j) /)
grdgr = inner(rnmi,rnmj)

return
end

end module surface_type

