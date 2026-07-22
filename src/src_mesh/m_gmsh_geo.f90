!Coded on March 9, 2016

module gmsh_geo
implicit none

type point_info
 !# info of node
 integer(4) :: npoint
 real(8),dimension(:,:),allocatable :: xyz
 real(8),dimension(:),  allocatable :: lc    ! size of element
end type point_info

!# info of line
type line_info
 integer(4) :: nline
 integer(4) :: max_npoint
 integer(4),dimension(:),  allocatable :: linetype ! 1 for line, 2 for spline
 integer(4),dimension(:),  allocatable :: npoint  ! # of point for each line
integer(4),dimension(:,:), allocatable :: node_id ! 2 nodes for line, 5 nodes for spline
end type line_info

!# info of lineloop
type lineloop_info
 integer(4) :: nlineloop ! # of lineloops
 integer(4) :: max_nline
 integer(4),dimension(:),  allocatable :: nline
 integer(4),dimension(:,:),allocatable :: line_id ! j-th line id of i-th lineloop
end type lineloop_info

!# info of plane surface
type plane_info
 integer(4) :: nplane ! # of plane surface
 integer(4) :: max_nlineloop
 integer(4),dimension(:),  allocatable :: nlineloop
 integer(4),dimension(:,:),allocatable :: lineloop_id ! j-th lineloop id of i-th
end type plane_info

!# info of plane loop
type planeloop_info
 integer(4) :: nplaneloop
 integer(4) :: max_nplane
 integer(4),dimension(:),  allocatable :: nplane
 integer(4),dimension(:,:),allocatable :: plane_id ! j-th lineloop id of i-th
end type planeloop_info

!# info of volume
type volume_info
 integer(4) :: nvolume
 integer(4) :: max_nplaneloop
 integer(4), dimension(:),  allocatable :: nplaneloop
 integer(4), dimension(:,:),allocatable :: planeloop_id
end type volume_info

type geo_info_3d
 type(point_info)     :: point
 type(line_info)      :: line
 type(lineloop_info)  :: lineloop
 type(plane_info)     :: plane
 type(planeloop_info) :: planeloop
 type(volume_info)    :: volume
end type

contains
!###################################################  INITPOINTS
subroutine INITPOINT(point1,npoint)
implicit none
type(point_info),intent(inout) :: point1
integer(4),      intent(in)    :: npoint

point1%npoint = npoint
allocate( point1%xyz(3,npoint) )
allocate( point1%lc(npoint))

return
end subroutine
!###################################################  INITLINES
subroutine INITLINE(line1,nline,max_npoint)
implicit none
integer(4),     intent(in)    :: nline,max_npoint
type(line_info),intent(inout) :: line1

 line1%nline      = nline
 line1%max_npoint = max_npoint
 allocate( line1%linetype(nline) )
 allocate( line1%npoint(nline)   )
 allocate( line1%node_id(max_npoint,nline))
 line1%node_id(:,:)=0

return
end subroutine
!###################################################  INITLINELOOP
subroutine INITLINELOOP(lineloop1,nlineloop,max_nline)
implicit none
type(lineloop_info),intent(inout) :: lineloop1
integer(4),         intent(in)    :: nlineloop,max_nline

 lineloop1%nlineloop = nlineloop
 lineloop1%max_nline = max_nline
 allocate( lineloop1%nline(nlineloop) )
 allocate( lineloop1%line_id(max_nline,nlineloop) )

return
end subroutine
!###################################################  INITPLANE
subroutine INITPLANE(plane1,nplane,max_nlineloop)
implicit none
integer(4),      intent(in)    :: nplane,max_nlineloop
type(plane_info),intent(inout) :: plane1

 plane1%nplane = nplane
 plane1%max_nlineloop = max_nlineloop
 allocate( plane1%nlineloop(nplane) )
 allocate( plane1%lineloop_id(max_nlineloop,nplane) )

return
end subroutine
!###################################################  INITPLANELOOP
subroutine INITPLANELOOP(planeloop1,nplaneloop,max_nplane)
implicit none
integer(4),          intent(in)    :: nplaneloop,max_nplane
type(planeloop_info),intent(inout) :: planeloop1

 planeloop1%nplaneloop = nplaneloop
 planeloop1%max_nplane = max_nplane
 allocate( planeloop1%nplane(nplaneloop) )
 allocate( planeloop1%plane_id(max_nplane,nplaneloop) )

return
end subroutine

!###################################################  INITVOLUME
subroutine INITVOLUME(volume1,nvolume,max_nplaneloop)
implicit none
integer(4),       intent(in)    :: nvolume,max_nplaneloop
type(volume_info),intent(inout) :: volume1

 volume1%nvolume        = nvolume
 volume1%max_nplaneloop = max_nplaneloop
 allocate( volume1%nplaneloop(nvolume) )
 allocate( volume1%planeloop_id(max_nplaneloop,nvolume) )

return
end subroutine

!###################################################  OUTGEOINFO3D(outgeo1,geo3d)
subroutine OUTGEOINFO3D(outgeo,geo3d)
implicit none
character(50),    intent(in) :: outgeo ! 2017.09.08
type(geo_info_3d),intent(in) :: geo3d
integer(4)                   :: i,j
integer(4)                   :: ifile=11
type(point_info)             :: point
type(line_info)              :: line
type(lineloop_info)          :: lineloop
type(plane_info)             :: plane
type(planeloop_info)         :: planeloop
type(volume_info)            :: volume
 
point   = geo3d%point ; plane     = geo3d%plane
line    = geo3d%line  ; planeloop = geo3d%planeloop
lineloop= geo3d%lineloop ; volume = geo3d%volume

open(ifile,file=outgeo)

!#[0]## copy

!---------------------------------------------------------------- output
!#[1]# point
do i=1,point%npoint
  write(ifile,*) "Point(",i,")={",&
  & point%xyz(1,i),",",&
  & point%xyz(2,i),",",&
  & point%xyz(3,i),",", point%lc(i),"};"
end do

!#[2]## Lines
do i=1,line%nline
 if (line%linetype(i) .eq. 1) then ! normal line
  write(ifile,*) "Line(",i,")={",line%node_id(1,i),",",line%node_id(2,i),"};"
!  write(*,*) "Line(",i,")={",line%node_id(1,i),",",line%node_id(2,i),"};"
 else if (line%linetype(i) .eq. 2) then ! spline
  write(ifile,*) "Spline(",i,")={",&
                 & (line%node_id(j,i),",",   j=1,line%npoint(i)-1),&
                 &  line%node_id(line%npoint(i),i),"};"
!  write(*,*) "Spline(",i,")={",&
!                 & (line%node_id(j,i),",",   j=1,line%npoint(i)-1),&
!                 &  line%node_id(line%npoint(i),i),"};"
 end if
end do

!#[3]## Lineloops and Plane Surface
do i=1,lineloop%nlineloop
  write(ifile,*)"Line Loop(",i,")={",&
                 & (lineloop%line_id(j,i),",",j=1,lineloop%nline(i)-1),&
                 &  lineloop%line_id(lineloop%nline(i),i),"};"
end do

!#[4]## Plane SUrface
do i=1,plane%nplane
  write(ifile,*)"Plane Surface(",i,")={",&
                 & (plane%lineloop_id(j,i),",",j=1,plane%nlineloop(i)-1),&
                 &  plane%lineloop_id(plane%nlineloop(i),i),"};"
end do

!#[4]## Surfaceloops
do i=1,planeloop%nplaneloop
  write(ifile,*)"Surface loop(",i,")={",&
                 & (planeloop%plane_id(j,i),",",j=1,planeloop%nplane(i)-1),&
                 &  planeloop%plane_id(planeloop%nplane(i),i),"};"
end do

!#[5]## Volumes
do i=1,volume%nvolume
  write(ifile,*)"Volume(",i,")={",&
                 & (volume%planeloop_id(j,i),",",j=1,volume%nplaneloop(i)-1),&
                 &  volume%planeloop_id(volume%nplaneloop(i),i),"};"
end do
!---------------------------------------------------------------- output
close(ifile)
return
end subroutine

end module gmsh_geo


