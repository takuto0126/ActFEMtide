! Coded on 2015.09.17
! to generate .geo file for simple 2D mesh with suitable .pos file
! from topo file around Nakadake
!
! Input file
! aso_caldera.dat
! topo127_134_29_36.xyz
!
! Output file
!  header.geo
!  header.pos
!
! Note x is eastward, y is northward, z is upward
program meshgen1
use param
use horizontalresolution
implicit none
type(param_forward)    :: g_param
type(param_source)     :: s_param
real(8),dimension(2)   :: xy1, xy2, xy, v12, v0, dv
integer(4)             :: ifile=1, ndat,iline
integer(4)             :: ipoint,ii,j,ishift, ispline, ilineloop
character(50)          :: posfile,header2d ! 2017.09.28
real(8)                :: si, sb  ! size of triangles [km]
real(8)                :: x1,x2,y1,y2,x01,x02,y01,y02, x0,y0, z0,lc
real(8),dimension(4,2) :: lonlat4
integer(4),allocatable,dimension(:,:) :: line12

![0]## read param
 call readparam(g_param,s_param) ! 2018.06.14

!#[1]# set
 si   = g_param%sizein
 sb   = g_param%sizebo
 x01 = g_param%xbound(1)
 x1  = g_param%xbound(2)
 x2  = g_param%xbound(3)
 x02 = g_param%xbound(4)
 y01 = g_param%ybound(1)
 y1  = g_param%ybound(2)
 y2  = g_param%ybound(3)
 y02 = g_param%ybound(4)
 lc  = g_param%sizebo
 header2d = g_param%header2d

![1]## generate .geo file
![1-1] Boundary points and lines
 open(ifile,file=header2d(1:len_trim(header2d))//".geo")
 write(ifile,*) "lc=",lc ,";"
 write(ifile,*) "Point(1)={",x02,",",y02,",0.0,lc} ;"! top right
 write(ifile,*) "Point(2)={",x02,",",y01,",0.0,lc} ;"! top left
 write(ifile,*) "Point(3)={",x01,",",y01,",0.0,lc} ;"! bottom left
 write(ifile,*) "Point(4)={",x01,",",y02,",0.0,lc} ;"! bottom right
 write(ifile,*) "Line(1)={1,2} ;"
 write(ifile,*) "Line(2)={2,3} ;"
 write(ifile,*) "Line(3)={3,4} ;"
 write(ifile,*) "Line(4)={4,1} ;"
 write(ifile,*) "Line Loop(1)={1,2,3,4} ;"
 write(ifile,*) "Plane Surface(1)={1};"
 close(ifile)

![2]## generate .pos file
 posfile=header2d(1:len_trim(header2d))//".pos"
 call calobsr(s_param,g_param) ! calculate sigma_r, A_r, nobsr
 call outbgmesh2d(posfile,g_param) ! outbgmesh2d.f90 2017.09.28
! call outbgfield2d(g_param) ! outbgfield2d.f90 2017.09.28

end program meshgen1
