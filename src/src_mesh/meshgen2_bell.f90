! Coded on 2015.09.17 by T. MINAMI
! to generate nakadake2d.geo reflecting real topography
!
program meshgen2_bell
use param
use bell  ! 2017.09.26
use mesh_type
implicit none
type(mesh)          :: h_mesh
type(mesh)          :: hz_mesh
type(param_forward) :: g_param
type(param_source)  :: s_param
type(param_cond)    :: g_cond  ! 2017.09.28
type(param_bell)    :: b_param
character(50)       :: mshfile2d, mshfile2d_z, mshfile2dz
character(50)       :: header2d, header3d, posfile
integer(4)          :: node, ecount, inib,lcount, n3dn
! n2(j,1) : the initial line number that j-th line belongs to
real(8),allocatable,dimension(:) ::  x3d,y3d,z3d
integer(4)                       :: itopoflag ! 2017.09.29
integer(4)                       :: nlinbry(4),linbry(4,100)
! nlinbry(i) is # of lines on i-th calculation boundary, linbry(i,j) is line number of j-th line of i-th boundary

!#[0]##
call READPARAM(g_param,s_param,g_cond) ! 2017.09.28
call READPARAMBELL(g_param,b_param)    ! 2017.09.28

!#[1]## set
 header2d    = g_param%header2d
 header3d    = g_param%header3d
 mshfile2d   = header2d(1:len_trim(header2d))//".msh"
 mshfile2d_z = header2d(1:len_trim(header2d))//"_z.msh"
 mshfile2dz  = header2d(1:len_trim(header2d))//"z.msh"
 posfile     = header3d(1:len_trim(header3d))//".pos"
 itopoflag   = g_param%itopoflag  ! 2017.09.29

!#[2]## read 2d mesh
 CALL READMESH_TOTAL(h_mesh,mshfile2d) ! see m_mesh_type.f90
 CALL GENXYZMINMAX(h_mesh,g_param)     ! 2017.09.28

 node = h_mesh%node
 lcount = h_mesh%nlin
 ecount = h_mesh%ntri

!#[3]## read topography
 if ( itopoflag .eq. 1 ) then ! 2017.09.29
  call calztopo2(h_mesh, g_param) ! z will be modified
 end if  ! 2017.09.29
 call polygonz3(mshfile2d,mshfile2dz,mshfile2d_z,h_mesh%xyz(3,:),node)
 !open(1,file=mshfile2d_z)
 !call MESHOUT(1,h_mesh)
 !close(1)

!#[4]## identify boundary lines
 call findbryline7(h_mesh,nlinbry,linbry)
 n3dn=node+8
 allocate(x3d(n3dn),y3d(n3dn),z3d(n3dn))

!#[4]## Prepare node
 call prepare7_5(h_mesh, x3d, y3d, z3d, n3dn, linbry, g_param)

!#[5]## generate 3d geometry file
 call outpregeo8(h_mesh, x3d, y3d, z3d, nlinbry, linbry, g_param,b_param)

!#[6]## out 3d posfile
 call calobsr(s_param,g_param) ! cal g_param%nobsr, xyz_r,sigma_r,A_r (m_param.f90)
 call calzobsr(h_mesh,g_param) ! calzobsr.f90 2017.09.28
! call outbgfield3d(g_param)    ! outbgfield.f90 2017.09.28
 call outbgmesh3d(posfile,g_param)

end program meshgen2_bell
!###################################################################
! modified for spherical on 2016.11.20
! iflag = 0 for xyz
! iflag = 1 for xyzspherical
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

!############################################ calztopo2
!Copied from extrude.f90 and modified on Sep. 17, 2015
!# modify to interpolate in lat and lon, then convert to UTM coordinate on Dec. 17, 2015
!# modified for multiple topofraphy files on 2017.09.26
!# this program assumes eariler files are narrow compared to later ones 2017.09.26
subroutine calztopo2(h_mesh,g_param)
use mesh_type
use param
implicit none
type(param_forward) :: g_param
type(mesh),intent(inout) :: h_mesh
real(8),dimension(h_mesh%node) :: x, y, z
integer(4)    :: node
real(8)       :: calz, lon_c, lat_c, x_c, y_c
real(8)       :: rlat,rlon,x_utm,y_utm
real(8), dimension(h_mesh%node) :: lon,lat
real(8), allocatable,dimension(:,:) :: lon1,lat1, z1
character(50),allocatable,dimension(:)   :: files           ! 2017.09.25
integer(4),   allocatable,dimension(:)   :: nt,nsouth,neast ! 2017.09.25
real(8),      allocatable,dimension(:,:) :: lonlatshift     ! 2017.09.26
real(8)       :: cx0_utm,cy0_utm
real(8)       :: lonorigin, latorigin
character(3)  :: UTM
character(120):: line
character(50) :: header2d,header3d,pregeo,topofile
integer(4)    :: ifile, ntmax
integer(4)    :: ntopo, i, j, ii, k, neast1,nsouth1,nfile
real(8)       :: zz1,zz2,zz3,zz4             ! 2017.09.26
real(8)       :: sbound,nbound,wbound,ebound ! 2017.09.26
real(8)       :: lonw,lone,latn,lats         ! 2017.09.26

!#[0]## set
nfile       = g_param%nfile           ! 2017.09.25
allocate(files(nfile))                ! 2017.09.25
allocate(nt(nfile),neast(nfile),nsouth(nfile)) ! 2017.09.25
allocate(lonlatshift(2,nfile) )       ! 2017.09.26
files       = g_param%topofile        ! 2017.09.25
lonlatshift = g_param%lonlatshift     ! 2017.09.26
node        = h_mesh%node
header2d    = g_param%header2d
header3d    = g_param%header3d
x           = h_mesh%xyz(1,:)
y           = h_mesh%xyz(2,:)
z           = h_mesh%xyz(3,:) ! will be modified and set in h_mesh again
UTM         = g_param%UTM   ! zone
lonorigin   = g_param%lonorigin
latorigin   = g_param%latorigin

!# [1] ### count lines of gebcofile and allocate
ntmax = 0     ! 2017.09.25
do ifile=1,nfile
  ntopo=0
  open(1,file=files(ifile))
  do while (ntopo .ge. 0)
    read(1,*,end=99)
    ntopo=ntopo+1 ! ntopo is # of lines in gebcofile
  end do
  99 nt(ifile)=ntopo
  close(1)
  write(*,*) "file=",files(ifile),"nt(ifile)=",nt(ifile)
  ntmax=max(ntmax, nt(ifile)) ! 2017.09.25
end do
allocate( lon1(ntmax,nfile),lat1(ntmax,nfile)) ! 2017.09.25
allocate( z1(ntmax,nfile) ) ! x1(i,j) is x coord of j-th node of i-th file

!# [2] ### read coordinates in gebcofile
do ifile=1,nfile
  open(1,file=files(ifile))
  write(*,*) "lonlatshift=",lonlatshift(1:2,ifile),"ifile",ifile
  do i=1,nt(ifile)
    read(1,*) lon1(i,ifile),lat1(i,ifile), z1(i,ifile)   ! 2017.09.25
    lon1(i,ifile) = lon1(i,ifile) + lonlatshift(1,ifile) ! 2017.09.25
    lat1(i,ifile) = lat1(i,ifile) + lonlatshift(2,ifile) ! 2017.09.25
  end do
  z1(:,ifile)=z1(:,ifile)/1000.d0 ! [m] -> [km]
  close(1)
end do

!# [3] ### change order and measure # of nodes in horizontal and vertical directions
!# reorder the values to ((west -> east), north -> south )
do i=1,nfile ! file loop
  write(*,*) "lon1(1:2,ifile)=",lon1(1:2,i),"ifile",i
  call changeorder(lon1(1:nt(i),i),lat1(1:nt(i),i),z1(1:nt(i),i), nt(i), nsouth(i), neast(i))
  write(*,*) "file=",files(i)
  write(*,*) "nsouth=",nsouth(i),"neast=",neast(i),"ntopo=",nt(i)
end do

!# [4] ### calculate depth "z" at all the horizontal nodes
rlon=lonorigin !copy because lonorigin and latorigin are parameters
rlat=latorigin
CALL UTMGMT(rlon,rlat,x_utm,y_utm,UTM,0) ! 0:LONLAT2UTM , 1:UTM2LONLAT
cx0_utm=x_utm
cy0_utm=y_utm
write(*,*) "(cx0_utm,cy0_utm)=(",cx0_utm,",",cy0_utm,")"

CALL UTMGMT_N(node,x(:)*1.d3+cx0_utm,y(:)*1.d3+cy0_utm,lon,lat,UTM,1) ! get LON and LAT for nodes
write(*,*) "### UTMGMT_N END!! ###"

z(:)=0.d0
do k=1,node ! node is total # of nodes of 2d mesh
  ! provided narrow area to broad area
  do ifile=1,nfile   ! broad topo file
    neast1=0; nsouth1=0  ! 2017.09.26
   !# lon1(2,1) is north west corner, while lon1(2,nt(2)) is south east corner
    wbound = lon1(1, ifile) ;    ebound = lon1(nt(ifile),ifile)
    nbound = lat1(1, ifile) ;    sbound = lat1(nt(ifile),ifile)
   if ( wbound .lt. lon(k) .and. lon(k) .lt. ebound .and. &
   &    sbound .lt. lat(k) .and. lat(k) .lt. nbound ) then   ! within the area for ifile
    do i=1,nsouth(ifile)-1 ! latitudial location
     if ( lat1( i *neast(ifile) + 1, ifile) .le. lat(k) .and. &
      & lat(k) .le. lat1((i-1)*neast(ifile) + 1, ifile) )           nsouth1=i ! south location
    end do
    do j=1,neast(ifile)-1 ! longitudinal location
     if ( lon1(j,ifile) .le. lon(k) .and. lon(k).le. lon1(j+1,ifile)) neast1=j !north south location
    end do
    if ( neast1 .eq. 0 .or. nsouth1 .eq. 0 ) then
     write(*,*) "GEGEGE stop! cannot find lon",lon(k),"lat",lat(k),"in file",ifile
     stop
    end if
    ii   = (nsouth1 - 1)*neast( ifile) + neast1 !       2017.09.26
    zz1  = z1(ii,               ifile) ! top left       2017.09.26
    zz2  = z1(ii+1,             ifile) ! top right      2017.09.26
    zz3  = z1(ii+neast(ifile),  ifile) ! bottom left    2017.09.26
    zz4  = z1(ii+neast(ifile)+1,ifile) ! bottom right   2017.09.26
    lonw = lon1(ii,             ifile) ! west lon       2017.09.26
    lone = lon1(ii+1,           ifile) ! east lon       2017.09.26
    latn = lat1(ii,             ifile) ! north lat      2017.09.26
    lats = lat1(ii+neast(ifile),ifile) ! south lat      2017.09.26
    z(k)=calz(lon(k),lat(k), lonw, lone, latn, lats, zz1,zz2,zz3,zz4) ! 2017.09.26
    if (z(k)  .lt. 0.d0 ) z(k)=0.d0
    if ( .not. ( z(k) .ge. 0.d0 ) ) then
     write(*,*) k,"lon(k),lat(k)",lon(k),lat(k),"z=",z(k),"ifile=",ifile
     write(*,*) "lonw,lone",lonw,lone,"lats,latn",lats,latn
     write(*,*) "zz1,zz2,zz3,zz4=",zz1,zz2,zz3,zz4
     stop
    end if
    goto 100  ! 2017.09.26
   end if     ! 2017.09.26  if within area
 end do       ! file loop end
 100 continue ! 2017.09.26
end do        ! node loop end

! output topo file
open(1,file=header2d(1:len_trim(header2d))//".msh")
open(2,file=header2d(1:len_trim(header2d))//"_lonlatz.msh")
!open(3,file="aso_caldera/aso_caldera_triwgs.dat")
open(3,file="lonlat.dat")
write(*,*) "check1"
do i=1,5
 read(1,'(a)') line
 write(2,'(a)') line
end do
write(*,*) "check1"
do k=1,node
 read(1,*)
 write(2,'(i10,3g15.7)') k,lon(k),lat(k),z(k)/10.D0
 write(3,'(3g15.7)') lon(k),lat(k),z(k)*1.d3
end do
write(*,*) "check2"
do
read(1,'(a)',end=98) line
write(2,'(a)') line
end do
98 continue
close(1)
close(2)
close(3)

!#[]## set output
h_mesh%xyz(3,:) = z

write(*,*) "### CALZTOPO2 END!! ###"

return
999 continue
write(*,*) "GEGEGE when node # is k=",k,"x(k),y(k)=",x(k),y(k)
write(*,*) "neast1=",neast1,"nsouth1=",nsouth1
stop
end subroutine calztopo2

!######################################################  function calz
function calz(x,y,x1,x2,y1,y2,z1,z2,z3,z4) ! note y1 > y > y2
!     x1<x<x2
!z1 ->  ----  y1 <- z2
!      |    | y
!z3 ->  ----  y2 <- z4
implicit none
real(8) :: calz,x,y,x1,x2,y1,y2,z1,z2,z3,z4,aa,a1,a2,a3,a4
aa=(x2-x1)*(y1-y2)
a1=(x2-x)*(y-y2)/aa !
a2=(x-x1)*(y-y2)/aa
a3=(x2-x)*(y1-y)/aa
a4=(x-x1)*(y1-y)/aa
calz=a1*z1+a2*z2+a3*z3+a4*z4
end function calz
!#######################################################
! Coded on Sep. 17, 2015
! change order and calculate n_south and n_east
! to ((west -> east), north -> south)
subroutine changeorder(x1,y1,z1, ntopo, n_south, n_east)
implicit none
integer(4), intent(in)    :: ntopo
real(8),    intent(inout) :: x1(ntopo),y1(ntopo),z1(ntopo)
integer(4), intent(out)   :: n_south, n_east
real(8)                   :: x2(ntopo),y2(ntopo),z2(ntopo)
integer(4)                :: i, j, ishift, jshift

if ( x1(1) .eq. x1(2) ) then ! [1] if the order is based on colmun
  write(*,*) "x1(1)=",x1(1),"x1(2)=",x1(2)
  if ( y1(1) .lt. y1(2) ) then  ! [1-1] upward
   i=1
   do while ( y1(i) .lt. y1(i+1) )
     i=i+1
   end do
!   open(1,file="check.dat")
   n_south=i ; n_east=ntopo/n_south ! ok
   do i=1,n_south
      ishift=(i-1)*n_east
      do j=1, n_east
	  jshift=(j-1)*n_south
        x2(ishift+j)=x1(jshift + n_south - (i-1) )
	  y2(ishift+j)=y1(jshift + n_south - (i-1) )
	  z2(ishift+j)=z1(jshift + n_south - (i-1) )
!	  write(1,*) x2(ishift+j), y2(ishift+j), z2(ishift+j)
	end do
    end do
!    close(1)
!
  else if ( y1(2) .lt. y1(1) ) then ! [1-2] downward
   i=1
   do while ( y1(i) .gt. y1(i+1) )
     i=i+1
   end do
   n_south=i ; n_east=ntopo/n_south
   do i=1,n_south
      ishift=(i-1)*n_east
      do j=1, n_east
	  jshift=(j-1)*n_south
        x2(ishift+j)=x1(jshift + i )
	  y2(ishift+j)=y1(jshift + i )
	  z2(ishift+j)=z1(jshift + i )
	end do
   end do
!
  else                             ! [1-3]
    write(*,*) "GEGEGE! x1(1) = x1(2) and y1(1)=y2(2)"
    stop
  end if
  x1=x2
  y1=y2
  z1=z2
  goto 100
end if !=======  if ( x1(1) = x1(2) ) end!!
    i=1
    do while ( x1(i) .lt. x1(i+1) )
     i=i+1
   end do
   n_east=i ; n_south=ntopo/n_east
100 return
end subroutine changeorder
!
!##################################################### polygonz3
subroutine polygonz3(mshfile,polyzfile,poly_zfile,z,nmax) ! output polyzfile (mshfile)
implicit none
real(8),dimension(nmax) :: z
integer(4) :: i,j,node,nmax
real(8) :: x1,y1,z1
character(50) ::  mshfile,polyzfile,poly_zfile
character(100) :: line
open(1,file=mshfile)
open(2,file=polyzfile)
open(3,file=poly_zfile)
do i=1,4
read(1,'(a100)') line
write(2,'(a100)') line
write(3,'(a100)') line
end do
read(1,*) node
write(2,*) node
write(3,*) node
do i=1,node
read(1,*) j,x1,y1,z1
write(2,'(i10,3g15.7)') j,x1,y1,z(i)*1.d0
write(3,'(i10,3g15.7)') j,x1,y1,z(i)*5.d0
end do
do while ( i .ge. 0)
read(1,'(a100)',end=99) line
write(2,'(a100)') line
write(3,'(a100)') line
end do
99 continue
close(1)
close(2)
close(3)
write(*,*) "### POLYGONZ3 END!! ###"
return
end subroutine polygonz3
!########################################## findbryline7
subroutine findbryline7(h_mesh,nlinbry,linbry)
use mesh_type
implicit none
type(mesh),intent(in)  :: h_mesh
integer(4),intent(out) :: nlinbry(4),linbry(4,100)
integer(4) :: nodek, nlink
integer(4) :: n2k(h_mesh%nlin,3)
! n2k(j,1) is the belonging old line # for line j, n2k(j,2:3) are the start and end node
real(8)   :: xyk(h_mesh%node,2),zk(h_mesh%node)
real(8),dimension(4) :: d=(/1.d0, -1.d0, -1.d0, -1.d0/)
! nlinbry(i) is # of lines on i-th calculation boundary
! linbry(i,j) is line number of j-th line of i-th boundary
! 1 - 4 calculation boundaries are east, south, west, and north boundaries, respectively.
integer(4) :: line,i,j,k,ii,jj,n1,n2,is,isnext,ilast

!#[1]## set
nodek    = h_mesh%node
nlink    = h_mesh%nlin
n2k(:,1) = h_mesh%n2flag(:,2)
n2k(:,2) = h_mesh%n2(:,1)
n2k(:,3) = h_mesh%n2(:,2)
xyk(:,1) = h_mesh%xyz(1,:)
xyk(:,2) = h_mesh%xyz(2,:)
zk(:)    = h_mesh%xyz(3,:)

!#
! make nlinbry and linebry
nlinbry(:)=0
ilast=0
is=1 ! starting node number, which is north east corner node
do i=1,4 ! 1: right; 2: bottom; 3: right; 4 : top boundary
write(*,*) "i=",i,"findbry"
  if ( i .eq. 1 .or. i .eq. 3 ) ii = 1 !   ( east, west)
  if ( i .eq. 2 .or. i .eq. 4 ) ii = 2 !   ( bottom, top)
  10 continue
  do j=1,nlink ! line loop  ; nlink is # of boundary lines in .msh file
    n1=n2k(j,2);n2=n2k(j,3) ! set the start and end node of j-th line
    if ( (n1 .eq. is .or. n2 .eq. is ) .and. abs(xyk(n1,ii)-xyk(n2,ii)) .lt. 1.d-4 .and. ilast .ne. j ) then
    ! search line start from is
    nlinbry(i)=nlinbry(i)+1
    linbry(i,nlinbry(i))=j
    isnext=n2
    if ( is .eq. n2)      linbry(i,nlinbry(i))=-j
    if ( is .eq. n2)      isnext=n1
    is=isnext
    ilast=j
    !write(*,*) "i=",i,"nlinbry(i)=",nlinbry(i)
    !if ( linbry(i,nlinbry(i)) .gt. 0) write(*,*) n1,"->",n2
    !if ( linbry(i,nlinbry(i)) .lt. 0) write(*,*) n2,"->",n1
    !write(*,*) "xyk(n1,ii),xyk(n2,ii), delt=",xyk(n1,ii),xyk(n2,ii),(xyk(n1,ii)-xyk(n2,ii)),1.d-5
    goto 10
    end if
    !if (n2k(j,1) .eq. 186)then
    !write(*,*) "n2k(j,2:3)",n2k(j,2),n2k(j,3)
    !write(*,*) "xyk(n1,ii),xyk(n2,ii), delt=",xyk(n1,ii),xyk(n2,ii),(xyk(n1,ii)-xyk(n2,ii))
    !end if
  end do
  ! if no lines selected for next linbry, move on next calculation boundary
end do

write(*,*) "pre3dbry.geo start!!"
!### write nodes
open(1,file="pre3dbry.geo")
write(1,*) "lc=500.0 ;"
do i=1,nodek
 write(1,*)"Point(",i,")={",xyk(i,1),",",xyk(i,2),",",zk(i),",lc} ;"
end do
!### write lines only on boundaries
jj=0
do i=1,4
 do j=1,nlinbry(i)
  jj=jj+1
  line=linbry(i,j)
  !  n1=n2k(abs(line),2);n2=n2k(abs(line),3)
  !  write(*,*) "i=",i,"j=",j,"line=",line,"x1,y1=",(xyk(n1,k),k=1,2),"x2,y2=",(xyk(n2,k),k=1,2)
  if ( line  .gt. 0) then
   write(1,*) "Line(",jj,")={",n2k(line,2),",",n2k(line,3),"} ;"
  else
   write(1,*) "Line(",jj,")={",n2k(-line,3),",",n2k(-line,2),"} ;"
  end if
 end do
end do
close(1)
write(*,*) "## pre3dbry.geo end! ##"

if ( is .ne. 1) goto 99
write(*,*) "### FINDBRYLINE7 END!! ###"
return
99 write(*,*) "GEGEGE boundary round loop cannot be achieved."
write(*,*) "The end of node is", isnext
write(*,*) "x,y=",xyk(isnext,1),xyk(isnext,2)
stop
end subroutine findbryline7
!########################################## prepare7_5
subroutine prepare7_5(h_mesh,x3d,y3d,z3d,n3dn,linbry,g_param)
use mesh_type
use param
implicit none
type(param_forward),intent(in) :: g_param
type(mesh),intent(in) :: h_mesh
integer(4),intent(in) :: n3dn, linbry(4,100)
real(8),dimension(n3dn),intent(out) :: x3d,y3d,z3d
integer(4) :: node,nlink, n2k(h_mesh%nlin,3)
real(8),dimension(h_mesh%node) :: x,y,z
real(8) :: zmax,zmin
integer(4) :: i,j,line,n1, n3d

!#[1]## set
node = h_mesh%node
x  = h_mesh%xyz(1,:)
y  = h_mesh%xyz(2,:)
z  = h_mesh%xyz(3,:)
nlink = h_mesh%nlin
n2k(:,1)   = h_mesh%n2flag(:,2)
n2k(:,2:3) = h_mesh%n2(:,1:2)
zmax = g_param%zmax
zmin = g_param%zmin

!! x3d, y3d, and z3d
n3d=node
write(*,*) "n3d=",n3d
do i=1,4
line=linbry(i,1)
if (line .gt. 0) n1=n2k(line,2)
if (line .lt. 0) n1=n2k(abs(line),3)
x3d(n3d+i)=x(n1)
y3d(n3d+i)=y(n1)
z3d(n3d+i)=zmax
x3d(n3d+4+i)=x(n1)
y3d(n3d+4+i)=y(n1)
z3d(n3d+4+i)=zmin
!write(*,*) "i=",i,"x3d-z3d=",x3d(i),y3d(i),z3d(i)
!write(*,*) "i=",n3d+4+i,"x3d-z3d=",x3d(n3d+4+i),y3d(n3d+4+i),z3d(n3d+4+i)
end do
x3d(1:n3d)=x(1:n3d) ! x,y,z are set in subroutine mknode4
y3d(1:n3d)=y(1:n3d)
z3d(1:n3d)=z(1:n3d)
!! x3d, y3d, and z3d is completed
write(*,*) "### PREPARE7_5 END!! ###"
return
end subroutine prepare7_5
!######################################### outpregeo8
! here assume in the calculation area, no ocean exist
subroutine outpregeo8(h_mesh,x3d,y3d,z3d,nlinbry,linbry,g_param,b_param)
use param
use bell  ! 2017.09.28
use mesh_type
implicit none
type (param_forward),intent(in) :: g_param
type (param_bell)   ,intent(in) :: b_param
type (mesh),intent(in)          :: h_mesh
integer(4),intent(in)           :: nlinbry(4),linbry(4,100)
real(8),dimension(h_mesh%node+8),intent(in) :: x3d,y3d,z3d
integer(4)                      :: n3k(h_mesh%ntri,3)
integer(4)                      :: nlink,ntrik,node,n2k(h_mesh%nlin,3)
! linbry(i,j) is positive/negative line # corresponding to j-th line
! along the i-th calculation boundary;
! 1 to 4-th boundaries are right, bottom, left, and top boundary
! negative # indicates opposite direction of the original line
integer(4) :: i,j,k,ii,jj,n3(3),line(3),linelist(h_mesh%ntri*6,3),toploop(3)
! linelist(i,j) is the starting node number (j=1), or, end node number (j=2), or,
!                       the belonging group (j=3), of the line of number i
integer(4) :: toplooplist(h_mesh%ntri),botlooplist(h_mesh%ntri)
integer(4) :: n23dlinbry(4,100,2),n23dlinbryb(4,100,2)
integer(4) :: lin,n1,n2,n0,n43(4,3),n4(4),n4b(4),iilast, n3d,n3dn
! n23dlinbry(i,j,k) is start node(k=1) or end node(k=2) of j-th line on the i-th boundary
integer(4),allocatable,dimension(:) :: linb
character(50) :: pregeo,header3d,bellgeofile
real(8)       :: x1,x2,y1,y2,lc,z0
integer(4)    :: nbell,ibell ! 2018.10.02

!#[1]## set
node        = h_mesh%node
header3d    = g_param%header3d
ntrik       = h_mesh%ntri
n3k         = h_mesh%n3
nlink       = h_mesh%nlin
n2k(:,1)    = h_mesh%n2flag(:,2)
n2k(:,2:3)  = h_mesh%n2(:,1:2)
lc          = g_param%sizebo
bellgeofile = b_param%bellgeofile


pregeo      = header3d(1:len_trim(header3d))//".geo"
n3d         = node
n3dn        = node+8
open(1,file=pregeo)
! Point Output
write(1,*) "lc=",lc,";"
do i=1,n3dn ! n3dn coordinates input, n3dn=n3d+8
    write(1,*)"Point(",i,")={",x3d(i),",",y3d(i),",",z3d(i),",lc};"
end do
! convert linbry to n23dlinbry
do i=1,4
  do j=1,nlinbry(i) ! # of lines in i-th boundary
  lin=linbry(i,j)
  if (lin .gt. 0) then
    n1=n2k(lin,2);  n2=n2k(lin,3)
  else
    n1=n2k(-lin,3); n2=n2k(-lin,2)
  end if
  n23dlinbry(i,j,1)=n1! land and ocean upper surface calculation boundary lines
  n23dlinbry(i,j,2)=n2
  n23dlinbryb(i,j,1)=n1 ! land and ocean bottom surface calculation boundary lines
  n23dlinbryb(i,j,2)=n2 ! 99999999 is a default number for ki22dsptr
end do
end do
write(*,*) "n23dlinbry and n23dlinbryb generated!"
! upper line, lineloop, and surface in air ####################################    
!#                           TOP TOP TOP  VOLUME
n4(1)=n3d+1; n4(2)=n3d+2
n4(3)=n3d+3; n4(4)=n3d+4
linelist(:,1:3)=0
write(1,*) "Line(1)={",n4(1),",",n4(2),"};";linelist(1,1:2)=(/n4(1),n4(2)/) ! top east
write(1,*) "Line(2)={",n4(2),",",n4(3),"};";linelist(2,1:2)=(/n4(2),n4(3)/) ! top south
write(1,*) "Line(3)={",n4(3),",",n4(4),"};";linelist(3,1:2)=(/n4(3),n4(4)/) ! top west
write(1,*) "Line(4)={",n4(4),",",n4(1),"};";linelist(4,1:2)=(/n4(4),n4(1)/) ! top north
write(1,*) "Line Loop( 1 )={1,2,3,4};"
write(1,*) "Plane Surface(1)={1};"
write(1,*) "Line(5)={",n4(1),",",n23dlinbry(1,1,1),"};"
                   linelist(5,1:2)=(/n4(1),n23dlinbry(1,1,1)/) ! start, end node, belonging
write(1,*) "Line(6)={",n4(2),",",n23dlinbry(2,1,1),"};"
                   linelist(6,1:2)=(/n4(2),n23dlinbry(2,1,1)/)
write(1,*) "Line(7)={",n4(3),",",n23dlinbry(3,1,1),"};"
                   linelist(7,1:2)=(/n4(3),n23dlinbry(3,1,1)/)
write(1,*) "Line(8)={",n4(4),",",n23dlinbry(4,1,1),"};"
                   linelist(8,1:2)=(/n4(4),n23dlinbry(4,1,1)/)
ii=8 ! line #
!#####################################################################
!### add lines for transmitter cables to linelist with belonging group of 1, while normal one is 0
do i=1, nlink
  if ( n2k(i,1) .eq. 9 ) then ! n2k(i,1) .eq. 90 means the transmitter cable (see meshgen1.f90)
    ii=ii+1
    write(1,*) "Line(",ii,")={",n2k(i,2),",",n2k(i,3),"};"
    linelist(1,1:2)=n2k(i,2:3) ! top east
    linelist(ii,3)=1
  end if
end do
write(1,*)  "Physical Line(1)={",(i,",",i=9,ii-1),ii,"};"
!####
! the node order is the same between nakadake2d.msh and nakadake3d.geo up to n3d
!#####################################################################
n43(1,1:3)=(/-6,-1,5/) ! top east surface
n43(2,1:3)=(/-7,-2,6/) ! top south surface
n43(3,1:3)=(/-8,-3,7/) ! top west surface
n43(4,1:3)=(/-5,-4,8/) ! top north surface
jj=1 ! line loop #
do i=1,4
  allocate(linb(nlinbry(i)))
  do j=1,nlinbry(i)
    call findline(n23dlinbry(i,j,1),n23dlinbry(i,j,2),linb(j),ii,linelist,ntrik,1)
  end do
  write(1,*) "Line Loop(",i+jj,")={",(n43(i,k),",",k=1,3) ! comma is the last
  write(1,*) (linb(k),",",k=1,nlinbry(i)-1),linb(nlinbry(i)),"};"
  write(1,*) "Plane Surface(",i+jj,")={",i+jj,"};"
  deallocate(linb)
end do
jj=5 ! line loop  & Plane Surface #
do i=1,ntrik
! top and bottom
n3(:)=n3k(i,:);line(:)=0
   call findline(n3(1),n3(2),line(1),ii,linelist,ntrik,1)
   call findline(n3(2),n3(3),line(2),ii,linelist,ntrik,1)
   call findline(n3(3),n3(1),line(3),ii,linelist,ntrik,1)
   toploop(1:3)=(/line(1),line(2),line(3)/)
   jj=jj+1
   write(1,*) "Line Loop(",jj,")={",line(1),",",line(2),",",line(3),"};"
   write(1,*) "Plane Surface(",jj,")={",jj,"};"
   toplooplist(i)=jj
   botlooplist(i)=jj
end do
write(1,*) "Surface loop(1)={1,2,3,4,5",(",",toplooplist(i),i=1,ntrik),"};"
write(1,*) "Volume(1)={1};"
write(1,*) "Physical Volume (1) = {1} ;"
!#############################################################
!                                  BOTTOM BOTTOM VOLUME
n4(1)=n3d+5; n4(2)=n3d+6
n4(3)=n3d+7; n4(4)=n3d+8
n4b(1)=n23dlinbryb(1,1,1); n4b(2)=n23dlinbryb(2,1,1) ! n4b(i) is the node number for ocean bottom node at the i-th corner
n4b(3)=n23dlinbryb(3,1,1); n4b(4)=n23dlinbryb(4,1,1)
! bottom line, lineloop, and surface in air
write(1,*) "Line(",ii+1,")={",n4b(1),",",n4(1),"};";linelist(ii+1,1:2)=(/n4b(1),n4(1)/)
write(1,*) "Line(",ii+2,")={",n4b(2),",",n4(2),"};";linelist(ii+2,1:2)=(/n4b(2),n4(2)/) ! top south
write(1,*) "Line(",ii+3,")={",n4b(3),",",n4(3),"};";linelist(ii+3,1:2)=(/n4b(3),n4(3)/) ! top west
write(1,*) "Line(",ii+4,")={",n4b(4),",",n4(4),"};";linelist(ii+4,1:2)=(/n4b(4),n4(4)/) ! top north
write(1,*) "Line(",ii+5,")={",n4(1),",",n4(2),"};";linelist(ii+5,1:2)=(/n4(1),n4(2)/)  ! bottom east
write(1,*) "Line(",ii+6,")={",n4(2),",",n4(3),"};";linelist(ii+6,1:2)=(/n4(2),n4(3)/)  ! bottom south
write(1,*) "Line(",ii+7,")={",n4(3),",",n4(4),"};";linelist(ii+7,1:2)=(/n4(3),n4(4)/)  ! bottom west
write(1,*) "Line(",ii+8,")={",n4(4),",",n4(1),"};";linelist(ii+8,1:2)=(/n4(4),n4(1)/)  ! bottom north
n43(1,1:3)=(/ii+2,-(ii+5),-(ii+1)/)
n43(2,1:3)=(/ii+3,-(ii+6),-(ii+2)/)
n43(3,1:3)=(/ii+4,-(ii+7),-(ii+3)/)
n43(4,1:3)=(/ii+1,-(ii+8),-(ii+4)/)
ii=ii+8
iilast=ii
do i=1,4
  allocate(linb(nlinbry(i)))
  do j=1,nlinbry(i)
    call findline(n23dlinbryb(i,j,1),n23dlinbryb(i,j,2),linb(j),ii,linelist,ntrik,1) ! now all the lines should be already defined
  end do
  jj=jj+1
  write(1,*) "Line Loop(",jj,")={",(n43(i,k),",",k=1,3) ! comma is the last
  write(1,*) (linb(k),",",k=1,nlinbry(i)-1),linb(nlinbry(i)),"};" !  FOUR BOTTOM SIDE LINE LOOP
  write(1,*) "Plane Surface(",jj,")={",jj,"};"
  deallocate(linb)
end do
if ( ii .ne. iilast ) goto 100
jj=jj+1
write(1,*) "Line Loop(",jj,")={",ii-3,",",ii-2,",",ii-1,",",ii,"};" !    LAST BOTTOM LINE LOOP
write(1,*) "Plane Surface(",jj,")={",jj,"};"
write(1,*) "Surface Loop(2)={",(botlooplist(i),",",i=1,ntrik) ! end with comma
write(1,*)  jj-4,",",jj-3,",",jj-2,",",jj-1,",",jj,"};"
!write(1,*) "Volume(2)={2};"
!write(1,*) "Physical Volume (1)={1,2};"

!# inserted on 2016.10.11
!goto 101
bellgeofile = b_param%bellgeofile
nbell = b_param%nbell ! 2018.10.02
write(1,*) "nbell=",b_param%nbell,";" ! 2018.10.02
do ibell = 1,nbell
call FINDZ_TRI(b_param%xyz_bell(1:2,ibell),z0,h_mesh)
write(1,20) "ztop["  ,ibell,"]=",b_param%ztop(ibell),";"
write(1,20) "zthic[" ,ibell,"]=",b_param%zthic(ibell),";"
write(1,20) "zbot["  ,ibell,"]=",b_param%zbot(ibell),";"
write(1,20) "xb["    ,ibell,"]=",b_param%xyz_bell(1,ibell),";"
write(1,20) "yb["    ,ibell,"]=",b_param%xyz_bell(2,ibell),";"
write(1,20) "z0["    ,ibell,"]=",z0,";"
write(1,20) "radius[",ibell,"]=",b_param%radius(ibell),";"
write(1,20) "width[" ,ibell,"]=",b_param%width(ibell),";"
write(1,20) "Lc2["   ,ibell,"]=",b_param%reso_bell(ibell),";"
write(1,*) " " ! Inserted on 2016.10.11

end do ! 2018.10.02
write(1,*) "Include '",bellgeofile(1:len_trim(bellgeofile)),"' ;" !2018.10.02
101 continue

close(1)
write(*,*) "### OUTPREGEO8 END!! ###"
return
100 write(*,*) "GEGEGE!! iilast=",iilast,"not equal to ii",ii
stop
20 format(a,i5,a,g15.7,a)
end subroutine outpregeo8

!########################################################################### FINDZ_TRI
! istart and iend are added on 10.13, 2016
! This subroutine finds corresponding elements that includes given coord
subroutine FINDZ_TRI(x2,z,h_mesh)
!use mesh_type
use outerinnerproduct
use mesh_type
implicit none
type(mesh),intent(in)  :: h_mesh
real(8),   intent(in)  :: x2(2) ! x,y[km]
real(8),   intent(out) :: z
integer(4) :: i, j,n(4),is,ie,ii
real(8) :: elm_xyz(3,4),x3_center(3),xyz(3,h_mesh%node)
real(8),dimension(3) :: o1,o2,o3,o4
real(8) :: a1,a2,a3
real(8),dimension(3) :: v12,v1c,v13,v23,v2c,x3c
integer(4) :: ntri,m(3),n3(h_mesh%ntri,3)

!#[1]## set
ntri = h_mesh%ntri
n3   = h_mesh%n3
xyz  = h_mesh%xyz

do i=1,ntri ! modified on May 17, 2016
 m(1:3) = n3(i,1:3)
 v12(1:2) = xyz(1:2,m(2)) - xyz(1:2,m(1)) ; v12(3) = 0.d0
 v13(1:2) = xyz(1:2,m(3)) - xyz(1:2,m(1)) ; v13(3) = 0.d0
 v23(1:2) = xyz(1:2,m(3)) - xyz(1:2,m(2)) ; v23(3) = 0.d0
 v1c(1:2) = x2(1:2) - xyz(1:2,m(1))       ; v1c(3) = 0.d0
 v2c(1:2) = x2(1:2) - xyz(1:2,m(2))       ; v2c(3) = 0.d0
 o1 = outer(v13,v1c)
 o2 = outer(v1c,v12)
 o3 = outer(v2c,v23)
 if ( o1(3) .ge. 0. .and. o2(3) .ge. 0. .and. o3(3) .ge. 0.) then
  o4  = outer(v13,v12)
  a1  = o3(3)/o4(3)
  a2  = o1(3)/o4(3)
  a3  = o2(3)/o4(3)
  x3c(:) = a3 * xyz(:,m(3)) + a2 * xyz(:,m(2)) + a1 * xyz(:,m(1))
  write(*,*) "x2=",x2(1:2)
  write(*,*) "x3c(1:3)=",x3c(1:3)
  write(*,'(i3,a5,3g15.7)') (j," xyz=",xyz(1:3,m(j)),j=1,3)
  goto 101
 end if
end do
goto 99

101 continue
z = x3c(3)

return

99 continue
write(*,*) "GEGEGE! x2(1:2)=",x2(1:2),"not found in triangles."
stop
end


!######################################
subroutine findline(n1,n2,line,ii,linelist,ntrik,ifile)
implicit none
integer(4),intent(in) :: n1,n2! n1 and n2 are the start and end node # of line of interest
integer(4),intent(in) :: ntrik ! nodek is the number of nodes on 2-D ocean/land  surface
integer(4),intent(out) :: line ! line is line # which conect node n1 and n2. if n2 -> n1 line exists, line should be negative value
integer(4),intent(inout) :: ii ! ii is the number of lines which has been already defined
integer(4),dimension(ntrik*6,3),intent(inout) :: linelist ! linelist(i,1) is the start node # and linelist(i,2) is the end node # regarding i-th line
integer(4),intent(in) :: ifile ! ifile is the device numver of .geo file where the Line data should be written
integer(4) :: i
! find n1 -> n2 line
line=0
do i=9,ii ! 1 - 8 line is upper unique lines
if ( n1 .eq. linelist(i,1) .and. n2 .eq. linelist(i,2)) goto 9
if ( n2 .eq. linelist(i,1) .and. n1 .eq. linelist(i,2)) goto 10
end do
! make new line
ii=ii+1
write(ifile,*) "Line(",ii,")={",n1,",",n2,"};"
!write(*,*) "ntrik*3=",ntrik*3,"ii=",ii
line=ii
linelist(ii,1:2)=(/n1,n2/)
return
9 line=i
return
10 line=-i
return
end subroutine findline
