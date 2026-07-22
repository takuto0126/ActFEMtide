  use mesh_type
  use param_mesh
  use constants
  
  implicit none

  real(8)       :: wb,eb,nb,sb,xmax,ymax,x,y,z,lat0,lon0
  real(8)       :: lonorigin,latorigin,colatorigin
  integer       :: i,noden
  character(70) :: mesh_ctl
  character(70) :: meshfile,oceanmeshfile
  real(8),allocatable,dimension(:) ::lon,lat,colat,alt
  type(meshpara):: g_meshpara
  type(mesh)    :: em3d_mesh
  type(mesh)    :: ocean_mesh
  type(mesh)    :: in_mesh

        
!mesh_ctl="mesh/mesh.ctl"
read(*,*) mesh_ctl ! 2024.03.07 Takuto Minam
write(*,*) "ctl file=",mesh_ctl ! 2024.03.07 Takuto Minami
write(6,*)"good"


call readmeshpara(g_meshpara,mesh_ctl)

 ! wb=138.5
 ! eb=149.5
 ! sb=34.5
 ! nb=43.5

wb= g_meshpara%wlon
eb= g_meshpara%elon
sb= g_meshpara%slat
nb= g_meshpara%nlat



!write(6,*)"good"

  latorigin=(nb+sb)/2
  lonorigin=(wb+eb)/2
  colatorigin=90.d0-latorigin
write(*,*) "latorigin, lonorigin",latorigin, lonorigin

!  meshfile="mesh/em3d.msh"   ! 2024.03.07 Takuto Minami
!  oceanmeshfile="mesh/ocean.msh"  ! 2024.03.07 Takuto Minami
read(*,*) oceanmeshfile ! 2024.03.07 Takuto Minami
read(*,*) meshfile       ! 2024.03.07 Takuto Minami
write(*,*) "oceanmesh",trim(oceanmeshfile) ! 2024.03.07 TM
write(*,*) "meshfile",trim(meshfile) ! 2024.03.07 TM
 

  CALL READMESH_TOTAL(em3d_mesh,meshfile)
  CALL READMESH_TOTAL(ocean_mesh,oceanmeshfile)

  call  mksphericalmesh(meshfile,mesh_ctl)


  allocate(lon(ocean_mesh%node),lat(ocean_mesh%node),colat(ocean_mesh%node))
 
!open(10,file="em3d.msh")
 ! open(11,file="polygon_ki.msh")!
 ! open(12,file="mesh_latlon.dat")

 ! read(10,*)noden,xmax,ymax   !x:lat   y:lon

 !   do i=1,5

    ! read(11,"()")

 ! end do

xmax=em3d_mesh%xyz(1,1)
ymax=em3d_mesh%xyz(2,1)

open(12,file="lat_lon_mesh")

  do i=1,ocean_mesh%node

 !   lat=x/xmax*(nb-lat0)
  !  lon=y/ymax*(eb-lon0)
!lat=ocean_mesh%xyz(2,i)/(ymax-g_meshpara%lenout)*(nb-lat0)+lat0
     !lon=ocean_mesh%xyz(1,i)/(xmax-g_meshpara%lenout)*(eb-lon0)+lon0
    

     
     lat(i) = ocean_mesh%xyz(2,i)/d2r/earthrad + latorigin
     colat(i)=90.0-lat(i)
     lon(i) = ocean_mesh%xyz(1,i)/d2r/earthrad/sin(colatorigin*d2r) + lonorigin

!write(6,*)lat-latorigin
write(12,*)lat(i),lon(i),ocean_mesh%xyz(3,i),ocean_mesh%xyz(1,i),ocean_mesh%xyz(2,i)
  end do

  
end program


subroutine mksphericalmesh(meshfile,mesh_ctl)

  use param_mesh
  use mesh_type
  
implicit none
    real(8),allocatable,dimension(:)::lon,lat,alt,colat
    real(8),allocatable,dimension(:)::x,y,z,r
  character(70) :: meshfile
  character(70) :: mesh_ctl
  real(8)::nnode,latorigin,lonorigin,colatorigin
  integer::j,info(9),nele,i
  type(mesh)::in_mesh
  type(meshpara)::g_meshpara

  write(6,*)"sphericalstart"

  
  latorigin=-20
  lonorigin=185
  colatorigin=90-latorigin

  call readmeshpara(g_meshpara,mesh_ctl)
  call readmesh_total(in_mesh,meshfile)
 
  allocate(lon(in_mesh%node),lat(in_mesh%node),colat(in_mesh%node),alt(in_mesh%node))
allocate(x(in_mesh%node),y(in_mesh%node),z(in_mesh%node),r(in_mesh%node))
   open(30,file=meshfile)
   open(31,file="em3dspherical.msh")


    
  lat(1:in_mesh%node) = in_mesh%xyz(2,1:in_mesh%node)/d2r/earthrad + latorigin
  colat(1:in_mesh%node)=90.0-lat(1:in_mesh%node) 
  lon(1:in_mesh%node) = in_mesh%xyz(1,1:in_mesh%node)/d2r/earthrad/sin(colatorigin*d2r) + lonorigin
  alt(1:in_mesh%node)=in_mesh%xyz(3,1:in_mesh%node)
  
  r(1:in_mesh%node)  = earthrad + alt(1:in_mesh%node) ! [km]
  x(1:in_mesh%node)  = r*sin(colat(1:in_mesh%node) *d2r)*cos(lon(1:in_mesh%node) *d2r) ! [km]
  y(1:in_mesh%node)  = r*sin(colat(1:in_mesh%node) *d2r)*sin(lon(1:in_mesh%node) *d2r) ! [km]
  z(1:in_mesh%node)  = r*cos(colat(1:in_mesh%node) *d2r)              ! [km]
  do i=1,in_mesh%node+7
     read(30,"()")
  end do
 
  write(31,*)"$MeshFormat"
  write(31,*)"2.2 0 8"
  write(31,*)"$EndMeshFormat"
  write(31,*)"$Nodes"
  write(31,*)in_mesh%node
  
 
  
  do i=1,in_mesh%node
     !write(6,*)in_mesh%xyz(1:3,i)
     !write(6,*)i,lon(i),lat(i),alt(i)
     !write(6,*)i,x(i),y(i),z(i)
     write(31,*)i,x(i),y(i),z(i)
  end do
  

  
  write(31,*)"$EndNodes"
  write(31,*)"$Elements"
  

read(30,*)nele
 
write(31,*)nele
 

 
  do i=1,in_mesh%npoi
  !write(6,*)i,info(1:6)
  read(30,*)info(1:6)
  write(31,*)info(1:6)

end do

do i=in_mesh%npoi+1,nele
  ! write(6,*)i,info(1:9)
  read(30,*)info(1:9)
  write(31,*)info(1:9)

end do
   write(31,*)"$EndElements"



end subroutine mksphericalmesh




  

  
 


