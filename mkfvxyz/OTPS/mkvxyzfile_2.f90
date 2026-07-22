  implicit none

  real(8)::lat,lon,alt,vamp,vphase,uamp,uphase,lat2,lon2,lat3,lon3,pi,x,y
  complex(8)::xyz(2),inum
  integer::i
pi=3.1415926535
  inum=(0,1)

  open(10,file="lat_lon_mesh")
  open(11,file="u_xyz_mesh")
  open(12,file="v_xyz_mesh")
  open(13,file="vxyz_mesh")!result
  open(14,file="checkvxyz")

  do i=1,3
     read(11,*)
     read(12,*)
  end do

  read(10,*)lat,lon,alt,x,y
  read(11,*)lat3,lon3,uamp,uphase
  read(12,*)lat3,lon3,vamp,vphase

  xyz(1)=uamp*(cos(uphase*pi/180)+inum*sin(uphase*pi/180))
  xyz(2)=vamp*(cos(vphase*pi/180)+inum*sin(vphase*pi/180))
  

  

  write(13,*)lat,lon,alt,xyz(1)/100,xyz(2)/100,1
 write(14,*)x,y,real(xyz(1))

  do i=2,1000000

     read(10,*,end=999)lat2,lon2,alt,x,y
     
     if(lat.eq.lat2.and.lon.eq.lon2)then
  
        write(13,*)lat,lon,alt,xyz(1)/100,xyz(2)/100,i


     elseif(lat.ne.lat2.or.lon.ne.lon2)then
    

        read(11,*,err=99)lat3,lon3,uamp,uphase
        read(12,*,err=99)lat3,lon3,vamp,vphase

         xyz(1)=uamp*(cos(uphase*pi/180)+inum*sin(uphase*pi/180))
         xyz(2)=vamp*(cos(vphase*pi/180)+inum*sin(vphase*pi/180))

         goto 100

99       xyz(1)=(0,0)
         xyz(2)=(0,0)

         

100      write(13,*)lat2,lon2,alt,xyz(1)/100,xyz(2)/100,i
         write(14,*)x,y,real(xyz(1))
         

         lat=lat2
         lon=lon2

      end if

   end do

999 close(13)
   close(14)

 end program



  

  
