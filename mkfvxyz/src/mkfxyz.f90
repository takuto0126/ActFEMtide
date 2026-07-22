implicit none
real(8)::date,alt,colat,elong,x,y,z,f,lat,lon
integer(4)::isv,itype

isv=0
itype=1

open(10,file="lat_lon_mesh")
open(11,file="fxyz_mesh")

date=2010
!alt=-2.362
!colat=109.8383
!elong=184.2947

do


read(10,*,end=999)lat,lon,alt

colat=90-lat
elong=lon
if(lon<0)then
   elong=360-lon
end if

call igrf13syn(isv,date,itype,alt,colat,elong,x,y,z,f)

write(11,*)lat,lon,alt,y,x,-z,f

end do



!write(*,*)"x",x,"nT"
!write(*,*)"y",y,"nT"
!write(*,*)"z",z,"nT"
!write(*,*)"f",f,"nT"

999 end program
