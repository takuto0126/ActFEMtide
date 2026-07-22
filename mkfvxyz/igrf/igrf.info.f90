program igrf
implicit none
real(8)::date,alt,colat,elong,x,y,z,f
integer(4)::isv,itype

isv=0
itype=1

date=2010
alt=-2.362
colat=109.8383
elong=184.2947

call igrf13syn(isv,date,itype,alt,colat,elong,x,y,z,f)

write(*,*)"x",x,"nT"
write(*,*)"y",y,"nT"
write(*,*)"z",z,"nT"
write(*,*)"f",f,"nT"

end program

