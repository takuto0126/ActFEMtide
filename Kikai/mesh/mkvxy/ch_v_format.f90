program chformat

 implicit none

 integer :: i
 integer, parameter :: nin = 195823
 real(8), parameter :: pi = 3.141592653589793d0
 real(8), dimension(nin) :: lon, lat, vea, vep, vna, vnp
 real(8), dimension(nin) :: ver, vei, vnr, vni, lonu, latv

 character*30 :: infile, outufile, outvfile

! u:ve, v:vn

!-----<  input file, output file >-----

 infile='./vfield.m2'
 outufile='./u_xyfield_m2.dat'
 outvfile='./v_xyfield_m2.dat'

!-----<  input data >-----

 open(20,file=infile,status='old')

 do i = 1,nin
    read(20,'(2f8.3,4f11.4)') lon(i), lat(i), vea(i), vep(i), vna(i), vnp(i)
 enddo

 close(20)

!-----<  main program >-----

 do i = 1,nin
    ver(i) = vea(i) * cos(vep(i)*pi/180.d0)
    vei(i) = vea(i) * sin(vep(i)*pi/180.d0)
    vnr(i) = vna(i) * cos(vnp(i)*pi/180.d0)
    vni(i) = vna(i) * sin(vnp(i)*pi/180.d0)
 enddo

 do i = 1,nin
    lonu(i) = lon(i) - 1.d0/24.d0
    latv(i) = lat(i) - 1.d0/24.d0
 enddo

!-----<  output data >-----

 open(31,file=outufile)
 
 do i =1,nin
    write(31,'(2f8.3,2f11.4)') lonu(i), lat(i), ver(i), vei(i)
 enddo

 close(31)

 open(32,file=outvfile)
 
 do i =1,nin
    write(32,'(2f8.3,2f11.4)') lon(i), latv(i), vnr(i), vni(i)
 enddo

 close(32)

end program chformat
