program blrv

 implicit none

 integer, parameter :: nv = 195823
 integer, parameter :: nlc = 19252
 integer :: i, sumundef, outlc
 real(8), parameter :: id = 12.d0
 real(8), parameter :: xmin = 109.97, xmax = 164.94, ymin = 20.01, ymax = 62.65
 real(8) :: dx, dy
 real(8), dimension(nlc) :: x, y, blvr, blvi

 character*30 :: inlfile, invfile, outfile
 
! v:vn
!-----< input file, output file >-----

! Input Output file

 inlfile='./lat_lon_mesh'
 invfile='./v_xyfield_m2.dat'
 outfile='./blu_xyfield_m2.dat'

!-----<  main program >-------------------------------------

! Set location

 open(21,file=inlfile,status='old')

 dx = 1.d0/id
 dy = 1.d0/id
 sumundef = 0
 outlc = 0
 
 do  i = 1,nlc
     read(21,'(2f8.3)') y(i), x(i)

     if (xmin <= x(i) .and. x(i) <= xmax .and. &
         ymin <= y(i) .and. y(i) <= ymax        ) then

!----------------------------------------------------------

!-----< Bilinear interpolation >-----

        open(22,file=invfile,status='old')

        call clblr(nlc, nv, i       , x   , y   ,&
                   dx , dy, sumundef, blvr, blvi  )

        close(22)
!----------------------------------------------------------

     else
        blvr(i) = 0.d0
        blvi(i) = 0.d0
        outlc = outlc +1
     endif

 enddo

 close(21)


 open(30,file=outfile)

 do i = 1,nlc
    write(30,'(2f8.3,2f11.4)') y(i), x(i), blvr(i), blvi(i)
 enddo

 close(30)

 write(6,*) sumundef 
 write(6,*) outlc

end program blrv

!clblr
!----------------------------------------------------------
subroutine clblr (nlc, nv, il      , x   , y   ,&
                  dx , dy, sumundef, blvr, blvi  )

!----------------------------------------------------------
 implicit none

 integer :: m, n, undef, sundef, ix, iy
 real(8) :: x1, y1, x2, y2,  s1, s2, t1, t2, subvnr, subvni
 real(8), dimension(2,2) :: wf, vnr, vni
 integer, intent(in) :: nlc, nv, il
 integer, intent(inout) :: sumundef
 real(8), intent(in) :: dx, dy
 real(8), intent(in), dimension(nlc) :: x, y
 real(8), intent(out), dimension(nlc) :: blvr, blvi


 blvr(il) = 0.d0
 blvi(il) = 0.d0

 m = int(x(il)/dx)
 n = int(y(il)/dy)

 x1 = dfloat(m)*dx
 y1 = dfloat(n)*dy - 1.d0/24.d0

 s1 = (x(il) - x1)/dx
 s2 = 1.d0 - s1
 t1 = (y(il) - y1)/dy
 t2 = 1.d0 -t1

 wf(1,1) = s2*t2
 wf(1,2) = s2*t1
 wf(2,1) = s1*t2
 wf(2,2) = s1*t1

 x2 = x1 + dx
 y2 = y1 + dy

 sundef=0
 
 call findlc(nv, x1, y1, subvnr, subvni, undef)
 
 if (undef==0) then
    vnr(1,1) = subvnr
    vni(1,1) = subvni

 elseif (undef==1) then
    sundef = sundef + undef

 endif

 call findlc(nv, x2, y1, subvnr, subvni, undef)
 
 if (undef==0) then
    vnr(2,1) = subvnr
    vni(2,1) = subvni

 elseif (undef==1) then
    sundef = sundef + undef

 endif

 call findlc(nv, x1, y2, subvnr, subvni, undef)
 
 if (undef==0) then
    vnr(1,2) = subvnr
    vni(1,2) = subvni

 elseif (undef==1) then
    sundef = sundef + undef

 endif
 
 call findlc(nv, x2, y2, subvnr, subvni, undef)
 
 if (undef==0) then
    vnr(2,2) = subvnr
    vni(2,2) = subvni

 elseif (undef==1) then
    sundef = sundef + undef

 endif


 if (sundef==0) then

    do iy = 1,2

       do ix = 1,2
          blvr(il) = blvr(il) + vnr(ix,iy)*wf(ix,iy)
          blvi(il) = blvi(il) + vni(ix,iy)*wf(ix,iy)
       enddo

    enddo

 else
    sumundef = sumundef + 1
 endif

end subroutine clblr

!findlc
!----------------------------------------------------------
subroutine findlc (nv, xl, yl, vnr, vni, undef)

!----------------------------------------------------------
 implicit none

 integer :: i
 real(8) :: xl3, yl3, lon3, lat3
 integer, intent(in) :: nv
 integer, intent(out) :: undef
 real(8), intent(in) :: xl, yl
 real(8), intent(out) :: vnr, vni
 real(8), dimension(nv) :: lon, lat, vr, vi


 xl3=xl*1000.d0
 yl3=yl*1000.d0
 vnr = 999.d0
 vni = 999.d0
 undef = 0

 do i = 1,nv
    read(22,'(2f8.3,2f11.4)') lon(i), lat(i), vr(i), vi(i)

    lon3 = lon(i)*1000.d0
    lat3 = lat(i)*1000.d0

    if (int(xl3)==int(lon3) .and. int(yl3)==int(lat3)) then
       vnr = vr(i)
       vni = vi(i)
    endif

 enddo

 if (int(vnr)==999 .or. int(vni)==999) then
    undef=1
 endif

end subroutine findlc
