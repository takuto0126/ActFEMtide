! coded on 2017.09.30
! This program assume the land is homogeneous without any topography
program calanacomp
use param
use m_param_ana
implicit none
type(param_forward) :: g_param
type(param_source)  :: s_param
type(param_cond)    :: g_cond
integer(4)          :: nsrc, nobs,nfreq
integer(4)          :: j,k,l
real(8), allocatable,dimension(:) :: freq
real(8)             :: obsxy(2), exy(2,2),bxyz(3,2)
real(8)             :: I
real(8), allocatable,dimension(:,:)     :: xyzobs
real(8), allocatable,dimension(:,:,:,:) ::bexyzobs ! (10,obs,nsrc,nfreq)
real(8), allocatable,dimension(:,:)     :: xs1,xs2
real(8)             :: ds ! source interval
character(50)       :: header, head,outfile
character(50),allocatable,dimension(:) :: obsname
character(1)        :: num

!#[1] read param
 call readparam(g_param,s_param,g_cond)
 nsrc      = s_param%nsource
 I         = s_param%I
 nfreq     = g_param%nfreq
 nobs      = g_param%nobs
 allocate(freq(nfreq))
 freq      = g_param%freq
 cond(1,0) = g_cond%sigma_air
 cond(1,1) = g_cond%sigma_land(1) ! assume homogeneous earth
 allocate( xs1(3,nsrc), xs2(3,nsrc))
 xs1       = s_param%xs1 * 1.d3  ! [km] -> [m]
 xs2       = s_param%xs2 * 1.d3  ! [km] -> [m]
 allocate(xyzobs(3,nobs))
 xyzobs    = g_param%xyzobs*1.d3 ! [km] -> [m]
 header    = g_param%outputfolder
 allocate(obsname(nobs))
 obsname   = g_param%obsname
 allocate( bexyzobs(10,nobs,nsrc,nfreq) )! 1:amp, 2 phase
 write(*,*) "set end!!"

!#[2] calculate
ds = 10.0 ! [m]
do k=1,nfreq
 write(*,*) "k=",k,"freq=",freq(k),"[Hz]"
 do j=1,nsrc
  do l=1,nobs
    obsxy = xyzobs(1:2,l)
   call dipole_l(xs1(1:2,j),xs2(1:2,j),I,ds,obsxy,exy,bxyz,freq(k))
   bexyzobs(1:2, l,j,k) = bxyz(1,1:2) ! bx amp, phase
   bexyzobs(3:4, l,j,k) = bxyz(2,1:2) ! by amp, phase
   bexyzobs(5:6, l,j,k) = bxyz(3,1:2) ! bz amp, phase
   bexyzobs(7:8, l,j,k) = exy( 1,1:2) ! ex amp, phase
   bexyzobs(9:10,l,j,k) = exy( 2,1:2) ! ey amp, phase
  end do
 end do
end do

!#[3] output
do l=1,nobs
  head=header(1:len_trim(header))//obsname(l)(1:len_trim(obsname(l)))
 do j=1,nsrc
  write(num,'(i1)') j
  outfile=head(1:len_trim(head))//"_S"//num//"_ana.dat"
  open(1,file=outfile)
   do k=1,nfreq
    write(1,'(11g15.7)') freq(k),bexyzobs(1:10,l,j,k)
   end do
  close(1)
 end do
end do

end program
