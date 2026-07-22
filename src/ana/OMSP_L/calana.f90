! Coded on 2017.02.17 by Takuto MINAMI
!
program calana
implicit none
integer(4) :: nobs=4
real(8),dimension(2,4) :: obsxy
real(8),dimension(2,2) :: xs
real(8) :: ds=1.0 ![m]
real(8) :: I=1.d0
integer(4),parameter :: nfreq = 19
integer(4) :: l,j,k,n
real(8),dimension(nfreq) :: freq
real(8)    :: exy(2,2),bxyz(3,2)
character(3),dimension(4) :: obsname
freq = (/ 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, &
&         20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0 /)
obsname(1) ="A02"
obsname(2) ="A04"
obsname(3) ="A01"
obsname(4) ="A03"

!#[-1]## site locations
 obsxy(1:2,1)=(/ -0.1322694   ,   0.1997259    /) ! A02 [km]
 obsxy(1:2,2)=(/ -0.2658363   ,  -0.1345957E-01/) ! A04 [km]
 obsxy(1:2,3)=(/ -0.1276669   ,  -0.2418020    /) ! A01 [km]
 obsxy(1:2,4)=(/  0.1995709   ,  -0.3178668    /) ! A03 [km]
 xs(1:2,1)   =(/ -0.606938230999978   ,    0.644824799999595 /)
 xs(1:2,2)   =(/ -0.324348940999946   ,    0.698154110000003 /)
 obsxy = obsxy*1.d3 ! [km] -> [m]
 xs    = xs*1.d3     ! [km] -> [m]

 write(*,*) "freq=",freq,"[Hz]"
 write(*,*) "I=",I,"[A]"

!#[2]# open files
 do l=1,nobs
 open(1,file=obsname(l)(1:3)//".dat")

!#[3]## calculation
 do j=1,nfreq
  call dipole_l(xs(:,1),xs(:,2),I,ds,obsxy(:,l),exy,bxyz,freq(j))
  write(1,'(11g15.7)') freq(j),((bxyz(k,n),n=1,2),k=1,3),&
  &                            ((exy(k,n),n=1,2),k=1,2)
 end do

!#[4]## output
 close(1)
 end do

end program calana
