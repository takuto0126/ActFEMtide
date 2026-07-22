program phase_shift
implicit none
integer(4) :: i,ndat
complex(8),allocatable,dimension(:) :: u,v
real(8),   allocatable,dimension(:,:) :: xyz
real(8) :: phase,time,dphi,omega,pi,period
real(8) :: time_mjd,dum=0.0,r2d
integer(4) :: ntime=1
integer :: yyyy,mm,dd,hh,mi,ss ! Date and time
complex(8) :: u_s,v_s,P,iunit=(0.d0,1.d0)
integer :: SecondsPerDay=86400
pi = 4.*atan(1.d0)
r2d=180./pi
!# read vxyz_mesh
open(1,file="vxyz_mesh")
  i=0
  do
     read(1,*,end=99)
     i=i+1
  end do
  99 ndat = i
  allocate(u(ndat),v(ndat),xyz(3,ndat))
  rewind(1)
  do i=1,ndat
      read(1,*) xyz(1:3,i),u(i),v(i)
  end do
close(1)

!# estimate phase
!# grinidge time
yyyy=2000
mm=1
dd=1
hh=0
mi=0
ss=0
open(1,file="tmp.dat")
write(1,*) dum,dum,yyyy,mm,dd,hh,mi,ss
close(1)
open(1,file="tmp.dat")
call read_time(ntime,1,1,yyyy,mm,dd,hh,mi,ss,time_mjd)
close(1)
! relatively Jan 1 1992 (48622mjd)       
         time=(time_mjd-dble(48622))*dble(SecondsPerDay) ! sec
         write(*,*) "1992 to 2000 time",time,"sec"
         write(*,*) "Days = time/86400",time/86400
         write(*,*) "years = time/86400/365.2",time/86400/365.2
!# from README of TOPEX
!#  height is h(t,x) = Re [ h(x) exp { i [w (t - t0) + V0(t0)] } ]
!# the velocity is similar
period = 12.42060122 *3600.d0 ! sec
omega = 2.*pi/period
dphi = time * omega ! rad
write(*,*) "dphi =",atan2(sin(dphi),cos(dphi))*r2d, "deg"

!# output phase shifted vyx_mesh 
open(1,file="vxyz_greenwich_mesh")
P = exp(iunit*dphi)
do i=1,ndat
 u_s = u(i)*P
 v_s = v(i)*P
 write(1,*) xyz(1:3,i),u_s,v_s,i
end do
close(1)

end program

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_time(ntime,iunit,iopt,yyyy,mm,dd,hh,&
                           mi,ss,time_mjd)
      implicit none
      integer k,iunit,iopt,julian,mjd,mm1,dd1,yyyy1
      integer ntime,yyyy(ntime),mm(ntime),dd(ntime)
      integer hh(ntime),mi(ntime),ss(ntime)
      real dum
      real*8 time_mjd(ntime)
      character*10 cdate,deblank
!
      do k=1,ntime
       if(iopt.eq.1)then
        read(iunit,*)dum,dum,yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       else
        read(iunit,*)yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       endif
! convert to mjd
       call date_mjd(mm(k),dd(k),yyyy(k),mjd)
! check if exists such a date
        julian=mjd+2400001        
        call CALDAT (JULIAN,MM1,DD1,YYYY1)
        if(mm(k).ne.mm1.or.dd(k).ne.dd1.or.yyyy(k).ne.yyyy1)then
         write(cdate,'(i2,a1,i2,a1,i4)')mm(k),'.',dd(k),'.',yyyy(k)
         cdate=deblank(cdate)
         write(*,*)'Wrong date in (lat_lon)_time file:',cdate
         stop
        endif
        time_mjd(k)=dble(mjd)+dble(hh(k))/24.D0+ &
                    dble(mi(k))/(24.D0*60.D0)  + &
                    dble(ss(k))/(24.D0*60.D0*60.D0)
      enddo
      return
      end
