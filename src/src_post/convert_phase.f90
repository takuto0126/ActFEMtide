!# TOPEX POSEIDON phase (Jan 1, 1992) to Grinidge Phase
!# OTPS Phase : phase at Jan 1, 1992 ()
program convert_phase
implicit none
character(70),allocatable,dimension(:) :: inputfile
integer(4) :: nfile,nfreq
real(8),allocatable,dimension(:,:,:) :: bxyz,exy

read(*,*) nfile
write(*,*) "nfile=",nfile
allocate(inputfile(nfile))
do i=1,nfile
 read(*,*) inputfile(i)
 write(*,*) "inputfile =",inputfile(i)
end do

do i=1,nfile
open(1,file=inputfile(i))

 10 read(1,*,end=11)
 nfreq=nfreq + 1
 goto 10
 11 rewind(1)
 
 write(*,*) "nfreq =",nfreq
 allocate(freq(nfreq),bxyz(2,3,nfreq),exy(2,2,nfreq))
 open(2,file=inputfile(i)//"2" )
 do i=1,nfreq
  read(1,*) freq(i),bxyz(:,:,i),exy(:,:,i)
  call change_phase_GN(bxyz(:,:,i),exy(:,:,i))
  write(2,*) freq(i),bxyz(:,:,i),exy(:,:,i)
 end do
 close(2)
 close(1)
end do

!##############################################################
subroutine change_phase_GN(bxyz,exy)
implicit none
real(8),intent(in) :: bxyz(2,3)
real(8),intent(in) :: exy(2,2)
real(8) :: pi,r2d,omega,period,dtime,second_per_day
real(8) :: mjd2,mjd2
real(8) :: dat(2,5),phase
complex(8) :: iunit=(0.,1,d0),c

dat(:,1:3) = bxyz
dat(:,4:5) = exy
pi = 4.*ata(1.d0)
r2d = 180./pi
d2r = pi/180.
period = 12.42060122 * 60*60 ! sec
omega = 2.*pi/period ! rad/s
second_per_day = 86400.0 
call date_mjd(1,1,2000,mjd2) ! 01.01.2000
call date_mjd(1,1,1992,mjd1) ! 01.01.1992
dtime = (mjd2 - mjd1 )*second_per_day

do i=1,5
 c = exp(iunit*(omega*dtime + dat(2,i)*d2r))
 phase = atan2(imag(c),real(c))*r2d
 dat(2,i) = phase
end do

!# update phase
bxyz(2,1:3) = dat(2,1:3)
exy(2,1:2)  = dat(2,4:5)

return
end

!##############################################################
      subroutine date_mjd(mm,id,iyyy,mjd)
      implicit none
! converts date to mjd
! INPUT:  id - day, mm - month, iyyy - year
! OUTPUT: mjd>0 - modified julian days
! date>=11.17.1858 corresponds to mjd=0
      integer dpm(12),days,i,nleap,k
      integer mm,id,iyyy,mjd
      data dpm/31,28,31,30,31,30,31,31,30,31,30,31/
      mjd=0
! NO earlier dates
      if(iyyy.lt.1858)iyyy=1858
      if(iyyy.eq.1858.and.mm.gt.11)mm=11
      if(iyyy.eq.1858.and.mm.eq.11.and.id.gt.17)id=17
!
      days=0 
      do i=1,mm-1
       days=days+dpm(i)
       if(i.eq.2.and.int(iyyy/4)*4.eq.iyyy)days=days+1         
      enddo
      days=days+id-321
! leap day correction
      do k=1900,iyyy,100
       if(iyyy.eq.k.and.mm.gt.2)days=days-1
      enddo
      do k=2000,iyyy,400
       if(iyyy.eq.k.and.mm.gt.2)days=days+1
      enddo
! EACH 4th year is leap year
      nleap=int((iyyy-1-1860)*0.25)            
      if(iyyy.gt.1860)nleap=nleap+1
! EXCEPT
      do k=1900,iyyy-1,100
       if(iyyy.gt.k)nleap=nleap-1
       if(iyyy.eq.k.and.mm.gt.2)days=days-1
      enddo
! BUT each in the row 2000:400:... IS LEAP year again
      do k=2000,iyyy-1,400
       if(iyyy.gt.k)nleap=nleap+1
       if(iyyy.eq.k.and.mm.gt.2)days=days+1
      enddo
      mjd=365*(iyyy-1858)+nleap &
         +days   
      return
      end