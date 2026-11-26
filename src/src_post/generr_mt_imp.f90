program generr_mt_imp
implicit none
real(8)    :: assq,freq,a,b,c,d,e,f,g,h,ratio
integer(4) :: i,j,k,nobs
complex(8) :: zssq,zxx,zxy,zyx,zyy,iunit=(0.,1.d0)
integer(4) :: nfreq
character(50) :: inputfolder,outputfolder
character(100),allocatable,dimension(:) :: file1,file2

!#[1]##
read(*,*) ! read header
read(*,'(20x,i10)') nobs
write(*,*) "nobs : ",nobs
read(*,'(20x,i10)') nfreq
write(*,*) "nfreq : ",nfreq
read(*,'(20x,a)') inputfolder
read(*,'(20x,a)') outputfolder
allocate(file1(nobs),file2(nobs))
do i=1,nobs
 read(*,'(20x,a)') file1(i)
 read(*,'(20x,a)') file2(i)
 write(*,*) i,"in  file ",trim(file1(i))
 write(*,*) i,"out file ",trim(file2(i))
end do
read(*,'(20x,g15.7)') ratio
write(*,*) "ratio",ratio

!#[2]## calculate impedance error based on ssq MT impedance

do k=1,nobs

open(1,file=trim(inputfolder)//file1(k))
open(2,file=trim(outputfolder)//file2(k))
open(3,file=trim(outputfolder)//"DUM_"//file2(k))
 do i=1,nfreq
  read(1,*) freq,a,b,c,d,e,f,g,h
  zxx = a + iunit*b
  zxy = c + iunit*d
  zyx = e + iunit*f
  zyy = g + iunit*h
!  write(*,*) freq,zxx,zxy,zyx,zyy
  zssq = sqrt((zxx**2. + zxy**2. + zyx**2. + zyy**2.)/2.) ! 2025.07.31
  assq = sqrt(real(zssq)**2. + imag(zssq)**2.)
  write(2,'(f15.7,a,8g15.7)') freq," 1 1 1 1 ",(assq*ratio,j=1,8)
  write(3,'(f15.7,a,8g15.7)') freq," 0 0 0 0 ",(assq*ratio,j=1,8)
 end do
close(1)
close(2)

end do ! obs loop

end program generr_mt_imp
