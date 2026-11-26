program generr_mt_tip
implicit none
real(8)    :: txabs,tyabs,freq,a,b,c,d,ratio
integer(4) :: i,j,k,nobs
complex(8) :: tx,ty,iunit=(0.,1.d0)
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
 do i=1,nfreq
  read(1,*) freq,a,b,c,d
  tx = a + iunit*b
  ty = c + iunit*d
  txabs = sqrt(real(tx)**2. + imag(tx)**2.)
  tyabs = sqrt(real(ty)**2. + imag(ty)**2.)
  write(2,'(f15.7,a,8g15.7)') freq," 1 1",&
  & a,b,c,d,txabs*ratio,txabs*ratio,tyabs*ratio,tyabs*ratio
 end do
close(1)
close(2)

end do ! obs loop

end program generr_mt_tip
