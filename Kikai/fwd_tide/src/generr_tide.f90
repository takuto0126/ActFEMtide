! coded on 2026.07.27
program generr_tide
implicit none
real(8)    :: ratio,bxa,bxp,bya,byp,bza,bzp,exa,exp,eya,eyp
real(8)    :: bxa_err,bya_err,bza_err,exa_err,eya_err
real(8)    :: bxp_err,byp_err,bzp_err,exp_err,eyp_err
real(8)    :: freq,pi,r2d
integer(4) :: i,j,k,nobs
complex(8) :: zssq,zxx,zxy,zyx,zyy,iunit=(0.,1.d0)
integer(4) :: nfreq
character(50) :: inputfolder,outputfolder
character(100),allocatable,dimension(:) :: file1,file2
pi = 4.d0*atan(1.d0)
r2d = 180.d0/pi

!#[1]## read parameter
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

!# output data files
!# ampfile
!# freq 1 0 1 0 0 amp_bx, amp_by, amp_bz, amp_ex, amp_ey, err_|bx|, err_|by|, err_|bz|, err_|ex|, err_|ey|
!# phase file
!# freq 1 0 1 0 0 phase_bx, phase_by, phase_bz, phase_ex, phase_ey, err_arg(bx), err_arg(by), err_arg(bz), err_arg(ex), err_arg(ey)
!#
do i=1,nobs
 open(unit=10,file=trim(inputfolder)//trim(file1(i)),status='old',form='formatted')
 open(unit=20,file=trim(outputfolder)//trim(file2(i))//"_amp.dat",status='replace',form='formatted')
 open(unit=21,file=trim(outputfolder)//trim(file2(i))//"_pha.dat",status='replace',form='formatted')
 do j=1,nfreq
  read(10,*) freq,bxa,bxp,bya,byp,bza,bzp,exa,exp,eya,eyp
  bxa_err = ratio*bxa
  bya_err = ratio*bya
  bza_err = ratio*bza
  exp_err = ratio*exp
  eya_err = ratio*eya
  exp_err = ratio*exp
  bxp_err = atan2(bxa*ratio, bxa)*r2d 
  byp_err = atan2(bya*ratio, bya)*r2d 
  bzp_err = atan2(bza*ratio, bza)*r2d
  exp_err = atan2(exp*ratio, exp)*r2d 
  eyp_err = atan2(eya*ratio, eya)*r2d 
  write(20,'(g15.7,a,10g15.7)') freq," 1 1 1 1 1 ",bxa,bya,bza,eya,eya,bxa_err,bya_err,bza_err,eya_err,eya_err 
  write(21,'(g15.7,1x,10g15.7)') freq," 1 1 1 1 1",bxp,byp,bzp,exp,eyp,bxp_err,byp_err,bzp_err,exp_err,eyp_err
 end do
 close(10)
 close(20)
 close(21)
end do


end program generr_tide
