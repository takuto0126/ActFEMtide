! moved to src/common/ to share with inversion
! Coded on 2016, Jan. 22
module outresp
use param
implicit none

type obsfiles
integer(4) :: nfile
character(50),allocatable,dimension(:) :: filename
integer(4),   allocatable,dimension(:) :: devicenumber
end type

type respdata
integer(4) :: nobs
complex(8),allocatable,dimension(:) :: ftobs  ! 2021.09.15
real(8),   allocatable,dimension(:) :: fpobsamp
real(8),   allocatable,dimension(:) :: fpobsphase
real(8),   allocatable,dimension(:) :: fsobsamp
real(8),   allocatable,dimension(:) :: fsobsphase
real(8),   allocatable,dimension(:) :: ftobsamp
real(8),   allocatable,dimension(:) :: ftobsphase
end type

! MT impedances are added 2021.09.15
type respmt ! 2021.09.14
integer(4) :: nobs ! 2021.09.15
complex(8),allocatable,dimension(:) :: zxx ! x north, y east  ! 2021.09.15
complex(8),allocatable,dimension(:) :: zxy   ! 2021.09.15
complex(8),allocatable,dimension(:) :: zyx   ! 2021.09.15
complex(8),allocatable,dimension(:) :: zyy   ! 2021.09.15
real(8),   allocatable,dimension(:) :: rhoxx ! 2021.09.15
real(8),   allocatable,dimension(:) :: rhoxy ! 2021.09.15
real(8),   allocatable,dimension(:) :: rhoyx ! 2021.09.15
real(8),   allocatable,dimension(:) :: rhoyy ! 2021.09.15
real(8),   allocatable,dimension(:) :: phaxx ! 2021.09.15
real(8),   allocatable,dimension(:) :: phaxy ! 2021.09.15
real(8),   allocatable,dimension(:) :: phayx ! 2021.09.15
real(8),   allocatable,dimension(:) :: phayy ! 2021.09.15
end type

contains
!######################################### ALLOCATERESPMT
! coded on 2021.09.14
subroutine ALLOCATERESPMT(nobs,resp_mt)
implicit none
integer(4),    intent(in)  :: nobs
type(respmt),  intent(out) :: resp_mt

resp_mt%nobs=nobs
allocate(resp_mt%zxx(  nobs))
allocate(resp_mt%zxy(  nobs))
allocate(resp_mt%zyx(  nobs))
allocate(resp_mt%zyy(  nobs))
allocate(resp_mt%rhoxx(nobs))
allocate(resp_mt%rhoxy(nobs))
allocate(resp_mt%rhoyx(nobs))
allocate(resp_mt%rhoyy(nobs))
allocate(resp_mt%phaxx(nobs))
allocate(resp_mt%phaxy(nobs))
allocate(resp_mt%phayx(nobs))
allocate(resp_mt%phayy(nobs))

return
end subroutine
!######################################### ALLOCATERESPDATA
subroutine ALLOCATERESPDATA(nobs,resp)
implicit none
integer(4),    intent(in) :: nobs
type(respdata),intent(out) :: resp

resp%nobs=nobs
allocate(resp%ftobs     (resp%nobs)) ! 2021.09.15
allocate(resp%fpobsamp  (resp%nobs))
allocate(resp%fpobsphase(resp%nobs))
allocate(resp%fsobsamp  (resp%nobs))
allocate(resp%fsobsphase(resp%nobs))
allocate(resp%ftobsamp  (resp%nobs))
allocate(resp%ftobsphase(resp%nobs))

return
end subroutine

!######################################################
!# 2017.10.12
!# copied from fluidity-4.1.10/femtools/Funits.F90
  function free_unit()
    !!< Find a free unit number. Start from unit 10 in order to ensure that
    !!< we skip any preconnected units which may not be correctly identified
    !!< on some compilers.
    integer :: free_unit
    
    logical :: connected

    do free_unit=10, 99

       inquire(unit=free_unit, opened=connected)

       if (.not.connected) return

    end do
    
    write(*,*) "No free unit numbers avalable"
    stop

  end function

!######################################## OUTFREQFILES
subroutine OUTFREQFILES(freq,resp,g_param,comp)
implicit none
real(8),intent(in) :: freq
type(respdata),intent(in) :: resp
type(param_forward),intent(in) :: g_param
character(2),intent(in) :: comp
character(50) :: filename
character(9) :: num
character(50) :: header
integer(4) :: i,j

header=g_param%outputfolder
write(num,'(e9.2)') freq
filename=header(1:len_trim(header))//"freqresp"//num(2:9)//comp(1:2)//".dat"

open(1,file=filename)
!write(*,*) "resp%nobs=",resp%nobs
do i=1,resp%nobs
 write(1,'(8e15.7)') (g_param%xyzobs(j,i),j=1,2),resp%ftobsamp(i),  resp%ftobsphase(i),&
   &  resp%fsobsamp(i),resp%fsobsphase(i),resp%fpobsamp(i),resp%fpobsphase(i)
end do
close(1)

return
end subroutine

!######################################### MAKEOBSFILES
subroutine MAKEOBSFILES(g_param,files)
implicit none
type(param_forward),intent(in) :: g_param
type(obsfiles),intent(out) :: files
integer(4) :: i
character(3) :: num
character(50) :: header,site

!#[1] set # of files
files%nfile=g_param%nobs
allocate( files%filename(files%nfile))
allocate( files%devicenumber(files%nfile))

!#[2] set filename
header=g_param%outputfolder
do i=1,files%nfile
write(num,'(i3.3)') i
 site=g_param%obsname(i)
 files%filename(i)=header(1:len_trim(header))//site(1:len_trim(site))//".dat"
end do

!#[3] set devicenumber
do i=1,files%nfile
 files%devicenumber(i)=i+30
end do

return
end subroutine

!######################################### OUTOBSFILES
subroutine OUTOBSFILES(freq,files,resp5)
implicit none
type(obsfiles),intent(in) :: files
type(respdata),dimension(5),intent(in) :: resp5
real(8),intent(in) :: freq
integer(4) :: i,j

do i=1,files%nfile
 write(files%devicenumber(i),'(11g15.7)') freq,(resp5(j)%ftobsamp(i),resp5(j)%ftobsphase(i),j=1,5)
end do

return
end subroutine

!######################################## OUTFREQFILES
!# 2017.10.12
subroutine OUTFREQFILES2(freq,nsr,resp,g_param,obs)
use obs_type
implicit none
real(8),            intent(in) :: freq
integer(4),         intent(in) :: nsr
type(respdata),     intent(in) :: resp(5,nsr)
type(param_forward),intent(in) :: g_param
type(obs_info),     intent(in) :: obs
character(50)                  :: filename
character(9)                   :: num
character(1)  :: nn
character(50) :: header,xyheader
integer(4)    :: i,j,k
integer(4)    :: len1,len2 ! 2017.10.12
integer(4)    :: idev      ! 2017.10.12

!#[1]## set
header   = g_param%outputfolder
xyheader = g_param%xyfilehead   ! 2017.10.12
write(num,'(e9.2)') freq
len1 = len_trim(header)
len2 = len_trim(xyheader)

do k=1,nsr
 write(nn,'(i1)') k
 filename=header(1:len1)//xyheader(1:len2)//num(2:9)//"_S"//nn//".dat"
 idev = free_unit() ! 2017.10.12
 open(idev,file=filename)
 do i=1,obs%nobs ! 2017.10.12
  write(idev,'(12e15.7)') (obs%xyz_obs(j,i),j=1,2),    &
  & (resp(j,k)%ftobsamp(i),resp(j,k)%ftobsphase(i),j=1,5)
 end do
 close(idev) ! 2017.10.12
end do

return
end subroutine

!######################################### OUTOBSFILESFWD 2017.07.11
!# coded on 2017.07.11
subroutine OUTOBSFILESFWD(g_param,sparam,nsr,tresp,nfreq)
implicit none
type(param_forward),intent(in) :: g_param
type(param_source), intent(in) :: sparam
integer(4),         intent(in) :: nsr
integer(4),         intent(in) :: nfreq
type(respdata),     intent(in) :: tresp(5,nsr,nfreq)
character(50)  :: head, sour, site
character(70)  :: filename
integer(4)     :: i,j,k,l,nhead,nsite,nsour,nobs
real(8)        :: freq

!#[1]## set
 nobs  = g_param%nobs
 head  = g_param%outputfolder
 nhead = len_trim(head)

!#[2]##
 do l=1,nobs
  site  = g_param%obsname(l)
  nsite = len_trim(site)

  do k=1,nsr
   sour     = sparam%sourcename(k)
   nsour    = len_trim(sour)
   filename = head(1:nhead)//site(1:nsite)//"_"//sour(1:nsour)//".dat"
   open(31,file=filename)

   do i=1,nfreq
    freq = g_param%freq(i)
    write(31,'(11g15.7)') freq,(tresp(j,k,i)%ftobsamp(l),&
    &                           tresp(j,k,i)%ftobsphase(l),j=1,5)
   end do

   close(31)
  end do
 end do

 write(*,*) "### OUTOBSFILESFWD END!! ###"

return
end subroutine
!######################################### OPENOBSFILES
subroutine OPENOBSFILES(files)
implicit none
type(obsfiles),intent(in) :: files
integer(4) :: i

do i=1,files%nfile
 open(files%devicenumber(i),file=files%filename(i))
end do

return
end subroutine

!######################################### CLOSEOBSFILES
subroutine CLOSEOBSFILES(files)
implicit none
type(obsfiles),intent(in) :: files
integer(4) :: i

do i=1,files%nfile
 close(files%devicenumber(i))
end do

return
end subroutine

end module outresp

