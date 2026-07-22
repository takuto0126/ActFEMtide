! coded on 2017.06.07
module freq_mpi
use param
implicit none

type freq_info
 integer(4) :: nfreq
 integer(4) :: nfreq_ip
 real(8),   allocatable,dimension(:) :: freq
 real(8),   allocatable,dimension(:) :: freq_ip
 integer(4),allocatable,dimension(:) :: ifreq_ip
end type

contains

!#############################################
!# Coded on 2017.06.05
subroutine SETFREQIP(g_param,ip,np,g_freq)
implicit none
integer(4),            intent(in)    :: ip,np
type(param_forward),   intent(in)    :: g_param
type(freq_info),       intent(inout) :: g_freq
real(8),   allocatable,dimension(:)  :: freq,freq_ip
integer(4),allocatable,dimension(:)  :: ifreq_ip
integer(4)                           :: i,ii
integer(4)                           :: nfreq,nfreq_ip

!#[1]## set

  nfreq    = g_param%nfreq
  nfreq_ip = g_param%nfreq/np
  if ( nfreq - nfreq_ip*np .gt. 0 ) nfreq_ip = nfreq_ip + 1
  allocate(  freq_ip(nfreq_ip) )
  allocate( ifreq_ip(nfreq_ip) )
  allocate(  freq   (nfreq)    )

!#[2]## cal freq, freq_ip, ifreq_ip
 freq       = g_param%freq
 freq_ip(:) = -1.d0 ! default
 do i=1,nfreq_ip
  ii = ip+1 + (i-1)*np
  if (ii .le. nfreq)  freq_ip(i) = freq(ii) ! freqency
  if (ii .le. nfreq)  ifreq_ip(i) = ii       ! frequency index
 end do
! do i=1,nfreq_ip
!  write(*,*) "ip=",ip,"i=",i,"freq_ip(i)=",freq_ip(i)
! end do

!#[3]## set output
  g_freq%nfreq    = nfreq
  g_freq%nfreq_ip = nfreq_ip
  allocate( g_freq%freq_ip(nfreq_ip) )
  allocate( g_freq%ifreq_ip(nfreq_ip) )
  allocate( g_freq%freq   (nfreq)    )
  g_freq%freq     = freq
  g_freq%freq_ip  = freq_ip
  g_freq%ifreq_ip = ifreq_ip

return
end subroutine



end module freq_mpi
