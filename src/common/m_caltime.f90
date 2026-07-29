!# coded on 2017.09.06
module caltime
implicit none ! 2025.06.24
integer, parameter :: ikind = selected_int_kind(9) ! 2025.06.24

type watch
 integer(ikind) :: t1
 integer(ikind) :: t2
 integer(ikind) :: t3
 integer(ikind) :: t_rate ! 2025.06.25
 integer(ikind) :: t_max
 real(8)    :: time  ! [min]
 integer(ikind) :: imin  ! 2025.09.17
 integer(ikind) :: ihour ! 2025.09.17
 real(8)        :: sec  ! 2025.09.17
end type

contains
!################ start
subroutine watchstart(t_watch)
implicit none
type(watch),intent(inout) :: t_watch

call system_clock(t_watch%t1,t_watch%t_rate,t_watch%t_max)

return
end subroutine watchstart

!################ watch sttop
subroutine watchstop(t_watch)
implicit none
type(watch),intent(inout) :: t_watch

call system_clock(t_watch%t2)

call calt(t_watch)

return
end subroutine watchstop

!############################################################# caltime
! 2017.07.21
subroutine calt(t_watch)
implicit none
type(watch),intent(inout) :: t_watch
integer(ikind) :: t2,t1,t_max,t_rate
real(8)      :: time ! [min]
real(8)      :: diff

!#[1]# set
  t2     = t_watch%t2
  t1     = t_watch%t1
  t_rate = t_watch%t_rate
  t_max   = t_watch%t_max

!#[2]# cal time
  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  time = diff/dble(t_rate)/60.d0 ! [min]

!#[3]# output
  t_watch%time = time !"[min]" 2025.06.24
  t_watch%imin = int(time) ! [min] 2025.09.17
  t_watch%ihour= int(time)/60  ! [hour] 2025.09.17
  t_watch%imin = mod(int(time),60) ! [min] 2025.09.17
  t_watch%sec = (time - int(time))*60  ! [sec] 2025.09
return
end subroutine calt

end module caltime
