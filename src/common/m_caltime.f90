!# coded on 2017.09.06
module caltime

type watch
 integer(4) :: t1
 integer(4) :: t2
 integer(4) :: t3
 integer(4) :: t_rate
 integer(4) :: t_max
 real(8)    :: time  ! [min]
end type

contains
!################ start
subroutine watchstart(t_watch)
implicit none
type(watch),intent(inout) :: t_watch

call system_clock(t_watch%t1)

return
end subroutine watchstart

!################ start
subroutine watchstop(t_watch)
implicit none
type(watch),intent(inout) :: t_watch

call system_clock(t_watch%t2,t_watch%t_rate,t_watch%t_max)
call calt(t_watch)

return
end subroutine watchstop

!############################################################# caltime
! 2017.07.21
subroutine calt(t_watch)
implicit none
type(watch),intent(inout) :: t_watch
integer(4)   :: t1,t2,t_rate,t_max
real(8)      :: time ! [min]
real(8)      :: diff

!#[1]# set
  t2     = t_watch%t2
  t1     = t_watch%t1
  t_max   = t_watch%t_max
  t_rate = t_watch%t_rate

!#[2]# cal time
  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  time = time + diff/dble(t_rate)/60.d0

!#[3]# output
  t_watch%time = time

return
end subroutine calt

end module caltime
