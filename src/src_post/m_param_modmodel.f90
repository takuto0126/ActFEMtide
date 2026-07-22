!###################################################
module param_modmodel

type param_mod
 integer(4)    :: icondflag
 character(50) :: condfile
 character(50) :: modelconn
 character(50) :: modelfile
end type

contains
!###################################################
subroutine readparammod(g_param_mod)
implicit none
type(param_mod),intent(out) :: g_param_mod
integer(4)                  :: input=5

read(input,10) g_param_mod%icondflag
if (    g_param_mod%icondflag .eq. 1) then
  write(*,*) "input path for reference conductivity structure file"! 2017.07.19
  read(input,20) g_param_mod%condfile ! ref cond file    2017.07.19
  write(*,*) "refcondfile=",g_param_mod%condfile              ! 2017.07.19

elseif ( g_param_mod%icondflag .eq. 2) then
  write(*,*) "Input modelconnect file and modelfile"               ! 2018.03.18
  read(input,'(20x,a)') g_param_mod%modelconn ! ref cond file   2018.03.18
  read(input,'(20x,a)') g_param_mod%modelfile ! ref cond file   2018.03.18
  write(*,*) "refmodelconnect",g_param_mod%modelconn          ! 2018.03.18
  write(*,*) "refmodelfile",   g_param_mod%modelfile          ! 2018.03.18else

else
 write(*,*) "GEGEGE"
 stop
end if

return
10 format(20x,i10)
20 format(20x,a)
end subroutine

end module param_modmodel

