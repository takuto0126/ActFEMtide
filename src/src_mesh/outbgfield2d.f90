!# Coded on 2016.09.28
subroutine outbgfield2d(g_param)
use param
implicit none
type(param_forward),intent(in)     :: g_param
character(50) :: header2d,fieldfile,pregeo
real(8)       :: sizein,sizebo
real(8),dimension(4) :: xbound,ybound, zbound
integer(4)                         :: nobsr
real(8),allocatable,dimension(:,:) :: xyz_r      ! 2017.09.08
real(8),allocatable,dimension(:)   :: sigma_r,A_r! 2017.09.08
integer(4) :: i

!#[1]set
sizein  = g_param%sizein
sizebo  = g_param%sizebo
xbound  = g_param%xbound
ybound  = g_param%ybound
zbound  = g_param%zbound
nobsr   = g_param%nobsr            ! 2017.09.08
allocate( xyz_r(3,nobsr)            ) ! 2017.09.08
allocate( sigma_r(nobsr),A_r(nobsr) ) ! 2017.09.08
A_r     = g_param%A_r ! [km]
sigma_r = g_param%sigma_r ! [km]
xyz_r   = g_param%xyz_r

!#[2] files
header2d   = g_param%header2d
fieldfile  = header2d(1:len_trim(header2d))//"_field.geo"
pregeo     = header2d(1:len_trim(header2d))//".geo"


!#[2]## output field file
open(1,file=fieldfile)

!write(1,'(a,a,a)') "Include '",pregeo(1:len_trim(pregeo)),"';"
write(1,'(a)') "Field[1]      = Box;"
write(1,10) "Field[1].VIn  =", sizein ,";"
write(1,10) "Field[1].VOut =", sizebo ,";"
write(1,10) "Field[1].XMax =", xbound(3),";"
write(1,10) "Field[1].XMin =", xbound(2),";"
write(1,10) "Field[1].YMax =", ybound(3),";"
write(1,10) "Field[1].YMin =", ybound(2),";"
write(1,10) "Field[1].ZMax =", zbound(3),";"
write(1,10) "Field[1].ZMin =", zbound(2),";"

write(1,'(a)') "Field[2]   = MathEval;"
write(1,'(a)') 'Field[2].F = "2.7^((x^2. + y^2. + z^2.)/(10.^2))";'

do i=1,nobsr
 write(1,*) "Field[",i+2,"]   = MathEval;"
 write(1,*) 'Field[',i+2,'].F ="',sizein,' - (',sizein,' -(', A_r(i),'))*'
 write(1,10) '2.718^(-((x-(',xyz_r(1,i),'))^2.  '
 write(1,10)       ' + (y-(',xyz_r(2,i),'))^2.)/'
 write(1,10)       '2./(',sigma_r(i),')^2.)";'
end do

write(1,11) "Field[",nobsr+3,"]   = Min;"
write(1,*) "Field[",nobsr+3,"].FieldsList = {1,2",(",",i+2,i=1,nobsr),"};"
write(1,11) "Background Field = ",nobsr+3,";"

close(1)

return
10 format(a,g15.7,a)
11 format(a,i10,a)

end subroutine outbgfield2d


