program hankel
  use omega_lib
  implicit none
  real(8) i,err
  real(8),external :: f
  integer k
  do k=1,10
     omega=dble(k)*0.1d0
     call intdeo(f, 0.d0, omega, 1.d-10, i, err)
     write(*,*) (i-exp(-omega**2*0.5d0))/i
  end do
end program hankel

function f(x)
  use omega_lib
  implicit none
  real(8),intent(in)::x
  real(8)f
  real(8),external::dbesj0
  f=exp(-x**2*0.5d0)*dbesj0(omega*x)*x
end function f
