! Modified for solution of long grounded wire on 2016.10.22
! See the detail of solution in p.237 of Ward and Hohann (1988)
!
! Coded on April 7, 2016
! to conduct Hankel transform with OMSP,
! Omura's Mathmatical Source Package
! Found that Ooura's package is much superior to Anderson's code!!!! 2016.04.07

!program hankel_takuto
subroutine dipole_l(xs1,xs2,I,ds,obsxy,exy,bxyz,freq)
use m_param_ana
implicit none
real(8),intent(in)  :: freq,xs1(2),xs2(2) ! source [m]
real(8),intent(in)  :: I  ! current [A]
real(8),intent(in)  :: ds ! source division scale [m]
real(8),intent(in)  :: obsxy(2) ! site location
real(8),intent(out) :: exy(2,2) ! amp and phase of electric field [uV/m, deg]
real(8),intent(out) :: bxyz(3,2)! amp and phase of magnetic field [nT, deg]
!
real(8) bz2(2),by0(2),bx2(2),err,ex2(2),ey2(2)
complex(8) :: bz,by,bxj1(2),bx,ex0,exj1(2),eyj1(2),ey0,BYJ1(2),ex_homo,k2,sqk2
real(8),external :: fz,fyj1,fyj0,fx,fexj0,fexj1,feyj1,funall
real(8) :: amp(3)   ! [nT]
real(8) :: phase(3) ! [deg]
integer :: k,nobs,j,m,nsource
integer(4),dimension(2) :: isign=(/1,-1/)
real(8) :: coef, coef1,coef2,coef3
real(8) :: D=0.01 ! 1 cm
integer(4) :: is
real(8) :: v(2),va,x0(2),x1(2),sin,cos,L
complex(8) :: ex_tmp,ey_tmp,bx_tmp,by_tmp
real(8),allocatable,dimension(:) :: x_source
complex(8) :: bz_total,by_total,bx_total,ex_total,ey_total
COMPLEX(8) :: FBZ,FBX,FBY0,FBY1,FEX0,FEX1,FEY1, ZH, ZHANKS,FEX1v1,FEX1v2
external FBZ,FBX,FBY0,FBY1,FEX0,FEX1,FEY1,FEX1v1,FEX1v2


!#[0]## set parameters
L=sqrt((xs1(1)-xs2(1))**2. + (xs1(2)-xs2(2))**2.)/2.d0 ! half length of dipole
omega = 2.*pi*freq
nsource = 2*L/ds ! number of sub dipole for numerical integration
allocate(x_source(nsource))
do is=1,nsource
 x_source(is) = -L + ds/2.d0 + ds*(is-1)
end do
!write(*,*) "length of dipole=",ds,"[m]"

!#[1]## make coefficient, I*ds*mu/4/pi
 y0hat= cond(0) + iunit*omega*epsiron  ! 2021.12.20
 z0hat= iunit*omega*mu

coef  = I*ds*mu/4./pi*1.d+9  ! [nT]
coef1 = I*mu/4./pi*1.d+9  ! [nT]
coef2 = I*ds/4.d0/pi *1.d+6  ! [mV/km]
coef3 = I/4.d0/pi *1.d+6  ! [mV/km]
!write(*,*) "coef=",coef

!#[2]## observation info and allocate
! v  = xs2(1:2) - xs1(1:2)
! x0 = (xs2(1:2) + xs1(1:2))/2.
! xobs = (  cos(theta)  sin(theta) ) (obsxy(1) - x0(1))
! yobs = ( -sin(theta)  cos(theta) ) (obsxy(2) - x0(2))
! sin(theta) = (1,0) times v/|v|
! cos(theta) = (1,0) dot   v/|v|

v(1:2)=xs2(1:2) -xs1(1:2)
va=sqrt(v(1)**2.d0 + v(2)**2.d0)
x0(1:2)=(xs2(1:2)+xs1(1:2))/2.d0
sin=v(2)/va
cos=v(1)/va
x1(1:2)=obsxy(1:2) - x0(1:2)
xobs   =  cos*x1(1) + sin*x1(2)
yobs   = -sin*x1(1) + cos*x1(2)

!#[3]## start calculation
 ! BZ
 FNAME="FBZ ";jflag=1
 if ( R .eq. 0.d0 ) R=0.001d0 ! 1 mm
 iflag=3 ! 1,3
    bz_total = 0.d0
    do is = 1, nsource
     R=dsqrt((xobs-x_source(is))**2.d0+yobs**2.d0) ![m]
     do iflag_realimag=1,2 ! 1 for real, 2 for imaginary
      call intdeo(funall, 0.d0, R, 1.d-10, bz2(iflag_realimag), err)
     end do
     bz=coef*yobs/R*(bz2(1) + iunit*bz2(2))
     bz_total = bz_total + bz
    end do

 ! BY
 iflag=3 ! 1,3
    FNAME="FBY1";jflag=1
    !# xobs + L, -L
    by_total = 0.d0
    do j=1,2
     R=dsqrt((xobs+isign(j)*L)**2.d0+yobs**2.d0) ! R(xobs+L) - R(xobs -L)
     do iflag_realimag=1,2 ! 1 for real, 2 for imaginary
      call intdeo(funall, 0.d0, R, 1.d-10, by0(iflag_realimag), err)
     end do
     BYJ1(j) = (xobs+isign(j)*L)/R*(by0(1)+iunit*by0(2)) ! 1 for first order Bessel
    end do
    by_total = - coef1*(BYJ1(1)-BYJ1(2)) ! [nT]
    !#[5-2] J0
    FNAME="FBY0";jflag=0
    do is = 1,nsource
     R=dsqrt((xobs-x_source(is))**2.d0+yobs**2.d0)
     if ( R .eq. 0.d0 ) R=0.001d0 ! 1 mm
     do iflag_realimag=1,2 ! 1 for real, 2 for imaginary
      call intdeo(funall, 0.d0, R, 1.d-10, by0(iflag_realimag), err)
     end do
     by_total = by_total - coef*(by0(1)+ iunit*by0(2))
    end do
    by_total= - by_total

 ! BX
 FNAME="FBX ";jflag=1
 iflag=3
    bx_total = 0.d0
    do j=1,2
     R=dsqrt((xobs+isign(j)*L)**2.d0+yobs**2.d0) ! D is replaced by L
     do iflag_realimag=1,2 ! 1 for real, 2 for imaginary
      call intdeo(funall, 0.d0, R, 1.d-10, bx2(iflag_realimag), err)
     end do
     bxj1(j)=coef1*yobs/R*(bx2(1) + iunit*bx2(2))
    end do
    bx_total= - (bxj1(1) - bxj1(2))

 !Ex
    iflag=3
    ex_total=0.d0
    !#[J1]
    FNAME="FEX1";jflag=1
    do j=1,2
     R=dsqrt((xobs+isign(j)*L)**2.d0+yobs**2.d0)
     do iflag_realimag=1,2 ! 1 for real, 2 for imaginary
      call intdeo(funall, 0.d0, R, 1.d-10, ex2(iflag_realimag), err)
     end do
     exj1(j)= (xobs+isign(j)*L)/R*(ex2(1) + iunit*ex2(2))
    end do
    ex_total = - coef3*(exj1(1) - exj1(2)) ! [uV/m] dQ/dx|^R1_R2
    !#[J0]
    FNAME="FEX0";jflag=0
    do is =1,nsource
     R=dsqrt((xobs-x_source(is))**2.d0+yobs**2.d0)
     if ( R .eq. 0.d0 ) R=0.001d0 ! 1 mm
     do iflag_realimag=1,2 ! 1 for real, 2 for imaginary
      call intdeo(funall, 0.d0, R, 1.d-10, ex2(iflag_realimag), err)
     end do
     ex_total = ex_total - z0hat*coef2*(ex2(1) + iunit*ex2(2)) !
    end do

  !
  k2=-z0hat*y0hat
  sqk2=cdsqrt(k2)
  r=dsqrt(xobs**2.d0+yobs**2.d0)
  ex_homo=coef2*cdexp(-iunit*sqk2*r)*(-z0hat/r -iunit/r/r*cdsqrt(-z0hat/y0hat)-1./r/r/r/y0hat)
  CALL amphase(ex_homo,amp(1),phase(1))
!  write(16,*) yobs,amp(1),phase(1)

 !Ey
  iflag=3
  ey_total = 0.d0
  !#[J1]
  FNAME="FEY1";jflag=1
  do j=1,2
   R=dsqrt((xobs+isign(j)*L)**2.d0+yobs**2.d0) ! [m]
   do iflag_realimag=1,2
    call intdeo(funall,0.d0,R,1.d-10,ey2(iflag_realimag),err)
   end do
   eyj1(j)=yobs/R*(ey2(1)+iunit*ey2(2))
  end do
  ey_total = - coef3*(eyj1(1) - eyj1(2))

!#[3] ## coordinate rotation
ex_tmp = cos*ex_total - sin*ey_total
ey_tmp = sin*ex_total + cos*ey_total
bx_tmp = cos*bx_total - sin*by_total
by_tmp = sin*bx_total + cos*by_total
CALL amphase(ex_tmp,exy(1,1),exy(1,2))
CALL amphase(ey_tmp,exy(2,1),exy(2,2))
CALL amphase(bx_tmp,bxyz(1,1),bxyz(1,2))
CALL amphase(by_tmp,bxyz(2,1),bxyz(2,2))
CALL amphase(bz_total,bxyz(3,1),bxyz(3,2))

return
end subroutine dipole_l
!=========================================================== amphase
subroutine amphase(cmp,amp,phase)
use m_param_ana
implicit none
complex(8),intent(in) :: cmp
real(8),intent(out) :: amp,phase
 amp=dsqrt(dreal(cmp)**2.d0 + dimag(cmp)**2.d0)
 phase=datan2(dimag(cmp), dreal(cmp))*r2d
 return
end

!################################### funall
function funall(G0)
 use m_param_ana
 implicit none
 real(8),intent(in) :: G0
 real(8) funall
 real(8),external :: dbesj0,dbesj1
 COMPLEX(8) :: FBZ,FBX,FBY0,FBY1,FEX0,FEX1,FEY1, FEX1v1,FEX1v2
 complex(8) FUN
 !#[1]# choose function and calculate
  if ( FNAME .eq. "FBZ   ") FUN=FBZ(G0)
  if ( FNAME .eq. "FBX   ") FUN=FBX(G0)
  if ( FNAME .eq. "FBY0  ") FUN=FBY0(G0)
  if ( FNAME .eq. "FBY1  ") FUN=FBY1(G0)
  if ( FNAME .eq. "FEX0  ") FUN=FEX0(G0)
  if ( FNAME .eq. "FEX1  ") FUN=FEX1(G0)
  if ( FNAME .eq. "FEY1  ") FUN=FEY1(G0)

 !#[2]# real or imag
  if (iflag_realimag .eq. 1)  funall=dreal(FUN) ! real
  if (iflag_realimag .eq. 2)  funall=dimag(FUN) ! imaginary

 !#[3]# J0 or J1
  if (jflag .eq. 0)   funall=funall*dbesj0(R*G0) ! 0th order
  if (jflag .eq. 1)   funall=funall*dbesj1(R*G0) ! 1st order

end function funall

!# Coded on Dec. 1, 2015
!========================================================== FBX
!# see Ward and Hohmann (1988), p.233
complex(8) function FBX(G0)
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0) ! get G, G2, RTE, RTM, u0 from G0
  if (iflag .eq. 1) FBX=0.d0
  if (iflag .ne. 1) FBX= (RTM + RTE)*cdexp(u0*zobs)
  FBX=FBX*expuz() ! see m_param_ana.f90
return
end
!========================================================== FBY1
complex(8) function FBY1(G0)
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0) ! get G, G2, RTE, RTM, u0 from G0
  if (iflag .eq. 1) FBY1=0.d0
  if (iflag .ne. 1) FBY1= (RTM + RTE) ! RTM + RTE
  FBY1=FBY1*expuz()
return
end
!========================================================== FBY0
complex(8) function FBY0(G0)
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0) ! get G, G2, RTE, RTM, u0 from G0
  if (iflag .eq. 1) FBY0 = 1.0       *G
  if (iflag .eq. 2) FBY0= - RTE      *G
  if (iflag .eq. 3) FBY0= (1.0 - RTE)*G
  FBY0=FBY0*expuz()
return
end
!========================================================== FBZ
!# see Ward and Hohmann (1988), p.233
complex(8) function FBZ(G0) ! G0 is lambda
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0)
  if (iflag .eq. 1)  FBZ=        1.d0 *G2/u0
  if (iflag .eq. 2)  FBZ=         RTE *G2/u0
  if (iflag .eq. 3)  FBZ= (1.d0 + RTE)*G2/u0
  FBZ=FBZ*expuz()
return
end
!========================================================== FEX1
!# see Ward and Hohmann (1988), p.233
complex(8) function FEX1(G0)
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0)
  if (iflag .eq. 1)  FEX1=  1.d0      *u0/y0hat -  1.d0     *z0hat/u0
  if (iflag .eq. 2)  FEX1=      - RTM *u0/y0hat -       RTE *z0hat/u0
  if (iflag .eq. 3)  FEX1= (1.d0- RTM)*u0/y0hat - (1.d0+RTE)*z0hat/u0
  FEX1=FEX1*expuz()
return
end
!========================================================== FEX0
!# see Ward and Hohmann (1988), p.233
complex(8) function FEX0(G0)
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0)
  if (iflag .eq. 1)  FEX0=  1.d0       *G/u0
  if (iflag .eq. 2)  FEX0=         RTE *G/u0
  if (iflag .eq. 3)  FEX0= (1.d0 + RTE)*G/u0
  FEX0=FEX0*expuz()
return
end
!========================================================== FEY1
!# see Ward and Hohmann (1988), p.233
complex(8) function FEY1(G0)
use m_param_ana
real(8),intent(in) :: G0
  CALL CALRTEM(G0)
  if (iflag .eq. 1)  FEY1=  1.d0      *u0/y0hat -  1.d0     *z0hat/u0
  if (iflag .eq. 2)  FEY1=      - RTM *u0/y0hat -       RTE *z0hat/u0
  if (iflag .eq. 3)  FEY1= (1.d0- RTM)*u0/y0hat - (1.d0+RTE)*z0hat/u0
  FEY1=FEY1*expuz()
return
end

!===========================================================  RTE and RTM
subroutine CALRTEM(G0)
use m_param_ana
implicit none
real(8),intent(in)  :: G0
complex(8) :: Y_hat(0:nlayer), Y(0:nlayer),u(0:nlayer)
complex(8) :: Z_hat(0:nlayer), Z(0:nlayer), th, k2
integer(4) :: j
  G=G0 ; G2=G*G

!#[1]## coefficients, G2, u, Y=u/i/omega/mu
 do j=0,nlayer
  k2=-iunit*mu*omega*(cond(j)+iunit*omega*epsiron) ! 2021.12.20
  u(j)=cdsqrt(G2-k2)
  Y(j)=u(j)/(iunit*omega*mu)
  Z(j)=u(j)/(cond(j)+iunit*omega*epsiron) ! 2021.12.20
 end do
!#[2]## RTE
 Y_hat(nlayer)=Y(nlayer)
 Z_hat(nlayer)=Z(nlayer)
 do j=nlayer-1, 1, -1
  Y_hat(j)=Y(j)*(Y_hat(j+1)+Y(j)      *th(u(j)*h(j)))&
  &            /(Y(j)      +Y_hat(j+1)*th(u(j)*h(j)))
  Z_hat(j)=Z(j)*(Z_hat(j+1)+Z(j)      *th(u(j)*h(j)))&
  &            /(Z(j)      +Z_hat(j+1)*th(u(j)*h(j)))
 end do
 RTE=(Y(0)-Y_hat(1))/(Y(0)+Y_hat(1))
 RTM=(Z(0)-Z_hat(1))/(Z(0)+Z_hat(1))
 u0=u(0)
 u1=u(1)
return
end

!########################################################### tanh
      function th(a)
      implicit real(selected_real_kind(8))(a-h,o-z)
      complex(8) :: th,a
      th=(1.d0 - cdexp(-2.*a))/(1.d0 + cdexp(-2.*a))
      return
      end
