!# Parameters for OMSP/hankel_takuto.f90
!# Created on Dec. 1, 2015
module m_param_ana
implicit none
!### Please set frep, nlayer, cond, h, I before running program ###
!character(50) ::  outfolder
!integer(4), parameter :: nlayer = 3  ! # of layers beneath the surface
!real(8)  :: cond(2,0:nlayer) ![S/m] conductivity of i-th layer
!real(8)  ::    h(0:nlayer) ! [m] thickness of i-th layer beneath the surface
!data h / 999999.d0,200.d0,800.d0, 999999.d0/ ! [m] h(0) and h(nlayer) wont be used
!data cond(1,0:nlayer) / 1.d-8, 0.01d0,0.01d0,0.01d0 /! cond(0) is air conductivity,
!data cond(2,0:nlayer) / 1.d-8, 1.d-8,1.d-8,1.d-8 /! cond(nlayer) is the lowermost
real(8),parameter :: pi=4.d0*datan(1.d0),r2d=180./pi
real(8),parameter :: mu=4.d0*pi*1.d-7
real(8),parameter :: epsiron=8.85418782*1.d-12 ! [m^-3*kg^-1*s^4*A^2] electric permittivity
real(8)           :: omega
integer(4)                       :: nlayer ! 2021.12.20
real(8),allocatable,dimension(:) :: h      ! h(0:nlayer) 2021.12.20
real(8),allocatable,dimension(:) :: cond   ! cond(0:nlayer) 2021.12.20

!###############################################################
! global variables
complex(8),parameter :: iunit=(0.d0,1.d0)
integer(4) :: iflag ! 1: primary, 2: secondary, 3:total
integer(4) :: iflag_realimag ! 1:real, 2:imaginary
complex(8) :: RTE,RTM,u0,u1,G,G2
real(8)    :: xobs      ! [m] Observation line
real(8)    :: yobs,R,R1,R2
real(8)    :: zobs=0.d0     ! [m] z should be equally larger than 0
complex(8) :: y0hat,z0hat ! calculated in CALRTEM
integer(4) :: istructure
!real(8)    :: L0=1.d0 ! Scale length [m]
character(6) :: FNAME ! added on April 20, 2016 for hankel_takuto.f90
integer(4) :: jflag ! 0 for J0, 1 for J1 on April 20, 2016, for hankel_takuto.f90

contains

!# NOTE!! ###
! This function assumes nlayer = 1 or
! z <= 0.d0 i.e. in the air
! if not, please modify hankel_takuto.f90 and this subroutine
complex(8) function expuz()
! z is downward positive in Ward and Hohmann (1988)

if (zobs .lt. 0.d0) expuz=cdexp( u0*zobs) ! above the surface
if (zobs .ge. 0.d0) expuz=cdexp(-u1*zobs)  ! below the ground surface

return
end function expuz

end module m_param_ana
