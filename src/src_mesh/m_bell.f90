!# coded on 2017.09.28
module bell
use param ! 2017.09.28
implicit none

type param_bell ! inserted on 2016.10.12 for Kusatsushirane
!# NOTE zthick should be less than width !!! 2017.09.28
 character(50) :: bellgeofile
 integer(2)    :: nbell       ! # of bells 2018.10.02
 integer(4),allocatable,dimension(:)   :: ilonlatflag ! 1 for lonlat, 2 for xy [km]
 real(8),   allocatable,dimension(:)   :: lon_bell, lat_bell
 real(8),   allocatable,dimension(:)   :: radius![km]
 real(8),   allocatable,dimension(:)   :: width ![km]
 real(8),   allocatable,dimension(:)   :: ztop  ! < 0 [km] how deep is the top of bell
 real(8),   allocatable,dimension(:)   :: zbot  ! < 0 [km] how deep is the bottom of bell
 real(8),   allocatable,dimension(:)   :: zthic ! [km] thickness of top of the bell
 real(8),   allocatable,dimension(:)   :: reso_bell ! [km] mesh reso for inside bell
 real(8),   allocatable,dimension(:,:) :: xyz_bell
end type

contains

!#################################################################
! Coded2017.09.27
subroutine readparambell(g_param,b_param)
implicit none
type(param_forward),intent(in)  :: g_param
type(param_bell),   intent(out) :: b_param
integer(4)   :: input = 5
character(3) :: UTM
real(8)      :: lonorigin, latorigin
real(8)      :: lonlat(2)
integer(4)   :: nbell,ibell   ! 2018.10.02

!#[1]## set
lonorigin = g_param%lonorigin ! 2017.09.28
latorigin = g_param%latorigin ! 2017.09.28
UTM       = g_param%UTM       ! 2017.09.28

!#[2]## read mesh param
read(input,*)   ! header
read(input,'(20x,a)')      b_param%bellgeofile
 write(*,*) "bell geofile",b_param%bellgeofile
read(input,'(20x,i10)')    b_param%nbell      ! 2018.10.02
!===
 nbell = b_param%nbell                   ! 2018.10.02
 allocate( b_param%ilonlatflag(nbell) )  ! 2018.10.02
 allocate( b_param%lon_bell(   nbell) )  ! 2018.10.02
 allocate( b_param%lat_bell(   nbell) )  ! 2018.10.02
 allocate( b_param%radius(     nbell) )  ! 2018.10.02
 allocate( b_param%width(      nbell) )  ! 2018.10.02
 allocate( b_param%ztop(       nbell) )  ! 2018.10.02
 allocate( b_param%zbot(       nbell) )  ! 2018.10.02
 allocate( b_param%zthic(      nbell) )  ! 2018.10.02
 allocate( b_param%reso_bell(  nbell) )  ! 2018.10.02
 allocate( b_param%xyz_bell(3, nbell) )  ! 2018.10.02
!==
do ibell=1,nbell ! 2018.10.02

read(input,'(20x,i10)')    b_param%ilonlatflag(ibell)
 write(*,*) "bell lonlatflag",b_param%ilonlatflag(ibell)
if ( b_param%ilonlatflag(ibell) .eq. 1) then
 read(input,'(20x,2g12.5)') b_param%lon_bell(ibell), b_param%lat_bell(ibell)!2018.10.02
 write(*,*) "lonlat=",      b_param%lon_bell(ibell), b_param%lat_bell(ibell)
 lonlat(1:2) =(/b_param%lon_bell(ibell),b_param%lat_bell(ibell) /)
 write(*,*) "UTM zone:",UTM
 call UTMXY(lonlat,lonorigin,latorigin,b_param%xyz_bell(1:2,ibell),UTM)
     write(*,*) "xyz_bell(1:2)=",b_param%xyz_bell(1:2,ibell)
else
 read(input,'(20x,2g12.5)') b_param%xyz_bell(1:2,ibell)
end if

read(input,'(20x,g15.7)') b_param%radius(ibell)![km]
read(input,'(20x,g15.7)') b_param%width(ibell) ![km]
read(input,'(20x,g15.7)') b_param%ztop(ibell)  ! < 0 [km] how deep is the top of bell
read(input,'(20x,g15.7)') b_param%zbot(ibell)  ! < 0 [km] how deep is the bottom of bell
read(input,'(20x,g15.7)') b_param%zthic(ibell) ! [km] thickness of top of the bell
read(input,'(20x,g15.7)') b_param%reso_bell(ibell) ! [km] mesh resolution for inside bell
!close(input)

end do ! 2018.10.02

return
end subroutine ! 2020.12.08

end module bell
