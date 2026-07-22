
! coded on May 20, 2016
!
module obs_type
use matrix
implicit none

type obs_info
 integer(4)                          :: nobs
 character(4)                        :: name
 real(8),allocatable,dimension(:,:)  :: xyz_obs
 type (real_crs_matrix)              :: coeff(2,3)
 integer(4),allocatable,dimension(:) :: devnum     ! 2017.10.12
 ! from tmtgem
 real(8),pointer,dimension(:,:)  :: xyzspherical ! 2016.11.20
 real(8),pointer,dimension(:,:)  :: lonlatalt    ! 2016.11.20
 type (real_crs_matrix)          :: coeff_vF
 integer(4),pointer,dimension(:) :: devnum_b(:)  ! 2021.07.26
 integer(4),pointer,dimension(:) :: devnum_e(:)  ! 2021.07.26


end type


end module obs_type
