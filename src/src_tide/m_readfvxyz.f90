module readFvxyz
  
implicit none
contains
  !##################################################

  subroutine readfxyz(nodes,node,fxyz,fxyzfile)! 2024.03.16
  use mesh_type

  implicit none

  integer(4),intent(in) :: nodes ! 2024.03.16 node for ocean_mesh
  integer(4),intent(in) :: node  ! 2024.03.16 node for em3d.msh
  real(8)               :: lon,lat,alt
  real(8)              ::fxyz(3,node)
  character(70)::fxyzfile    ! 70 2024.03.14
 
  integer(8)::i


  fxyz(:,:) = 0.d0

    open(10,file=fxyzfile)
    do i=1,nodes
       
       read(10,*)lat,lon,alt,fxyz(1,i),fxyz(2,i),fxyz(3,i)
    end do
    close(10)
    return
  end subroutine readfxyz

  !##################################################

  subroutine readvxyz(nodes,node,vxyz,vxyzfile)
  use mesh_type

  implicit none

  integer(4),intent(in) :: nodes ! 2024.03.16 node for ocean_mesh
  integer(4),intent(in) :: node  ! 2024.03.16 node for em3d.msh
  real(8)::lon,lat,alt
  complex(8)::vxyz(3,node)
  character(70)::vxyzfile!,oceanmeshfile ! 70 2024.03.14
  !character(70)::oceanmeshfile
  integer(8)::i

  vxyz(:,:) = 0.d0 ! 2024.03.16

   open(11,file=vxyzfile)

    do i=1,nodes

       read(11,*)lat,lon,alt,vxyz(1,i),vxyz(2,i)

       vxyz(3,i)=(0,0)

    end do

    close(11)

  return
  end subroutine readvxyz

  !###################################################

end module readfvxyz



  

    
  
