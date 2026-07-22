
!# coded on 2022.10.24
program changecond
use param
use mesh_type
use cond_type ! see m_cond_type.f90
implicit none
type(param_cmod_surf) :: g_param_surf ! see m_cond_type.f90 ! 2022/.10.24
character(100)        :: inmshfile,topomshfile
type(mesh)            :: h_mesh,topo_mesh
type(param_cond)      :: g_cond

!#[1]## read control
 call readcmodparam_surf(g_param_surf) ! see m_cond_type.f90

!#[2]## read mesh
 inmshfile       = g_param_surf%mshfile
 topomshfile     = g_param_surf%topo_2d_mshfile
 CALL READMESH_TOTAL(h_mesh,inmshfile)
 CALL READMESH_TOTAL(topo_mesh,topomshfile)

!#[3]## read input cond
 g_cond%condfile = g_param_surf%inputcondfile
 CALL READCOND(g_cond)          ! read conductivity structure



end program
