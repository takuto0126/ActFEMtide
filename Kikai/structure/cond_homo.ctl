!------10!-------20!----
input 3d mshfile   !../mesh/em3d.msh
0homo,1cond,2model !0
homo resistivity   !100.0
output cond        !./cond_homo.msh
0:elevation,1:depth!0
# of cuboid        !1
1  xminmax [km]    !  -50.0          50.0
1  yminmax [km]    !  -50.2          50.2
1 minmax z [km]    !  -50.0          -5.0
1  rho    [Ohm.m]  ! 100.0
