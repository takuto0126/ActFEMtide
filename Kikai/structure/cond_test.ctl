!------10!-------20!----
input 3d mshfile   !../mesh/em3d.msh
0homo,1cond,2model !0
homo resistivity   !100.0
output cond        !./cond_test.msh
0:elevation,1:depth!0
# of cuboid        !2
1  xminmax [km]    !  -50.0          -10.0
1  yminmax [km]    !  -50.2          50.2
1 minmax z [km]    !  -50.0          -5.0
2  rho    [Ohm.m]  ! 10.0
2  xminmax [km]    !  10.0           50.0
2  yminmax [km]    !  -50.2          50.2
2 minmax z [km]    !  -50.0          -5.0
2 rho              ! 1000.0
