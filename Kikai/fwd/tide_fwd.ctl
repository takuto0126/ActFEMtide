## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is forward control file
##-----10!-------20!
itopofile 0 or 1   !1
# of topofile      !1
topofile           !../topo/etopo_kikai-l.xyz
lon lat shift      !0.0         0.0
mesh file          !../mesh/em3d.msh
2d triangle z file !../mesh/polygonz.msh
local line file    !../mesh/lineinfo.dat
Ocean mesh file.   !../mesh/ocean.msh
mesh ctl file      !../mesh/mesh.ctl
Vxyzfile           !../fvxyz/vxyz_mesh
Fxyzfile           !../fvxyz/fxyz_mesh
output folder      !./result/
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !1
Frequency [Hz]     !0.00002237
west bound         !-200
east bound         !200
south bound        !-200
north bound        !200
lenout [km]        !100.0
upz in [km]  (>0)  !500
downz in [km](<0)  !-500
zmax   [km]        !1000
zmin   [km]        !-1000
sizein [km]        !0.10
sizebo [km]        !2.0
sigma_obs [km]     !0.4
A_obs     [km]     !0.01
dlen_source [km]   !0.1
sigma_src [km]     !0.3
A_src     [km]     !0.005
# of observatory   !14
lonlat(1),xyz (2)  !0
UTM ZONE           |52N
lonlatorigin       !130.500000     30.500000
12  Name           !P5
12  xyz            !130.3583       30.6500      -0.001
11  Name           !P6
11  xyz            !130.1833   30.8400  -0.001
11  Name           !P7
11  xyz            !130.5500   30.6666  -0.001
11  Name           !P8
11  xyz            !130.3150   30.8500  -0.001
11  Name           !P9
11  xyz            !130.2833   31.0833  -0.001
12  Name           !P10
12  xyz            !130.5333   31.0166  -0.001
11  Name           !P11
11  xyz            !130.0833   30.8833  -0.001
11  Name           !P12
11  xyz            !130.0833   30.7150  -0.001
11  Name           !P13
11  xyz            !130.0866   30.3783  -0.001
11  Name           !P14
11  xyz            !129.6220   31.0871  -0.001
12  Name           !P15
12  xyz            !129.7143   31.0531  -0.001
12  Name           !P16
12  xyz            !129.8270   31.0030  -0.001
12  Name           !P17
12  xyz            !129.9753   30.9331  -0.001
12  Name           !P18
12  xyz            !130.6743   30.6226  -0.001
0:no,1:nxny,2:ntrik!1
file head          !bxyz_xy2D
nx, ny             !200        200
1:same z,2:seafloor!2
ki23dfile          !../mesh/polygonki.msh
nx,ny->ki23dptr fl !../mesh/ki23dptr.dat
# of sources       !1
Source Name        !S1
source start point !189.354722   -24.33944  -0.001
source end   point !189.254722   -24.33944  -0.001
Elcetric current[A]!0.0
sigma_air    [S/m] !1.e-8
condflag 0:home,1: !1
##nvolume            !1
##cond               !0.01
condfile           !../structure/cond_test.msh

