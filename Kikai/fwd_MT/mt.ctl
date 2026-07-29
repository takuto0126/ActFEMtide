## lines starting with "##" work as lines for comments ! 2021.06.08
## this file is forward control file
##-----10!-------20!
## itopofile controls whether topo file is provided or not
## itopofile = 0   : no topo file and z=0 surface is assumed
##                   obs and src coordinates should be provided by cartesian position
## itopofile = 1   : topofile(s) is(are) provided, iflag_map should be provided
   itopofile       !1
## iflag_map controls the type of Map projection
## iflag_map = 1 (ECP) : lon lat topo + Equidistant Cylindrical Projection  (ECP)
##                 x=planetrad*cos(latorigin)*(lon-lonorigin), y=planetrad*(lat-latorigin)
##                 one additional line: lonlat origin (2f15.7) 
## iflag_map = 2 (UTM) : lon lat topo + Universal Transvers Mercatol (UTM) projection
##                 two additional line: UTM zone like 52S (a3) and lonlat origin (2f15..7)
iflag_map          !1
lonlatorigin       !131.084782     32.884882
# of topofile      !1
topofile           !../mesh/topo.dat
lon lat shift      !0.0         0.0
mesh file          !../mesh/em3d.msh
2d triangle z file !../mesh/polygonz.msh
local line file    !../mesh/lineinfo.dat
angle              !0.0
output folder      !./result/
header2d  (a50)    !nakadake2d
header3d  (a50)    !nakadake3d
# of frequency     !3
Frequency [Hz]     !0.001d0
Frequency [Hz]     |0.03d0
Frequency [Hz]     !0.01d0
Frequency [Hz]     !0.03d0
Frequency [Hz]     !0.1d0
west bound         !-1.7
east bound         !1.7
south bound        !-1.5
north bound        !1.5
lenout [km]        !25.0
upz in [km]  (>0)  !1.3
downz in [km](<0)  !-1.1
zmax   [km]        !50.0
zmin   [km]        !-50.0
sizein [km]        !0.10
sizebo [km]        !2.0
sigma_obs [km]     !0.4
A_obs     [km]     !0.01
dlen_source [km]   !0.1
sigma_src [km]     !0.3
A_src     [km]     !0.005
# of observatory   !14
lonlat(1),xyz (2)  !1
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
ixyflg 0:no,1:surfv!0
sigma_air    [S/m] !1.e-8
condflag 0:homo,1: !1
##Volume             !1
##cond               !0.01
condfile           !../structure/cond_test.msh

