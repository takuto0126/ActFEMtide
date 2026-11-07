# modified on 2020.12.08
# coded on 2018.02.27
# to draw a slice
####################################################################################
## Explanation of slice.exe: (ActFEMv1.0/src/src_post/slice.f90)
## slice.exe can generate two types of slices
## one is vertical slice, the other is horizontal slice
## You can choose either of them by setting three coefficients, A, B, and D
##
## [1] For vertical slices
##     slice plane is designed by Ax + By + D = 0
##     x is eastward, y is northward, z is upward, sea level is z=0, x,y,z are in [km]
##     For example, A=1,B=0,D=-1 generate x = 1 km plane.
##     Note that A should be > 0 in slice.exe
##
## [2] For horizontal slices
##     You shold set A=0, B=0, and D not equal to 0.
##     Then slice is generated for the plane of z + D = 0 (i.e., z = - D [km] plane)
##     z is upward, sea level is z=0
##
##  In both above two cases, four boudary values should be set after setting of A,B,D
##  1. left   boundary (lb) [km]
##  2. right  boundary (rb) [km]
##  3. bottom boundary (bb) [km]
##  4. top    boundary (tb) [km]
##
##  lb [km]     rb [km]
##  v           v
##  |-----------| <- tb [km] (elevation)
##  |  cross    |
##  |  section  |
##  |-----------| <- bb [km] (elevation, 0 means sea level)
##
##  These four parameters control the width and height of rectangles cut out of the slice plane
##  The origin in the slice plane is assumed to be the point
##  where the origin of the 3-D mesh is normally prohjected to the slice plane
##  The vertical axis in the slice plane is upward, the same as z axis in 3-D mesh
##  The horizontal axis in the slice plane is rightward when you see the plane
##  from the origin of the original 3-D mesh
##
######################################################################################

#!/bin/bash

folder=./


ID=model$1
echo $ID
modelfile=result_inv/${ID}.dat
z_elevation=-5.0 # for horizontal slice [km], write with decimal point # 2026.11.07 CC
y_slice=40.0     # WE vertical slice y = y_slice [km] # 2026.11.07 CC
x_slice=30.0     # SN vertical slice x = x_slice [km] # 2026.11.07 CC
xmin=-100.0 ; xmax=100.0 # W and E boundaries [km] # 2026.11.07 CC
ymin=-100.0 ; ymax=100.0 # S and N boundaries [km] # 2026.11.07 CC
zmin=-120.0 ; zmax=0.0  # bottom and top of vertical slices [km] # 2026.11.07 CC
igendat=1 # 0: no dat file generation, 1: dat file generation

# CODE

#[1]## generate control file
cat << EOF > tmp0.ctl
##--------------- control file for slice.exe starts here!! ----------------------
## 2020.12.10 Takuto Minami
##
## Note that lines starting with "##" work for comments ---
##
##             20->|
## 2026.11.07 CC mesh is changed to Kikai mesh
mshfile            !../mesh/em3d.msh
## You should chose either of two types of resistivity structure file format,
## cond.msh type or model**.dat type
## 0 for cond.msh type, 1 for model**.dat type
##  case for cond.msh type
##0:cond,1:model     !0
##ncondfile          !1
##polygonhead        !$ID
##condfile           !./cond_G3.msh
## case for model**.dat type
0:cond,1:model     !1
model connect      !result_inv/model_connect.dat
nmodelfile         !1
## polygonhead: header of output file
polygonhead        !$ID
modelfile          !$modelfile
## ibound 0: no bound, 1:lower, 2: upper bound
## If you want not to show some resistivity values larger / smaller than some threshold
## you can set by ibound =1 for lower bound, and 2 for upper bound of resistivity
## no bound case
ibound             !0
## lower bound case
##iboud              !1
##lower bo.log10[O.m]!1.0
## upper bound case
##iboud              !2
##upper bo.log10[O.m]!2.8
## nslice: number of slices you want to generate
nslice             !3
##---- Details for the first slice ------
## First slice file will be ${ID}_01.dat
## slice y = y_slice (WE vertical slice)
## 2026.11.07 CC D and boundaries are changed for Kikai, sitefile line is removed
A                  !0.0
B                  !1.0
D                  !$(echo "$y_slice" | awk '{printf "%.3f",-$1}')
boundary left      !$xmin
boundary right     !$xmax
boundary bottom    !$zmin
boundary top       !$zmax
##--- Details for the second slice -----
## second slice file will be ${ID}_02.dat
## slice x = x_slice (SN vertical slice)
## 2026.11.07 CC D and boundaries are changed for Kikai, sitefile line is removed
A                  !1.0
B                  !0.0
D                  !$(echo "$x_slice" | awk '{printf "%.3f",-$1}')
boundary left      !$ymin
boundary right     !$ymax
boundary bottom    !$zmin
boundary top       !$zmax
##--- Details for the third slice -----
## horizontal slice z = z_elevation
## 2026.11.07 CC D and boundaries are changed for Kikai, sitefile line is removed
A                  !0.0
B                  !0.0
D                  !$(echo "$z_elevation" | awk '{printf "%.3f",-$1}')
boundary left      !$xmin
boundary right     !$xmax
boundary bottom    !$ymin
boundary top       !$ymax
##------------------ control file for slice.exe ends here !! ----------------------
EOF

#[2]## generate tmp.dat
if [ $igendat -eq 1 ];then
SRC="../../src/src_post"
cd $SRC
make clean -f Makefile_gfort
make -f Makefile_gfort # for the use of gfortran
cd -

${SRC}/slice.exe <<EOF1
tmp0.ctl
EOF1
rm tmp0.ctl
fi

#[3]## Draw picture
SL1="WE"
SL2="SN"

CPT="slice.cpt"
range=$xmin/$xmax/$zmin/$zmax # 2026.11.07 CC
range_sn=$ymin/$ymax/$zmin/$zmax # 2026.11.07 CC
#bb=a0.25/a0.25:"Distance\(km\)":WeSn
bb=a0.25/a0.25wesn
scl=8/5
range2=0/8/0/5

xshift=2.5   ; yshift=2.5 # x=2.5 is default

gmt begin $ID pdf
gmt makecpt -Crainbow -T1.0/2.5/0.01 -I -Dred,purple

scl2=20/20
range3=0/20/0/20

yy=`echo " scale=1; 2.5 + 5.3" | bc`
gmt basemap -JX$scl -R$range -Bxa30 -Bya20+l"Elevation[km]" -BWeSn  -X${xshift} -Y${yy} # 2026.11.07 CC
if [ -e "./${ID}_01.dat" ];then
gmt plot "${ID}_01.dat" -C -L 
fi
gmt text -JX$scl -R$range2  -F+f12,Helvetica+jLM -G255 <<EOF 
0.9   4.6 $ID
0.15  4.6 W
7.00  4.6 E
EOF
gmt basemap -JX$scl -R$range_sn -Bxa30+l"Distance[km]" -Bya20+l"Elevation[km]" -BWeSn -Y-6 # 2026.11.07 CC
if [ -e  "${ID}_02.dat" ];then
 gmt plot "${ID}_02.dat" -C -L 
fi
gmt text -JX$scl -R$range2  -F+f12,Helvetica+jLM -G255 <<EOF 
0.9  4.6   $ID
0.15  4.6  S
7.0  4.6  N
EOF

# horizontal slice
gmt basemap -JX8/8 -R$xmin/$xmax/$ymin/$ymax -Bxa30+l"Distance[km]" -Bya30+l"Distance[km]" -BWeSn -X10 # 2026.11.07 CC
if [ -e  "${ID}_03.dat" ];then
 gmt plot "${ID}_03.dat" -C -L 
fi
# site locations
fldr=result_inv/
awk '{print($1,$2)}' ${fldr}site_act.dat | gmt plot -Sc0.2 -W1,black -Gred #active sites
awk '{print($1,$2)}' ${fldr}S1_site.dat  | gmt plot -W2p,red # source S1
if [ -e "${fldr}S2_site.dat" ];then # 2026.11.07 CC S2 does not exist for Kikai
awk '{print($1,$2)}' ${fldr}S2_site.dat  | gmt plot -W2p,red # source S2
fi # 2026.11.07 CC
awk '{print($1,$2)}' ${fldr}site_mt.dat  | gmt plot -Sc0.2 -W1,black -Ggreen # MT site

# ID and W E
gmt text -JX8/10 -R0/8/0/10  -F+f12,Helvetica+jLM -G255 <<EOF 
0.9   7.5  $ID
0.15  4.0  W
7.5   4.0  E
0.9   8.5  z = $z_elevation km
EOF

gmt colorbar -Dx9.1/0+w9/0.3 -B0.2  -C 

gmt text -R0/20/0/20 -JX20/20 -F+f13p,Helvetica+jCM << EOF 
9.5 11.0 log@-10@- @%30%@~r@~@%%
9.5 10.3  (@~W@~m)
EOF

gmt end
open ${ID}.pdf &
rm -f tmp.ctl ${ID}_0[1-3].dat # 2026.11.07 CC
