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
z_elevation=0.75 # for horizontal slice
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
mshfile            !../mesh_joint/nakadake3d.msh
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
## slice y = 0 (WE vertical slice)
A                  !0.0
B                  !1.0
D                  !0.0
boundary left      !-1.0
boundary right     !1.0
boundary bottom    !0.25
boundary top       !1.5
sitefile           !dummy
##--- Details for the second slice -----
## second slice file will be ${ID}_02.dat
## slice x = 0 (SN vertical slice)
A                  !1.0
B                  !0.0
D                  !0.0
boundary left      !-1.0
boundary right     !1.0
boundary bottom    !0.25
boundary top       !1.5
sitefile           !dummy
##--- Details for the third slice -----
## horizontal slie z = 1.0
A                  !0.0
B                  !0.0
D                  !-$z_elevation
boundary left      !-1.0
boundary right     !1.0
boundary bottom    !-1.0
boundary top       !1.0
sitefile           !dummy
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
range=-1.0/1.0/0.25/1.5
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
gmt basemap -JX$scl -R$range -Bxa0.25 -Bya0.25+l"Elevation" -BWeSn  -X${xshift} -Y${yy} 
if [ -e "./${ID}_01.dat" ];then
gmt plot "${ID}_01.dat" -C -L 
fi
gmt text -JX$scl -R$range2  -F+f12,Helvetica+jLM -G255 <<EOF 
0.9   4.6 $ID
0.15  4.6 W
7.00  4.6 E
EOF
gmt basemap -JX$scl -R$range -Bxa0.25+l"Distance[km]" -Bya0.25+l"Elevation[km]" -BWeSn -Y-6
if [ -e  "${ID}_02.dat" ];then
 gmt plot "${ID}_02.dat" -C -L 
fi
gmt text -JX$scl -R$range2  -F+f12,Helvetica+jLM -G255 <<EOF 
0.9  4.6   $ID
0.15  4.6  S
7.0  4.6  N
EOF

# horizontal slice
gmt basemap -JX8/8 -R-1/1/-1/1 -Bxa0.25+l"Distance[km]" -Bya0.25+l"Distance[km]" -BWeSn -X10
if [ -e  "${ID}_03.dat" ];then
 gmt plot "${ID}_03.dat" -C -L 
fi
# site locations
fldr=result_inv/
awk '{print($1,$2)}' ${fldr}site_act.dat | gmt plot -Sc0.2 -W1,black -Gred #active sites
awk '{print($1,$2)}' ${fldr}S1_site.dat  | gmt plot -W2p,red # source S1
awk '{print($1,$2)}' ${fldr}S2_site.dat  | gmt plot -W2p,red # source S2
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
rm test*.ctl *.dat
