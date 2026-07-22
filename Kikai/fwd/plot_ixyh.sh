#!/bin/bash

infile=result_xy2/ixyh_xy2D0.22E-04.dat

topo=../../Kikai/mesh_light/topo.xyz

wesn=-150/150/-150/150
di=2


gmt begin ixyh pdf
gmt surface $topo -R$wesn -I1/1 -Gtopo.grd

#---------------vh [m2/sec] -------------
gmt makecpt -Cjet -T-0.01/0.18/0.001
# real
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"|vh| real 10@+3@+[m2/s]" -Y55
awk '{print($1,$2,sqrt($16 ^2 + $18 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$16,$18,0,0,0)}' $infile | gmt psvelo -Se5/0.0/0 -A+a30+ea+jc+n30 -W0.7
gmt grdcontour topo.grd -C50 -L-40/40 -W1

# imag
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"|vh| imag 10@+3@+[m2/s]" -X11
awk '{print($1,$2,sqrt($17 ^2 + $19 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$17,$19,0,0,0)}' $infile | gmt psvelo -Se5/0.0/0 -A+a30+ea+jc+n30 -W0.7

gmt colorbar -Dx11/0+w10/0.3 -B0.02+l"Transport (v*depth) 10@+3@+[m2/sec]"

#---------------jext = sigma (v by F) * h [A/m] -------------
gmt makecpt -Cjet -T-4/18/0.001
# real
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"J@-ext@-h=@~s@~(v*F)*h real [A/km]" -Y-13 -X-11
awk '{print($1,$2,sqrt($8 ^2 + $10 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$8,$10,0,0,0)}' $infile | gmt psvelo -Se0.05/0.0/0 -A+a30+ea+jc+n30 -W0.7
gmt grdcontour topo.grd -C50 -L-40/40 -W1

# imag
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"J@-ext@-h=@~s@~(v*F)*h imag [A/km]" -X11
awk '{print($1,$2,sqrt($9 ^2 + $11 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$9,$11,0,0,0)}' $infile | gmt psvelo -Se0.05/0.0/0 -A+a30+ea+jc+n30 -W0.7

gmt colorbar -Dx11/0+w10/0.3 -B2+l"current density * depth [A/km]"  -G0/18

#---------------jtotal = sigma (E + v by F) * h [A/m] -------------
gmt makecpt -Cjet -T-4/18/0.001
# real
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"J@-total@-h=@~s@~(E+v*F)*h real [A/km]" -Y-13 -X-11
awk '{print($1,$2,sqrt($12 ^2 + $14 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$12,$14,0,0,0)}' $infile | gmt psvelo -Se0.09/0/0 -A+a30+ea+jc+n30 -W0.7
gmt grdcontour topo.grd -C50 -L-40/40 -W1

# imag
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"J@-total@-h=@~s@~(E+v*F)*h imag [A/km]" -X11
awk '{print($1,$2,sqrt($13 ^2 + $15 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$13,$15,0,0,0)}' $infile | gmt psvelo -Se0.09/0/0 -A+a30+ea+jc+n30 -W0.7

gmt colorbar -Dx11/0+w10/0.3 -B2+l"current density * depth [A/km]" -G0/18

#---------------bh [nT] -------------
gmt makecpt -Cjet -T-1/5/0.001
# real
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"bh seafloor real [nT]" -Y-13 -X-11
awk '{print($1,$2,sqrt($20 ^2 + $22 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$20,$22,0,0,0)}' $infile | gmt psvelo -Se0.2/0/0 -A+a30+ea+jc+n30 -W0.7
gmt grdcontour topo.grd -C50 -L-40/40 -W1

# imag
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"bh seafloor imag [nT]" -X11
awk '{print($1,$2,sqrt($21 ^2 + $23 ^2))}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1
awk '(NR%200)%4 == 0 && (NR - NR%200)%800 == 0 {print($1,$2,$21,$23,0,0,0)}' $infile | gmt psvelo -Se0.2/0/0 -A+a30+ea+jc+n30 -W0.7

gmt colorbar -Dx11/0+w10/0.3 -B1+l"Bh [nT]" -G0/5

#---------------bz [nT] -------------
gmt makecpt -Cpolar -T-4/4/0.001
# real
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"b upward seafloor real [nT]" -Y-13 -X-11
awk '{print($1,$2,$24)}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

# imag
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"b upward seafloor imag [nT]" -X11
awk '{print($1,$2,$25)}' $infile | gmt surface -R$wesn -I$di -T0 -Ggrd.grd
gmt grdimage grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

gmt colorbar -Dx11/0+w10/0.3 -B1+l"Bz (upward) [nT]"

gmt end show