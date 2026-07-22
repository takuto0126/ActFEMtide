#!/bin/bash

infile=./result_gwch_xy2/bxyz_xy2D0.22E-04_S1.dat
topo=../../Kikai/mesh_light/topo.xyz

wesn=-150/150/-150/150
di=2

gmt begin bxyz pdf
gmt surface $topo -R$wesn -I1/1 -Gtopo.grd

#---------------------------   Amplitude   -----------------------
gmt makecpt -Cjet -T-1/5/0.01 

#[bz] East
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"|B| east [nT]" -Y16
awk '{print($1,$2,$4)}' $infile | gmt surface -R$wesn -I$di -Ggrd.grd -T0
gmt grdimage  grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

#[bx] North
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"|B| north [nT]" -X11
awk '{print($1,$2,$6)}' $infile | gmt surface -R$wesn -I$di -Ggrd.grd -T0
gmt grdimage  grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

#[by] UP
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"|B| up [nT]" -X11
awk '{print($1,$2,$8)}' $infile | gmt surface -R$wesn -I$di -Ggrd.grd -T0
gmt grdimage  grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

gmt colorbar -Dx11/0+w10/0.3 -B1+l"Amp nT" -G0/5

#-----------------------------   phase   ---------------------------
gmt makecpt -Ccyclic -T-180/180/0.01 

#[bz] East
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50+l"Northing[km]" -BWeSn+t"arg(B east) [nT]" -X-22 -Y-13
awk '{print($1,$2,$5)}' $infile | gmt surface -R$wesn -I$di -Ggrd.grd -T0
gmt grdimage  grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

#[bx] North
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"arg(B north) [nT]"  -X11
awk '{print($1,$2,$7)}' $infile | gmt surface -R$wesn -I$di -Ggrd.grd -T0
gmt grdimage  grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

#[by] UP
gmt basemap -R$wesn -JX10/10 -Bxa50+l"Easting[km]" -Bya50 -BWeSn+t"arg(B up) [nT]" -X11
awk '{print($1,$2,$9)}' $infile | gmt surface -R$wesn -I$di -Ggrd.grd -T0
gmt grdimage  grd.grd -C -R$wesn
gmt grdcontour topo.grd -C50 -L-40/40 -W1

gmt colorbar -Dx11/0+w10/0.3 -B90+l"Phase deg" -G-180/180
gmt end show