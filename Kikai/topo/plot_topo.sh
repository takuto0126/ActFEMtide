#!/bin/bash

file=etopo_kikai-l.xyz

gmt begin topo pdf
gmt makecpt -Cglobe -T-5000/5000/1 -I
gmt basemap -JM20 -R95/165/0/60 -Ba10
gmt xyz2grd $file -R95/165/0/60 -I0.2 -Ggrd.grd
gmt grdimage grd.grd -C
gmt colorbar -C -B1000 -Dx22/1+w10/0.4
gmt end show