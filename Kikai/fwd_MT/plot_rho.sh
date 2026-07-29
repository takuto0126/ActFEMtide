#!/bin/bash

infile=result/P10_MT.dat

gmt begin rho_${ext} pdf
gmt gmtset FONT_ANNOT_PRIMARY="10p,Helvetica,black"
# 1: freq, 2:rhoxx, 3:phaxx,4rhoxy,5phaxy,6rhoyx,7phayx,8rhoyy,9phayy
gmt basemap -JX7l/7l -R0.8/110/0.0001/1000 -Bxa1f3g3+l"Frequency [Hz]" -Bya1f3g3+l"App. Resistivity [Ohm.m]" -BWeSn -Y12

awk '{print($1,$2)}' $infile | gmt plot  -Gred    -Sc0.3 -l"xx"
awk '{print($1,$4)}' $infile | gmt plot  -Gblue   -Sc0.3 -l"xy"
awk '{print($1,$6)}' $infile | gmt plot  -Ggreen  -Sc0.3 -l"yx"
awk '{print($1,$8)}' $infile | gmt plot  -Gpurple -Sc0.3 -l"yy"
gmt legend -DjBR+w1.5/2+o0.2 -F+gwhite+p1p


gmt basemap -JX7l/7 -R0.8/110/90/110 -Bxa1f3g3+l"Frequency [Hz]" -Bya2f2g2+l"Apparent Resistivity [Ohm.m]" -BWeSn -Y-10
awk '{print($1,$4)}' $infile | gmt plot  -Gblue -Sc0.3  -l"xy"
awk '{print($1,$6)}' $infile | gmt plot  -Ggreen -Sc0.3  -l"yx"
gmt legend -DjBR+w1.5/2+o0.2 -F+gwhite+p1p

# Phase
gmt basemap -JX7l/7 -R0.8/110/0/90 -Bxa1f3g3+l"Frequency [Hz]" -Bya15f15g15+l"Phase [deg]" -BWeSn -X10 -Y10 
awk '{print($1,$3)}'       $infile | gmt plot  -Gred -Sc0.3  -l"xx"
awk '{print($1,$5 + 180)}' $infile | gmt plot  -Gblue -Sc0.3  -l"xy"
awk '{print($1,$7      )}' $infile | gmt plot  -Ggreen -Sc0.3  -l"yx"
awk '{print($1,$9 + 180)}' $infile | gmt plot  -Gpurple -Sc0.3  -l"yy"
gmt legend -DjBR+w1.5/2+o0.2 -F+gwhite+p1p

# Phase for 45 deg
gmt basemap -JX7l/7 -R0.8/110/40/50 -Bxa1f3g3+l"Frequency [Hz]" -Bya1f1g1+l"Phase [deg]" -BWeSn -Y-10 
awk '{print($1,$5 + 180)}'       $infile | gmt plot  -Gblue -Sc0.3  -l"xy"
awk '{print($1,$7      )}'       $infile | gmt plot  -Ggreen -Sc0.3 -l"yx"
gmt legend -DjBR+w1.5/2+o0.2 -F+gwhite+p1p

gmt end
open rho_${ext}.pdf &
