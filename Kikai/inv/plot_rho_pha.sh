#!/bin/bash

infile=result_fwd/A01_MT.dat
outps=resp.ps
outpdf=resp.pdf
gmt gmtset FONT_ANNOT_PRIMARY="10p,Helvetica,black"
# 1: freq, 2:rhoxx, 3:phaxx,4rhoxy,5phaxy,6rhoyx,7phayx,8rhoyy,9phayy
gmt psbasemap -JX7l/7l -R0.8/110/0.1/1000 -Ba1f3g3:"Frequency [Hz]":/a1f3g3:"App. Resistivity [Ohm.m]":WeSn -P -K -Y12 > $outps

cat $infile | awk '{print($1,$2)}' | gmt psxy -J -R -Gred -Sc0.3 -K -O >> $outps

cat $infile | awk '{print($1,$4)}' | gmt psxy -J -R -Gblue -Sc0.3 -K -O >> $outps

cat $infile | awk '{print($1,$6)}' | gmt psxy -J -R -Ggreen -Sc0.3 -K -O >> $outps

cat $infile | awk '{print($1,$8)}' | gmt psxy -J -R -Gpurple -Sc0.3 -K -O >> $outps


gmt psbasemap -JX7l/7 -R0.8/110/90/110 -Ba1f3g3:"Frequency [Hz]":/a2f2g2:"Apparent Resistivity [Ohm.m]":WeSn -K -O -Y-10 >> $outps

cat $infile | awk '{print($1,$4)}' | gmt psxy -J -R -Gblue -Sc0.3 -K -O >> $outps

cat $infile | awk '{print($1,$6)}' | gmt psxy -J -R -Ggreen -Sc0.3 -K -O >> $outps



# Phase
gmt psbasemap -JX7l/7 -R0.8/110/0/90 -Ba1f3g3:"Frequency [Hz]":/a15f15g15:"Phase [deg]":WeSn -X10 -Y10 -K -O >> $outps
cat $infile | awk '{print($1,$3)}' | gmt psxy -J -R -Gred -Sc0.3 -K -O >> $outps
cat $infile | awk '{print($1,$5)}' | gmt psxy -J -R -Gblue -Sc0.3 -K -O >> $outps
cat $infile | awk '{print($1,$7 +180 )}' | gmt psxy -J -R -Ggreen -Sc0.3 -K -O >> $outps
cat $infile | awk '{print($1,$9 + 180)}' | gmt psxy -J -R -Gpurple -Sc0.3 -K -O >> $outps

# Phase for 45 deg
gmt psbasemap -JX7l/7 -R0.8/110/40/50 -Ba1f3g3:"Frequency [Hz]":/a1f1g1:"Phase [deg]":WeSn -Y-10 -K -O >> $outps
cat $infile | awk '{print($1,$5)}' | gmt psxy -J -R -Gblue -Sc0.3 -K -O >> $outps
cat $infile | awk '{print($1,$7 +180 )}' | gmt psxy -J -R -Ggreen -Sc0.3 -O >> $outps


convert -density 300 $outps $outpdf
open $outpdf &
