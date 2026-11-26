#!/bin/bash

it=04 # iteration number

#=========================================================== function plot_ap
plot_ap(){
it=$1    # iteration
site=$2  # site number
X=$3     # x shift
Y=$4     # y shift
echo $site
infile1=result_inv/MT${site}_MT${it}.dat
Da="cat result_inv/MT${site}_MTIN.dat " # amp   data file
Dp="cat result_inv/MT${site}_MTIN.dat " # phase data file
S1="cat $infile1 "
S="-Sc0.2"
FF="-JX6/8 -R0/6/0/8 -F+f12p,Helvetica,black+jLM -Gwhite -W1,black"

gmt basemap -JX7l/7l -R0.8/110/0.001/10000 -Bxa1pf3g1 -Bya1pf3g1+l"App. Res. [Ohm.~m]" -X$X -Y$Y
$Da | awk '{print($1,$2,$10)}' | gmt plot $S -Ey+p0.6,blue+w0.4 -W1,blue
#$Da | awk '{print($1,$2)}'     | gmt plot $S -W1,blue -Gwhite
#$Da | awk '{print($1,$2+$10)}' | gmt plot -Si0.2 -W1,blue -Gblue
#$Da | awk '{print($1,$2-$10)}' | gmt plot -St0.2 -W1,blue -Gblue
$Da | awk '{print($1,$4,$12)}' | gmt plot $S -Ey+p0.6,red+w0.4 -W1,red
#$Da | awk '{print($1,$4)}'     | gmt plot $S -W1,red -Gwhite
#$Da | awk '{print($1,$4+$12)}' | gmt plot -Si0.2  -W1,red  -Gred
#$Da | awk '{print($1,$4-$12)}' | gmt plot -St0.2  -W1,red  -Gred
$Da | awk '{print($1,$6,$14)}' | gmt plot $S -Ey+p0.6,purple+w0.4 -W1,purple
#$Da | awk '{print($1,$6)}'     | gmt plot $S -W1,purple -Gwhite
#$Da | awk '{print($1,$6+$14)}' | gmt plot -Si0.2 -W1,purple -Gpurple
#$Da | awk '{print($1,$6-$14)}' | gmt plot -St0.2 -W1,purple -Gpurple
$Da | awk '{print($1,$8,$16)}' | gmt plot $S -Ey+p0.6,green+w0.4 -W1,green
#$Da | awk '{print($1,$8)}'      | gmt plot $S -W1,green -Gwhite
#$Da | awk '{print($1,$8+$16)}' | gmt plot -Si0.2 -W1,green -Ggreen
#$Da | awk '{print($1,$8-$16)}' | gmt plot -St0.2 -W1,green -Ggreen
# Calculation
$S1 | awk '{print($1,$8)}' | gmt plot -W1,green,-  # Rho yy
$S1 | awk '{print($1,$2)}' | gmt plot -W1,blue,-   # Rho xx 
$S1 | awk '{print($1,$4)}' | gmt plot -W1,red    # Rho xy
$S1 | awk '{print($1,$6)}' | gmt plot -W1,purple # Rho yx 

echo "2.5 6.5 $site" | gmt text $FF

# phase
gmt basemap -JX7l/7 -R0.8/110/-180/180 -Bxa1pf3g1+l"Frequency [Hz]" -Bya90g90+l"Phase [deg]" -Y-8
$Dp | awk '{print($1,$3,$11)}' | gmt plot $S -Ey+p0.6,blue+w0.4 -W1,blue
#$Dp | awk '{print($1,$3)}'      | gmt plot $S -W1,blue -Gwhite
#$Dp | awk '{print($1,$3+$11)}' | gmt plot -Si0.2 -W1,blue -Gblue
#$Dp | awk '{print($1,$3-$11)}' | gmt plot -St0.2 -W1,blue -Gblue
$Dp | awk '{print($1,$5,$13)}' | gmt plot $S -Ey+p0.6,red+w0.4 -W1,red
#$Dp | awk '{print($1,$5)}'      | gmt plot $S -W1,red -Gwhite
#$Dp | awk '{print($1,$5+$13)}' | gmt plot -Si0.2 -W1,red -Gred
#$Dp | awk '{print($1,$5-$13)}' | gmt plot -St0.2 -W1,red -Gred
$Dp | awk '{print($1,$7,$15)}' | gmt plot $S -Ey+p0.6,purple+w0.4 -W1,purple
#$Dp | awk '{print($1,$7)}'      | gmt plot $S -W1,purple -Gwhite
#$Dp | awk '{print($1,$7+$15)}' | gmt plot -Si0.2 -W1,purple -Gpurple
#$Dp | awk '{print($1,$7-$15)}' | gmt plot -St0.2 -W1,purple -Gpurple
$Dp | awk '{print($1,$9,$17)}' | gmt plot $S -Ey+p0.6,green+w0.4 -W1,green
#$Dp | awk '{print($1,$9)}'      | gmt plot $S -W1,green -Gwhite
#$Dp | awk '{print($1,$9+$17)}' | gmt plot -Si0.2 -W1,green -Ggreen
#$Dp | awk '{print($1,$9-$17)}' | gmt plot -St0.2 -W1,green -Ggreen
# calculation
$S1 | awk '{print($1,$3)}'     | gmt plot -W1,blue,-   # phase xx 
$S1 | awk '{print($1,$5)}'     | gmt plot -W1,red    # phase xy
$S1 | awk '{print($1,$7)}'     | gmt plot -W1,purple # phase yx 
$S1 | awk '{print($1,$9)}'     | gmt plot -W1,green,-  # phase yy 

echo "2.5 6.5 $site" | gmt text $FF
}
#============================================================= function plot_ap end


gmt begin MT_fit_${it} pdf

#       it  site    X  Y
plot_ap $it   1    2  14
plot_ap $it   2    9  8
plot_ap $it   3    9  8
plot_ap $it   4    9  8

gmt end

open MT_fit_${it}.pdf