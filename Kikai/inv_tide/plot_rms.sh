#!/bin/bash
# plot rms for tide-only inversion (rms_mt.dat is not generated)

inact=result_inv/rms.dat
nite=`wc $inact | awk -F" " '{print($1 -1)}'`
echo $nite
gmt begin rms pdf

gmt basemap -JX5/5 -R0/${nite}/0/15 -Bxa1+l"Iteration" -Bya1+l"Normalized RMS misfit"
tail -$nite $inact | awk -F" " '{print($1, $2)}' | gmt plot -W1,red -l"nrms_tide"
gmt plot -W1,gray,- <<EOS
0 1
${nite} 1
EOS
gmt legend -DjTR+w2.5c -F+gwhite+p1p

gmt basemap -JX5/3l -R0/${nite}/0.07/1000 -Bxa1+l"Iteration" -Bya5+l"Hyper parameter (black)" -Y6
tail -$nite $inact | awk -F" " '{print($1, $4)}' | gmt plot -W1,black
tail -$nite $inact | awk -F" " '{print($1, $5)}' | gmt plot -W1,red -JX5/3 -R0/${nite}/0/5000 -Bya1000+l"Roughness (red)" -BE

gmt end
open rms.pdf &
