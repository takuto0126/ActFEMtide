#!/bin/bash
# Create etopo_kikai-l.xyz from etopo-l.tiff (GeoTIFF, GDAL-readable by GMT)
# etopo-l.tiff stores values as depth-positive (land negative, sea positive),
# so the z column is negated to give the standard elevation-positive-up convention.

set -e

tiff=etopo-l.tiff
xyz=etopo_kikai-l.xyz

gmt grd2xyz $tiff | awk '{printf "%.10f %.10f %g\n", $1, $2, -$3}' > $xyz

rm -f gmt.history
