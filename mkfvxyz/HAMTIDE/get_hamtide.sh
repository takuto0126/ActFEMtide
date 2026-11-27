#!/bin/bash
# Download HAMTIDE11a tidal current files (NetCDF) from ICDC, University of Hamburg   2026.11.27 CC
#   https://www.geo.uni-hamburg.de/en/service/icdc/data/ocean/hamtide.html             2026.11.27 CC
#   Data citation: Taguchi, E., Stammer, D., and Zahel, W. (2014), J. Geophys. Res.     2026.11.27 CC
#   Oceans, 119, 4573-4592, doi:10.1002/2013JC009766. "Distributed in netCDF format by  2026.11.27 CC
#   the Integrated Climate Data Center (ICDC), CEN, University of Hamburg, Germany."    2026.11.27 CC
# Usage: ./get_hamtide.sh [constituents ...]   (default: m2)                           2026.11.27 CC
#   available: m2 s2 n2 k2 2n k1 o1 p1 q1                                              2026.11.27 CC
# Output: HAMcurrent11a_<const>.nc (UAMP, VAMP [cm/s], UPHA, VPHA [deg, Greenwich])    2026.11.27 CC
set -e
url=http://icdc.cen.uni-hamburg.de/thredds/fileServer/ftpthredds/hamtide
consts=${@:-m2}
cd "$(dirname "$0")"
for c in $consts; do
  f=HAMcurrent11a_${c}.nc
  if [ -s $f ]; then echo "$f exists"; continue; fi
  echo "downloading $f ..."
  curl -f -L -sS -o $f.part $url/$f
  mv $f.part $f
done
ls -l HAMcurrent11a_*.nc
