#!/bin/bash
# Coded on 11 Nov 2018
# shebang is moved to the first line 2026.11.26 CC
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1 # 2023.12.07

head="polygon"
sdir="../../tmtgem/mesh"
ctl="mesh.ctl"

OPT="-setnumber General.RandomSeed 12345" # 2026.07.23 for the same mesh generation, random seed is fixed

#===== error check functions ===== 2026.11.26 CC
# Fortran "stop" returns exit code 0, so each step is judged by
#  (1) exit code, (2) error messages in its log ("GEGEGE" for Fortran, "Error" for gmsh),
#  (3) existence of the output files (non-empty) required by the next step.
# Log of each step is saved as log_<step>.log                                 2026.11.26 CC
fail() {                                                                     # 2026.11.26 CC
  echo ""                                                                    # 2026.11.26 CC
  echo "##### tetmeshgen.sh STOPPED at step [$1] #####"                      # 2026.11.26 CC
  echo "##### $2"                                                            # 2026.11.26 CC
  [ -n "$3" ] && echo "##### see $3"                                         # 2026.11.26 CC
  exit 1                                                                     # 2026.11.26 CC
}                                                                            # 2026.11.26 CC
# run_step <step name> <log file> <stdin file or "-"> <command ...>          # 2026.11.26 CC
run_step() {                                                                 # 2026.11.26 CC
  local step=$1 log=$2 input=$3; shift 3                                     # 2026.11.26 CC
  echo "===== [${step}] $* ====="                                            # 2026.11.26 CC
  if [ "$input" = "-" ]; then                                                # 2026.11.26 CC
    "$@" 2>&1 | tee "$log"                                                   # 2026.11.26 CC
  else                                                                       # 2026.11.26 CC
    "$@" < "$input" 2>&1 | tee "$log"                                        # 2026.11.26 CC
  fi                                                                         # 2026.11.26 CC
  local status=${PIPESTATUS[0]}                                              # 2026.11.26 CC
  [ "$status" -ne 0 ] && fail "$step" "exit code = $status" "$log"           # 2026.11.26 CC
  grep -q "GEGEGE" "$log" && fail "$step" "$(grep -m1 GEGEGE "$log")" "$log" # 2026.11.26 CC
  grep -q "^Error" "$log" && fail "$step" "$(grep -m1 '^Error' "$log")" "$log" # gmsh 2026.11.26 CC
  return 0                                                                   # 2026.11.26 CC
}                                                                            # 2026.11.26 CC
# check_files <step name> <file ...> : stop if a file is missing or empty    # 2026.11.26 CC
check_files() {                                                              # 2026.11.26 CC
  local step=$1; shift                                                       # 2026.11.26 CC
  for f in "$@"; do                                                          # 2026.11.26 CC
    [ -s "$f" ] || fail "$step" "output file '$f' was not created (or is empty)" # 2026.11.26 CC
  done                                                                       # 2026.11.26 CC
}                                                                            # 2026.11.26 CC

# [0] coastline.exe
cd $sdir || fail "make" "cannot cd to $sdir"           # 2026.11.26 CC
make clean
make || fail "make" "compilation in $sdir failed"      # 2026.11.26 CC
cd - > /dev/null                                       # 2026.11.26 CC

./clean.sh # erase existing file

#${sdir}/msh2spherical.exe < ${ctl}
#exit

#[1]## coastline
run_step "1 coastline" log_1_coastline.log ${ctl} ${sdir}/coastline.exe   # 2026.11.26 CC
check_files "1 coastline" ${head}ki.geo bgmesh.pos                       # 2026.11.26 CC

#[2]## gmsh polygonki.geo
run_step "2 gmsh 2D" log_2_gmsh2d.log - gmsh ${head}ki.geo $OPT -2 -format msh2 -bgm bgmesh.pos # 2026.11.26 CC
check_files "2 gmsh 2D" ${head}ki.msh                                    # 2026.11.26 CC
rm bgmesh.pos
#gmsh ${head}ki.msh  >& /dev/null &

#[3]## mshki2ocean.f90, to extract horizontal ocean mesh
run_step "3 mshki2ocean" log_3_mshki2ocean.log ${ctl} ${sdir}/mshki2ocean.exe # 2026.11.26 CC
check_files "3 mshki2ocean" ${head}_ki.msh ki2oceanptr.dat               # 2026.11.26 CC
#gmsh ${head}_ki.msh  >& /dev/null &

#[4]## extrude.f90, to generate ocean.msh
run_step "4 extrude" log_4_extrude.log ${ctl} ${sdir}/extrude.exe        # 2026.11.26 CC
check_files "4 extrude" ocean.msh ${head}z.msh                           # 2026.11.26 CC
#gmsh ocean.msh  >& /dev/null &


#[5]## mk3dgeo.f90
run_step "5 mk3dgeo" log_5_mk3dgeo.log ${ctl} ${sdir}/mk3dgeo.exe        # 2026.11.26 CC
check_files "5 mk3dgeo" pre3d.geo bgmesh3d.pos                           # 2026.11.26 CC

#gmsh pre3d.geo $OPT -3 -format msh2
# gmsh 3D output is written only to mesh3d.log (not to the screen) as before 2026.11.26 CC
# Mesh.AngleToleranceFacetOverlap (default 0.1 deg) is decreased to 0.01 deg, because very thin ocean wedges
# (e.g. shallow shelf with coarse outer mesh when lenout is large) make sea-surface and seafloor facets
# almost coplanar and gmsh judges them as "overlapping facets"                          2026.11.26 CC
OPT3D="-setnumber Mesh.AngleToleranceFacetOverlap 0.01"                                # 2026.11.26 CC
run_step "5 gmsh 3D" mesh3d.log - gmsh pre3d.geo $OPT $OPT3D -3 -bgm bgmesh3d.pos -format msh2 > /dev/null # 2026.11.26 CC
check_files "5 gmsh 3D" pre3d.msh                                        # 2026.11.26 CC

#[6]## combine3d.f90
run_step "6 combine3d" log_6_combine3d.log ${ctl} ${sdir}/combine3d.exe  # 2026.11.26 CC
check_files "6 combine3d" em3d.msh                                       # 2026.11.26 CC

#[7]##
run_step "7 mkline" log_7_mkline.log - ${sdir}/mkline.exe                # 2026.11.26 CC
check_files "7 mkline" lineinfo.dat                                      # 2026.11.26 CC

echo ""                                                  # 2026.11.26 CC
echo "##### tetmeshgen.sh finished successfully #####"   # 2026.11.26 CC
