#!/bin/bash
# Generate fxyz_mesh (background field, IGRF) and vxyz_mesh_<tidemodel> (tidal velocity) at the ocean-mesh nodes
# Output of the tidal velocity: vxyz_mesh_HAMTIDE or vxyz_mesh_OTPS (set Vxyzfile in tide_fwd.ctl)  2026.11.27 CC
#  tidemodel=HAMTIDE : HAMTIDE11a (free; ../../mkfvxyz/HAMTIDE/get_hamtide.sh downloads it)  2026.11.27 CC
#  tidemodel=OTPS    : TPXO9 atlas via OTPS extract_HC (data are provided by the authors)     2026.11.27 CC
# In both cases, ocean nodes that fall on land or outside the tidal model grid are filled    2026.11.27 CC
# with the velocity of the nearest valid node (fill_nearest_uv.py).                           2026.11.27 CC
set -e                                    # 2026.11.27 CC
tidemodel=${1:-HAMTIDE}                   # HAMTIDE (default) or OTPS; e.g. ./mkfvxyz_mesh.sh OTPS 2026.11.27 CC
const=m2                                  # tidal constituent 2026.11.27 CC

src="../../mkfvxyz/src"  # 2024.03.07 Tatkuo Minami
OTPS="../../mkfvxyz/OTPS"
HAMTIDE="../../mkfvxyz/HAMTIDE"           # 2026.11.27 CC
igrf="../../mkfvxyz/igrf"

cd $src  # 2024.03.07 Tatkuo Minami
make clean
make
cd -
if [ "$tidemodel" = "OTPS" ]; then        # compile OTPS only when it is used 2026.11.27 CC
cd $OTPS # 2024.03.07 Tatkuo Minami
make clean
make
cd -
fi                                        # 2026.11.27 CC
# ok up to here 2026.07.22

fldr=../mesh # 2024.03.07 Tatkuo Minami
meshctl=${fldr}/mesh.ctl            # 2024.03.07 Tatkuo Minami
oceanmesh=${fldr}/ocean.msh               # 2024.03.07 Tatkuo Minami
em3dmesh=${fldr}/em3d.msh                 # 2024.03.07 Tatkuo Minami

#[1]## generate lat_lon_mesh
${src}/mkmshlonlat.exe <<EOF
"${meshctl}"
"${oceanmesh}"
"${em3dmesh}"
EOF
# ok up to here 2026.07.22

#[2]## generate fxyz_mesh
${src}/meshfxyz.exe # input is lat_lon_mesh

#[3]## generate u_xyz_mesh and v_xyz_mesh (lat, lon, amplitude [cm/s], Greenwich phase [deg])
if [ "$tidemodel" = "OTPS" ]; then        # 2026.11.27 CC
#[3-1]## generate tide model list
cat <<EOF > Model_tpxo9_atlas_v3
${OTPS}/DATA/TPXO9_atlas_v3_test/h_*_tpxo9_atlas_30_v3
${OTPS}/DATA/TPXO9_atlas_v3_test/u_*_tpxo9_atlas_30_v3
${OTPS}/DATA/TPXO9_atlas_v3_test/grid_tpxo9_atlas_30_v3
EOF
# ok up to here 2026.07.22

#[3-2]## ge
${OTPS}/extract_HC << finish
Model_tpxo9_atlas_v3
lat_lon_mesh
u
${const}
AP
oce
1
u_xyz_mesh
finish

${OTPS}/extract_HC << finish
Model_tpxo9_atlas_v3
lat_lon_mesh
v
${const}
AP
oce
1
v_xyz_mesh
finish
# ok up to here 2026.07.22

elif [ "$tidemodel" = "HAMTIDE" ]; then   # 2026.11.27 CC
  nc=${HAMTIDE}/HAMcurrent11a_${const}.nc                                        # 2026.11.27 CC
  [ -s $nc ] || { echo "$nc not found. Run ${HAMTIDE}/get_hamtide.sh ${const}"; exit 1; } # 2026.11.27 CC
  python3 ${src}/hamtide2uv.py lat_lon_mesh $nc ${const} u_xyz_mesh v_xyz_mesh   # 2026.11.27 CC
else                                                                             # 2026.11.27 CC
  echo "GEGEGE unknown tidemodel=$tidemodel (HAMTIDE or OTPS)"; exit 1           # 2026.11.27 CC
fi                                                                               # 2026.11.27 CC

#[3-3]## fill nodes on land / out of the model grid with the nearest valid node 2026.11.27 CC
python3 ${src}/fill_nearest_uv.py u_xyz_mesh v_xyz_mesh                          # 2026.11.27 CC

${src}/mkvxyzfile.exe # input is u_vxyz_mesh and v_vxyz_mesh, generate vxyz_mesh
mv vxyz_mesh vxyz_mesh_${tidemodel}   # the file name shows the tidal model 2026.11.27 CC
echo "output: vxyz_mesh_${tidemodel} (tidal model: ${tidemodel}, constituent: ${const})" # 2026.11.27 CC

# ok up to here 2026.07.22
exit

gfortran m_readfvxyz.f90
./a.out



