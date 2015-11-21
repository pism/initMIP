#!/bin/bash
set -x -e
pw=$1
N=4
INGRID=5

pism_grid=4500
smbfile=dsmb_0${INGRID}B13_ISMIP6.nc
wget -nc --ftp-user=searise --ftp-password=$pw -nc ftp://cryoftp1.gsfc.nasa.gov/ISMIP6/initMIP/dSMB/$smbfile

griddes=g${GRID}km.nc
create_greenland_bamber_grid.py -g ${INGRID}000 $griddes

# FIX SMB
ncks -A -v x,y,mapping $griddes $smbfile


pism_grid_file=g${pism_grid}m.nc
create_greenland_epsg3413_grid.py -g ${pism_grid} $pism_grid_file
nc2cdo.py $pism_grid_file

pismsmbfile=smb_$pism_grid
EXTRAPOLATE=on cdo -P $N remapbil,$pism_grid_file $smbfile dsmb_$pism_grid_file
ncks -A -v x,y,mapping $pism_grid_file dsmb_$pism_grid_file
ncap2 -O -s "climatic_mass_balance=DSMB*910;" dsmb_$pism_grid_file dsmb_$pism_grid_file

ncatted -O -a standard_name,climatic_mass_balance,o,c,"land_ice_surface_specific_mass_balance" dsmb_$pism_grid_file
ncatted -O -a grid_mapping,climatic_mass_balance,o,c,"mapping" -a units,climatic_mass_balance,o,c,"kg m-2 year-1" dsmb_$pism_grid_file

mpiexec -n $N fill_missing_petsc.py -v climatic_mass_balance dsmb_$pism_grid_file tmp_dsmb_$pism_grid_file
ncks -A -v climatic_mass_balance tmp_dsmb_$pism_grid_file dsmb_$pism_grid_file

#ncatted -a _FillValue,climatic_mass_balance,d,, -a missing_value,climatic_mass_balance,d,, dsmb_$pism_grid_file
#ncap2 -O -s "where(climatic_mass_balance<=-10000) climatic_mass_balance=-10000;" dsmb_$pism_grid_file dsmb_$pism_grid_file

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc

echo "fetching master file ... "
wget -nc ${DATAURL}${DATANAME}   # -nc is "no clobber"
echo "  ... done."
echo


PISMVERSION=pism_$DATANAME
echo -n "creating bootstrapable $PISMVERSION from $DATANAME ... "
# copy the vars we want, and preserve history and global attrs
ncks -O -v mapping,lat,lon,bheatflx,topg,thk,presprcp,smb,airtemp2m $DATANAME $PISMVERSION

# convert from water equiv to ice thickness change rate; assumes ice density 910.0 kg m-3
ncap2 -O -s "precipitation=presprcp*(1000.0/910.0)" $PISMVERSION $PISMVERSION

ncatted -O -a units,precipitation,c,c,"m/year" $PISMVERSION
ncatted -O -a long_name,precipitation,c,c,"ice-equivalent mean annual precipitation rate" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
# use pism-recognized name for 2m air temp
ncrename -O -v airtemp2m,ice_surface_temp  $PISMVERSION
ncatted -O -a units,ice_surface_temp,c,c,"Celsius" $PISMVERSION
# use pism-recognized name and standard_name for surface mass balance, after
# converting from liquid water equivalent thickness per year to [kg m-2 year-1]
ncap2  -O -s "climatic_mass_balance=1000.0*smb" $PISMVERSION $PISMVERSION
# Note: The RACMO field smb has value 0 as a missing value, unfortunately,
# everywhere the ice thickness is zero. Here we replace with 100 m a-1 ablation.
# This is a *choice* of the model of surface mass balance in thk==0 areas.
ncap2 -O -s "where(thk <= 0.0){climatic_mass_balance=-10000.0;}" $PISMVERSION $PISMVERSION
ncatted -O -a standard_name,climatic_mass_balance,m,c,"land_ice_surface_specific_mass_balance" $PISMVERSION
ncatted -O -a units,climatic_mass_balance,m,c,"kg m-2 year-1" $PISMVERSION
# de-clutter by only keeping vars we want
ncks -O -v mapping,lat,lon,bheatflx,topg,thk,precipitation,ice_surface_temp,climatic_mass_balance \
  $PISMVERSION $PISMVERSION
# straighten dimension names
ncrename -O -d x1,x -d y1,y -v x1,x -v y1,y $PISMVERSION $PISMVERSION
nc2cdo.py $PISMVERSION
echo "done."

EXTRAPOLATE=on cdo -P $N remapbil,$pism_grid_file $PISMVERSION smb_Greenland_${pism_grid}m.nc
mpiexec -n $N fill_missing_petsc.py -v climatic_mass_balance,ice_surface_temp smb_Greenland_${pism_grid}m.nc tmp_smb_Greenland_${pism_grid}m.nc
ncks -A -v climatic_mass_balance,ice_surface_temp tmp_smb_Greenland_${pism_grid}m.nc smb_Greenland_${pism_grid}m.nc

python create_anomalies.py -a dsmb_$pism_grid_file -b smb_Greenland_${pism_grid}m.nc initMIP_climate_forcing_100a_dsmb.nc
ncks -A -v mapping ${pism_grid_file} initMIP_climate_forcing_100a_dsmb.nc
ncatted -a grid_mapping,climatic_mass_balance,o,c,"mapping" -a grid_mapping,ice_surface_temp,o,c,"mapping" initMIP_climate_forcing_100a_dsmb.nc
ncks -d time,0 initMIP_climate_forcing_100a_dsmb.nc initMIP_climate_forcing_100a_ctrl.nc
