#!/bin/bash
set -x -e
targetres=1000
bagrid=bamber${targetres}.nc
create_greenland_bamber_grid.py -g $targetres $bagrid
nc2cdo.py $bagrid
N=4
method=con
rmweights=remapweigths.nc
sgrid=sgrid.nc

inprefix=ex_gris_g3600m_relax_v2_ppq_0.6_tefo_0.02_calving_ocean_kill_forcing_type_ctrl_hydro_diffuse_100a
# extract thickness to create a source grid
ncks --64 -O -v thk ${inprefix}_ctrl.nc $sgrid
nc2cdo.py $sgrid
cdo -P $N gen$method,$bagrid $sgrid $rmweights

for exp in "ctrl" "asmb"; do
    infile=${inprefix}_${exp}.nc
    tmpfile1=tmp_$infile
    ncks --64 -O -v mapping,thk,topg,usurf,climatic_mass_balance,bmelt,uvelsurf,vvelsurf,discharge_flux_cumulative,tempsurf,taub_mag,mask $infile $tmpfile1
    ncap2 -O -s '*sz_idt=time.size(); *discharge_flux[$time,$y,$x]= 0.f; *dlithkdt[$time,$y,$x]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {discharge_flux(idt,:,:)=-(discharge_flux_cumulative(idt,:,:)-discharge_flux_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*1e12; dlithkdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1));} discharge_flux.ram_write(); dlithkdt.ram_write(); uvelsurf=uvelsurf/31556925.9747; vvelsurf=vvelsurf/31556925.9747; climatic_mass_balance=climatic_mass_balance/31556925.9747; sftgif[$time,$y,$x]=0b; where(mask==2 || mask==3) sftgif=1; sftflf[$time,$y,$x]=0b; where(mask==3) sftflf=1; sfrgrf[$time,$y,$x]=0b; where(mask==2) sfrgrf=1;' $tmpfile1 $tmpfile1
    ncatted -a units,discharge_flux,o,c,"kg s-1" -a units,dlithkdt,o,c,"m s-1" -a units,uvelsurf,o,c,"m s-1" -a units,vvelsurf,o,c,"m s-1" -a units,climatic_mass_balance,o,c,"kg m-2 s-1" $tmpfile1

    ncks -O -x -v mask,discharge_flux_cumulative $tmpfile1 $tmpfile1
    nc2cdo.py $tmpfile1
    tmpfile2=tmp2_$infile
    cdo timselavg,5 $tmpfile1 $tmpfile2
    ncrename -v thk,lithk -v usurf,orog -v climatic_mass_balance,acabf -v bmelt,libmassbf  -v tempsurf,litempsnic -v taub_mag,strbasemag -v discharge_flux,licalvf $tmpfile2
    ncatted -a standard_name,acabf,o,c,"land_ice_surface_specific_mass_balance_flux"  -a standard_name,dlithkdt,o,c,"tendency_of_land_ice_thickness" -a standard_name,libmassbf,o,c,"land_ice_basal_specific_mass_balance_flux" -a standard_name,licalvf,o,c,"land_ice_specific_mass_flux_due_to_calving" -a standard_name,sftgif,o,c,"land_ice_area_fraction" -a standard_name,sfrgrf,o,c,"grounded_ice_sheet_area_fraction" -a standard_name,sftflf,o,c,"floating_ice_sheet_area_fraction" -a standard_name,strbasemag,o,c,"magnitude_of_land_ice_basal_drag" -a standard_name,litempsnic,o,c,"temperature_at_ground_level_in_snow_or_firn" $tmpfile2
    expu=$(echo $exp | tr '[:lower:]' '[:upper:]')
    obase=GIS_UAF_PISM_${expu}
    outfile=${obase}.nc
    cdo -P $N remap,$bagrid,$rmweights $tmpfile2 $outfile
    ncks -A -v x,y,mapping $bagrid $outfile
    cdo splitname,swap $outfile _$obase
    for file in *_$obase.nc; do
        ncks -A -v x,y,mapping $bagrid $file
        ncks -A -v run_stats,pism_config $tmpfile1 $file
    done
done
