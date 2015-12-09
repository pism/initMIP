#!/bin/bash
set -x -e

inprefix=ts_greenland_g3600m_const_ctrl_v2_sia_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_ssa_e_1.0_phi_min_5.0_phi_max_40.0_topg_min_-700_topg_max_700_hydro_100a

for exp in "ctrl" "asmb"; do
    infile=${inprefix}_${exp}.nc
    # expu=${exp^^}
    expu=$(echo $exp | tr '[:lower:]' '[:upper:]')
    outfile=scalar_GIS_UAF_PISM_${expu}.nc
    cdo yearmean -selvar,imass,discharge_flux,grounded_basal_ice_flux,surface_ice_flux,nonneg_flux,iareag,iareaf $infile $outfile
    ncrename -v imass,lim -v discharge_flux,tendacabf -v surface_ice_flux,tendlibmassbf $outfile
    ncatted -a standard_name,lim,o,c,"land_ice_mass" -a standard_name,iareag,o,c,"grounded_land_ice_area" -a standard_name,iareaf,o,c,"floating_ice_shelf_area" -a standard_name,tendacabf,o,c,"tendency_of_land_ice_mass_due_to_surface_mass_balance" -a standard_name,tendlibmassbf,o,c,"tendency_of_land_ice_mass_due_to_basal_mass_balance" $outfile
done
