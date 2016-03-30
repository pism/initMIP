#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden

import os
import time
import numpy as np
from pyproj import Proj
import cf_units
import sh
try:
    import subprocess32 as sub
except:
    import subprocess as sub
    
from argparse import ArgumentParser
try:
    from pypismtools import unit_converter
except:
    from pypismtools.pypismtools import unit_converter

from netCDF4 import Dataset as CDF
# Set up the option parser
parser = ArgumentParser()
parser.description = "Script to make ISMIP6-conforming 2D time series."
#parser.add_argument("INIT_FILE", nargs=1)
parser.add_argument("EXP_FILE", nargs=1)
parser.add_argument("-n", '--n_procs', dest="n_procs", type=int,
                    help='''number of cores/processors. default=4.''', default=4)

parser.add_argument("-t", "--target_resolution", dest="target_resolution", type=int,
                    choices=[1000, 5000],
                    help="Horizontal grid resolution", default=1000)

parser.add_argument("-e", "--experiment", dest="experiment",
                    choices=['ctrl', 'asmb'],
                    help="Output size type", default='ctrl')

options = parser.parse_args()
experiment = options.experiment
# initfile = options.INIT_FILE[0]
infile = options.EXP_FILE[0]
n_procs = options.n_procs
target_resolution = options.target_resolution
target_grid_filename = 'searise_grid_{}m.nc'.format(target_resolution)
IS = 'GIS'
GROUP = 'UAF'
MODEL = 'PISM'
EXP = experiment
project = '{IS}_{GROUP}_{MODEL}'.format(IS=IS, GROUP=GROUP, MODEL=MODEL)

offset = 6
ismip6_reporting_inverval = 5
ismip6_request_vars = (
    'lithk',
    'orog',
    'topg',
    'hfgeoubed',
    'acabf',
    'libmassbf',
    'dlithkdt',
    'uvelsurf',
    'vvelsurf',
    'wvelsurf',
    'uvelbase',
    'vvelbase',
    'wvelbase',
    'uvelmean',
    'vvelmean',
    'litempsnic',
    'litempbot',
    'strbasemag',
    'licalvf'
    )

pism_copy_vars = [
    'beta',
    'climatic_mass_balance',
    'discharge_flux',
#    'grounded_basal_flux',
    'taub_mag',
    'thk',
    'usurf',
    'topg',
    'hfgeoubed',
    'uvelsurf',
    'vvelsurf',
    'wvelsurf',
    'uvelbase',
    'vvelbase',
    'wvelbase',
    'sftflf',
    'sftgif',
    'sftgrf',
    'mapping',
    'pism_config',
    'run_stats'
]

pism_diag_vars = [
    'taub_mag',
    'thk',
    'usurf',
    'topg',
    'hfgeoubed',
    'uvelsurf',
    'vvelsurf',
    'wvelsurf',
    'uvelbase',
    'vvelbase',
    'wvelbase',
    'sftflf',
    'sftgif',
    'sftgrf',
    'mapping',
    'pism_config',
    'run_stats'
]

pism_flux_vars = [
    'climatic_mass_balance',
    'discharge_flux',
#    'grounded_basal_flux'
    ]


class ISMIP6Var(object):
    pism_name = None
    units = None
    standard_name = None
    def __init__(self, ismip6_name, pism_name, units, standard_name):
        self.ismip6_name = ismip6_name
        self.pism_name = pism_name
        self.units = units
        self.standard_name

    def __repr__(self):
        return "ISMIP6 Variable"
    
    
ismip6_vars = [
    ISMIP6Var("acabf"         , "climatic_mass_balance" , "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux"),
    ISMIP6Var("dlithkdt"      , "dHdt"                  , "m s-1", "tendency_of_land_ice_thickness"),
    ISMIP6Var("hfgeoubed"     , "hfgeoubed"             , "W m-2", "upward_geothermal_heat_flux_at_ground_level"),
    ISMIP6Var("iareaf"        , "iareaf"                , "m2", "floating_ice_shelf_area"),
    ISMIP6Var("iareag"        , "iareag"                , "m2", "grounded_land_ice_area"),
    ISMIP6Var("libmassbf"     , ""                      , "kg m-2 s-1", "land_ice_basal_specific_mass_balance_flux"),
    ISMIP6Var("licalvf"       , "discharge_flux"        , "kg m-2 s-1", "land_ice_specific_mass_flux_due_to_calving"),
    ISMIP6Var("lim"           , "imass"                 , "kg", "land_ice_mass"),
    ISMIP6Var("limnsw"        , "limnsw"                , "kg", "land_ice_mass_not_displacing_sea_water"),
    ISMIP6Var("litempbot"     , "tempbase"              , "K", "land_ice_basal_temperature"),
    ISMIP6Var("litempsnic"    , "tempsurf"              , "K", "temperature_at_ground_level_in_snow_or_firn"),
    ISMIP6Var("lithk"         , "thk"                   , "m", "land_ice_thickness"),
    ISMIP6Var("orog"          , "usurf"                 , "m", "surface_altitude"),
    ISMIP6Var("sfrgrf"        , "sfrgrf"                , "1", "grounded_ice_sheet_area_fraction"),
    ISMIP6Var("sftflf"        , "sftflf"                , "1", "floating_ice_sheet_area_fraction"),
    ISMIP6Var("sftgif"        , "sftgif"                , "1", "land_ice_area_fraction"),
    ISMIP6Var("strbasemag"    , "taub"                  , "Pa", "magnitude_of_land_ice_basal_drag"),
    ISMIP6Var("tendacabf"     , ""                      , "kg s-1", "tendency_of_land_ice_mass_due_to_surface_mass_balance"),
    ISMIP6Var("tendlibmassbf" , ""                      , "kg s-1", "tendency_of_land_ice_mass_due_to_basal_mass_balance"),
    ISMIP6Var("tendlicalvf"   , ""                      , "kg s-1", "tendency_of_land_ice_mass_due_to_calving"),
    ISMIP6Var("topg"          , "topg"                  , "m", "bedrock_altitude"),
    ISMIP6Var("uvelbase"      , "uvelbase"              , "m s-1", "land_ice_basal_x_velocity"),
    ISMIP6Var("uvelmean"      , "uvelbar"               , "m s-1", "land_ice_vertical_mean_x_velocity"),
    ISMIP6Var("uvelsurf"      , "uvelsurf"              , "m s-1", "land_ice_surface_x_velocity"),
    ISMIP6Var("vvelbase"      , "vvelbase"              , "m s-1", "land_ice_basal_y_velocity"),
    ISMIP6Var("vvelmean"      , "vvelbar"               , "m s-1", "land_ice_vertical_mean_y_velocity"),
    ISMIP6Var("vvelsurf"      , "vvelsurf"              , "m s-1", "land_ice_surface_y_velocity"),
    ISMIP6Var("wvelbase"      , "wvelbase"              , "m s-1", "land_ice_basal_upward_velocity"),
    ISMIP6Var("wvelsurf"      , "wvelsurf"              , "m s-1", "land_ice_surface_upward_velocity")
]

ismip6_to_pism_dict = dict((v.ismip6_name, v.pism_name) for v in ismip6_vars)
pism_to_ismip6_dict = dict((v.pism_name, v.ismip6_name) for v in ismip6_vars)


def create_searise_grid(filename, grid_spacing, **kwargs):
    '''
    Create dummy grid description
    '''

    if 'fileformat' not in kwargs.keys():
        fileformat = 'NETCDF4'
    else:
        fileformat = str.upper(kwargs['fileformat'])
        
    
    xdim = 'x'
    ydim = 'y'

    # define output grid, these are the extents of the Bamber domain
    e0 = -800000.0
    n0 = -3400000.0
    e1 = 700000.0
    n1 = -600000.0

    # Shift to cell centers
    e0 += grid_spacing / 2
    n0 += grid_spacing / 2
    e1 -= grid_spacing / 2
    n1 -= grid_spacing / 2

    de = dn = grid_spacing  # m
    M = int((e1 - e0) / de) + 1
    N = int((n1 - n0) / dn) + 1

    easting = np.linspace(e0, e1, M)
    northing = np.linspace(n0, n1, N)
    ee, nn = np.meshgrid(easting, northing)

    # Set up SeaRISE Projection
    projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=-39 +lat_0=90 +lat_ts=71 +units=m"
    proj = Proj(projection)

    lon, lat = proj(ee, nn, inverse=True)

    # number of grid corners
    grid_corners = 4
    # grid corner dimension name
    grid_corner_dim_name = "nv4"

    # array holding x-component of grid corners
    gc_easting = np.zeros((M, grid_corners))
    # array holding y-component of grid corners
    gc_northing = np.zeros((N, grid_corners))
    # array holding the offsets from the cell centers
    # in x-direction (counter-clockwise)
    de_vec = np.array([-de / 2, de / 2, de / 2, -de / 2])
    # array holding the offsets from the cell centers
    # in y-direction (counter-clockwise)
    dn_vec = np.array([-dn / 2, -dn / 2, dn / 2, dn / 2])
    # array holding lat-component of grid corners
    gc_lat = np.zeros((N, M, grid_corners))
    # array holding lon-component of grid corners
    gc_lon = np.zeros((N, M, grid_corners))
    
    for corner in range(0, grid_corners):
        ## grid_corners in x-direction
        gc_easting[:, corner] = easting + de_vec[corner]
        # grid corners in y-direction
        gc_northing[:, corner] = northing + dn_vec[corner]
        # meshgrid of grid corners in x-y space
        gc_ee, gc_nn = np.meshgrid(
            gc_easting[:, corner], gc_northing[:, corner])
        # project grid corners from x-y to lat-lon space
        gc_lon[:, :, corner], gc_lat[:, :, corner] = proj(
            gc_ee, gc_nn, inverse=True)


    nc = CDF(filename, 'w', format=fileformat)

    nc.createDimension(xdim, size=easting.shape[0])
    nc.createDimension(ydim, size=northing.shape[0])
    
    var = xdim
    var_out = nc.createVariable(var, 'f', dimensions=(xdim))
    var_out.axis = xdim
    var_out.long_name = "X-coordinate in Cartesian system"
    var_out.standard_name = "projection_x_coordinate"
    var_out.units = "meters"
    var_out[:] = easting

    var = ydim
    var_out = nc.createVariable(var, 'f', dimensions=(ydim))
    var_out.axis = ydim
    var_out.long_name = "Y-coordinate in Cartesian system"
    var_out.standard_name = "projection_y_coordinate"
    var_out.units = "meters"
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
    var_out.units = "degrees_east"
    var_out.valid_range = -180., 180.
    var_out.standard_name = "longitude"
    var_out.bounds = "lon_bnds"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
    var_out.units = "degrees_north"
    var_out.valid_range = -90., 90.
    var_out.standard_name = "latitude"
    var_out.bounds = "lat_bnds"
    var_out[:] = lat

    nc.createDimension(grid_corner_dim_name, size=grid_corners)

    var = 'lon_bnds'
    # Create variable 'lon_bnds'
    var_out = nc.createVariable(
        var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lon_bnds'
    var_out.units = "degreesE"
    # Assign values to variable 'lon_nds'
    var_out[:] = gc_lon
        
    var = 'lat_bnds'
    # Create variable 'lat_bnds'
    var_out = nc.createVariable(
        var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lat_bnds'
    var_out.units = "degreesN"
    # Assign values to variable 'lat_bnds'
    var_out[:] = gc_lat

    var = 'dummy'
    var_out = nc.createVariable(
        var,
        'f',
        dimensions=(
            "y",
            "x"),
        fill_value=-2e9)
    var_out.units = "meters"
    var_out.long_name = "Just A Dummy"
    var_out.comment = "This is just a dummy variable for CDO."
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = 0.

    mapping = nc.createVariable("mapping", 'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = 90.
    mapping.standard_parallel = 71.
    mapping.straight_vertical_longitude_from_pole = -39.
    
    from time import asctime
    historystr = 'Created ' + asctime() + '\n'
    nc.history = historystr
    nc.proj4 = projection
    nc.Conventions = 'CF-1.5'
    nc.close()

    
if __name__ == "__main__":


    project_dir = project
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)
    
    
    tmp_filename = 'tmp_{}.nc'.format(EXP)
    tmp_file = os.path.join(project_dir, tmp_filename)
    try:
        os.remove(tmp_file)
    except OSError:
        pass
    cmd = ['ncks', '-O',
           '-v', '{}'.format(','.join(pism_copy_vars)),
           infile, tmp_file]
    sub.call(cmd)

    # source_grid_filename = 'source_grid.nc'
    # source_grid_file = os.path.join(project_dir, source_grid_filename)
    # ncks_cmd = ['ncks', '-O', '-v', 'thk,mapping', infile, source_grid_file]
    # sub.call(ncks_cmd)
    # nc2cdo_cmd = ['nc2cdo.py', source_grid_file]
    # sub.call(nc2cdo_cmd)
    
    # # If exist, remove target grid description file
    # target_grid_file = os.path.join(project_dir, target_grid_filename)
    # try:
    #     os.remove(target_grid_file)
    # except OSError:
    #     pass

    # # Create target grid description file
    # create_searise_grid(target_grid_file, target_resolution)
    # tmp_filename = 'tmp_{}.nc'.format(EXP)
    # tmp_file = os.path.join(project_dir, tmp_filename)
    # try:
    #     os.remove(tmp_file)
    # except OSError:
    #     pass

    # ncks_cmd = ['ncks', '-O', '-v',
    #             ','.join(map(str, pism_copy_vars)),
    #             infile,
    #             tmp_file]
    # sub.call(ncks_cmd)

    # # Open file
    # nc = CDF(tmp_file, 'a')


    
    # cdo_weights_filename = 'searise_grid_{}m_weights.nc'.format(target_resolution)
    # cdo_weights_file = os.path.join(project_dir, cdo_weights_filename)
    # if not os.path.isfile(cdo_weights_filename):
    #     cdo_cmd = ['cdo', '-P', '{}'.format(n_procs),
    #         'gencon,{}'.format(target_grid_file),
    #         source_grid_file,
    #         cdo_weights_file]
    #     sub.call(cdo_cmd)
    # cp_cmd = ['cp', infile, tmp_file]
    # sub.call(cp_cmd)
    # nc = CDF(tmp_file, 'a')
    # # for m_var in ['velsurf_mag']:
    # #     print m_var
    # #     nc.renameVariable(m_var, 'velsurfmag')
    # #     nc.sync()
    # #     nc_var = nc.variables['velsurfmag']
    # #     i_unit = nc_var.units
    # #     i_f = cf_units.Unit(i_unit)
    # #     o_unit = 'm s-1'
    # #     o_f = cf_units.Unit(o_unit)
    # #     nc_var[:] = i_f.convert(nc_var[:], o_f)
    # #     #nc_var[:] = unit_converter(nc_var[:], i_unit, o_unit) 
    # #     nc_var.units = o_unit
    # # nc.close()
    # outfilename = '{}.nc'.format(project)
    # outfile = os.path.join(project_dir, outfilename)
    # try:
    #     os.remove(outfile)
    # except OSError:
    #     pass
    
    
# targetres=1000
# bagrid=bamber${targetres}.nc
# create_greenland_bamber_grid.py -g $targetres $bagrid
# nc2cdo.py $bagrid
# N=4
# method=con
# rmweights=remapweigths.nc
# sgrid=sgrid.nc

# inprefix=ex_gris_g3600m_relax_v2_ppq_0.6_tefo_0.02_calving_ocean_kill_forcing_type_ctrl_hydro_diffuse_100a
# # extract thickness to create a source grid
# ncks --64 -O -v thk ${inprefix}_ctrl.nc $sgrid
# nc2cdo.py $sgrid
# cdo -P $N gen$method,$bagrid $sgrid $rmweights

# for exp in "ctrl"; do
#     infile=${inprefix}_${exp}.nc
#     tmpfile1=tmp_$infile
#     ncks --64 -O -v mapping,thk,topg,usurf,climatic_mass_balance,bmelt,uvelsurf,vvelsurf,discharge_flux_cumulative,tempsurf,taub_mag,mask $infile $tmpfile1
#     ncap2 -O -s '*sz_idt=time.size(); *discharge_flux[$time,$y,$x]= 0.f; *dlithkdt[$time,$y,$x]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {discharge_flux(idt,:,:)=-(discharge_flux_cumulative(idt,:,:)-discharge_flux_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*1e12; dlithkdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1));} discharge_flux.ram_write(); dlithkdt.ram_write(); uvelsurf=uvelsurf/31556925.9747; vvelsurf=vvelsurf/31556925.9747; climatic_mass_balance=climatic_mass_balance/31556925.9747; sftgif[$time,$y,$x]=0b; where(mask==2 || mask==3) sftgif=1; sftflf[$time,$y,$x]=0b; where(mask==3) sftflf=1; sfrgrf[$time,$y,$x]=0b; where(mask==2) sfrgrf=1;' $tmpfile1 $tmpfile1
#     ncatted -a units,discharge_flux,o,c,"kg s-1" -a units,dlithkdt,o,c,"m s-1" -a units,uvelsurf,o,c,"m s-1" -a units,vvelsurf,o,c,"m s-1" -a units,climatic_mass_balance,o,c,"kg m-2 s-1" $tmpfile1

#     ncks -O -x -v mask,discharge_flux_cumulative $tmpfile1 $tmpfile1
#     nc2cdo.py $tmpfile1
#     tmpfile2=tmp2_$infile
#     cdo timselavg,5 $tmpfile1 $tmpfile2
#     ncrename -v thk,lithk -v usurf,orog -v climatic_mass_balance,acabf -v bmelt,libmassbf  -v tempsurf,litempsnic -v taub_mag,strbasemag -v discharge_flux,licalvf $tmpfile2
#     ncatted -a standard_name,acabf,o,c,"land_ice_surface_specific_mass_balance_flux"  -a standard_name,dlithkdt,o,c,"tendency_of_land_ice_thickness" -a standard_name,libmassbf,o,c,"land_ice_basal_specific_mass_balance_flux" -a standard_name,licalvf,o,c,"land_ice_specific_mass_flux_due_to_calving" -a standard_name,sftgif,o,c,"land_ice_area_fraction" -a standard_name,sfrgrf,o,c,"grounded_ice_sheet_area_fraction" -a standard_name,sftflf,o,c,"floating_ice_sheet_area_fraction" -a standard_name,strbasemag,o,c,"magnitude_of_land_ice_basal_drag" -a standard_name,litempsnic,o,c,"temperature_at_ground_level_in_snow_or_firn" $tmpfile2
#     expu=$(echo $exp | tr '[:lower:]' '[:upper:]')
#     obase=GIS_UAF_PISM_${expu}
#     outfile=${obase}.nc
#     cdo -P $N remap,$bagrid,$rmweights $tmpfile2 $outfile
#     ncks -A -v x,y,mapping $bagrid $outfile
#     cdo splitname,swap $outfile _$obase
#     for file in *_$obase.nc; do
#         ncks -A -v x,y,mapping $bagrid $file
#         ncks -A -v run_stats,pism_config $tmpfile1 $file
#     done
# done
