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
parser.add_argument("-e", "--experiment", dest="experiment",
                    choices=['ctrl', 'asmb'],
                    help="Output size type", default='ctrl')
parser.add_argument("-r", "--remap_method", dest="remap_method",
                    choices=['con', 'bil'],
                    help="Remapping method", default='con')
parser.add_argument("-t", "--target_resolution", dest="target_resolution", type=int,
                    choices=[1000, 5000],
                    help="Horizontal grid resolution", default=1000)

parser.add_argument("-w", "--override_weights_file",
                    dest="override_weights_file", action="store_true",
                    help="Override weights file", default=False)

options = parser.parse_args()
experiment = options.experiment
infile = options.EXP_FILE[0]
n_procs = options.n_procs
override_weights_file = options.override_weights_file
remap_method = options.remap_method
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
    'cell_area',
    # 'surface_mass_balance_average',
    # 'basal_mass_balance_average',
    'dHdt',
    'discharge_flux',
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
    'ubar',
    'vbar',
    'sftflf',
    'sftgif',
    'sftgrf',
    'mapping',
    'pism_config',
    'run_stats'
]

pism_mass_to_vol_vars = [
    'discharge_flux'
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
    'basal_mass_balance_average',
    'discharge_flux',
    'surface_mass_balance_averge'
    ]

pism_stats_vars = ('pism_config', 'run_stats')


class ISMIP6Var(object):
    ismip6_name = None
    pism_name = None
    units = None
    standard_name = None
    def __init__(self, ismip6_name, pism_name, units, standard_name):
        self.ismip6_name = ismip6_name
        self.pism_name = pism_name
        self.units = units
        self.standard_name = standard_name

    def __repr__(self):
        return "ISMIP6 Variable"
    

# ismip6_vars = {
#     "acabf": ISMIP6Var("acabf"         , "surface_mass_balance_average" , "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux"),
#     "dlithkdt": ISMIP6Var("dlithkdt"      , "dHdt"                  , "m s-1", "tendency_of_land_ice_thickness"),
#     "hfgeoubed": ISMIP6Var("hfgeoubed"     , "hfgeoubed"             , "W m-2", "upward_geothermal_heat_flux_at_ground_level"),
#     "iareaf": ISMIP6Var("iareaf"        , "iareaf"                , "m2", "floating_ice_shelf_area"),
#     "iareag": ISMIP6Var("iareag"        , "iareag"                , "m2", "grounded_land_ice_area"),
#     "libmassbf": ISMIP6Var("libmassbf"     , "basal_mass_balance_average"                      , "kg m-2 s-1", "land_ice_basal_specific_mass_balance_flux"),
#     "licalvf": ISMIP6Var("licalvf"       , "discharge_flux"        , "kg m-2 s-1", "land_ice_specific_mass_flux_due_to_calving"),
#     "lim": ISMIP6Var("lim"           , "imass"                 , "kg", "land_ice_mass"),
#     "limnsw": ISMIP6Var("limnsw"        , "limnsw"                , "kg", "land_ice_mass_not_displacing_sea_water"),
#     "litempbot": ISMIP6Var("litempbot"     , "tempbase"              , "K", "land_ice_basal_temperature"),
#     "litempsnic": ISMIP6Var("litempsnic"    , "tempsurf"              , "K", "temperature_at_ground_level_in_snow_or_firn"),
#     "lithk": ISMIP6Var("lithk"         , "thk"                   , "m", "land_ice_thickness"),
#     "orog": ISMIP6Var("orog"          , "usurf"                 , "m", "surface_altitude"),
#     "sftgrf": ISMIP6Var("sftgrf"        , "sftgrf"                , "1", "grounded_ice_sheet_area_fraction"),
#     "sftflf": ISMIP6Var("sftflf"        , "sftflf"                , "1", "floating_ice_sheet_area_fraction"),
#     "sftgif": ISMIP6Var("sftgif"        , "sftgif"                , "1", "land_ice_area_fraction"),
#     "strbasemag": ISMIP6Var("strbasemag"    , "taub_mag"                  , "Pa", "magnitude_of_land_ice_basal_drag"),
#     "tendacabf": ISMIP6Var("tendacabf"     , ""                      , "kg s-1", "tendency_of_land_ice_mass_due_to_surface_mass_balance"),
#     "tendlibmassbf": ISMIP6Var("tendlibmassbf" , ""                      , "kg s-1", "tendency_of_land_ice_mass_due_to_basal_mass_balance"),
#     "tendlicalvf": ISMIP6Var("tendlicalvf"   , ""                      , "kg s-1", "tendency_of_land_ice_mass_due_to_calving"),
#     "topg": ISMIP6Var("topg"          , "topg"                  , "m", "bedrock_altitude"),
#     "uvelbase": ISMIP6Var("uvelbase"      , "uvelbase"              , "m s-1", "land_ice_basal_x_velocity"),
#     "uvelmean": ISMIP6Var("uvelmean"      , "ubar"               , "m s-1", "land_ice_vertical_mean_x_velocity"),
#     "uvelsurf": ISMIP6Var("uvelsurf"      , "uvelsurf"              , "m s-1", "land_ice_surface_x_velocity"),
#     "vvelbase": ISMIP6Var("vvelbase"      , "vvelbase"              , "m s-1", "land_ice_basal_y_velocity"),
#     "vvelmean": ISMIP6Var("vvelmean"      , "vbar"               , "m s-1", "land_ice_vertical_mean_y_velocity"),
#     "vvelsurf": ISMIP6Var("vvelsurf"      , "vvelsurf"              , "m s-1", "land_ice_surface_y_velocity"),
#     "wvelbase": ISMIP6Var("wvelbase"      , "wvelbase"              , "m s-1", "land_ice_basal_upward_velocity"),
#     "wvelsurf": ISMIP6Var("wvelsurf"      , "wvelsurf"              , "m s-1", "land_ice_surface_upward_velocity")
# }

ismip6_vars = {
    "acabf": ISMIP6Var("acabf"         , "surface_mass_balance_average" , "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux"),
    "dlithkdt": ISMIP6Var("dlithkdt"      , "dHdt"                  , "m s-1", "tendency_of_land_ice_thickness"),
    "hfgeoubed": ISMIP6Var("hfgeoubed"     , "hfgeoubed"             , "W m-2", "upward_geothermal_heat_flux_at_ground_level"),
    "libmassbf": ISMIP6Var("libmassbf"     , "basal_mass_balance_average"                      , "kg m-2 s-1", "land_ice_basal_specific_mass_balance_flux"),
    "licalvf": ISMIP6Var("licalvf"       , "discharge_flux"        , "kg m-2 s-1", "land_ice_specific_mass_flux_due_to_calving"),
    "litempbot": ISMIP6Var("litempbot"     , "tempbase"              , "K", "land_ice_basal_temperature"),
    "litempsnic": ISMIP6Var("litempsnic"    , "tempsurf"              , "K", "temperature_at_ground_level_in_snow_or_firn"),
    "lithk": ISMIP6Var("lithk"         , "thk"                   , "m", "land_ice_thickness"),
    "orog": ISMIP6Var("orog"          , "usurf"                 , "m", "surface_altitude"),
    "sftgrf": ISMIP6Var("sftgrf"        , "sftgrf"                , "1", "grounded_ice_sheet_area_fraction"),
    "sftflf": ISMIP6Var("sftflf"        , "sftflf"                , "1", "floating_ice_sheet_area_fraction"),
    "sftgif": ISMIP6Var("sftgif"        , "sftgif"                , "1", "land_ice_area_fraction"),
    "strbasemag": ISMIP6Var("strbasemag"    , "taub_mag"                  , "Pa", "magnitude_of_land_ice_basal_drag"),
    "topg": ISMIP6Var("topg"          , "topg"                  , "m", "bedrock_altitude"),
    "uvelbase": ISMIP6Var("uvelbase"      , "uvelbase"              , "m s-1", "land_ice_basal_x_velocity"),
    "uvelmean": ISMIP6Var("uvelmean"      , "ubar"               , "m s-1", "land_ice_vertical_mean_x_velocity"),
    "uvelsurf": ISMIP6Var("uvelsurf"      , "uvelsurf"              , "m s-1", "land_ice_surface_x_velocity"),
    "vvelbase": ISMIP6Var("vvelbase"      , "vvelbase"              , "m s-1", "land_ice_basal_y_velocity"),
    "vvelmean": ISMIP6Var("vvelmean"      , "vbar"               , "m s-1", "land_ice_vertical_mean_y_velocity"),
    "vvelsurf": ISMIP6Var("vvelsurf"      , "vvelsurf"              , "m s-1", "land_ice_surface_y_velocity"),
    "wvelbase": ISMIP6Var("wvelbase"      , "wvelbase"              , "m s-1", "land_ice_basal_upward_velocity"),
    "wvelsurf": ISMIP6Var("wvelsurf"      , "wvelsurf"              , "m s-1", "land_ice_surface_upward_velocity")
}

ismip6_to_pism_dict = dict((k, v.pism_name) for k, v in ismip6_vars.iteritems())
pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars.iteritems())

def adjust_time_axis(filename):
    '''
    Adjusts the time axis
    '''
    
    nc = CDF(filename, 'a')
    time = nc.variables['time']
    time_bnds_var = time.bounds
    time_bnds = nc.variables[time_bnds_var]
    nt = len(time[:])
    new_timeline = np.linspace(0, 100, nt, endpoint=True)
    time[:] = new_timeline
    time_bnds[:,0] = new_timeline - 5
    time_bnds[:,1] = new_timeline
    time.units = 'years'
    time_bnds.units = 'years'
    nc.close()

def make_ismip6_conforming(filename):
    '''
    Make file ISMIP6 conforming
    '''
    
    # Open file
    nc = CDF(filename, 'a')

    cell_area_var = nc.variables['cell_area']
    cell_area_units = cell_area_var.units
    cell_area = cell_area_var[:]

    for pism_var in nc.variables:
        nc_var = nc.variables[pism_var]
        if pism_var in pism_to_ismip6_dict.keys():
            ismip6_var = pism_to_ismip6_dict[pism_var]
            print('Processing {} / {}'.format(pism_var, ismip6_var))
            if not pism_var == ismip6_var:
                print('  Renaming {pism_var} to {ismip6_var}'.format(pism_var=pism_var, ismip6_var=ismip6_var))
                nc.renameVariable(pism_var, ismip6_var)
                nc.sync()
            if pism_var in pism_mass_to_vol_vars:
                print('  Converting {pism_var} from mass to volume'.format(pism_var=pism_var))
                nc_var[:] /= cell_area
                i_units = nc_var.units
                o_units = cf_units.Unit(i_units) / cf_units.Unit(cell_area_units)
                nc_var.units = o_units.format()
            if not nc_var.units == ismip6_vars[ismip6_var].units:
                o_units = ismip6_vars[ismip6_var].units            
                i_units = nc_var.units
                print('  Converting {pism_var} from {i_units} to {o_units}'.format(pism_var=pism_var, i_units=i_units, o_units=o_units))    
                i_f = cf_units.Unit(i_units)
                o_f = cf_units.Unit(o_units)
                nc_var[:] = i_f.convert(nc_var[:], o_f)
                nc_var.units = o_units
                nc_var.standard_name = ismip6_vars[ismip6_var].standard_name
    nc.close()


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
    
    
    tmp_filename = 'tmp_{}.nc'.format(EXP.upper())
    tmp_file = os.path.join(project_dir, tmp_filename)
    try:
        os.remove(tmp_file)
    except OSError:
        pass
    # Check if request variables are present
    nc = CDF(infile, 'r')
    for m_var in pism_copy_vars:
        if m_var not in nc.variables:
            print("Requested variable '{}' missing".format(m_var))
    nc.close()
    cmd = ['ncks', '-O',
           '-v', '{}'.format(','.join(pism_copy_vars)),
           infile, tmp_file]
    sub.call(cmd)
    
    # Make the file ISMIP6 conforming
    make_ismip6_conforming(tmp_file)
    # Should be temporary until new runs
    nc2cdo_cmd = ['nc2cdo.py', tmp_file]
    sub.call(nc2cdo_cmd)
                
    # Create source grid definition file
    source_grid_filename = 'source_grid.nc'
    source_grid_file = os.path.join(project_dir, source_grid_filename)
    ncks_cmd = ['ncks', '-O', '-v', 'thk,mapping', infile, source_grid_file]
    sub.call(ncks_cmd)
    nc2cdo_cmd = ['nc2cdo.py', source_grid_file]
    sub.call(nc2cdo_cmd)

    # If exist, remove target grid description file
    target_grid_file = os.path.join(project_dir, target_grid_filename)
    try:
        os.remove(target_grid_file)
    except OSError:
        pass

    # Create target grid description file
    create_searise_grid(target_grid_file, target_resolution)
    
    # Generate weights if weights file does not exist yet
    cdo_weights_filename = 'searise_grid_{}m_weights.nc'.format(target_resolution)
    cdo_weights_file = os.path.join(project_dir, cdo_weights_filename)
    if (not os.path.isfile(cdo_weights_filename)) or (override_weights_file is True):
        print('Generating CDO weights file {}'.format(cdo_weights_file))
        cdo_cmd = ['cdo', '-P', '{}'.format(n_procs),
                   'gen{method},{grid}'.format(method=remap_method, grid=target_grid_file),
            source_grid_file,
            cdo_weights_file]
        sub.call(cdo_cmd)

    # Remap to SeaRISE grid
    
    out_filename = '{project}_{exp}.nc'.format(project=project, exp=EXP.upper())
    out_file = os.path.join(project_dir, out_filename)
    try:
        os.remove(out_file)
    except OSError:
        pass
    cdo_cmd = ['cdo', '-P', '{}'.format(n_procs),
               'remap,{},{}'.format(target_grid_file, cdo_weights_file),
               tmp_file,
               out_file]
    sub.call(cdo_cmd)

    # Adjust the time axis
    adjust_time_axis(out_file)

    vars_dir = os.path.join(project_dir, EXP)
    if not os.path.exists(vars_dir):
        os.mkdir(vars_dir)

    for m_var in ismip6_request_vars:
        final_file = '{}/{}_{}.nc'.format(vars_dir, m_var, project)
        print('Finalizing file {}'.format(final_file))
        # Generate file
        ncks_cmd = ['ncks', '-O',
                    '-v', m_var,
                    out_file,
                    final_file]
        sub.call(ncks_cmd)
        # Add stats vars
        ncks_cmd = ['ncks', '-A',
                    '-v', ','.join(pism_stats_vars),
                    tmp_file,
                    final_file]
        sub.call(ncks_cmd)
        # Add coordinate vars and mapping
        ncks_cmd = ['ncks', '-A', '-v', 'x,y,mapping',
                    target_grid_file,
                    final_file]
        sub.call(ncks_cmd)
        # Update attributes
        nc = CDF(final_file, 'a')
        try:
            nc_var = nc.variables[m_var]
            nc_var.mapping = 'mapping'
            nc_var.units = ismip6_vars[m_var].units
            nc_var.standard_name = ismip6_vars[m_var].standard_name
        except:
            pass
        nc.close()

    
