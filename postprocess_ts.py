#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden

import os
import numpy as np
import csv
import cf_units
try:
    import subprocess32 as sub
except:
    import subprocess as sub
    
from argparse import ArgumentParser
from netCDF4 import Dataset as CDF
from resources_ismip6 import *

# Set up the option parser
parser = ArgumentParser()
parser.description = "Script to make ISMIP6-conforming scalar time series."
#parser.add_argument("INIT_FILE", nargs=1)
parser.add_argument("EXP_FILE", nargs=1)
parser.add_argument("-e", "--experiment", dest="experiment",
                    choices=['ctrl', 'asmb'],
                    help="Output size type", default='ctrl')
parser.add_argument("-t", "--target_resolution", dest="target_resolution", type=int,
                    choices=[1000, 5000],
                    help="Horizontal grid resolution", default=1000)

options = parser.parse_args()
experiment = options.experiment
infile = options.EXP_FILE[0]
target_resolution = options.target_resolution

# Need to get grid resolution from file
nc = CDF(infile, 'r')
pism_grid_dx = int(round(nc.variables['run_stats'].grid_dx_meters))
nc.close()
PISM_GRID_RES_ID = str(pism_grid_dx / 100)
TARGET_GRID_RES_ID = str(target_resolution / 1000)

IS = 'GIS'
GROUP = 'UAF'
MODEL = 'PISM' + PISM_GRID_RES_ID
EXP = experiment
TYPE = '_'.join([EXP, '0' + TARGET_GRID_RES_ID])
project = '{IS}_{GROUP}_{MODEL}'.format(IS=IS, GROUP=GROUP, MODEL=MODEL)

pism_stats_vars = ['pism_config',
                   'run_stats']

ismip6_vars_dict = get_ismip6_vars_dict('ismip6vars.csv', 1)
ismip6_to_pism_dict = dict((k, v.pism_name) for k, v in ismip6_vars_dict.iteritems())
pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars_dict.iteritems())

pism_copy_vars = [x for x in (ismip6_to_pism_dict.values() + pism_stats_vars)]


if __name__ == "__main__":


    project_dir = os.path.join(GROUP, MODEL, TYPE)
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
    
    out_filename = 'scalars_{project}_{exp}.nc'.format(project=project, exp=EXP)
    out_file = os.path.join(project_dir, out_filename)
    try:
        os.remove(out_file)
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
           infile, out_file]
    sub.call(cmd)
    
    # Adjust the time axis
    print('Adjusting time axis')
    adjust_time_axis(out_file)

#     vars_dir = os.path.join(project_dir, EXP)
#     if not os.path.exists(vars_dir):
#         os.mkdir(vars_dir)

#     for m_var in ismip6_vars_dict.keys():
#         final_file = '{}/{}_{}.nc'.format(vars_dir, m_var, project)
#         print('Finalizing variable {}'.format(m_var))
#         # Generate file
#         print('  Copying to file {}'.format(final_file))
#         ncks_cmd = ['ncks', '-O',
#                     '-v', m_var,
#                     out_file,
#                     final_file]
#         sub.call(ncks_cmd)
#         # Add stats vars
#         print('  Adding config/stats variables')
#         ncks_cmd = ['ncks', '-A',
#                     '-v', ','.join(pism_stats_vars),
#                     tmp_file,
#                     final_file]
#         sub.call(ncks_cmd)
#         # Add coordinate vars and mapping
#         print('  Adding coordinte and mapping variables')
#         ncks_cmd = ['ncks', '-A', '-v', 'x,y,mapping',
#                     target_grid_file,
#                     final_file]
#         sub.call(ncks_cmd)
#         # Update attributes
#         print('  Adjusting attributes')
#         nc = CDF(final_file, 'a')
#         try:
#             nc_var = nc.variables[m_var]
#             nc_var.mapping = 'mapping'
#             nc_var.units = ismip6_vars_dict[m_var].units
#             nc_var.standard_name = ismip6_vars_dict[m_var].standard_name
#         except:
#             pass
# #        nc.variables['pism_config'].delncattr('calendar')
#         nc.Conventions = 'CF-1.6'
#         nc.close()
#         print('  Done finalizing variable {}'.format(m_var))

    
