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


def make_scalar_vars_ismip6_conforming(filename, ismip6_vars_dict):
    '''
    Make file ISMIP6 conforming
    '''
    
    # Open file
    nc = CDF(filename, 'a')

    pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars_dict.iteritems())
    
    for pism_var in nc.variables:
        nc_var = nc.variables[pism_var]
        if pism_var in pism_to_ismip6_dict.keys():
            ismip6_var = pism_to_ismip6_dict[pism_var]
            print('Processing {} / {}'.format(pism_var, ismip6_var))
            if not pism_var == ismip6_var:
                print('  Renaming {pism_var} to {ismip6_var}'.format(pism_var=pism_var, ismip6_var=ismip6_var))
                nc.renameVariable(pism_var, ismip6_var)
                nc.sync()
            if not nc_var.units == ismip6_vars_dict[ismip6_var].units:
                o_units = ismip6_vars_dict[ismip6_var].units            
                i_units = nc_var.units
                print('  Converting {pism_var} from {i_units} to {o_units}'.format(pism_var=pism_var, i_units=i_units, o_units=o_units))    
                i_f = cf_units.Unit(i_units)
                o_f = cf_units.Unit(o_units)
                nc_var[:] = i_f.convert(nc_var[:], o_f)
                nc_var.units = o_units
                nc_var.standard_name = ismip6_vars_dict[ismip6_var].standard_name
    nc.close()

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
    make_scalar_vars_ismip6_conforming(out_file, ismip6_vars_dict)

    # Update attributes
    print('Adjusting attributes')
    nc = CDF(out_file, 'a')
    nc.Conventions = 'CF-1.6'
    nc.close()
    print('Finished processing scalars file {}'.format(out_file))

    
