#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden

import os
import numpy as np
from argparse import ArgumentParser
from netCDF4 import Dataset as NC

try:
    import pypismtools.pypismtools as ppt
except:
    import pypismtools as ppt

# Set up the option parser
parser = ArgumentParser()
parser.description = "Analyze flux gates. Used for 'Complex Greenland Outlet Galcier Flow Captured'."
parser.add_argument("-a", "--anomaly_file", dest="anomaly_file",
                    help='''Anomaly smb file''')
parser.add_argument("-b", "--background_file", dest="background_file",
                    help='''Background smb file''')
parser.add_argument("OUTFILE", nargs=1)
options = parser.parse_args()
anomaly_file = options.anomaly_file
background_file = options.background_file

outfile = options.OUTFILE[0]

nc_a = NC(anomaly_file, 'r')
nc_b = NC(background_file, 'r')

try:
    os.remove(outfile)
except OSError:
    pass
nc = NC(outfile, 'w')

xdim_a, ydim_a, zdim_a, tdim_a = ppt.get_dims(nc_a)
xdim_b, ydim_b, zdim_b, tdim_b = ppt.get_dims(nc_b)

assert xdim_a == xdim_b
assert ydim_a == ydim_b

xdim = xdim_a
ydim = ydim_a
tdim = 'time'

nx = len(nc_a.dimensions[xdim_a])
ny = len(nc_a.dimensions[ydim_b])
nt = 100

nc.createDimension(xdim, size = (nx))
nc.createDimension(ydim, size = (ny))
nc.createDimension(tdim)

time_var = nc.createVariable(tdim, 'float64', dimensions=(tdim))
time_var.units = 'years'
time_var.axis = 'T'
time_var[:] = range(nt)

smb_background = nc_b.variables['climatic_mass_balance']
smb_background_untis = smb_background.units
smb_background_standard_name = smb_background.standard_name

temp_background = nc_b.variables['ice_surface_temp']
temp_background_untis = temp_background.units
temp_background_standard_name = temp_background.standard_name

smb_anomaly = nc_a.variables['climatic_mass_balance']
smb_anomaly_units = smb_anomaly.units

smb_var = nc.createVariable('climatic_mass_balance', 'float64', dimensions=(tdim, ydim, xdim))
smb_var.units = smb_background_untis
smb_var.standard_name = smb_background.standard_name

temp_var = nc.createVariable('ice_surface_temp', 'float64', dimensions=(tdim, ydim, xdim))
temp_var.units = temp_background_untis
temp_var.standard_name = temp_background.standard_name

for t in range(nt):
    temp_var[t,::] = np.squeeze(temp_background[:])
    if t <= 40:
        smb_var[t,::] = np.squeeze(smb_background[:]) + np.squeeze(smb_anomaly[:]) * np.floor(t) / 40
    else:
        smb_var[t,::] = np.squeeze(smb_background[:]) + np.squeeze(smb_anomaly[:])


nc.close()
nc_a.close()
nc_b.close()
