#!/usr/bin/env python
# Copyright (C) 2014-2015 Andy Aschwanden

import itertools
import os
from argparse import ArgumentParser

# Set up the option parser
parser = ArgumentParser()
parser.description = "Generating scripts for parameter study."
parser.add_argument("regridfile", nargs=1)
parser.add_argument("-n", '--n_procs', dest="n", type=int,
                    help='''number of cores/processors. default=64.''', default=64)
parser.add_argument("-w", '--wall_time', dest="walltime",
                    help='''walltime. default: 12:00:00.''', default="12:00:00")
parser.add_argument("-q", '--queue', dest="queue", choices=['standard_4', 'standard_16', 'standard', 'gpu', 'gpu_long', 'long', 'normal'],
                    help='''queue. default=standard_4.''', default='standard_4')
parser.add_argument("-c", "--climate", dest="climate",
                    choices=['const'],
                    help="climate", default='const')
parser.add_argument("-d", "--domain", dest="domain",
                    choices=['greenland', 'jakobshavn'],
                    help="sets the modeling domain", default='greenland')
parser.add_argument("-f", "--o_format", dest="oformat",
                    choices=['netcdf3', 'netcdf4_parallel', 'pnetcdf'],
                    help="output format", default='netcdf4_parallel')
parser.add_argument("-g", "--grid", dest="grid", type=int,
                    choices=[18000, 9000, 4500, 3600, 1800, 1500, 1200, 900, 600, 450],
                    help="horizontal grid resolution", default=1500)
parser.add_argument("--o_size", dest="osize",
                    choices=['small', 'medium', 'big', '2dbig'],
                    help="output size type", default='2dbig')
parser.add_argument("-s", "--system", dest="system",
                    choices=['pleiades', 'fish', 'pacman', 'debug'],
                    help="computer system to use.", default='pacman')
parser.add_argument("-t", "--etype", dest="etype",
                    choices=['ctrl', 'old_bed', 'ba01_bed', '970mw_hs', 'jak_1985', 'cresis'],
                    help="subglacial topography type", default='ctrl')
parser.add_argument("--dataset_version", dest="version",
                    choices=['1.1', '1.2', '2'],
                    help="Input data set version", default='2')

options = parser.parse_args()
args = options.regridfile

nn = options.n
oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system

climate = options.climate
grid = options.grid
etype = options.etype
version = options.version

domain = options.domain
if domain.lower() in ('greenland'):
    pism_exec = 'pismr'
elif domain.lower() in ('jakobshavn'):
    x_min = -280000
    x_max = 320000
    y_min = -2410000
    y_max = -2020000
    pism_exec = '''\'pismo -x_range {x_min},{x_max} -y_range {y_min},{y_max} -bootstrap\''''.format(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
else:
    print('Domain {} not recognized, exiting'.format(domain))
    import sys
    sys.exit(0)


def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result



def uniquify_list(seq, idfun=None):
    '''
    Remove duplicates from a list, order preserving.
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    '''

    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def make_pbs_header(system, cores, walltime, queue):
    systems = {}
    systems['debug'] = {}
    systems['fish'] = {'gpu' : 16,
                       'gpu_long' : 16,
                       'standard' : 12}
    systems['pacman'] = {'standard_4' : 4,
                        'standard_16' : 16}
    systems['pleiades'] = {'long' : 20,
                           'normal': 20}

    assert system in systems.keys()
    if system not in 'debug':
        assert queue in systems[system].keys()
        assert cores > 0

        ppn = systems[system][queue]
        nodes = cores / ppn

    if system in ('debug'):

        header = ''
        
    elif system in ('pleiades'):
        
        header = """
#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=ivy
#PBS -j oe

cd $PBS_O_WORKDIR

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn)
    else:
        header = """
#!/bin/bash
#PBS -q {queue}
#PBS -l walltime={walltime}
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -j oe

cd $PBS_O_WORKDIR

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn)

    return header

    
regridfile=args[0]
infile = ''
pism_dataname = 'pism_Greenland_{}m_mcb_jpl_v{}_{}.nc'.format(grid, version, etype)
etype = '{}_v{}'.format(etype, version)
dura = 100


# ########################################################
# set up parameter sensitivity study: tillphi
# ########################################################

hydro = 'null'

sia_e = (1.25)
ppq = (0.6)
tefo = (0.02)
ssa_n = (3.25)
ssa_e = (1.0)

phi_min_values = [5.0]
phi_max_values = [40.]
topg_min_values = [-700]
topg_max_values = [700]
combinations = list(itertools.product(phi_min_values, phi_max_values, topg_min_values, topg_max_values))

tsstep = 'daily'
exstep = 'yearly'
regridvars = 'litho_temp,enthalpy,tillwat,bmelt,Href,thk'

scripts = []
posts = []
for n, combination in enumerate(combinations):

    phi_min, phi_max, topg_min, topg_max = combination

    ttphi = '{},{},{},{}'.format(phi_min, phi_max, topg_min, topg_max)

    experiment='{}_{}_sia_e_{}_ppq_{}_tefo_{}_ssa_n_{}_ssa_e_{}_phi_min_{}_phi_max_{}_topg_min_{}_topg_max_{}_hydro'.format(climate, etype, sia_e, ppq, tefo, ssa_n, ssa_e, phi_min, phi_max, topg_min, topg_max, hydro)
    script = 'initMIP_{}_g{}m_{}.sh'.format(domain.lower(), grid, experiment)
    scripts.append(script)
    post = 'initMIPq_{}_g{}m_{}_post.sh'.format(domain.lower(), grid, experiment)
    posts.append(post)
    
    for filename in (script, post):
        try:
            os.remove(filename)
        except OSError:
            pass

    pbs_header = make_pbs_header(system, nn, walltime, queue)
        
    
    os.environ['PISM_EXPERIMENT'] = experiment
    os.environ['PISM_TITLE'] = 'PISM initMIP'
    
    with open(script, 'w') as f:

        f.write(pbs_header)

        exp_type = 'ctrl'
        
        outfile = '{domain}_g{grid}m_{experiment}_{dura}a_{exp_type}.nc'.format(domain=domain.lower(),grid=grid, experiment=experiment, dura=dura, exp_type=exp_type)
            
        params_dict = dict()
        params_dict['PISM_DO'] = ''
        params_dict['PISM_OFORMAT'] = oformat
        params_dict['PISM_OSIZE'] = osize
        params_dict['PISM_EXEC'] = pism_exec
        params_dict['PISM_DATANAME'] = pism_dataname
        params_dict['PISM_SURFACE_BCFILE'] = 'initMIP_climate_forcing_{dura}a_{exp_type}.nc'.format(dura=dura, exp_type=exp_type)
        params_dict['REGRIDFILE'] = regridfile
        params_dict['TSSTEP'] = tsstep
        params_dict['EXSTEP'] = exstep
        params_dict['REGRIDVARS'] = regridvars
        params_dict['SIA_E'] = sia_e
        params_dict['SSA_E'] = ssa_e
        params_dict['SSA_N'] = ssa_n
        params_dict['PARAM_NOAGE'] = 'foo'
        params_dict['PARAM_PPQ'] = ppq
        params_dict['PARAM_TEFO'] = tefo
        params_dict['PARAM_TTPHI'] = ttphi
        
        
        params = ' '.join(['='.join([k, str(v)]) for k, v in params_dict.items()])
        cmd = ' '.join([params, './run.sh', str(nn), climate, str(dura), str(grid), 'hybrid', hydro, outfile, infile, '2>&1 | tee job.${PBS_JOBID}'])

        f.write(cmd)
        f.write('\n\n')


        exp_type = 'asmb'
        
        outfile = '{domain}_g{grid}m_{experiment}_{dura}a_{exp_type}.nc'.format(domain=domain.lower(),grid=grid, experiment=experiment, dura=dura, exp_type=exp_type)
            
        params_dict['PISM_SURFACE_BCFILE'] = 'initMIP_climate_forcing_{dura}a_{exp_type}.nc'.format(dura=dura, exp_type=exp_type)
                
        params = ' '.join(['='.join([k, str(v)]) for k, v in params_dict.items()])
        cmd = ' '.join([params, './run.sh', str(nn), climate, str(dura), str(grid), 'hybrid', hydro, outfile, infile, '2>&1 | tee job.${PBS_JOBID}'])

        f.write(cmd)
        f.write('\n')

        if etype in 'ctrl_v2':
            mytype = "MO14 2015-04-27"
        elif etype in 'cresis_v2':
            mytype = "MO14+CReSIS 2015-04-27"
        elif etype in 'ctrl_v1.2':
            mytype = "MO14 2014-11-19"
        elif etype in ('ctrl', 'ctrl_v1.1'):
            mytype = "MO14 2014-06-26"
        elif etype in ('old_bed', 'old_bed_v1.1', 'old_bed_v1.2', 'old_bed_v2'):
            mytype = "BA01"
        elif etype in 'searise':
            mytype = "SR13"
        else:
            import sys
            print('etype {} not recognized, exiting'.format(etype))
            sys.exit(0)



        
scripts = uniquify_list(scripts)

submit = 'submit_{domain}_g{grid}m_{climate}_{etype}_tillphi.sh'.format(domain=domain.lower(), grid=grid, climate=climate, etype=etype)
try:
    os.remove(submit)
except OSError:
    pass

with open(submit, 'w') as f:

    f.write('#!/bin/bash\n')

    for k in range(len(scripts)):
        f.write('JOBID=$(qsub {script})\n'.format(script=scripts[k]))
        #f.write('qsub -W depend=afterok:${{JOBID}} {post}\n'.format(post=posts[k]))

print("\nRun {} to submit all jobs to the scheduler\n".format(submit))

