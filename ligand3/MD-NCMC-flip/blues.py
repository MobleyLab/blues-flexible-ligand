"""
example.py: Provides an example script to run BLUES and
benchmark the run on a given platform

Authors: Samuel C. Gill
Contributors: Nathan M. Lim, David L. Mobley

* Benchmarking related code adapted from:
https://github.com/pandegroup/openmm/blob/master/examples/benchmark.py
(Authors: Peter Eastman)
"""

from __future__ import print_function
from blues.moves import FlipMove 
from blues.engine import MoveEngine
from blues import utils
from blues.simulation import Simulation, SimulationFactory
import parmed
from parmed.openmm import RestartReporter
from mdtraj.reporters import NetCDFReporter
from simtk import openmm
from optparse import OptionParser
import sys
import logging


def runNCMC(platform_name, nstepsNC, nprop, outfname):

    #Generate the ParmEd Structure
    prmtop = '../complex_wat.prmtop'
    inpcrd = '../ligand3.equi.rst'
    struct = parmed.load_file(prmtop, xyz=inpcrd)
    print('Structure: %s' % struct.topology)

    nstepsNC = 1000

    #Define some options
    opt = { 'temperature' : 300.0, 'friction' : 1, 'dt' : 0.002,
            'nIter' : 25000, 'nstepsNC' : nstepsNC, 'nstepsMD' : 1000, 'nprop' : nprop,
            'nonbondedMethod' : 'PME', 'nonbondedCutoff': 10,
            'constraints': 'HBonds', 'freeze_distance' : 5.0,
            'trajectory_interval' : 1000, 'reporter_interval' : 1000,
            'ncmc_traj' : False, 'write_move' : False,
            'platform' : platform_name,
            'outfname' : 'gmx', 'NPT' : False,
            'verbose' : False}



    #Define the 'model' object we are perturbing here.
    # Calculate particle masses of object to be moved
    resname = "LIG"
    dihedral_atoms = ["C10", "C9", "C3", "C2" ]
    alch_list = ['C9', 'H92', 'H93', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'H1','H4','H6', 'C7', 'H71', 'H72', 'H73', 'C8', 'H81', 'H82', 'H83']
    ligand = FlipMove(struct, prmtop, inpcrd, dihedral_atoms, alch_list, resname)

    # Initialize object that proposes moves.
    ligand_mover = MoveEngine(ligand)

    # Generate the MD, NCMC, ALCHEMICAL Simulation objects
    simulations = SimulationFactory(struct, ligand_mover, **opt)
    simulations.createSimulationSet()

    # Add reporters to MD simulation.
    #traj_reporter = openmm.app.DCDReporter(outfname+'-nc{}.dcd'.format(nstepsNC), opt['trajectory_interval'])
    progress_reporter = openmm.app.StateDataReporter(sys.stdout, separator="\t",
                                reportInterval=opt['reporter_interval'],
                                step=True, totalSteps=opt['nIter']*opt['nstepsMD'],
                                time=True, speed=True, progress=True, remainingTime=True)
    #simulations.md.reporters.append(traj_reporter)
    simulations.md.reporters.append(progress_reporter)
    netcdf_reporter = NetCDFReporter('gmx.nc', 500)
    simulations.md.reporters.append(netcdf_reporter)

    # Run BLUES Simulation
    blues = Simulation(simulations, ligand_mover, **opt)
    blues.run(opt['nIter'])
    # saving last restart
    restrt = RestartReporter('blues.rst', 1, struct.ptr('natom') );
    state = simulations.md.context.getState(getPositions=True, getEnergy=True, getVelocities=True, enforcePeriodicBox=True)
    restrt.report(simulations.md, state)


parser = OptionParser()
parser.add_option('-f', '--force', action='store_true', default=False,
                  help='run BLUES example without GPU platform')
parser.add_option('-n','--ncmc', dest='nstepsNC', type='int', default=5000,
                  help='number of NCMC steps')
parser.add_option('-p','--nprop', dest='nprop', type='int', default=5,
                  help='number of propgation steps')
parser.add_option('-o','--output', dest='outfname', type='str', default="blues",
                  help='Filename for output DCD')
(options, args) = parser.parse_args()



platformNames = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
if 'CUDA' in platformNames:
    runNCMC('CUDA', options.nstepsNC, options.nprop, options.outfname)
elif 'OpenCL' in platformNames:
    runNCMC('OpenCL',options.nstepsNC, options.nprop, options.outfname)
else:
    if options.force:
        runNCMC('CPU', options.nstepsNC, options.outfname)
    else:
        print('WARNING: Could not find a valid CUDA/OpenCL platform. BLUES is not recommended on CPUs.')
        print("To run on CPU: 'python blues/example.py -f'")


