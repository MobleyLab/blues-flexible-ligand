#!/usr/bin/env python
from __future__ import division, print_function
# Usage : python3 step8.py

import sys

# OpenMM Imports
import simtk.openmm as mm
import simtk.openmm.app as app

# ParmEd Imports
from parmed import load_file, unit as u
from parmed.openmm import StateDataReporter, NetCDFReporter, RestartReporter

# Load the Amber files
print('Loading AMBER files...')
inpcrdName='rst/step7.rst.20000';
gmx_solv = load_file('complex_wat.prmtop', xyz=inpcrdName );
inpcrd=app.AmberInpcrdFile(inpcrdName);

#minimization info
minimizationTrue = False;
minStep = 4000;

#restart info
restartSim = True;

#simulation info
tempK = 300;
simTime = 0.5*1000; #units of picoseconds, total 5 ns simulation
simType = 'npt'; pressureAtm = 1;

#filename
outFname='step8';
mdoutInterval = 100 #for mdout file
netcdfInterval = 2500 #for mdout file
rstInterval = 5000000 #for rst file

# Time step in ps
stepLen = 0.002 * u.picoseconds
#HMR option
constraints = app.HBonds
steps=int(round(simTime/(stepLen.in_units_of(u.picoseconds)/u.picoseconds)));

# Create the OpenMM system
print('Creating OpenMM System')
system = gmx_solv.createSystem(nonbondedMethod=app.PME,
                                nonbondedCutoff=10.0*u.angstroms,
                                constraints=app.HBonds, removeCMMotion=True
)

#box = gmx_solv.box;
#print("Starting with new velocities ....");
printfile = sys.stdout;

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        tempK*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
#                        2.0*u.femtoseconds, # Time step
                        stepLen # Time step
)

# Add Force Barostat to the system
if simType == 'npt':
   system.addForce(mm.MonteCarloBarostat(pressureAtm*u.atmospheres, tempK*u.kelvin, 25))

# Define the platform to use; OpenCL, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('OpenCL')
#prop = dict(CudaPrecision='mixed') # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(gmx_solv.topology, system, integrator, platform)
#sim = app.Simulation(gmx_solv.topology, system, integrator, platform, prop)

# Set the particle positions
sim.context.setPositions(gmx_solv.positions)

# Set Box dimensions
if inpcrd.boxVectors is not None:
    sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
#sim.context.setPeriodicBoxVectors((box[0],0,0), (0, box[1],0), (0,0, box[2]));
#Restart reporter
restrt = RestartReporter('rst/'+outFname+'.rst', rstInterval, gmx_solv.ptr('natom') );
sim.reporters.append(restrt);

# If the velocities are not present in the Parmed structure
# new velocity vectors are generated otherwise the system is
# restarted from the previous State
if restartSim:
    print('RESTARTING simulation from a previous State..........%s' %inpcrdName)
    velocities = gmx_solv.velocities;
    sim.context.setVelocities(inpcrd.velocities)
else:
    # Set the velocities drawing from the Boltzmann distribution at the selected temperature
    print('GENERATING a new starting State.........')
    sim.context.setVelocitiesToTemperature(tempK*u.kelvin)

####################################################################
########          Minimize the energy   ############################
####################################################################
if minimizationTrue:
   print("#"*30);
   print('Minimizing energy.....')
   print('Minimization steps: {}'.format(minStep))
   
   state = sim.context.getState(getEnergy=True)
   
   print('Initial energy = {}'.format(state.getPotentialEnergy()), file=printfile)
   
   sim.minimizeEnergy(maxIterations=minStep)
   
   state = sim.context.getState(getPositions=True, getEnergy=True, getVelocities=True, enforcePeriodicBox=True)
   restrt.report(sim, state)
   
   print('Minimized energy = {}'.format(state.getPotentialEnergy()), file=printfile)
   print("#"*30);

####################################################################


# Set up the reporters to report energies and coordinates every * steps
sim.reporters.append(
        StateDataReporter('mdout/'+outFname+'.mdout', separator="\t",
                                               reportInterval=mdoutInterval,step=True,
                                               potentialEnergy=True, totalEnergy=True,
                                               volume=True, temperature=True, density=True)
)

# progress_reporter: prints simulation progress to 'sys.stdout'
progress_reporter = app.StateDataReporter('mdinfo/'+ outFname+'.mdinfo', separator="\t",
                                                  reportInterval=mdoutInterval,
                                                  step=True, totalSteps=steps,
                                                  time=True, speed=True, progress=True,
                                                  elapsedTime=True, remainingTime=True)
sim.reporters.append(progress_reporter)

# writes trajectory to '.nc' file. AMBER NetCDF(3.0)
sim.reporters.append(NetCDFReporter('netcdf/'+outFname+'.nc', netcdfInterval, crds=True))

# OpenMM platform information
mmver = mm.version.version
mmplat = sim.context.getPlatform()

# Host information
from platform import uname
for k,v in uname()._asdict().items():
    print(k, ':', v, file=printfile)
# Platform properties
for prop in mmplat.getPropertyNames():
    val = mmplat.getPropertyValue(sim.context, prop)
    print(prop, ':', val, file=printfile)

# Run dynamics
print('Running %d ps = %d steps of %s at %d K' %(simTime, steps, simType, tempK))
print('Running dynamics')

###### Real MD !! ###########
sim.step(steps)

# saving last restart
state = sim.context.getState(getPositions=True, getEnergy=True, getVelocities=True, enforcePeriodicBox=True)
restrt.report(sim, state)
