## This is a very quick skeleton driver script for umbrella sampling, created by David Mobley, Kalistyn Burley
## This script  sets the values of a specified dihedral angle to a given value and adds harmonic restraint for rotating the dihedral bond. The simulation box is then simulated for 10 ns.
## Usage: python3 simulate_umbrella.py #dihedralAngle #torsionConstant #indexWindow 

import numpy as np
import os
import parmed
from parmed.openmm import topsystem, RestartReporter
from parmed import unit as u
import sys

import simtk.openmm as mm
from simtk.openmm import Platform
from simtk.openmm import app
from simtk.unit import *
# There need to be more imports here relating to OpenMM

from openeye.oechem import *


def runOpenMM(parm, inpcrdFile, system, rad, K, Indices, solvate, out_dcd, out_csv, out_rst ):
    """

    Minimize molecule with specified topology, system, and positions
       using OpenMM. Return the positions of the optimized moelcule.

    Parameters
    ----------
    parm:      ParmEd object for system (parameters/corodinates)
    inpcrd:    Inpcrd file name from the last restart point
    system:    OpenMM system for this mol's Prmtop and Inpcrd files
    rad:       float of reference angle around which to restrain
    K:         float of force constant for harmonic restraint
    Indices:   list of integers of zero-based atom indices for restraint
    out_dcd:   filename for output dcd file for constant pressure equilibration/production
    out_csv:   filename for output csv file for constant pressure equilibration/production

    Returns
    -------
    topology.positions: OpenMM topology positions for minimized mol

    """

    
    def newIntegrator():
        integrator = mm.LangevinIntegrator(
                300.0 * u.kelvin,
                10.0 / u.picosecond,
                2.0 * u.femtosecond)
        return integrator

    def pmdStructureToOEMol(parm, resname):

        from oeommtools.utils import openmmTop_to_oemol
        mask = "!(:%s)" %resname
        structure_LIG = parmed.load_file( '../2gmx_wat.prmtop', xyz = '../equilibration/rst/step8.rst.125000' )
        structure_LIG.strip(mask)
        pos = structure_LIG.positions
        top = structure_LIG.topology
        molecule = openmmTop_to_oemol(top, pos, verbose=False)
        OEPerceiveBondOrders(molecule)
        OEAssignAromaticFlags(molecule)
        OEFindRingAtomsAndBonds(molecule)

        return molecule
  
    def getAtomIndices( structure, resname ):
        """
        Get atom indices of a ligand from ParmEd Structure.
        Arguments
        ---------
        resname : str
            String specifying the resiue name of the ligand.
        structure: parmed.Structure
            ParmEd Structure object of the atoms to be moved.
        Returns
        -------
        atom_indices : list of ints
            list of atoms in the coordinate file matching lig_resname
        """
        atom_indices_ligand = []
        topology = structure.topology
        for atom in topology.atoms():
           if str(resname) in atom.residue.name:
              atom_indices_ligand.append(atom.index)

        return atom_indices_ligand


    """
    Rotate the torsion to an angle rad using openeye toolkits
    """   
    molecule = pmdStructureToOEMol( parm, "LIG" )
    atom_indices_ligand = getAtomIndices( parm, "LIG" )


    dihedral_atoms = ["C10", "C9", "C3", "C2" ]
    atom1 = molecule.GetAtom(OEHasAtomName(dihedral_atoms[0]))
    atom2 = molecule.GetAtom(OEHasAtomName(dihedral_atoms[1]))
    atom3 = molecule.GetAtom(OEHasAtomName(dihedral_atoms[2]))
    atom4 = molecule.GetAtom(OEHasAtomName(dihedral_atoms[3]))
    if OESetTorsion(molecule, atom1, atom2, atom3, atom4, rad ) == False :
       print("Torsional bond couldn't be rotated. Please enter correct atoms!"); 
       exit()

    # Update ligand positions in nc_sim
    updated_pos = molecule.GetCoords()

    for index, atomidx in enumerate(atom_indices_ligand): 
        parm.positions[atomidx] = np.array(updated_pos[index])*u.nanometers

    """
    harmonically restrain dihedral angle
    see units, http://docs.openmm.org/6.3.0/userguide/theory.html
    """
    pi = np.pi
    harmonic = mm.CustomTorsionForce("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = %.5f" % pi);
    harmonic.addPerTorsionParameter("theta0");
    harmonic.addPerTorsionParameter("k");
    system.addForce(harmonic)
    harmonic.addTorsion(Indices[0], Indices[1], Indices[2], Indices[3], (rad, K))

    # Restraint non-moving part of the ligand
    restraintWt = 200 #kcal/mol/A2
    # define the custom force to restrain atoms to their starting positions
    force_restr = mm.CustomExternalForce('k_restr*periodicdistance(x, y, z, x0, y0, z0)^2')
    # Add the restraint weight as a global parameter in kcal/mol/A^2
    force_restr.addGlobalParameter("k_restr", restraintWt*u.kilocalories_per_mole/u.angstroms**2)
    # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
    force_restr.addPerParticleParameter("x0")
    force_restr.addPerParticleParameter("y0")
    force_restr.addPerParticleParameter("z0")
    alch_list = ['C9', 'H92', 'H93', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'H1', 'H2', 'H4', 'H5', 'H6']
    for idx, atom_crd in enumerate( parm.positions ):
            name=parm.atoms[idx].name;
            resname=parm.atoms[idx].residue.name;
            if resname == "LIG":
              if not name in alch_list:
                xyz = parm.positions[idx].in_units_of(u.nanometers)/u.nanometers
                force_restr.addParticle(idx, xyz)
    system.addForce( force_restr )

    # build simulaion
    platform = mm.Platform.getPlatformByName('CUDA')
    integ1 = newIntegrator()
    simulation = app.Simulation(parm.topology, system, integ1)
    simulation.context.setPositions( parm.positions )

    # Set Box dimensions
    inpcrd = app.AmberInpcrdFile( inpcrdFile );
    if inpcrd.boxVectors is not None:
       simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    print('RESTARTING simulation from a previous State..........%s' %inpcrdFile)
    velocities = parm.velocities 
    simulation.context.setVelocities( inpcrd.velocities ) 

    # perform minimization
    print('Minimizing...')
    simulation.minimizeEnergy( tolerance = 0.5 * kilojoule/mole  )
    
    # adding simulation reporters
    simulation.context.setVelocitiesToTemperature(300*u.kelvin)
    simulation.reporters.append(app.DCDReporter(out_dcd, 1000))
    simulation.reporters.append(app.StateDataReporter(csv_file, 1000, step=True, potentialEnergy=True, totalEnergy=True, volume=True,temperature=True, separator='\t'))
    restrt = RestartReporter( out_rst, 10000000, parm.ptr('natom') );
    state = simulation.context.getState(getPositions=True, getEnergy=True, getVelocities=True, enforcePeriodicBox=True)
    restrt.report(simulation, state)


    print('Production run at NVT...')
    simulation.step(5000000) # 10 ns
    
    # saving last restart
    state = simulation.context.getState(getPositions=True, getEnergy=True, getVelocities=True, enforcePeriodicBox=True)
    restrt.report(simulation, state)
    return 


## Load umbrella centers and force constants
umbr_num = int(sys.argv[1])            # index for window
center = int(sys.argv[2]) * np.pi/180. #Convert to radians
fc = int(sys.argv[3])                  #force constant

# Let's assume we already have input input files
prmFile = '../../complex_wat.prmtop'
# Load inputs; restarting from the endpoint of previous umbrella-sukanya 
inpcrd = '../../ligand1.equi.rst'

if not os.path.isfile(prmFile):
    raise IOError('Missing input file %s' %prmFile)
if not os.path.isfile(inpcrd):
    raise IOError('Missing input file %s' %inpcrd)

# specifying index numbers of atoms describing the dihedral angle
atomlist=[ 5955, 5954, 5950, 5949 ]

# Note that if you want to make this run faster, you can trivially parallelize by
# just making this script run one particular center and having a separate
# shell script run each umbrella

print( "Umbrella number " , umbr_num )
#Retrieve force constants, file names, etc.
prefix = 'umbrella%d' % umbr_num
dcd_file = 'dcd/'+prefix+'.dcd'
csv_file = 'csv/'+prefix+'.csv'
rst_file = 'rst/'+prefix+'.rst'

parm = parmed.load_file(prmFile, inpcrd)

# Create OpenMM system
system = parm.createSystem(nonbondedMethod=app.PME,constraints=app.HBonds,rigidWater=True,removeCMMotion=True)
# I would probably add a check here that the atoms actually have the names you expect them to have.

# Run simulation
runOpenMM(parm, inpcrd, system, center, fc, atomlist, False, dcd_file, csv_file, rst_file )

