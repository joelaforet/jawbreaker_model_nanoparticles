'''

plumed_OpenMM.py
by Joe Laforet Jr.
jrl78@duke.edu

Inspired by CHARMM-GUI OpenMM python scripts
(http://www.charmm-gui.org)


usage: plumed_OpenMM.py [-h] -p TOPOLOGY -c COORDINATE -i INPFILE
         [-t TOPPAR] [-b SYSINFO] -ff CHARMM [-icrst RSTFILE]
         [-irst RSTFILE] [-ichk CHKFILE] [--savechk]
         [--loadchk] [-c {0,1,2,3}] [-p] [-i]
         [--solvent_type {1,2,3,4,5}] [-T TEMPERATURE]
         [--timestep TIMESTEP] [--interval INTERVAL]
         [-l SIMULATION_LENGTH] [--timeit]

    optional arguments:
      -h, --help            show this help message and exit
      -p TOPOLOGY
                Input topology file
      -c COORDINATE
                Input coordinate file
      -i INPFILE
                Input paramter file
      -t TOPPAR
                Input CHARMM-GUI toppar stream file (optional)
      -b SYSINFO
                Input CHARMM-GUI sysinfo stream file (optional)
      -ff FORCE_FIELD
                Input force field type (default: CHARMM)
      -icrst RSTFILE
                Input CHARMM RST file (optional)
      -irst RSTFILE
                Input AMBER RST file (optional)
      -ichk CHKFILE
                Input checkpoint file (optional)
      -opdb PDBFILE
                Output PDB file (optional)
      -orst RSTFILE
                Output restart file (optional)
      -ochk CHKFILE
                Output checkpoint file (optional)
      -odcd DCDFILE
                Output trajectory file (optional)
      -rewrap
                Re-wrap the coordinates in a molecular basis. Default = False. (optional)



      -o OUTPUT, --output OUTPUT
                Output pdb filename
      --savechk             If activated, save a checkpoint state at every 10
                intervals
      --loadchk             If activated, load a checkpoint state with .chk
                extension of the output filename
      -c {0,1,2,3}, --constraint {0,1,2,3}
                0 = None (default); 1 = HBonds ; 2 = AllBonds ; 3 =
                HAngles
      -p, --periodic        If activated, runs the simulation in periodic box with
                PME method used for long range interaction (default =
                NoCutoff)
      -i, --implicit        If activated, runs the simulation in implicit water
                solvent (default = vacuum, unless explicitly solvated)
      --solvent_type {1,2,3,4,5}
                1 = HCT (default); 2 = OBC1 ; 3 = OBC2 ; 4 = GBn ; 5 =
                GBn2
      -T TEMPERATURE, --temperature TEMPERATURE
                Set simulation temperature (default = 300K)
      --timestep TIMESTEP   Set simulation time step in units of picosecond
                (default = 0.002 picosecond)
      --interval INTERVAL   Set interval of saved frames in the unit of picosecond
                (default = 10 ps)
      -l SIMULATION_LENGTH, --simulation_length SIMULATION_LENGTH
                Set duration of simulation time in units of nanosecond
                (default = 20ns)
      --timeit              If activated, creates time_log file
'''

from openmm.app import *
from openmm import *
from openmm.unit import *
from __future__ import print_function
from sys import stdout
import sys
import time
import os
import argparse

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *
from omm_rewrap import *

from openmmplumed import PlumedForce

CUR_DIR = os.getcwd()
start_time = time.time()

parser = argparse.ArgumentParser()
# System description input files
parser.add_argument('-p', dest='topfile', help='Input topology file', required=True)
parser.add_argument('-c', dest='crdfile', help='Input coordinate file', required=True)
parser.add_argument('-i', dest='inpfile', help='Input parameter file', required=True)

#Stuff for CHARMM FF
parser.add_argument('-t', dest='toppar', help='Input CHARMM-GUI toppar stream file (optional)')

#Force Field to use
parser.add_argument('-ff', dest='fftype', help='Input force field type. (CHARMM or AMBER). (default: CHARMM)', default='CHARMM')

#Loading input/output files
parser.add_argument("--savechk", action="store_true", help="If activated, save a checkpoint state at every 10 intervals")
parser.add_argument("--loadchk", action="store_true", help="If activated, load a checkpoint state with .chk extension of the output filename")
parser.add_argument('-icrst', metavar='RSTFILE', dest='icrst', help='Input CHARMM RST file (optional)')
parser.add_argument('-irst', metavar='RSTFILE', dest='irst', help='Input restart file (optional)')
parser.add_argument('-opdb', metavar='PDBFILE', dest='opdb', help='Output PDB file. Useful for PyMol. (optional)')
parser.add_argument('-orst', metavar='RSTFILE', dest='orst', help='Output restart file (optional)')
parser.add_argument('-odcd', metavar='DCDFILE', dest='odcd', help='Output DCD trajectory file. Useful for VMD/Chimera. (optional)')
parser.add_argument('-rewrap', dest='rewrap', help='Re-wrap the coordinates in a molecular basis (optional)', action='store_true', default=False)

#Parameters for Production Simulation. All sims use the same 6-step equilibration scheme.
#parser.add_argument('-c','--constraint', type=int, default=0, choices=range(0, 4), help="0 = None (default); 1 = HBonds ; 2 = AllBonds ; 3 = HAngles")
#parser.add_argument('-p','--periodic', action="store_true", help="If activated, runs the simulation in periodic box with PME method used for long range interaction (default = NoCutoff)")
#parser.add_argument('--implicit', action="store_true", help="If activated, runs the simulation in implicit water solvent (default = vacuum, unless explicitly solvated)")
#parser.add_argument('--solvent_type', type=int, default=1, choices=range(1, 6), help="Activated if '--implicit' is specified. 1 = HCT ; 2 = OBC1 ; 3 = OBC2 (default) ; 4 = GBn ; 5 = GBn2 ")
#parser.add_argument('-T','--temperature', type=int, default=300, help="Set simulation temperature (default = 300K)")
#parser.add_argument('--timestep', type=float, default=0.002, help="Set Production simulation time step in units of picosecond (default = 0.002 picosecond)")
#parser.add_argument('--interval', type=int, default=10, help="Set interval of saved frames in the unit of picosecond (default = 10 ps)")
#parser.add_argument('-l','--simulation_length', type=int, default=20, help="Set duration of Production simulation time in units of nanosecond (default = 20ns)")

#Plumed Force for Enhanced Sampling
parser.add_argument('--script', type = str, help="(Optional) Pass in a filepath adding custom force restraints to simulation system")

#Time log for computation time
parser.add_argument('--timeit', action="store_true", help="If activated, creates time_log file")
args = parser.parse_args()

# Load parameters
print("Loading parameters")
inputs = read_inputs(args.inpfile) # Needs omm_readinputs.py for this command

constraint_level=None
if args.constraint==1:
    constraint_level="HBonds"
elif args.constraint==2:
    constraint_level="AllBonds"
elif args.constraint==3:
    constraint_level="HAngles"

if args.periodic:
    periodicity = PME
else:
    periodicity = NoCutoff

if not args.implicit:
    solvent_type=None #0
else:
    if args.solvent_type==1:
        solvent_type=HCT
    if args.solvent_type==2:
        solvent_type=OBC1
    if args.solvent_type==3:
        solvent_type=OBC2
    if args.solvent_type==4:
        solvent_type=GBn
    if args.solvent_type==5:
        solvent_type=GBn2        

interval = int(args.interval / args.timestep)
simlength = int(args.simulation_length * 1000 / args.timestep)

script = """"""
if args.script != None :
    text = open(args.script, 'r')
    script = text.read()
    text.close()

### Actual code
try:
    top = read_top(args.topfile, args.fftype.upper()) # from omm_readparams.py
    crd = read_crd(args.crdfile, args.fftype.upper()) # from omm_readparams.py

    #Load input params based on FF type
    if args.fftype.upper() == 'CHARMM':
        params = read_params(args.toppar) # from omm_readparams.py
        top = read_box(top, args.sysinfo) if args.sysinfo else gen_box(top, crd) # from omm_readparams.py

    # Build system
    nboptions = dict(nonbondedMethod=inputs.coulomb,
                 nonbondedCutoff=inputs.r_off*nanometers,
                 constraints=inputs.cons,
                 ewaldErrorTolerance=inputs.ewald_Tol)

    if inputs.vdw == 'Switch': nboptions['switchDistance'] = inputs.r_on*nanometers
    if inputs.vdw == 'LJPME':  nboptions['nonbondedMethod'] = LJPME

    if inputs.implicitSolvent:
        nboptions['implicitSolvent'] = inputs.implicitSolvent
        nboptions['implicitSolventSaltConc'] = inputs.implicit_salt*(moles/liter)
        nboptions['temperature'] = inputs.temp*kelvin
        nboptions['soluteDielectric'] = inputs.solut_diele
        nboptions['solventDielectric'] = inputs.solve_diele
        nboptions['gbsaModel'] = inputs.gbsamodel

    if   args.fftype.upper() == 'CHARMM': system = top.createSystem(params, **nboptions)
    elif args.fftype.upper() == 'AMBER':  system = top.createSystem(**nboptions)

    if inputs.vdw == 'Force-switch': system = vfswitch(system, top, inputs) # from omm_vfswitch.py
    if inputs.lj_lrc == 'yes':
        for force in system.getForces():
            if isinstance(force, NonbondedForce): force.setUseDispersionCorrection(True)
            if isinstance(force, CustomNonbondedForce) and force.getNumTabulatedFunctions() != 1:
                force.setUseLongRangeCorrection(True)
    if inputs.e14scale != 1.0:
        for force in system.getForces():
            if isinstance(force, NonbondedForce): nonbonded = force; break
        for i in range(nonbonded.getNumExceptions()):
            atom1, atom2, chg, sig, eps = nonbonded.getExceptionParameters(i)
            nonbonded.setExceptionParameters(i, atom1, atom2, chg*inputs.e14scale, sig, eps)

    if inputs.pcouple == 'yes':      system = barostat(system, inputs) # from omm_barostat.py
    if inputs.rest == 'yes':         system = restraints(system, crd, inputs) # from omm_restraints.py
    integrator = LangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.dt*picoseconds)

    #Add custom forces for Enhanced Sampling with Plumed
    if script != """""":
        force = PlumedForce(script)
        system.addForce(force)

    # Set platform
    DEFAULT_PLATFORMS = 'CUDA', 'OpenCL', 'CPU'
    enabled_platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    if args.platform:
        if not args.platform[0] in enabled_platforms:
            print("Unable to find OpenMM platform '{}'; exiting".format(args.platform[0]), file=sys.stderr)
            sys.exit(1)

        platform = Platform.getPlatformByName(args.platform[0])
    else:
        for platform in DEFAULT_PLATFORMS:
            if platform in enabled_platforms:
                platform = Platform.getPlatformByName(platform)
                break
        if isinstance(platform, str):
            print("Unable to find any OpenMM platform; exiting".format(args.platform[0]), file=sys.stderr)
            sys.exit(1)

    print("Using platform:", platform.getName())
    prop = dict(CudaPrecision='single') if platform.getName() == 'CUDA' else dict()

    # Build simulation context
    simulation = Simulation(top.topology, system, integrator, platform, prop)
    simulation.context.setPositions(crd.positions)

    #Load from Checkpoint/Restart if specified
    if args.icrst:
        charmm_rst = read_charmm_rst(args.icrst) # from omm_readparams.py
        simulation.context.setPositions(charmm_rst.positions)
        simulation.context.setVelocities(charmm_rst.velocities)
        simulation.context.setPeriodicBoxVectors(charmm_rst.box[0], charmm_rst.box[1], charmm_rst.box[2])
    if args.irst:
        with open(args.irst, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))
    if args.ichk:
        with open(args.ichk, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())

    # Re-wrap
    if args.rewrap:
        simulation = rewrap(simulation) # from omm_rewrap.py

    #Makes the system log file
    simulation.reporters.append(StateDataReporter(args.output+'.log', interval, step=True,
            time = True, kineticEnergy = True, potentialEnergy=True, totalEnergy = True, temperature=True))

    # create time_log files
    if args.timeit:
        log_filename = CUR_DIR + args.output + ".time_log"
        file_handle = open(log_filename, 'w')
        file_handle.writelines("--- %s seconds ---" % (time.time() - start_time))
        file_handle.close()

    ########################################################################
                                # PRODUCTION STUFF#
    ########################################################################

    # Calculate initial system energy
    print("\nInitial system energy")
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    # Energy minimization
    if inputs.mini_nstep > 0:
        print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
        simulation.minimizeEnergy(tolerance=inputs.mini_Tol*kilojoule/mole, maxIterations=inputs.mini_nstep)
        print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    # Generate initial velocities
    if inputs.gen_vel == 'yes':
        print("\nGenerate initial velocities")
        if inputs.gen_seed:
            simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
        else:
            simulation.context.setVelocitiesToTemperature(inputs.gen_temp) 

    # Production
    if inputs.nstep > 0:
        print("\nMD run: %s steps" % inputs.nstep)
        if inputs.nstdcd > 0:
            if not args.odcd: args.odcd = 'output.dcd'
            simulation.reporters.append(DCDReporter(args.odcd, inputs.nstdcd))
        simulation.reporters.append(
            StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, progress=True,
                          remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
        )
        # Simulated annealing?
        if inputs.annealing == 'yes':
            interval = inputs.interval
            temp = inputs.temp_init
            for i in range(inputs.nstep):
                integrator.setTemperature(temp*kelvin)
                simulation.step(1)
                temp += interval
        else:
            simulation.step(inputs.nstep)

    # Write restart file
    if not (args.orst or args.ochk): args.orst = 'output.rst'
    if args.orst:
        state = simulation.context.getState( getPositions=True, getVelocities=True )
        with open(args.orst, 'w') as f:
            f.write(XmlSerializer.serialize(state))
    if args.ochk:
        with open(args.ochk, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
    if args.opdb:
        crd = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(top.topology, crd, open(args.opdb, 'w'))

except Exception as e:
    print(e)
    failfile = open("failed.txt", "a")
    failfile.writelines(args.output + '\n')
    failfile.close()
