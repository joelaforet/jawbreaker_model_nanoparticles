
import argparse
import sys
from openmm import *
from openmm.app import *
from openmm.unit import *


parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-s", "--sys", type=str, required=True, help="Input name of system to simulate")
parser.add_argument("--frames", type=int, default = 200, required=True, help="Frames to output")
parser.add_argument("--time", type=int, required=True, help="Input time to simulate for in nanoseconds.")
parser.add_argument("-r", "--run", type=int, required=True, help="Input replicate number for specific simulation.")
args = parser.parse_args()


system_Name = args.sys
frames_to_write = args.frames
simulation_time = args.time # nanoseconds

run_Number = args.run

# Input Files

prmtop = AmberPrmtopFile('{}.prmtop'.format(system_Name))
inpcrd = AmberInpcrdFile('{}.inpcrd'.format(system_Name))

# System Configuration

nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
ewaldErrorTolerance = 0.0005
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001

# Integration Options

dt = 0.002*picoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25

# Simulation Options


steps = (simulation_time * 1000)/ 0.002
#steps = 5000000
equilibrationSteps = 5000
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'single'}
#dcdReporter = DCDReporter('openmm_NPT_{}ns_{}.dcd'.format(str(simulation_time), system_Name), steps/frames_to_write)
pdbReporter = PDBReporter('openmm_NPT_{}ns_{}_run{}.pdb'.format(str(simulation_time), system_Name, str(run_Number)), steps/frames_to_write)
dataReporter = StateDataReporter('log_NPT_{}ns_{}_run{}.txt'.format(str(simulation_time),system_Name, str(run_Number)), 500, totalSteps=steps,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')
checkpointReporter = CheckpointReporter('NPT_{}ns_{}_run{}.chk'.format(str(simulation_time),system_Name, str(run_Number)), 1000000)

# Prepare the Simulation

print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions
system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Minimize and Equilibrate

print('Performing energy minimization...')
simulation.minimizeEnergy()
print('Equilibrating...')
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.reporters.append(pdbReporter)
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# Simulate

print('Simulating...')
#simulation.reporters.append(dcdReporter)
#simulation.reporters.append(dataReporter)
#simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
simulation.step(steps)

print('Simulation Complete!')
# Write file with final simulation state

state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())
