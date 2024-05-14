"""
Langevin Dynamics in OpenMM
By Joe Laforet Jr.
jrl78
"""


import simtk.openmm as openmm
import simtk.unit as unit

import numpy as np
from scipy.spatial.distance import cdist

from ipywidgets import IntProgress
from IPython.display import display
import time

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

import os
import argparse

def sample_points_inside_cube(num_points, size, min_distance):
    """
    Sample points with minimum distance inside a cube.

    Parameters:
        num_points (int): Number of points to sample.
        size (float): Side length of the cube.
        min_distance (float): Minimum distance between points.

    Returns:
        points (list of tuples): List of sampled points.
    """
    points = []
    while len(points) < num_points:
        new_point = np.random.rand(3) * size - size / 2
        if not any(np.linalg.norm(new_point - point) < min_distance for point in points):
            points.append(new_point)

    return points

parser = argparse.ArgumentParser()
parser.add_argument("-num1", help="Number of first particle type.", required=True)
parser.add_argument("-num2", help="Number of second particle type.", required=True)

parser.add_argument("-sig1", help="Sigma of first particle type. (angstroms)", required=True)
parser.add_argument("-sig2", help="Sigma of second particle type. (angstroms)", required=True)

parser.add_argument("-eps1", help="Epsilon of first particle type. (kcal/mol)", required=True)
parser.add_argument("-eps2", help="Epsilon of second particle type. (kcal/mol)", required=True)

parser.add_argument("-eps12", help="Epsilon value for interaction between type 1 and type 2 particles. (kcal/mol)", required = True)

parser.add_argument("-size", help="Side length of cubic box (nm). Default = 12nm", default=12)
parser.add_argument("-temp", help="Temperature of system. Default = 300K", default=300)
parser.add_argument("-time", help="Time to run simulation for. Default = 10ns", default=10)
parser.add_argument("-cutoff", help="Cutoff distance for force calculations. Default = 12 angstroms", default=12)

parser.add_argument("-npt", action="store_true", help="If activated, perform a simulation using NPT ensemble.")

args = parser.parse_args()


output = f"A_{args.num1}_B_{args.num2}_vol={args.size}_{args.time}ns_e1={args.eps1}_e2={args.eps2}_e12={args.eps12}_s1={args.sig1}_s2={args.sig2}"

if args.npt:
    output = 'NPT_'+output
else:
    output = 'NVT_'+output

# Define system parameters
num_particles_type1 = args.num1
num_particles_type2 = args.num2
box_size =args.size
temperature = args.temp * unit.kelvin
friction_coefficient = 1.0 / unit.picosecond
timestep = 1.0 * unit.femtosecond
simulation_time = args.time * unit.nanosecond  # in picoseconds

# For NPT
pressure = 1.0*unit.atmospheres
barostatInterval = 25

# Lennard-Jones parameters for particle types
lj_params = {
    'type1': {'epsilon': args.eps1 * unit.kilocalories_per_mole, 'sigma': args.sig1 * unit.angstroms},
    'type2': {'epsilon': args.eps2 * unit.kilocalories_per_mole, 'sigma': args.sig2 * unit.angstroms}
}

# make eps_12 an arbitrary

epsilon_12 = args.eps12 * unit.kilocalories_per_mole

# Create an OpenMM System
system = openmm.System()

# Add particles to the system with given mass
for _ in range(num_particles_type1):
    system.addParticle(39.948 * unit.amu)  # Argon mass

for _ in range(num_particles_type2):
    system.addParticle(39.948 * unit.amu)  # Argon mass

# Set up the simulation box
system.setDefaultPeriodicBoxVectors(
    unit.Quantity((box_size, 0, 0), unit.nanometer),
    unit.Quantity((0, box_size, 0), unit.nanometer),
    unit.Quantity((0, 0, box_size), unit.nanometer)
)


# Define the Lennard-Jones potential and add it to the system
lj_force = openmm.NonbondedForce()
lj_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
lj_force.setCutoffDistance(args.cutoff * unit.angstroms)

# Add particle interactions to the Lennard-Jones force
for particle_index in range(num_particles_type1):
    particle_1 = lj_force.addParticle(charge=0.0, sigma=lj_params['type1']['sigma'], epsilon = lj_params['type1']['epsilon'])

for particle_index in range(num_particles_type2):
    particle_2 = lj_force.addParticle(charge=0.0, sigma=lj_params['type2']['sigma'], epsilon = lj_params['type2']['epsilon'])

particle_12_Exception = lj_force.addException(particle1=particle_1, particle2=particle_2, chargeProd = 0.0, sigma= 0.5 * (lj_params['type2']['sigma'] +lj_params['type1']['sigma']),
                                              epsilon = epsilon_12, replace = True)

system.addForce(lj_force)

"""NPT HERE"""
# Add Barostat to sample NPT
if args.npt:
    system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))

# Create a Topology object representing the particles
num_particles = num_particles_type1 + num_particles_type2
topology = openmm.app.Topology()
chain = topology.addChain()


for i in range(num_particles_type1):
    residue = topology.addResidue(f"1", chain)
    atom = topology.addAtom(f"1", None, residue)

for i in range(num_particles_type1, num_particles):
    residue = topology.addResidue(f"2", chain)
    atom = topology.addAtom(f"2", None, residue)


# Create a Langevin integrator
platform = openmm.Platform.getPlatformByName('CUDA')  # Change to 'CUDA' for GPU acceleration if available
integrator = openmm.LangevinIntegrator(temperature, friction_coefficient, timestep)
simulation = openmm.app.Simulation(topology, system, integrator, platform=platform)

# Create a Context and set initial particle positions

#context = openmm.Context(system, integrator, platform)


"""BAD WAY OF INITIALIZING, MAKE LATTICE"""

# Randomly distribute particles within the box for each type
initial_positions = sample_points_inside_cube(num_particles, box_size, min_distance=.5)

# Set initial particle positions
simulation.context.setPositions(initial_positions)

# Set up a PDB file reporter to write particle positions
pdb_reporter = openmm.app.PDBReporter(f'{output}.pdb', 1000)  # Writes every 100 steps

simulation.reporters.append(openmm.app.StateDataReporter(output+'.log', 100, step=True,
            time = True, kineticEnergy = True, potentialEnergy=True, totalEnergy = True, temperature=True))

# Run Langevin dynamics simulation
simulation_steps = int(simulation_time / timestep)

for step in range(simulation_steps):
    simulation.step(1)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    position = state.getPositions()
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    #volume = state.getPeriodicBoxVolume()
    

    if (step + 1) % 1000 == 0:  # Write positions to PDB every 100 steps
        pdb_reporter.report(simulation = simulation, state=state)


# Clean up
pdb_reporter.close()

