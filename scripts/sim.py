from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, argv
import os
from utils import *
import time

input_file = argv[1]
output_file = argv[2]
steps = int(argv[3])
output_interval = int(argv[4])
minimize_first = len(argv) > 5 and argv[5] == "y"
implicit_solvent = len(argv) > 6 and argv[6] == "y"

input_path = os.path.join("inputs", input_file)
output_path = os.path.join("outputs", output_file)

pdb = PDBFile(input_path)
platform = get_best_platform()
modeller = Modeller(pdb.topology, pdb.positions)

if implicit_solvent:
    print("Using implicit solvent")
    forcefield = ForceField('amber99sb.xml', 'amber99_obc.xml')
    modeller.addHydrogens(forcefield)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
else:
    print("Using real solvent")
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield,
                        model='tip3p',         # Water model (must match forcefield)
                        padding=1.0*nanometers, # Box padding on all sides
                        neutralize=True)        # Add counterions if net charge exists
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.reporters.append(PDBReporter(output_path, output_interval))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

if minimize_first:
    # Run energy minimization
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    print("Minimization complete.")

start = time.time()
# Write time=0 structure
PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open(output_path, 'w'))
simulation.step(steps)
end = time.time()
print(f"Execution time: {end - start:.2f} seconds")