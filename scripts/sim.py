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

input_path = os.path.join("inputs", input_file)
output_path = os.path.join("outputs", output_file)

pdb = PDBFile(input_path)
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

# Nalu
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield,
                    model='tip3p',         # Water model (must match forcefield)
                    padding=1.0*nanometers, # Box padding on all sides
                    neutralize=True)        # Add counterions if net charge exists

platform = get_best_platform()
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.reporters.append(PDBReporter(output_path, output_interval))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

start = time.time()
simulation.step(steps)
end = time.time()
print(f"Execution time: {end - start:.2f} seconds")