from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, argv
import os

input_file = argv[1]
steps = int(argv[2])
output_interval = int(argv[3])

input_path = os.path.join("inputs", input_file)

pdb = PDBFile(input_path)
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

# Nalu
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield,
                    model='tip3p',         # Water model (must match forcefield)
                    padding=1.0*nanometers, # Box padding on all sides
                    neutralize=True)        # Add counterions if net charge exists

system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.reporters.append(PDBReporter('output.pdb', output_interval))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(input_path)