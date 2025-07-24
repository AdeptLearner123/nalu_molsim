from openmm.app import *
from openmm import *
from openmm.unit import *
import os
from sys import argv

# Input and output file names
input_pdb = argv[1]
output_pdb = argv[2]
input_path = os.path.join("inputs", input_pdb)
output_path = os.path.join("outputs", output_pdb)
                           
# Load the structure)
pdb = PDBFile(input_path)
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

# Add hydrogens (and solvent if needed)
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

# Create system with constraints on hydrogen bonds only
system = forcefield.createSystem(modeller.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer,
                                 constraints=HBonds)

# Set up integrator (not used for minimization, but required by Simulation)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Create simulation
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Run energy minimization
print("Minimizing energy...")
simulation.minimizeEnergy()
print("Minimization complete.")

# Save minimized structure
positions = simulation.context.getState(getPositions=True).getPositions()
with open(output_path, 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)

print(f"Minimized structure saved to {output_pdb}")