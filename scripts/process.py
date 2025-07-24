from openmm.app import PDBFile, Modeller, ForceField
from openmm.unit import *
from sys import argv
import os
import subprocess

FOLDX5_PATH_DIR = "/Users/naluzou/Projects/MolSim/foldx5_1Mac/"
FOLDX5_PATH = "/Users/naluzou/Projects/MolSim/foldx5_1Mac/foldx_20251231"

input_file = argv[1]
input_path = os.path.join("drafts", input_file)

print("Entering")
subprocess.run(["cp", input_path, "./"])
result = subprocess.run([FOLDX5_PATH, "--command=RepairPDB", f"--pdb={input_file}"], capture_output=True, text=True)
print("STDOUT:", result.stdout)
print("STDERR:", result.stderr)

fixed_file = input_file.removesuffix('.pdb') + "_Repair.pdb"
# Load repaired PDB
pdb = PDBFile(fixed_file)

# Load force field (adjust if you're using a different one)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create modeller object
modeller = Modeller(pdb.topology, pdb.positions)

# Add hydrogens at physiological pH (~7.0)
modeller.addHydrogens(forcefield, pH=7.0)

# Save to output file
output_path = os.path.join("inputs", input_file)
with open(output_path, 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

subprocess.run(["rm", input_file])
subprocess.run(["rm", fixed_file])
subprocess.run(["rm", fixed_file.removesuffix(".pdb") + ".fxout"])