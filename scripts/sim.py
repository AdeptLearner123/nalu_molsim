from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import os
from utils import *
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_pdb', type=str, help='Input PDB file')
parser.add_argument('output_pdb', type=str, help='Output PDB file')
parser.add_argument('nsteps', type=int, help='Number of MD steps')
parser.add_argument('report_interval', type=int, help='Reporting interval')
parser.add_argument('-f', '--freeze_ranges', type=str, default=None)
parser.add_argument('-m', '--minimize', action='store_true')
parser.add_argument('-i', '--implicit-solvent', action='store_true')
args = parser.parse_args()

input_file = args.input_pdb
output_file = args.output_pdb
steps = args.nsteps
output_interval = args.report_interval
minimize_first = args.minimize
implicit_solvent = args.implicit_solvent
freeze_ranges = args.freeze_ranges

input_path = os.path.join("inputs", input_file)
output_path = os.path.join("outputs", output_file)

pdb = PDBFile(input_path)
platform = get_best_platform()
modeller = Modeller(pdb.topology, pdb.positions)

if implicit_solvent:
    print("Using implicit solvent")
    forcefield = ForceField('amber99sb.xml', 'amber99_obc.xml')
    modeller.addHydrogens(forcefield)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
else:
    print("Using real solvent")
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield,
                        model='tip3p',         # Water model (must match forcefield)
                        padding=1.0*nanometers, # Box padding on all sides
                        neutralize=True)        # Add counterions if net charge exists
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)


if freeze_ranges is not None:
    restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    system.addForce(restraint)
    restraint.addGlobalParameter('k', 100.0*kilojoules_per_mole/nanometer)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    freeze_ranges = freeze_ranges.split(",")
    freeze_ranges = [value.split(":") for value in freeze_ranges]
    freeze_ranges = [(chain, int(start), int(end)) for chain, start, end in freeze_ranges]
    freeze_chains = sum([[chain] * (end - start + 1) for chain, start, end in freeze_ranges])
    freeze_idxs = sum([list(range(start, end + 1)) for _, start, end in freeze_ranges])
    freeze_idxs = list(zip(freeze_chains, freeze_idxs))


    # Add the atoms you want to freeze (e.g., residue 10)
    for atom in pdb.topology.atoms():
        res = atom.residue
        chain_id = res.chain.id
        res_index = res.index  # 0-based index

        if (chain_id, res_index) in freeze_idxs:
            print("Freezing: ", res_index, chain_id)
            #system.setParticleMass(atom.index, 0.0*dalton)
            if atom.name == 'CA':
                restraint.addParticle(atom.index, pdb.positions[atom.index])


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