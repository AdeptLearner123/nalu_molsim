from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import os
from utils import *
import time
import argparse

def restrain_residues(system,
                      modeller,                 # use modeller.topology & modeller.positions
                      freeze_keys,              # {('A','10'), ('B','25')} or {('A',10), ...}
                      k=1e4 * kilojoule_per_mole / nanometer**2,
                      ca_only=False,
                      use_pbc=True,
                      force_name="PositionalRestraint"):
    """
    Add harmonic positional restraints for atoms in residues matching freeze_keys,
    using coordinates from modeller.positions.
    """
    top = modeller.topology
    pos = modeller.positions  # Quantity (n_atoms, 3)
    if pos is None:
        raise ValueError("modeller.positions is None; set positions before calling.")

    n_top = top.getNumAtoms()
    if len(pos) != n_top or system.getNumParticles() != n_top:
        raise ValueError(
            f"Atom count mismatch: top={n_top}, pos={len(pos)}, system={system.getNumParticles()}. "
            "Ensure the System was created from modeller.topology before calling."
        )

    expr = ("k*periodicdistance(x,y,z,x0,y0,z0)^2"
            if use_pbc else
            "k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")

    # Reuse an existing named restraint if present, else create
    restraint = None
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        try:
            name = f.getName()
        except Exception:
            name = ""
        if isinstance(f, CustomExternalForce) and name == force_name:
            restraint = f
            break
    if restraint is None:
        restraint = CustomExternalForce(expr)
        try:
            restraint.setName(force_name)
        except Exception:
            pass
        restraint.addGlobalParameter("k", k)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        system.addForce(restraint)

    def in_freeze_set(chain_id, res):
        if (chain_id, res.id) in freeze_keys:
            return True
        try:
            return (chain_id, int(res.id)) in freeze_keys
        except ValueError:
            return False

    n_added = 0
    for atom in top.atoms():
        res = atom.residue
        if not in_freeze_set(res.chain.id, res):
            continue
        if ca_only and atom.name != "CA":
            continue
        i = atom.index
        p = pos[i]  # Quantity Vec3
        restraint.addParticle(i, [p[0], p[1], p[2]])
        n_added += 1

    print(f"[restrain_residues] Added restraints to {n_added} atom(s) "
          f"({'CA-only' if ca_only else 'all atoms'}), k={k}.")

    return restraint


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
print("Before:", [c.id for c in pdb.topology.chains()])
modeller = Modeller(pdb.topology, pdb.positions)
print("After:", [c.id for c in modeller.topology.chains()])

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

if freeze_ranges is not None:    
    freeze_ranges = freeze_ranges.split(",")
    freeze_ranges = [value.split(":") for value in freeze_ranges]
    freeze_ranges = [(chain, int(start), int(end)) for chain, start, end in freeze_ranges]
    freeze_chains = sum([[chain] * (end - start + 1) for chain, start, end in freeze_ranges])
    freeze_idxs = sum([list(range(start, end + 1)) for _, start, end in freeze_ranges])
    freeze_idxs = list(zip(freeze_chains, freeze_idxs))
    restrain_residues(system, modeller, freeze_idxs)
    simulation.context.reinitialize(preserveState=True)

start = time.time()
# Write time=0 structure
PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open(output_path, 'w'))
simulation.step(steps)
end = time.time()
print(f"Execution time: {end - start:.2f} seconds")