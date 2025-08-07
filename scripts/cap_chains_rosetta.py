from pyrosetta import *
from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue
from pyrosetta.rosetta.core.chemical import VariantType
from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue, remove_variant_type_from_pose_residue

from sys import argv

# === 1. Initialize PyRosetta ===
init("-ignore_unrecognized_res true")

# === 2. Load the PDB ===
file = argv[1]
pose = pose_from_pdb(f"pyrosetta/{file}.pdb")

# === 3. Get chain information ===
chain_to_residues = {}
for i in range(1, pose.total_residue() + 1):
    chain = pose.pdb_info().chain(i)
    if chain not in chain_to_residues:
        chain_to_residues[chain] = []
    chain_to_residues[chain].append(i)

# === 4. Add ACE/NME to each chain ===
for chain, res_indices in chain_to_residues.items():
    n_term = res_indices[0]
    c_term = res_indices[-1]

    print(f"Capping chain {chain}: N-term = {n_term}, C-term = {c_term}")

    # === 1. Remove conflicting terminal variants ===
    #remove_variant_type_from_pose_residue(pose, VariantType.name_to_variant_type("NtermProteinFull"), n_term)
    #remove_variant_type_from_pose_residue(pose, VariantType.name_to_variant_type("CtermProteinFull"), c_term)

    # === 2. Add ACE/NME caps ===
    add_variant_type_to_pose_residue(pose, VariantType.N_ACETYLATION, n_term)
    add_variant_type_to_pose_residue(pose, VariantType.C_METHYLAMIDATION, c_term)

# === 5. Save the capped structure ===
pose.dump_pdb("pyrosetta/{file}_capped.pdb")