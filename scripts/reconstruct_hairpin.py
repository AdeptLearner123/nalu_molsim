from pyrosetta import *
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.chemical import ChemicalManager

# === INIT PYROSETTA ===
init("-mute all")

# === USER INPUTS ===
input_pdb = "drafts/hairpin_incomplete.pdb"   # Replace with your file
output_pdb = "drafts/hairpin_p1.pdb"
missing_seq = "DVSLAFSEISVGAE"    # Human C9 CMH2 missing β-hairpin
insert_after_res = 30           # Where to insert the hairpin

# === LOAD INCOMPLETE C9 STRUCTURE ===
pose = pose_from_pdb(input_pdb)

# === SPLIT POSE AT INSERTION POINT ===
pose_pre = Pose()
pose_pre.assign(pose)
pose_pre.delete_residue_range_slow(insert_after_res + 1, pose.total_residue())

pose_post = Pose()
pose_post.assign(pose)
pose_post.delete_residue_range_slow(1, insert_after_res)

# === COMBINE PIECES ===
combined = Pose()
combined.assign(pose_pre)

def create_residue(code):
    residue_set = ChemicalManager.get_instance().residue_type_set('fa_standard')
    res_type = residue_set.name_map(code)
    return ResidueFactory.create_residue(res_type)

insert_start = combined.total_residue()
residues = ["ASP", "VAL", "SER", "LEU", "ALA", "PHE", "SER", "GLU", "ILE", "SER", "VAL", "GLY", "ALA", "GLU"]
angles = [(-135, 135), (-135, 135), (-135, 135), (-135, 135), (-135, 135), (-135, 135), (-60, -30), (-90, 0),
          (-110, -90), (90, 20), (-135, 135), (-135, 135), (-135, 135), (-135, 135), (-135, 135)]

for resname in residues:
    combined.append_residue_by_bond(create_residue(resname), build_ideal_geometry=True)

combined.set_omega(insert_start, 180.0)
combined.set_psi(insert_start, 110)

for i in range(len(residues)):
    phi, psi = angles[i]
    print(f"Setting angle of {insert_start + i + 1}: {phi}, {psi}")
    combined.set_phi(insert_start + i + 1, phi)
    combined.set_psi(insert_start + i + 1, psi)
    combined.set_omega(insert_start + i + 1, 180.0)

# === SAVE OUTPUT ===
combined.dump_pdb(output_pdb)
pose_post.dump_pdb("drafts/hairpin_p2.pdb")
print(f"✅ Hairpin inserted and relaxed. Saved as: {output_pdb}")