from modeller import *
from modeller.automodel import *

# Guide: https://salilab.org/modeller/wiki/Missing_residues

log.verbose()
env = environ()

# Load alignment and atom files
env.io.atom_files_directory = ['drafts']

class MyModel(LoopModel):
    def select_loop_atoms(self):
        return Selection(self.residue_range('31:A', '44:A'))
    
    def select_atoms(self):
        return Selection(self.residue_range('29:A', '46:A'))

    def special_restraints(self, aln):
        rsr = self.restraints

        atom1 = self.atoms['N:31:A']  # Atom name : Residue number : Chain ID
        atom2 = self.atoms['O:44:A']  # Opposite strand atom (anti-parallel -> reverse order)

        sheet_restraint = secondary_structure.Sheet(atom1, atom2, sheet_h_bonds=-4)  # -4 = 4 anti-parallel H-bonds
        rsr.add(sheet_restraint)

        # Example: enforce a beta strand on residues 45–49 and 54–58 (adjust as needed)
        #rsr.add(secondary_structure.strand(self.residue_range('45:A', '49:A')))
        #rsr.add(secondary_structure.strand(self.residue_range('54:A', '58:A')))

        # Or a beta-hairpin:
        #rsr.add(secondary_structure.Strand(self.residue_range('1:A', '35:A')))
        #rsr.add(secondary_structure.Strand(self.residue_range('40:A', '74:A')))

        # Or force a helix:
        # rsr.add(secondary_structure.alpha(self.residue_range('43:A', '60:A')))

# Initialize the model
a = MyModel(env, alnfile = 'drafts/hairpin_model.ali', knowns = 'hairpin_incomplete', sequence = 'hairpin_complete')

# Only needed for loop modeling
a.loop.starting_model = 1
a.loop.ending_model   = 10
a.loop.md_level       = refine.slow  # more accurate

a.starting_model = 1
a.ending_model   = 1

a.make()