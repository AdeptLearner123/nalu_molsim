from modeller import *
from modeller.automodel import *
import argparse

# Guide: https://salilab.org/modeller/wiki/Missing_residues

parser = argparse.ArgumentParser()
parser.add_argument('alignment', type=str, help='Input PDB file')
parser.add_argument('templates', type=str, help='Output PDB file')
parser.add_argument('-a', '--alpha-helices', type=str, default=None)
parser.add_argument('-b', '--beta-sheets', type=str, default=None)
parser.add_argument('-l', '--loops', type=str, default=None)
args = parser.parse_args()

log.verbose()
env = environ()
env.io.atom_files_directory = ['modeler']

file_name = args.alignment
templates = args.templates
helices = args.alpha_helices
sheets = args.beta_sheets
loops = args.loops


class MyAutoModel(AutoModel):
    def special_restraints(self, aln):
        rsr = self.restraints

        if helices is not None:
            for helix_range in helices.split(","):
                chain, start, end = helix_range.split(":")
                start_res = f"{start}:{chain}"
                end_res = f"{end}:{chain}"
                rsr.add(secondary_structure.alpha(self.residue_range(start_res, end_res)))
        
        if sheets is not None:
            for sheet_range in sheets.split(","):
                chain, start, end, count = sheet_range.split(":")
                atom1 = self.atoms[f'N:{start}:{chain}']  # Atom name : Residue number : Chain ID
                atom2 = self.atoms[f'O:{end}:{chain}']  # Opposite strand atom (anti-parallel -> reverse order)
                sheet_restraint = secondary_structure.Sheet(atom1, atom2, sheet_h_bonds=int(count))  # -4 = 4 anti-parallel H-bonds
                rsr.add(sheet_restraint)


class MyLoopModel(LoopModel):
    def _get_loop_atoms(self):
        loops_list = loops.split(",")
        loop_ranges = [loop.split(":") for loop in loops_list]
        selections = [Selection(self.residue_range(f'{start}:{chain}', f'{end}:{chain}')) for (chain, start, end) in loop_ranges]
        return Selection(selections)

    def select_loop_atoms(self):
        if loops is not None:
            return self._get_loop_atoms()
        else:
            return Selection(self.atoms) # dummy selection

    def select_atoms(self):
        if loops is not None:
            return self._get_loop_atoms()
        else:
            return Selection(self.atoms) # dummy selection

use_loop_model = loops is not None
ModelClass = MyLoopModel if use_loop_model else MyAutoModel

a = ModelClass(env,
              alnfile=f"modeler/{file_name}.ali",
              knowns=templates,
              sequence="target")

if loops is not None:
    # Only needed for loop modeling
    a.loop.starting_model = 1
    a.loop.ending_model   = 5
    a.loop.md_level       = refine.slow  # more accurate
    a.starting_model = 1
    a.ending_model   = 1
else:
    a.starting_model = 1
    a.ending_model   = 5

a.make()