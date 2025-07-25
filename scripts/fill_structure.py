from modeller import *
from modeller.automodel import *
from sys import argv

# Guide: https://salilab.org/modeller/wiki/Missing_residues

log.verbose()
env = environ()
env.io.atom_files_directory = ['modeler']

file_name = argv[1]
loop_ranges = argv[2].split(",")
loop_ranges = [(loop_range.split(":")[0], loop_range.split(":")[1]) for loop_range in loop_ranges]

class MyModel(LoopModel):
    def select_loop_atoms(self):
        selections = [Selection(self.residue_range(f"{start}:A", f"{end}:A")) for start, end in loop_ranges]
        return Selection(*selections)
    
    def select_atoms(self):
        return self.select_loop_atoms()

a = MyModel(env, alnfile = f"modeler/{file_name}.ali", knowns = f"{file_name}_incomplete", sequence = f"{file_name}_complete")

# Only needed for loop modeling
a.loop.starting_model = 1
a.loop.ending_model   = 10
a.loop.md_level       = refine.slow  # more accurate

a.starting_model = 1
a.ending_model   = 1

a.make()