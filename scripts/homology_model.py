from modeller import *
from modeller.automodel import *
from sys import argv

# Guide: https://salilab.org/modeller/wiki/Missing_residues

log.verbose()
env = environ()
env.io.atom_files_directory = ['modeler']

file_name = argv[1]
templates = argv[2].split(",")
helices = argv[3] if len(argv) > 3 else None


class MyModel(automodel):
    def special_restraints(self, aln):
        if helices is not None:
            rsr = self.restraints

            for helix_range in helices.split(","):
                start, end = helix_range.split(":")
                start_res = f"{start}:A"
                end_res = f"{end}:A"
                rsr.add(secondary_structure.alpha(self.residue_range(start_res, end_res)))


a = MyModel(env,
              alnfile=f"modeler/{file_name}.ali",
              knowns=templates,
              sequence="target")

a.starting_model = 1
a.ending_model   = 5

a.make()