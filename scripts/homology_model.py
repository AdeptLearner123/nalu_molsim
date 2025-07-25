from modeller import *
from modeller.automodel import *
from sys import argv

# Guide: https://salilab.org/modeller/wiki/Missing_residues

log.verbose()
env = environ()
env.io.atom_files_directory = ['modeler']

file_name = argv[1]
templates = argv[2].split(",")

a = automodel(env,
              alnfile=f"modeler/{file_name}.ali",
              knowns=templates,
              sequence="target")

a.starting_model = 1
a.ending_model   = 5

a.make()