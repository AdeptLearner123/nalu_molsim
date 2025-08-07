from modeller import environ, model
from modeller.automodel import automodel
from sys import argv

filename = argv[1]

# Initialize the environment
env = environ()
env.io.atom_files_directory = ['capper']  # current directory

# https://newreef.csb.vanderbilt.edu/comp/soft/modeller/manual/node16.html
env.patch_default = False

# Define a subclass of automodel to customize termini
class CappedModel(automodel):
    def special_patches(self, aln):
        from collections import defaultdict
        chains = defaultdict(list)

        for r in self.residues:
            chains[r.chain].append(r)

        for chain_id, residues in chains.items():
            start = residues[0]
            end = residues[-1]
            print(f"Patching chain {chain_id}: start = {start.num}, end = {end.num}")
            try:
                self.patch(residue_type='ACE', residues=(start,))
                self.patch(residue_type='CT3', residues=(end,))
            except Exception as e:
                print(f"⚠️ Could not patch chain {chain_id} ({start.num}–{end.num}): {e}")

# Create alignment
capped_file = f'capper/{filename}_capped.ali'

# Build the model
a = CappedModel(env,
                alnfile=capped_file,
                knowns=filename,
                sequence=filename,
                assess_methods=())  # no scoring needed

a.starting_model = 1
a.ending_model = 1

a.make()