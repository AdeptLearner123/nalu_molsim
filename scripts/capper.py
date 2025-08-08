import os
import sys
from modeller import *
from modeller.automodel import *

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    # Include unusual residues as needed
}


class NoMoveModel(AutoModel):
    def select_atoms(self):
        last_residue = self.residues[-1]
        return Selection(last_residue)


def extract_sequence(pdb_file):
    sequence = []
    seen_residues = set()

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15] == "CA":
                resname = line[17:20].strip()
                chain = line[21].strip()
                resnum = line[22:26].strip()
                unique_id = (chain, resnum)

                if unique_id in seen_residues:
                    continue

                seen_residues.add(unique_id)
                one_letter = three_to_one.get(resname, 'X')  # X for unknown
                sequence.append(one_letter)

    return ''.join(sequence)


def write_alignment_file(seq, template_name, target_name):
    """Write PIR-format alignment file with ACE/NME caps."""
    ali_text = f""">P1;{template_name}
structure:{template_name}:FIRST:@:LAST:@::::
{seq}*

>P1;{target_name}
sequence:{target_name}:::::::0.00: 0.00
{seq}*
"""
    with open(f"capper/{template_name}.ali", "w") as f:
        f.write(ali_text)


def run_modeller(template_name, target_name):
    """Run Modeller automodel using the generated alignment."""
    env = environ()
    env.io.atom_files_directory = ['.', './capper']
    a = NoMoveModel(env,
                  alnfile=f'capper/{template_name}.ali',
                  knowns=template_name,
                  sequence=target_name)
    a.starting_model = 1
    a.ending_model = 1
    a.make()


def main(pdb_file):
    pdb_path = os.path.join("capper", pdb_file)
    sequence = extract_sequence(pdb_path)
    template_name = pdb_file.replace(".pdb", "")
    target_name = template_name + "_capped"

    print("[+] Writing alignment file...")
    write_alignment_file(sequence, template_name, target_name)

    print("[+] Running Modeller...")
    run_modeller(template_name, target_name)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python add_ace_nme_modeller.py input.pdb")
        sys.exit(1)
    pdb_file = sys.argv[1]
    main(pdb_file)