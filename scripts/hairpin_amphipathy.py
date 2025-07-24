import matplotlib.pyplot as plt
import numpy as np

def moving_average_edge_friendly(values, window):
    n = len(values)
    result = []
    for i in range(n):
        # Determine the window bounds around index i
        start = max(0, i - window // 2)
        end = min(n, i + window // 2 + 1)
        avg = sum(values[start:end]) / (end - start)
        result.append(avg)
    return result


def get_stacked_bar2(bar1, bar2):
    stacked_mask = (bar1 >= 0) & (bar2 >= 0) | (bar1 <= 0) & (bar2 <= 0)
    bar2_stacked = np.where(stacked_mask, bar2, 0)
    bar2_unstacked = np.where(~stacked_mask, bar2, 0)
    return bar2_stacked, bar2_unstacked


hydropathy_scores = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# K is pointing into the pore
#STRAND1 = "KASMKRKGVELKDIKRCLGYHLDVSLAFS"
#STRAND2 = "LNESTINVARGEGRKVCDDKNFEAGVSIE"
STRAND1 = "KASMKEKGVDLNDVKHCLGFNMDLRIPLQDD"
STRAND2 = "INDRTIDVTKGNDTKICGDANVSATVSADKL"

residue_labels = [f"{res1}\n{res2}" for res1, res2 in zip(STRAND1, STRAND2)]
x = np.arange(len(STRAND1))

strand1_hydropathy = np.array([hydropathy_scores[res] for res in STRAND1])
strand2_hydropathy = np.array([hydropathy_scores[res] for res in STRAND2])
strand2_hydropathy_stacked, strand2_hydropathy_unstacked = get_stacked_bar2(strand1_hydropathy, strand2_hydropathy)

colors1 = ["red" if i % 2 == 0 else "blue" for i in x]
colors2 = ["pink" if i % 2 == 0 else "cyan" for i in x]

total_hydropathy = [h1 + h2 for h1, h2 in zip(strand1_hydropathy, strand2_hydropathy)]
amphipathy = [-h if i % 2 == 0 else h for i, h in enumerate(total_hydropathy)]
amphipathy_smoothed = moving_average_edge_friendly(amphipathy, 7)

plt.figure(figsize=(len(STRAND1) / 2, 4))  # Wider figure for readability
plt.bar(x, strand1_hydropathy, label="Strand 1", color=colors1)
plt.bar(x, strand2_hydropathy_stacked, bottom=strand1_hydropathy, label="Strand 2", color=colors2)
plt.bar(x, strand2_hydropathy_unstacked, label="Strand 2", color=colors2)
plt.plot(x, amphipathy_smoothed)
plt.title("Hydropathy Plot (Kyte-Doolittle)")
plt.xlabel("Residue")
plt.ylabel("Hydropathy")
plt.axhline(0, color='gray', linestyle='--')
plt.xticks(ticks=x, labels=residue_labels)
plt.tight_layout()
plt.show()


# Net side-chain charges at pH ~7.4
amino_acid_charges = {
    'A':  0,   # Alanine
    'R': +1,   # Arginine
    'N':  0,   # Asparagine
    'D': -1,   # Aspartic Acid
    'C':  0,   # Cysteine (neutral at pH 7.4)
    'E': -1,   # Glutamic Acid
    'Q':  0,   # Glutamine
    'G':  0,   # Glycine
    'H': +0.1, # Histidine (partially protonated)
    'I':  0,   # Isoleucine
    'L':  0,   # Leucine
    'K': +1,   # Lysine
    'M':  0,   # Methionine
    'F':  0,   # Phenylalanine
    'P':  0,   # Proline
    'S':  0,   # Serine
    'T':  0,   # Threonine
    'W':  0,   # Tryptophan
    'Y':  0,   # Tyrosine (neutral at pH 7.4)
    'V':  0    # Valine
}

aromatic_residues = set(['F', 'Y', "W"])

strand1_charge = np.array([amino_acid_charges[res] for res in STRAND1])
strand2_charge = np.array([amino_acid_charges[res] for res in STRAND2])
strand2_charge_stacked, strand2_charge_unstacked = get_stacked_bar2(strand1_charge, strand2_charge)
total_charge = [c1 + c2 for c1, c2 in zip(strand1_charge, strand2_charge)]
charge_smoothed = moving_average_edge_friendly(total_charge, 3)

aromatic_bars1 = np.array([0.5 if res in aromatic_residues else 0 for res in STRAND1])
aromatic_bars2 = np.array([0.5 if res in aromatic_residues else 0 for res in STRAND2])
aromatic_totals = aromatic_bars1 + aromatic_bars2

plt.figure(figsize=(len(STRAND1) / 2, 4))  # Wider figure for readability
plt.bar(x, strand1_charge, label="Strand 1", color="blue")
plt.bar(x, strand2_charge_stacked, bottom=strand1_charge, label="Strand 2", color="cyan")
plt.bar(x, strand2_charge_unstacked, label="Strand 2", color="cyan")
plt.bar(x, aromatic_totals, color="yellow")
plt.plot(x, charge_smoothed)
plt.title("Charge Plot")
plt.xlabel("Residue")
plt.ylabel("Net Charge")
plt.axhline(0, color='gray', linestyle='--')
plt.xticks(ticks=x, labels=residue_labels)
plt.tight_layout()
plt.show()