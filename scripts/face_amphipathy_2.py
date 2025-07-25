import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import math

def moving_moment_edge_friendly(values, angles, window):
    coss = [np.cos(np.deg2rad(angle)) for angle in angles]
    sins = [np.sin(np.deg2rad(angle)) for angle in angles]

    n = len(values)
    result = []
    for i in range(n):
        # Determine the window bounds around index i
        start = max(0, i - window // 2)
        end = min(n, i + window // 2 + 1)

        window_vals = values[start:end]
        window_cos = coss[start:end]
        window_sin = sins[start:end]
        hx = sum([score * cos for score, cos in zip(window_vals, window_cos)]) / len(window_vals)
        hy = sum([score * sin for score, sin in zip(window_vals, window_sin)]) / len(window_vals)
        magnitude = math.sqrt(hx**2 + hy**2)
        result.append(magnitude)
    return result


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


# Kyte-Doolittle hydropathy scores
hydropathy = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

#SEQUENCE = "VTNSLTTFFNQVSKLKNKLEPILKGEVVGAAISYSIILGFP"
#SEQUENCE = "NASGKLLALDSQLTNDFSEKSSYFQSQVDKIRKEAYAGAAA"
#SEQUENCE = "SLEKITEDFTQWPIVQDLYKNYLDLAGDATEIANKVVEVTK"
#SEQUENCE = "KASMKRKGVELKDIKRCLGYHLDVSLAFSEISVGAEFNKDDCVKRGEGRAVNITSENLI"

SEQUENCE = "DVYFRTTETETKIEGIAVIETTLKLKAADIDKNAQKVTNSLTTFFNQVSKLKNKLEPILKGEVVGAAISYSIILGFP" # red
#SEQUENCE = "KASAQKDILIKVLDDGITKLNEAQKSLLVSSQSFNNASGKLLALDSQLTNDFSEKSSYFQSQVDKIRKEAYAGAAA" # blue
#SEQUENCE = "FLLIYAALLQTAVGCWEYVTQTAEFYKDQSDMLLTKIDGVLVSAAQSYEQKFRSLEKITEDFTQWPIVQDLYKNYLDLAGDATEIANKVVEVTK" # green
ANGLE_PER_RESIDUE = 100

# Set up color normalization and colormap
min_score = -4.5
max_score = 4.5
norm = mcolors.Normalize(vmin=min_score, vmax=max_score)
cmap = cm.get_cmap("RdBu_r")  # Red for hydrophobic, Blue for hydrophilic

hydropathies = [hydropathy[res] for res in SEQUENCE]
angles = [((i * ANGLE_PER_RESIDUE) % 360) for i in range(len(SEQUENCE))]
moments = moving_moment_edge_friendly(hydropathies, angles, 16)

fig, ax = plt.subplots()

xticks = range(len(SEQUENCE))
ax.set_xticks(xticks)
ax.set_xticklabels(SEQUENCE)

for label, angle, hydropathy in zip(ax.get_xticklabels(), angles, hydropathies):
    label.set_y(-angle / 360 * 0.1)  # Shift each label vertically (downward = negative)
    color = cmap(norm(hydropathy))
    label.set_color(color)

plt.plot(xticks, moments)
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

charges = np.array([amino_acid_charges[res] for res in SEQUENCE])
charge_smoothed = moving_average_edge_friendly(charges, 3)

aromatic_bars = np.array([0.5 if res in aromatic_residues else 0 for res in SEQUENCE])

plt.figure(figsize=(len(SEQUENCE) / 2, 4))  # Wider figure for readability
plt.bar(xticks, charges, label="Strand 1", color="blue")
plt.bar(xticks, aromatic_bars, color="yellow")
plt.plot(xticks, charge_smoothed)
plt.title("Charge Plot")
plt.xlabel("Residue")
plt.ylabel("Net Charge")
plt.axhline(0, color='gray', linestyle='--')
plt.xticks(ticks=xticks, labels=SEQUENCE)
plt.tight_layout()
plt.show()