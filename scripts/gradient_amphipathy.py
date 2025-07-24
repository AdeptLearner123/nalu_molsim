import matplotlib.pyplot as plt

hydropathy_scores = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

sequence = "KASMKRKGVELKDIKRCLGYHLDVSLAFSEISVGAEFNKDDCVKRGEGRAVNITSENLI"
values = [hydropathy_scores[aa] for aa in sequence]

def moving_average(values, window):
    return [sum(values[i:i+window]) / window for i in range(len(values) - window + 1)]

window_size = 9
smoothed = moving_average(values, window_size)

positions = range(window_size // 2, len(sequence) - window_size // 2)
residue_labels = sequence[window_size // 2 : len(sequence) - window_size // 2]

plt.figure(figsize=(len(residue_labels) / 2, 4))  # Wider figure for readability
plt.plot(positions, smoothed, marker='o')
plt.title("Hydropathy Plot (Kyte-Doolittle)")
plt.xlabel("Residue")
plt.ylabel("Hydropathy")
plt.axhline(0, color='gray', linestyle='--')
plt.xticks(ticks=positions, labels=residue_labels, rotation=90)
plt.tight_layout()
plt.show()