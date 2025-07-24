import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Kyte-Doolittle hydropathy scores
hydropathy = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Example sequence
sequence = "KASMKEKGVDLNDVKHCLGFNMDLRIPLQDDLKDASVTASVNADGCIKTDNGKTVDITRDNI"

# Parameters
angle_per_residue = 100  # degrees per residue
radius = 1.5

# Set up color normalization and colormap
min_score = -4.5
max_score = 4.5
norm = mcolors.Normalize(vmin=min_score, vmax=max_score)
cmap = cm.get_cmap("RdBu_r")  # Red for hydrophobic, Blue for hydrophilic

# Setup figure
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
ax.axis('off')

# Store coordinates
coords = []
hx, hy = 0, 0  # Components of the hydrophobic moment

base_radius = 1.5
radius_increment = 0.02  # Small spiral expansion

for i, aa in enumerate(sequence):
    angle_deg = i * angle_per_residue
    angle_rad = np.deg2rad(angle_deg)
    radius = base_radius + i * radius_increment  # Increase radius slightly each step
    x = radius * np.cos(angle_rad)
    y = radius * np.sin(angle_rad)
    coords.append((x, y))

    score = hydropathy.get(aa, 0)
    color = cmap(norm(score))

    # Draw residue circle
    ax.plot(x, y, 'o', color=color, markersize=20)
    ax.text(x, y, aa, color='white', ha='center', va='center', fontsize=12)

    # Add vector component for hydrophobic moment
    hx += score * np.cos(angle_rad)
    hy += score * np.sin(angle_rad)

max_radius = base_radius + len(sequence) * radius_increment
ax.set_xlim(-max_radius - 0.5, max_radius + 0.5)
ax.set_ylim(-max_radius - 0.5, max_radius + 0.5)

# Draw connecting lines
for i in range(len(coords) - 1):
    x1, y1 = coords[i]
    x2, y2 = coords[i + 1]
    ax.plot([x1, x2], [y1, y2], color='gray', linewidth=1)

# Draw hydrophobic moment arrow
ax.arrow(0, 0, hx/20, hy/20, color='black', width=0.01, head_width=0.1, length_includes_head=True)
ax.text(hx/20, hy/20, "← Hydrophobic face", ha='left', va='bottom', fontsize=10, color='black')

plt.tight_layout()
plt.title("Helical Wheel Projection with Hydropathy Gradient")
plt.show()
