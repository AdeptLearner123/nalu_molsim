from sys import argv
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os

for i in range(1, len(argv)):
    label, input_files = argv[i].split(":")
    all_alpha_counts = []

    for input_file in input_files.split(","):
        print("Calculating " + input_file)

        # Load your trajectory or structure file
        # Replace with your filenames
        input_path = os.path.join("outputs", input_file)
        traj = md.load(input_path)

        # Compute secondary structure using DSSP
        ss = md.compute_dssp(traj)

        # ss is a 2D array: (n_frames, n_residues)
        # Each element is a DSSP code like 'H' for alpha helix

        # Count alpha helix ('H') residues per frame
        alpha_counts = np.sum(ss == 'H', axis=1)
        all_alpha_counts.append(alpha_counts)

        # Print average and plot
        print(f"Average number of alpha-helix residues: {np.mean(alpha_counts):.2f}")
    
    all_alpha_counts = np.array(all_alpha_counts)
    elementwise_avg = np.mean(all_alpha_counts, axis=0)
    plt.plot(elementwise_avg, label=label)

plt.xlabel("Frame")
plt.ylabel("# Alpha-Helix Residues")
plt.title("Alpha Helix Content Over Time")
plt.grid(True)
plt.legend()
plt.show()