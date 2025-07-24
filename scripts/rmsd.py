import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import os

for i in range (1, len(argv), 2):
    traj_file = argv[i]       # Trajectory file (e.g., traj.dcd or traj.pdb)
    ref_file = argv[i+1]        # Reference structure (e.g., reference.pdb)
    traj_path = os.path.join("outputs", traj_file)
    ref_path = os.path.join("inputs", ref_file)

    # ---- LOAD FILES ----
    print("Loading trajectory...")
    traj = md.load(traj_path)
    reference = md.load(ref_path)

    print(f"# backbone atoms in traj: {len(traj.topology.select('backbone'))}")
    print(f"# backbone atoms in reference: {len(reference.topology.select('backbone'))}")

    # ---- SELECT ATOMS ----
    # Here we use backbone atoms (N, CA, C) for RMSD
    atom_indices = traj.topology.select("backbone")

    # ---- ALIGN & CALCULATE RMSD ----
    print("Calculating RMSD...")
    traj.superpose(reference, atom_indices=atom_indices)
    rmsd = md.rmsd(traj, reference, atom_indices=atom_indices)

    # ---- SAVE TO FILE ----
    np.savetxt("rmsd.txt", np.column_stack([traj.time, rmsd]), header="Time(ps) RMSD(nm)")

    # ---- PLOT ----
    plt.plot(traj.time, rmsd, label=traj_file)

plt.xlabel("Frame")
plt.ylabel("RMSD (nm)")
plt.title("RMSD vs Time")
plt.grid(True)
plt.savefig("rmsd_plot.png")
plt.legend()
plt.show()