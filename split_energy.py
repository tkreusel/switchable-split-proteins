import pyrosetta; pyrosetta.init()
from pyrosetta import *
init()
from pyrosetta.teaching import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyrosetta.toolbox import pose_from_rcsb
import os
def split_energy(pdb:str, window = 3):
    pose = pose_from_rcsb(pdb) # Load the protein structure from PDB
    n_residues = pose.total_residue() # Get the total number of residues
    scorefxn = get_score_function(True)  # Get a scoring function (standard Rosetta energy function)
    E_nat = scorefxn(pose) # Calculate the total energy of the full protein

     # Check if the file exists and delete it
    if os.path.exists(f"{pdb}.pdb"):
        os.remove(f"{pdb}.pdb")
    if os.path.exists(f"{pdb}.clean.pdb"):
        os.remove(f"{pdb}.clean.pdb")

    split_energies = []  
    for i in range(1, n_residues-1):
       # Create Pose A (residues 1 to i)
        pose_A = Pose()
        pose_A.assign(pose) 
        pose_A.delete_residue_range_slow(i+1, n_residues)
        # Create Pose B (residues i+1 to n)
        pose_B = Pose()
        pose_B.assign(pose)
        pose_B.delete_residue_range_slow(1, i)
        # Calculate energies for fragments A and B, and SE
        energy_A = scorefxn(pose_A)
        energy_B = scorefxn(pose_B)
        E_split = E_nat - (energy_A + energy_B)
        
        split_energies.append((i, E_split))

    se_profile = pd.DataFrame(split_energies, columns=['Residue', 'Split Energy'])

    smoothed_values = se_profile['Split Energy'].copy()

    for _ in range(5): # Apply the smoothing 5 times
        smoothed_values = smoothed_values.rolling(window=window, center=True, min_periods=1).mean()
    se_profile['Smoothed Split Energy'] = smoothed_values

    # Plot both the original and smoothed split energy profiles
    fig, ax = plt.subplots()
    ax.plot(se_profile['Residue'], se_profile['Split Energy'], label='Original Split Energy', alpha=0.5)
    ax.plot(se_profile['Residue'], se_profile['Smoothed Split Energy'], label='Smoothed Split Energy (5x)', color='red', linewidth=2)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Split Energy")
    ax.set_title("Split Energy Profile (Smoothed)")
    ax.legend()
    se_plot = fig
    
    return se_profile, se_plot 
