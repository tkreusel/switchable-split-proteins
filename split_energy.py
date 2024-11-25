import pyrosetta; pyrosetta.init()
from pyrosetta import *
init()
from pyrosetta.teaching import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from pyrosetta.toolbox import pose_from_rcsb
from scipy.signal import argrelextrema
import os
def split_energy(pdb:str, window = 3):
    pose = pose_from_pdb(pdb) # Load the protein structure from PDB
    n_residues = pose.total_residue() # Get the total number of residues
    scorefxn = get_score_function(True)  # Get a scoring function (standard Rosetta energy function)
    scorefxn.set_weight(dslf_fa13, 0) #remove disulfide bond energy term
    E_nat = scorefxn(pose) # Calculate the total energy of the full protein

     # Check if the file exists and delete it
    if os.path.exists(f"{pdb}.pdb"):
        os.remove(f"{pdb}.pdb")
    if os.path.exists(f"{pdb}.clean.pdb"):
        os.remove(f"{pdb}.clean.pdb")

    split_energies = []  
    for i in range(1, n_residues):
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

    # Smooth the split energy profile five times with rolling window of 3
    smoothed_values = se_profile['Split Energy'].copy()
    for _ in range(5):
        smoothed_values = smoothed_values.rolling(window=window, center=True, min_periods=1).mean()
    se_profile['Smoothed Split Energy'] = smoothed_values

    return pose, se_profile

def energy_intervals(se_profile:pd.DataFrame, order:int = 20):
    # Calculate numeric derivative and identify low gradient residues
    se_profile['Derivative'] = np.gradient(se_profile['Smoothed Split Energy'].values, se_profile['Residue'].values)
    selected_residues = se_profile[np.abs(se_profile['Derivative']) < 0.5]['Residue'].values  # Select residues where the derivative is less than the threshold

    # Group consecutive residues into intervals
    intervals = []
    if len(selected_residues) > 0:
        start = selected_residues[0]
        for i in range(1, len(selected_residues)):
            if selected_residues[i] != selected_residues[i-1] + 1: # If the current residue is not consecutive with the previous one, create a new interval
                end = selected_residues[i-1]
                intervals.append((start, end))
                start = selected_residues[i]
        intervals.append((start, selected_residues[-1])) # Append the last interval
    
    smoothed_energy = se_profile.set_index('Residue')['Smoothed Split Energy']
    extended_intervals = []

    # Extend the intervals to the left and right based on condition: |E(i) - E(start)| < |min(E(i))|/20
    for start, end in intervals:
        left = start
        right = end
        while left > 1: # Extend to the left
            if np.abs(smoothed_energy[left - 1] - smoothed_energy[start]) < (np.abs(se_profile['Smoothed Split Energy'].min()))/20:
                left -= 1
            else:
                break
        while right < se_profile['Residue'].max(): # Extend to the right
            if np.abs(smoothed_energy[right + 1] - smoothed_energy[end]) < (np.abs(se_profile['Smoothed Split Energy'].min()))/20:
                right += 1
            else:
                break
        extended_intervals.append((left, right))

    # Merge overlapping intervals
    extended_intervals.sort()
    merged_intervals = []
    current_start, current_end = extended_intervals[0]
    for start, end in extended_intervals[1:]:
        if start <= current_end + 1:  # Overlapping or adjacent intervals
            current_end = max(current_end, end) # Merge intervals
        else:
            merged_intervals.append((current_start, current_end))
            current_start, current_end = start, end
    merged_intervals.append((current_start, current_end)) # Append the last interval

    # Remove intervals that contain local minima
    energy_profile = se_profile['Smoothed Split Energy'].values
    residues = se_profile['Residue'].values
    filtered_intervals = []
    for start, end in merged_intervals:
        contains_minimum = any((start <= r <= end) for r in residues[argrelextrema(energy_profile, np.less, order=order)[0]]) # Check if the interval contains any local minima
        is_boundary_interval = (start == residues[0]) or (end == residues[-1]) # Check if the first or last intervals include the first or last residues
        if not contains_minimum and not is_boundary_interval:  #If the interval contains a minimum or is a boundary interval, discard it
            filtered_intervals.append((start, end))

    return np.array(merged_intervals), np.array(filtered_intervals)

def plot(se_profile:pd.DataFrame, merged_intervals:np.ndarray=np.empty((0,2)), filtered_intervals:np.ndarray=np.empty((0,2)), type:str = 'profile'):
    # Plot the original and smoothed split energy profiles without intervals
    if type == 'profile':
        fig, ax = plt.subplots()
        ax.plot(se_profile['Residue'], se_profile['Split Energy'], label='Original Split Energy', alpha=0.5)
        ax.plot(se_profile['Residue'], se_profile['Smoothed Split Energy'], label='Smoothed Split Energy (5x)', color='red', linewidth=2)
        ax.set_xlabel("Residue")
        ax.set_ylabel("Split Energy")
        ax.set_title("Split Energy Profile (Smoothed)")
        ax.legend()
        plt.show()

    elif type == 'intervals':
        fig, ax = plt.subplots(1, 2, figsize=(15, 5))
        # Plot the smoothed split energy profile with the extended intervals
        ax[0].plot(se_profile['Residue'], se_profile['Smoothed Split Energy'], label='Smoothed Split Energy', color='blue')
        for start, end in merged_intervals:
            ax[0].axvspan(start, end, color='red', alpha=0.3, label='Extended Interval' if start == merged_intervals[0][0] else "")
        ax[0].set_xlabel('Residue')
        ax[0].set_ylabel('Smoothed Split Energy')
        ax[0].set_title('Split Energy Profile with Extended Intervals')
        ax[0].xaxis.set_major_locator(MultipleLocator(10))
        ax[0].tick_params(axis='x', rotation=90)
        for tick in ax[0].xaxis.get_majorticklocs():
            ax[0].axvline(x=tick, color='grey', linestyle='--', alpha=0.5)
        for tick in [7, 42, 47, 86, 133]:
            ax[0].axvline(x=tick, color='black', linestyle='--', alpha=0.5)



        handles, labels = ax[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax[0].legend(by_label.values(), by_label.keys())
        # Plot the smoothed split energy profile with the filtered intervals
        ax[1].plot(se_profile['Residue'], se_profile['Smoothed Split Energy'], label='Smoothed Split Energy', color='blue')
        for start, end in filtered_intervals:
            ax[1].axvspan(start, end, color='red', alpha=0.3, label='Extended Interval' if start == filtered_intervals[0][0] else "")
        ax[1].set_xlabel('Residue')
        ax[1].set_ylabel('Smoothed Split Energy')
        ax[1].set_title('Split Energy Profile with Extendend Intervals (without Local Minima)')
        ax[1].xaxis.set_major_locator(MultipleLocator(10))
        ax[1].tick_params(axis='x', rotation=90)
        for tick in ax[1].xaxis.get_majorticklocs():
            ax[1].axvline(x=tick, color='grey', linestyle='--', alpha=0.5)
        for tick in [7, 42, 47, 86, 133]:
            ax[1].axvline(x=tick, color='black', linestyle='--', alpha=0.5)
            
        handles, labels = ax[1].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax[1].legend(by_label.values(), by_label.keys())
        plt.tight_layout()
        plt.show()
