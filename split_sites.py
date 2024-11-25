
import pyrosetta; pyrosetta.init()
from pyrosetta import *
init()
from pyrosetta.toolbox import pose_from_rcsb
from pyrosetta.toolbox import cleanATOM
from pyrosetta.teaching import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from split_energy import energy_intervals
from split_energy import split_energy
from stride import stride
from stride import no_secondary

def split_sites_func(pdb):
    pdb = f'./data/{pdb}.pdb'
    temp = stride(pdb, system='mac', stride_dir='./stride', output_dir='/data/strideout')
    temp[['resnum', 'saa']] = temp[['resnum', 'saa']].astype('float64')
    pose, se_pose = split_energy(pdb)
    merged_pose, filtered_pose = energy_intervals(se_pose)
    split_sites = []
    for i in range(1,pose.total_residue()):
        sites = (i, i+1)
        if (temp.loc[temp['resnum']==sites[0], 'saa'].values[0] > 30) and (temp.loc[temp['resnum']==sites[1], 'saa'].values[0] > 30):
        ### If-condition for conservation score: < 2
            if (temp.loc[temp['resnum']==sites[0], 'code'].isin(['C','T']).iloc[0]) | (temp.loc[temp['resnum']==sites[1], 'code'].isin(['C','T']).iloc[0]):
                split_sites.append(sites)
    loops = []
    start = None
    for i in range(1,pose.total_residue()+1):
        if (temp.loc[temp['resnum']==i, 'code'].isin(['C','T']).iloc[0]):
            if start is None:
                start = i
        else:
            if start is not None:
                loops.append((start, i-1))
                start = None
    if start is not None:
        loops.append((start, i))
    loops = []
    start = None
    for i in range(1,pose.total_residue()+1):
        if (temp.loc[temp['resnum']==i, 'code'].isin(['C','T']).iloc[0]):
            if start is None:
                start = i
        else:
            if start is not None:
                loops.append((start, i-1))
                start = None
    if start is not None:
        loops.append((start, i))
    dump_list = []
    sites_loop = []
    for start, end in loops:
        for i in split_sites:
            if start <= i[0]+1 <= end+1:
                dump_list.append(i)
            else:
                sites_loop.append(dump_list)
                dump_list = []
    sites_loop.append(dump_list)
    sites_loop = [x for x in sites_loop if x != []]
    sites_loop_filtered = []
    for i in sites_loop:
        list_sum = []
        if len(i) > 1:
            for site in i:
                sum = temp.loc[temp['resnum']==site[0], 'saa'].values[0] + temp.loc[temp['resnum']==site[1], 'saa'].values[0]
                list_sum.append(sum)
            sorted = list(np.argsort(list_sum)[::-1][:2])
            if np.abs(i[sorted[0]][0] - i[sorted[1]][0]) > 5:
                sites_loop_filtered.append([i[j] for j in sorted])
            else:
                sites_loop_filtered.append([i[sorted[0]]])
        else:
            sites_loop_filtered.append(i)
    for i in sites_loop_filtered[:]:
        bool = []
        if len(i) < 2:
            for start, end in filtered_pose:
                if not (start <= i[0][0] + 1 <= end + 1):
                    bool.append(False)
                else:
                    bool.append(True)
            if not any(bool):
                sites_loop_filtered.remove(i)
    return sites_loop_filtered, pose