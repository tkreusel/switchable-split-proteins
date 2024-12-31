import pyrosetta; pyrosetta.init("-corrections::beta_nov16")
from pyrosetta import *
from pyrosetta.rosetta import *
from rosetta.protocols.rigid import *
import sys
from stride import stride

energy_calc_dict = {}



def split_pose_by_residue(pose, split_site, to_pdb=False, name = "Protein"):
    
    """
    Splits a Pose object into two parts at the given residue position.

    Args:
        pose (Pose): The PyRosetta Pose object to be split.
        split_site (int): The residue index after which to split the Pose.

    Returns:
        pose_N (Pose): N-terminal Pose (residues 1 to split_site).
        pose_C (Pose): C-terminal Pose (residues split_site+1 to end).
    """

    if split_site < 1 or split_site >= pose.total_residue():
        raise ValueError(
            f"Invalid split position: {split_site}. "
            f"Valid range is 1 to {pose.total_residue() - 1}."
        )
    
    # Create N-terminal Pose (residues 1 to split_site)
    pose_N = Pose()
    pose_N.assign(pose)
    pose_N.delete_residue_range_slow(split_site + 1, pose_N.total_residue())

    
    # Create C-terminal Pose (residues split_site+1 to end)
    pose_C = Pose()
    pose_C.assign(pose)
    pose_C.delete_residue_range_slow(1, split_site)

    if to_pdb is True:
        pose_N.dump_pdb(f"{name}_split_N.pdb")
        pose_C.dump_pdb(f"{name}_split_C.pdb")

    return pose_N, pose_C



def fast_relax(pose):
    
    from pyrosetta.rosetta.protocols.relax import FastRelax
    
    # Initialize Fastrelax
    scorefxn = get_fa_scorefxn()
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    
    relax.apply(pose)
    return pose



def protein_protein_docking(pose_1, pose_2, output_prefix = "dock", to_pdb = False, name = "protein", **kwargs):  
    
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover
    from pyrosetta.rosetta.protocols.cluster import ClusterPose


    jobs = kwargs.get("jobs", 70)
    translation = kwargs.get("translation", 3)
    rotation = kwargs.get("rotation", 8.0)
    scorefxn_high = kwargs.get("scorefxn_high", create_score_function("docking"))
    flexible_residues = kwargs.get("flexible_residues", None)
    min_dist = kwargs.get("min_dist", 0)
    max_dist = kwargs.get("max_dist", 20)
        
    
    #apply chain ID for poses
    for i in range(1, pose_1.total_residue() + 1):
        pose_1.pdb_info().chain(i, "A")
        
    for i in range(1, pose_2.total_residue() + 1):
        pose_2.pdb_info().chain(i, "B")
    
    #combined pose for docking
    combined_pose = Pose()
    combined_pose.assign(pose_1)
    combined_pose.append_pose_by_jump(pose_2, jump_anchor_residue = pose_1.total_residue())
    #set distance constraint
    #not yet implemented

    # set up docking FoldTree
    dock_jump = 1
    protocols.docking.setup_foldtree(combined_pose, "A_B", Vector1([dock_jump]))

    # centroid and fa movers
    to_centroid = SwitchResidueTypeSetMover("centroid")
    to_fullatom = SwitchResidueTypeSetMover("fa_standard")

    # original side chain for recovery
    recover_sidechains = protocols.simple_moves.ReturnSidechainMover(combined_pose)

    # convert to centroid mode
    to_centroid.apply(combined_pose)

    # create a (centroid) test pose
    test_pose = Pose()
    test_pose.assign(combined_pose)

    # create scoring function
    scorefxn_low = create_score_function("interchain_cen")
    scorefxn_high = create_score_function("interchain_docking")
    scorefxn_high_min = create_score_function("docking", "docking_min")

    #pertubation movers
    randomize_upstream = RigidBodyRandomizeMover(combined_pose, dock_jump, partner_upstream)
    randomize_downstream = RigidBodyRandomizeMover(combined_pose, dock_jump, partner_downstream)
    dock_perturb = RigidBodyPerturbMover(dock_jump, translation, rotation)
    spin = RigidBodySpinMover(dock_jump)
    slide_into_contact = protocols.docking.DockingSlideIntoContact(dock_jump)

    #setup MinMover
    movemap = MoveMap()
    movemap.set_jump(dock_jump, True)
    if flexible_residues is not None:
        for res in flexible_residues:
            movemap.set_bb(res, True)

    minmover = protocols.minimization_packing.MinMover()
    minmover.movemap(movemap)
    minmover.score_function(scorefxn_high_min)

    # SequenceMover fpr pertubation step
    perturb = protocols.moves.SequenceMover()
    perturb.add_mover(randomize_upstream)
    perturb.add_mover(randomize_downstream)
    perturb.add_mover(dock_perturb)
    perturb.add_mover(spin)
    perturb.add_mover(slide_into_contact)
    perturb.add_mover(to_fullatom)
    perturb.add_mover(recover_sidechains)
    perturb.add_mover(minmover)

    
    #docking protocol
    dock_prot = protocols.docking.DockingProtocol()
    dock_prot.set_movable_jumps(Vector1([1]))
    dock_prot.set_lowres_scorefxn(scorefxn_low)
    dock_prot.set_highres_scorefxn(scorefxn_high_min)

    
    # job distributor for multiple trajectories
    job_distributor = PyJobDistributor(output_prefix, jobs, scorefxn_high)
    job_distributor.native_pose = combined_pose

    #perform docking
    best_pose = None
    
    while not job_distributor.job_complete:
        test_pose.assign(combined_pose)

        # Perturb the pose
        perturb.apply(test_pose)

        # Perform docking
        dock_prot.apply(test_pose)

        # Output the docked structure
        to_fullatom.apply(test_pose)
        recover_sidechains.apply(test_pose)
        job_distributor.output_decoy(test_pose)

        # Store the best pose (lowest energy)
        if best_pose is None or scorefxn_high(test_pose) < scorefxn_high(best_pose):
            best_pose = Pose()
            best_pose.assign(test_pose)


    print("Docking finished")
        
    return best_pose



def calculate_energy(pose, name, scorefxn = None):
    
    global energy_calc_dict
    
    if scorefxn is None:
        scorefxn = create_score_function("beta_nov16")

    energy = scorefxn(pose)
    energy_calc_dict[name] = energy

    return energy_calc_dict

    
def find_interface_residues(path_to_pdb, split_site:int, to_pdb=False, **kwargs):
    from stride import stride
    
    stride_dir = kwargs.get('stride_dir', './data/stride/stride')

    pose = pose_from_pdb(path_to_pdb)
    calculate_energy(pose, "native_protein")

    relaxed_native_pose = fast_relax(pose)
    relaxed_native_pose.dump_pdb("relaxed_native_pose")
    calculate_energy(relaxed_native_pose, "relaxed_native_pose")

    #secondary structure information
    sec_structure_stride_df = stride(pdb_file = path_to_pdb, stride_dir = stride_dir)
    flexible_residues = sec_structure_stride_df.loc[sec_structure_stride_df["code"].isin(["C", "T"]), "pdbpos"].astype("int").tolist()


    #split protein at residue
    pose_N, pose_C = split_pose_by_residue(pose, split_site)

    calculate_energy(pose_N, "split_pose-N")
    calculate_energy(pose_C, "split-pose-C")
    
    #relax split proteins
    relaxed_pose_N = fast_relax(pose_N)
    relaxed_pose_C = fast_relax(pose_C)

    relaxed_pose_N.dump_pdb("relaxed_pose_N")
    relaxed_pose_C.dump_pdb("relaxed_pose_C")
    calculate_energy(relaxed_pose_N, "relaxed_split_pose-N")
    calculate_energy(relaxed_pose_C, "relaxed_split-pose-C")
    
    #perform docking
    docked_pose = protein_protein_docking(relaxed_pose_N, relaxed_pose_C, to_pdb = True, flexible_residues = flexible_residues)
    calculate_energy(docked_pose, "docked_split-pose")

    docked_pose.dump_pdb(f"best_docked_pose.pdb")

    print (energy_calc_dict)
    return relaxed_pose_N, relaxed_pose_C, docked_pose 



def main():
    if len(sys.argv) != 3:
        print("Usage: python find_interface_residues.py <path_to_pdb> <split_site>")
        sys.exit(1)

    # Get inputs from command-line arguments
    path_to_pdb = sys.argv[1]
    split_site = int(sys.argv[2])  # Convert split site to integer

    # Run find_interface_residues
    relaxed_pose_N, relaxed_pose_C, docked_pose = find_interface_residues(path_to_pdb, split_site, to_pdb=True)

    print(f"Relaxed N-terminal pose and C-terminal pose saved.")
    print(f"Docked pose saved.")

if __name__ == "__main__":
    main()
    
    