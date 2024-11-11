def find_nearest(pdb, n_nearest):
    import numpy as np
    from scipy.spatial import distance
    """
    This function finds the n_nearest atoms to a given atom in a pdb file.
    """
    with open(pdb, 'r') as f:
        lines = [line for line in f if line.startswith('ATOM') and line[13:15] == 'CA']
    coords = np.array([line.split()[6:9] for line in lines], dtype=float)
    distance_matrix = distance.cdist(coords, coords)
    nearest_ind = np.array([np.argsort(distance_matrix[i])[1:n_nearest+1] for i in range(len(distance_matrix))])
    nearest_dist = np.array([np.sort(distance_matrix[i])[1:n_nearest+1] for i in range(len(distance_matrix))])
    return nearest_ind, nearest_dist

def jsd(p, q):
    """
    This function calculates the Jensen-Shannon divergence between two probability distributions.
    """
    from scipy.special import rel_entr
    import numpy as np
    m = 0.5 * (p + q)
    jsd = 0.5 * np.sum(rel_entr(p, m)) + 0.5 * np.sum(rel_entr(q, m))
    return jsd

def mmseqs_cons(input_fasta : str, database_fasta : str, input_pdb : str = None, max_seqs:int = 1000, include_neighbors:bool=False, n_nearest:int=5, decay:float=0.5, alpha:float=0.5, output_path = './data/mmseqs_data/alignment_msa.fasta', method = 'jsd'):
    """
    This function calculates a conservation score for each residue in an input protein (fasta) using the provided database (fasta). If use_neighbors is set to True it also requires an input_pdb to calculate conservation using n_nearest spatial neighbors weighted by distance with a decay. Alpha determines how much the neighbors affect the scoring, proposed by https://academic.oup.com/bioinformatics/article/23/15/1875/203579 . Score is calculated either by Jensen-Shannon divergence (jsd) or Kullback-Leibler divergence (kl) (method parameter)."""
    import subprocess
    import numpy as np
    import os
    from Bio import AlignIO
    from scipy.special import rel_entr
    if not os.path.exists('./data/mmseqs_data'):
        os.makedirs('./data/mmseqs_data')
        subprocess.run(["mmseqs","createdb",database_fasta,"./data/mmseqs_data/database_db"])                #create mmseq db file from fasta (database)
    subprocess.run(["mmseqs","createdb", input_fasta, "./data/mmseqs_data/protein_db"])   
    subprocess.run(["mmseqs","search","./data/mmseqs_data/protein_db","./data/mmseqs_data/database_db","./data/mmseqs_data/result","./data/mmseqs_data/tmp","--max-seqs",str(max_seqs)])
    subprocess.run(["mmseqs","result2msa","./data/mmseqs_data/protein_db","./data/mmseqs_data/database_db","./data/mmseqs_data/result",output_path,"--msa-format-mode", "2"])     # generate multiple sequence alignment from search result
    with open(output_path, "rb+") as file:     # remove NUL character in fasta
        file.seek(-1, 2)       # move pointer to last byte and go back one more
        file.truncate()

        aa_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    blosum62_background = np.array([0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072])
    alignment = AlignIO.read('./data/mmseqs_data/alignment_msa.fasta', 'fasta')   # read msa using biopython
    alignment_array = np.array([list(rec.seq) for rec in alignment], dtype = 'U1')
    aa_prop_array = np.zeros((20, alignment_array.shape[1]))   # initiate array to store amino acid proportions
    for n_col, col in enumerate(alignment_array.T):            # calculate amino acid proportions for each column
        col_no_gaps = np.char.upper(col[col != '-'])
        aa_percs = np.array([np.sum(col_no_gaps == aa)/len(col_no_gaps) for aa in aa_order])    
        aa_prop_array[:, n_col] = aa_percs                              
    conservation_scores = np.zeros(len(alignment_array.T))
    for n_col, col in enumerate(aa_prop_array.T):
        if method == 'jsd':
            conservation_scores[n_col] = jsd(col, blosum62_background)  # compute Jensen-Shannon divergence for each column
        elif method == 'rel_entr':
            conservation_scores[n_col] = rel_entr(col,blosum62_background).sum()
    if include_neighbors:
        ind,dist = find_nearest(input_pdb, n_nearest)
        weights = np.array([np.exp(-decay*row) for row in dist])
        weightsums = np.sum(weights, axis = 1)
        neighbors = conservation_scores[ind]
        neighbors_weighted = np.sum(neighbors*weights, axis = 1)/weightsums
        conservation_scores_neighbors = alpha * conservation_scores + (1-alpha) * neighbors_weighted
        return conservation_scores, conservation_scores_neighbors
    else:
        return conservation_scores