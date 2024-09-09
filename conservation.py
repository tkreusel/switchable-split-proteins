def pdb_to_seq(pdb):
    """
    This function returns the amino acid sequence from a PDB file.
    """
    aacodes = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    with open(pdb, 'r') as f:
        seqdict = {}
        lines = f.readlines()
        for line in lines:
            if line.startswith('ATOM'):
                seqdict[int(line.split()[5])] = line.split()[3]
        seq3letter = [seqdict[i] for i in seqdict.keys()]
        seq = [aacodes[i] for i in seq3letter]
        return ''.join(seq)

def pdb_to_fasta(pdb, fasta_path = None):
    import os
    """
    This function returns the amino acid sequence from a PDB file in FASTA format.
    """
    if fasta_path == None:
        fasta_path = pdb.replace('.pdb', '.fasta')
    with open(fasta_path, 'w') as f:
        f.write('>' + os.path.basename(pdb).replace('.pdb', '') + '\n')
        f.write(pdb_to_seq(pdb))
    return None

def jsd(p, q):
    """
    This function calculates the Jensen-Shannon divergence between two probability distributions.
    """
    from scipy.special import rel_entr
    import numpy as np
    m = 0.5 * (p + q)
    jsd = 0.5 * np.sum(rel_entr(p, m)) + 0.5 * np.sum(rel_entr(q, m))
    return jsd

def conservation_score(pdb :str, Pfam_A : str = './data/Pfam/Pfam-A.hmm', hmmscanout : str = './data/result.tbl', family = None):
    """
    This function calculates the conservation score of each residue in a PDB file. 
    """
    import numpy as np
    import subprocess
    from Bio import AlignIO
    # background distribution of amino acids from BLOSUM62 dataset
    # amino acid order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
    aa_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    blosum62_background = np.array([0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072])
    pdb_to_fasta(pdb)
    fasta_path = pdb.replace('.pdb', '.fasta')
    if family == None:
        subprocess.run(['wsl', '--tblout', hmmscanout, 'hmmscan', 'Pfam-A.hmm', fasta_path], stdout=subprocess.PIPE, text=True)
        with open('./data/result.tbl', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line[0] != '#': 
                    if float(line.split()[4]) < 0.01:
                        family = line.split()[1].split('.')[0]
                        print(family)
                        break
                    else:
                        print('No significant Pfam family found.')
                        break
    ######## missing: use found family msa to align protein fasta to
    if family == 'PF00107':
        family = './data/Pfam/PF00107.alignment.full'
    # align protein sequence to Pfam family MSA using clustal omega
    subprocess.run(['wsl', 'clustalo', '--p1', family, '--p2', fasta_path, '-o', './data/output_alignment.fasta', '--threads=1'], stdout = subprocess.PIPE, text = True)
    alignment = AlignIO.read('./data/Pfam/output_alignment.fasta', 'fasta')   # read msa using biopython
    alignment_array = np.array([list(rec.seq) for rec in alignment], dtype = 'U1')
    prot_seq = alignment[-1].seq    # last sequence is protein of interest
    columns = np.where(np.array(prot_seq) != '-')[0]
    alignment_array_filtered = alignment_array[:, columns]
    aa_prop_array = np.zeros((20, alignment_array_filtered.shape[1]))   # initiate array to store amino acid proportions
    for n_col, col in enumerate(alignment_array_filtered.T):            # calculate amino acid proportions for each column
        col_no_gaps = np.char.upper(col[col != '-'])
        aa_percs = np.array([np.sum(col_no_gaps == aa)/len(col_no_gaps) for aa in aa_order])    
        aa_prop_array[:, n_col] = aa_percs
    conservation_scores = np.zeros(len(columns))
    for n_col, col in enumerate(aa_prop_array.T):
        conservation_scores[n_col] = jsd(col, blosum62_background)  # compute Jensen-Shannon divergence for each column
    return conservation_scores