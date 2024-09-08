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

def conservation_score(pdb, Pfam_A = './data/Pfam/Pfam-A.hmm'):
    """
    This function calculates the conservation score of each residue in a PDB file. 
    """
    import numpy as np
    import subprocess
    # background distribution of amino acids from BLOSUM62 dataset
    # amino acid order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
    aa_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    blosum62_background = np.array([0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072])
    pdb_to_fasta(pdb)
    ##### to do: parse hmmscan output, align msa with fasta (clustalo), amino acid composition in relevant columns, calculate conservation score
    fasta_path = pdb.replace('.pdb', '.fasta')
    result_hmmscan = subprocess.run(['wsl', 'hmmscan', 'Pfam-A.hmm', fasta_path], stdout=subprocess.PIPE, text=True)