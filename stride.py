def stride(pdb_file, stride_dir = './stride/stride', output_dir = './data/strideout', save = False, modify_pdb = False, system = 'mac'):
    """
    This function runs the Stride program on a PDB file and returns the secondary structure information and solvent accessable area of each residue. The output is a pandas dataframe. If save = True, the output is saved as a csv file in the specified output directory.
    """
    import subprocess
    import pandas as pd
    import os
    if system == 'windows':
        result = subprocess.run(['wsl', stride_dir, pdb_file], stdout=subprocess.PIPE, text=True)
    elif system == 'mac':
        result = subprocess.run(['stride', pdb_file], stdout=subprocess.PIPE, text=True)
    secondary = [line for line in result.stdout.splitlines() if line.startswith('ASG')]
    data_str = "n".join(secondary)
    if len(secondary[0].split()) == 10:
        secondary_df = pd.DataFrame([line.split() for line in secondary], columns=['ASG', 'aa', 'chain', 'pdbpos', 'resnum', 'code', 'name', 'phi', 'psi', 'saa'])
    elif len(secondary[0].split()) == 11:
        secondary_df = pd.DataFrame([line.split() for line in secondary], columns=['ASG', 'aa', 'chain', 'pdbpos', 'resnum', 'code', 'name', 'phi', 'psi', 'saa', 'temp']).drop(columns='temp')
    if save:
        secondary_df.to_csv(os.path.join(output_dir, os.path.basename(pdb_file).replace('.pdb', '_stride.csv')))
    if modify_pdb:
        saa_dict = {int(i) + 1 : code for i, code in enumerate(secondary_df['saa'])}
        amino_acid_dict = {int(i) + 1 : code for i, code in enumerate(secondary_df['code'])}
        with open(pdb_file, "r+") as f:
            lines = f.readlines()
            if lines[0].split()[-1] != 'MODIFIED':
                lines[0] = lines[0].rstrip() + ' ' + 'MODIFIED\n'
                for i,line in enumerate(lines):
                    if line.startswith('ATOM'):
                        lines[i] = line.rstrip() + '    ' + amino_acid_dict[int(line.split()[5])] + '    ' + saa_dict[int(line.split()[5])] + '\n'
            else:
                print('PDB file already modified.')
            f.seek(0)
            f.writelines(lines)
            f.truncate()
    return secondary_df

def no_secondary(pdb_file, stride_dir = './stride/stride', output_dir = './data/strideout', save = False, modify_pdb = False, system = 'windows'):
    """
    This function returns the residues that are not part of an alpha helix or beta sheet. The output is a pandas dataframe of only those amino acids.
    """
    from stride import stride
    secondary_df = stride(pdb_file, stride_dir=stride_dir, output_dir=output_dir, save=save, modify_pdb=modify_pdb, system = system)
    no_sec = secondary_df[secondary_df['code'].isin(['C', 'T'])]
    return no_sec

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