import os
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import prody

prody.confProDy(verbosity='none')

def read_npz(sample_id, npz_data_dir):
    npz_file = os.path.join(npz_data_dir, sample_id + '.npz')
    loaded = np.load(npz_file)
    features = dict(loaded)
    return features

def calc_rmsd(struct_dir, mask, align_larger_chain=True, rmsd_ca_atoms=True, save_to_pdbs=False):
    """Align using CA atoms and calculate RMSD of the peptide.

    Args:
        struct_dir (str): Directory with native and predicted structures
        mask (1D array): Binary mask applied to filter out residues.
        align_larger_chain (bool, optional): Whether to align using CA atoms 
            of larger chain (True) or all chains (False). Defaults to True.
        rmsd_ca_atoms (bool, optional): Whether to calculate RMSD from positions 
            of CA atoms of the peptide (True) or all heavy atoms (False). Defaults to True.
        save_to_pdbs (bool, optional): Whether to save structures after filtering 
            and aligning. Defaults to False.

    Returns:
        int: value of RMSD
    """
    true_pdb_path = str(next(Path(struct_dir).glob('true*.pdb')))
    pred_pdb_path = str(next(Path(struct_dir).glob('val*.pdb')))
    
    # Read in PDBs
    true = prody.parsePDB(true_pdb_path)
    pred = prody.parsePDB(pred_pdb_path)
    
    # Filter out residues using mask
    select_str = " ".join([str(idx) for idx in np.where(mask)[0]])
    true = true.select(f'resindex {select_str}').toAtomGroup()
    pred = pred.select(f'resindex {select_str}').toAtomGroup()

    # Align structures
    if align_larger_chain:
        prody.superpose(true.select('chain B and calpha'), pred.select('chain B and calpha'))
    else:
        prody.superpose(true.select('calpha'), pred.select('calpha'))

    # Select peptide chain
    if rmsd_ca_atoms:
        chain_true = true.select(f'chain C and calpha')
        chain_pred = pred.select(f'chain C and calpha')
    else:
        chain_true = true.select(f'chain C')
        chain_pred = pred.select(f'chain C')

    # Calculate RMSD
    rmsd = prody.calcRMSD(chain_true, chain_pred)

    # Save filtered and aligned complexes
    if save_to_pdbs:
        prody.writePDB(os.path.join(struct_dir, 'aligned_filtered_true.pdb'), true)
        prody.writePDB(os.path.join(struct_dir, 'aligned_filtered_val.pdb'), pred)

    return rmsd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Add columns rmsd1, rmsd2, rmsd3 '
                    'to metrics.csv and save the new '
                    'table to metrics_with_rmsds.csv. '
                    'Save aligned structures to pdb files'
    )
    parser.add_argument('--data_dir', type=str, required=True)
    parser.add_argument('--npz_data_dir', type=str, required=True)
    args = parser.parse_args()

    data_dir = args.data_dir
    metrics_path = os.path.join(data_dir, 'metrics.csv')
    structures_dir = os.path.join(data_dir, 'structures')
    npz_data_dir = args.npz_data_dir

    df = pd.read_csv(metrics_path)

    df['rmsd1'] = -1.0
    df['rmsd2'] = -1.0
    df['rmsd3'] = -1.0
    for idx, row in df.iterrows():
        struct_dir = os.path.join(structures_dir, row['sample_name'])

        features = read_npz(row['sample_name'], npz_data_dir)
        mask = features['all_atom_mask'].any(axis=1)

        rmsd = calc_rmsd(struct_dir, mask, save_to_pdbs=True)
        df.loc[idx, 'rmsd1'] = rmsd

        rmsd = calc_rmsd(struct_dir, mask, rmsd_ca_atoms=False)
        df.loc[idx, 'rmsd2'] = rmsd

        rmsd = calc_rmsd(struct_dir, mask, align_larger_chain=False, rmsd_ca_atoms=False)
        df.loc[idx, 'rmsd3'] = rmsd

    df.to_csv(os.path.join(data_dir, 'metrics_with_rmsds.csv'), index=False)
