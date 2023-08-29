import os
import subprocess
import argparse

import pandas as pd
import prody


def add_text_to_video(input_path: str, output_path: str, text: str):
    command = [
        './add_text.sh',
        text,
        input_path,
        output_path
    ]

    ffmpeg_process = subprocess.Popen(command)
    ffmpeg_process.wait()

def calc_rmsd(true_pdb, pred_pdb, align_larger_chains=True, rmsd_ca_atoms=True):
    """Align using CA atoms and calculate RMSD of the peptide.

    Args:
        true_pdb (str): Native
        pred_pdb (str): Predicted
        align_larger_chains (bool, optional): 
            Whether to align using CA atoms 
            of larger chains or all chains. Defaults to True.
        rmsd_ca_atoms (bool, optional):
            Whether to calculate RMSD from 
            positions of CA atoms of the peptide or all heavy atoms. Defaults to True.

    Returns:
        int: value of RMSD
    """
    true = prody.parsePDB(true_pdb)
    pred = prody.parsePDB(pred_pdb)

    if align_larger_chains:
        prody.superpose(true.select('chain B and calpha'), pred.select('chain B and calpha'))
    else:
        prody.superpose(true.select('calpha'), pred.select('calpha'))

    if rmsd_ca_atoms:
        chain_true = true.select(f'chain C and calpha')
        chain_pred = pred.select(f'chain C and calpha')
    else:
        chain_true = true.select(f'chain C')
        chain_pred = pred.select(f'chain C')

    rmsd = prody.calcRMSD(chain_true, chain_pred)

    return rmsd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Overlay text to each video in videos_dir')
    parser.add_argument('--metrics_path', type=str, required=True, 
                        help="Path to .csv with metrics and RMSDs")
    parser.add_argument('--videos_dir', type=str, required=True, 
                        help="Path to directory with videos")
    parser.add_argument('--videos_with_text_dir', type=str, required=True, 
                        help="Path to directory where to put the videos with text overlays")
    args = parser.parse_args()

    metrics_path = args.metrics_path
    videos_dir = args.videos_dir
    videos_with_text_dir = args.videos_with_text_dir

    df = pd.read_csv(metrics_path)

    for _, row in df.iterrows():
        # Input path
        input_path = os.path.join(videos_dir, row['sample_name'] + '.mp4')

        # Output path
        os.makedirs(videos_with_text_dir, exist_ok=True)
        output_path = os.path.join(videos_with_text_dir, row['sample_name'] + '.mp4')

        # Text
        text = (
            f"   sample:\t{row['sample_name'].split('_', 1)[0]}  \t\tRMSD1: {row['rmsd1']:.2f}\n"
            f"   pLDDT:\t{row['masked_plddt']:.0f} \t\t\tRMSD2: {row['rmsd2']:.2f}\n"
            f"   LDDT:\t{row['masked_lddt']:.0f} \t\t\tRMSD3: {row['rmsd3']:.2f}"
        )

        filepath = '/tmp/add_text_all_1803571'
        with open(filepath, 'w') as f:
            f.write(text)

        # Subprocess
        add_text_to_video(input_path, output_path, filepath)
