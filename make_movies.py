import time
import os
from glob import glob
import argparse

import pymol
from pymol import cmd

pymol.finish_launching()
time.sleep(1) #TODO: pyMOL doesn't open without delay

def create_aligned_movie(true_pdb_path, val_pdb_path, movie_path):
    cmd.reinitialize()

    cmd.load(true_pdb_path, 'true')
    cmd.load(val_pdb_path, 'val')
    cmd.color('tv_green', 'true')
    cmd.color('cyan', 'val')
    # res = cmd.cealign('val', 'true')
    cmd.orient('true')
    cmd.zoom('true', buffer=-2)

    nframes = 60
    cmd.mset(specification=f'1 x{nframes}')

    cmd.mview(first=1, power=1)
    cmd.turn(axis='y', angle=120)
    cmd.mview(first=round(nframes/3+1), power=1)
    cmd.turn(axis='y', angle=120)
    cmd.mview(first=round(2*nframes/3+1), power=1)

    cmd.movie.produce(movie_path, encoder='ffmpeg')
    time.sleep(5) #TODO: it doesn't work without delay

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make movies')
    parser.add_argument('--structures_dir', type=str, required=True, 
                        help="Path to directory containing directories with structures")
    parser.add_argument('--videos_dir', type=str, required=True,
                        help="Directory where to put generated videos")
    args = parser.parse_args()

    structures_dir = args.structures_dir
    videos_dir = args.videos_dir
    
    os.makedirs(videos_dir, exist_ok=True)

    for entry in os.listdir(structures_dir):
        entry_path = os.path.join(structures_dir, entry)

        true_pdb_path = glob('aligned_filtered_true*.pdb', root_dir=entry_path)[0]
        val_pdb_path = glob('aligned_filtered_val*.pdb', root_dir=entry_path)[0]

        true_pdb_path, val_pdb_path = map(lambda x: os.path.join(entry_path, x), (true_pdb_path, val_pdb_path))
        
        create_aligned_movie(true_pdb_path, val_pdb_path, 
                             movie_path=os.path.join(videos_dir, entry + '.mp4'))

cmd.quit()
