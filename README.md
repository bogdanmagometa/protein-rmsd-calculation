# protein-rmsd-calculation

## Project structure

- `calc_rmsd.py` - Add `rmsd1`, `rmsd2`, `rmsd3` columns to the `metrics.csv` table, save aligned filtered PDBs.
- `make_movies.py` - Create movies from aligned filtered PDBs.
- `add_text_all.py` - Overlay text on all movies.
- `arrange_grid.sh` - Arrange all movies into one grid of specified grid width and height.

## Prerequisites

- Linux OS
- Maybe macOS (at least [finish_launching](https://pymolwiki.org/index.php/Launching_From_a_Script) might not work)
- You need to have conda installed
- To generate videos without watermarks, license might be needed. Educational license is free: https://pymol.org/edu/

## Environment setup

```bash
$ conda create -n rmsd_calculation python=3.10
$ conda activate rmsd_calculation
$ conda install numpy pandas
$ conda install -c conda-forge prody    # used in calc_rmsd.py only
$ conda install -c conda-forge -c schrodinger pymol-bundle  # used in make_movies.py only
$ conda install -c conda-forge ffmpeg   # used in add_text_all.py, add_text.sh
```

## Usage example

Assuming the following directory tree structure:
```
./
    - output/
        - data/
            - structures/
                ...directories each containing true and predicted structures...
            - metrics.csv
    - datasets/
        -npz_data
            ...npz-files...
```

```bash
$ python3 calc_rmsd.py --data_dir ./output/data --npz_data ./datasets/npz_data
$ python3 make_movies.py --structures_dir ./output/data/structures --videos_dir ./videos
$ python3 add_text_all.py --metrics_path ./output/data/metrics_with_rmsds.csv --videos_dir ./videos --videos_with_text_dir ./videos_with_text
$ ./arrange_grid.sh ./videos_with_text ./all_structs_video.mp4 5 8
```
