#!/bin/bash
# Arrange videos into a grid of specified size
#   ./arrange_grid.sh *path_to_dir_with_videos* *output_file* *grid_width* *grid_height*
#

# Arguments
videos_with_text_dir=$1
output_fpath=$2
grid_width=$3
grid_height=$4

# Calculate width and height of any input video
one_inpt_file=$videos_with_text_dir/$(ls $videos_with_text_dir | tail -n 1)
width=$(ffprobe -v error -select_streams v:0 -show_entries stream=width -of default=noprint_wrappers=1:nokey=1 $one_inpt_file)
height=$(ffprobe -v error -select_streams v:0 -show_entries stream=height -of default=noprint_wrappers=1:nokey=1 $one_inpt_file)

# Generate input streams for empty grid cells
empty_stream="-f lavfi -i color=size=${width}x${height}:rate=30:color=black:d=2"
num_empty_streams=$(($grid_width * $grid_height - $(ls $videos_with_text_dir | wc -l)))
if [ $num_empty_streams -lt 0 ]; then
    echo "Error: Too many videofiles for grid of size ${grid_width}x${grid_height}" >&2
    exit 1
else
    for i in $(seq 1 $num_empty_streams); do
        empty_streams="${empty_streams} ${empty_stream}"
    done
    echo $empty_streams
fi

# Arrange videos into the grid
ffmpeg $(echo ${videos_with_text_dir}/*.mp4 | xargs printf " -i %s ") $empty_streams -filter_complex "xstack=grid=${grid_width}x${grid_height}" $output_fpath
