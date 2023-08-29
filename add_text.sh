#!/bin/bash

text_file=$1
input_video=$2
output_video=$3

width=$(ffprobe -v error -select_streams v:0 -show_entries stream=width -of default=noprint_wrappers=1:nokey=1 $input_video)
height=$(ffprobe -v error -select_streams v:0 -show_entries stream=height -of default=noprint_wrappers=1:nokey=1 $input_video)

ffmpeg -y -i "$input_video" -vf "drawtext=textfile=${text_file}:fontsize=24:x=-14:fontcolor=white:fontfile=/home/bohdan/Downloads/Arial.ttf:fontsize=42, pad=width=iw+4:height=ih+4:color=red:x=2:y=2" -c:v libx264 -preset veryslow -qp 0 "$output_video"
# ffmpeg -y -i $2 -vf "drawtext=textfile=${1}:fontsize=24:x=0:fontcolor=white:fontfile=/home/bohdan/Downloads/Arial.ttf:fontsize=38, pad=width=$((width+2)):height=$((height+2)):color=red" -c:v libx264 -preset veryslow -qp 0 $3
# ffmpeg -y -i $2 -vf "drawtext=textfile=${1}:fontsize=24:x=(w-text_w)/2:fontcolor=white:fontfile=/home/bohdan/Downloads/Arial.ttf:fontsize=48" -c:v libx264 -preset veryslow -qp 0 $3
