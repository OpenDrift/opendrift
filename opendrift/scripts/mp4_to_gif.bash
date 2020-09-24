#!/usr/bin/env bash

# Convert a mp4 animation to gif, with minimum filesize
# Usage: mp42gif <mp4file>


mp4file=$1
ext=${mp4file##*.}
fname=`basename $mp4file $ext`

ffmpeg -y -i ${fname}mp4 -vf palettegen palette.png

ffmpeg -y -i ${fname}mp4 -i palette.png -filter_complex paletteuse -r 10 ${fname}gif

rm palette.png
