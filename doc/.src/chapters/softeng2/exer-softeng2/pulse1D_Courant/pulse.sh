#!/bin/sh
set -x

# Simulation code
pulse_types="gaussian cosinehat plug half-cosinehat"
C_values="1.0 0.9 0.75"
for pulse_tp in $pulse_types; do
    if [ ! -d "$pulse_tp" ]; then
	mkdir $pulse_tp
    fi
    cd $pulse_tp
    for C in $C_values; do
	#python ../pulse.py $pulse_tp $C
	echo
    done
    python ../../../src-softeng2/animate_archives.py -0.5 1.3 ${pulse_tp}*/.*.npz
    sleep 2
    ffmpeg -i tmp_%04d.png -r 25 -vcodec libtheora movie.ogg
    ffmpeg -i tmp_%04d.png -r 25 -vcodec libx264 movie.mp4
    ffmpeg -i tmp_%04d.png -r 25 -vcodec libvpx movie.webm
    cd ..
done
