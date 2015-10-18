#!/bin/bash

function make_movies {
  dir=pulse${1}_in_two_media
  if [ -d $dir ]; then
    rm -rf $dir
  fi
  mkdir $dir
  mv -f frame*.png $dir
  cd $dir
  sh ../../../../make_moviefiles.sh frame_%04d.png 10
  doconce combine_images frame_0084.png frame_0108.png plot.png
  cd ..
}

# Make two representative demos of how a smooth pulse and a less
# smooth pulse is numerically treated when it goes from one
# medium into another (focus: reflected and transmitted waves,
# numerical noise).
# The end result is a movie for each case.

rm -f frame*.png
python -c "import wave1D_dn_vc as w; w.pulse(loc='left', pulse_tp='half-cosinehat', Nx=80, every_frame=1, slowness_factor=4, sigma=0.05)"
make_movies 1

rm -f frame*.png
python -c "import wave1D_dn_vc as w; w.pulse(loc='left', pulse_tp='gaussian', Nx=80, every_frame=1, slowness_factor=4, sigma=0.075)"
h
make_movies 2
