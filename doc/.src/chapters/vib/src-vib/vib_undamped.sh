#!/bin/sh
prog=vib_undamped.py
movieprog=ffmpeg
opt="--SCITOOLS_easyviz_backend matplotlib"

# Compare two time steps for 5 periods
python $prog --dt 0.1 --num_periods 5 $opt
mv vib1.png tmp1_vib.png
mv vib1.pdf tmp1_vib.pdf
python $prog --dt 0.05 --num_periods 5 $opt
mv vib1.png tmp2_vib.png
mv vib1.pdf tmp2_vib.pdf

doconce combine_images tmp1_vib.png tmp2_vib.png vib_phase_err1.png
doconce combine_images tmp1_vib.pdf tmp2_vib.png vib_phase_err1.pdf

cp vib_phase_err1.* ../fig-vib

# Make movies

dt=0.05
python $prog --dt $dt --num_periods 39 --savefig
scitools movie output_file=vib.html fps=4 tmp_vib*.png
dir=vib_undamped_movie_dt$dt
rm -rf $dir ../mov-vib/$dir
mkdir $dir
mv tmp_vib* $dir
mv vib.html $dir/index.html
cd $dir
$movieprog -r 12 -i tmp_vib%04d.png -c:v flv movie.flv
$movieprog -r 12 -i tmp_vib%04d.png -c:v libx264 movie.mp4
$movieprog -r 12 -i tmp_vib%04d.png -c:v libvpx movie.webm
$movieprog -r 12 -i tmp_vib%04d.png -c:v libtheora movie.ogg
cd ..
cp -r $dir ../mov-vib
git add ../mov-vib/$dir
