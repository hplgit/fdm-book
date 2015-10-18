#!/bin/sh
prog=wave1D_u0_s.py

grep 'if n == 90:' $prog
if [ $? -ne 0 ]; then
  echo "insert if n == 90: st.savefig('frame_C%s.pdf' % C) in $prog"
  exit
fi

C_values="1.0 0.95 0.2 1.0015"
for C in $C_values; do
python $prog $C
scitools movie output_file=index.html fps=2 frame*.png
scitools movie encoder=convert output_file=movie.gif fps=4 frame*.png
dir=guitar_C$C
rm -rf $dir
mkdir $dir
mv movie.gif index.html frame*.png $dir
done
scitools rename frame_C wave1D_guitar_C frame_C*.pdf
