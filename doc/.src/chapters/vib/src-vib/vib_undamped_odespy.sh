#!/bin/sh
prog=vib_undamped_odespy.py
plot="--SCITOOLS_easyviz_backend matplotlib"

python $prog 20 theta 1 $plot
python $prog 40 theta 1 $plot
doconce combine_images vib_20_1_u.png vib_40_1_u.png vib_theta_1_u.png
doconce combine_images vib_20_1_pp.png vib_40_1_pp.png vib_theta_1_pp.png
doconce combine_images vib_20_1_u.pdf vib_40_1_u.pdf vib_theta_1_u.pdf
doconce combine_images vib_20_1_pp.pdf vib_40_1_pp.pdf vib_theta_1_pp.pdf

python $prog 10 RK 1 $plot
python $prog 20 RK 1 $plot
doconce combine_images vib_10_1_u.png vib_20_1_u.png vib_RK_1_u.png
doconce combine_images vib_10_1_pp.png vib_20_1_pp.png vib_RK_1_pp.png
doconce combine_images vib_10_1_u.pdf vib_20_1_u.pdf vib_RK_1_u.pdf
doconce combine_images vib_10_1_pp.pdf vib_20_1_pp.pdf vib_RK_1_pp.pdf

python $prog 10 CN 10 $plot
python $prog 20 CN 10 $plot
doconce combine_images vib_10_10_u.png vib_20_10_u.png vib_CN_10_u.png
doconce combine_images vib_10_10_pp.png vib_20_10_pp.png vib_CN_10_pp.png
doconce combine_images vib_10_10_u.pdf vib_20_10_u.pdf vib_CN_10_u.pdf
doconce combine_images vib_10_10_pp.pdf vib_20_10_pp.pdf vib_CN_10_pp.pdf

python $prog 10 RK 10 $plot
python $prog 20 RK 10 $plot
doconce combine_images vib_10_10_u.png vib_20_10_u.png vib_RK_10_u.png
doconce combine_images vib_10_10_pp.png vib_20_10_pp.png vib_RK_10_pp.png
doconce combine_images vib_10_10_u.pdf vib_20_10_u.pdf vib_RK_10_u.pdf
doconce combine_images vib_10_10_pp.pdf vib_20_10_pp.pdf vib_RK_10_pp.pdf

mv -f vib_RK*.p* vib_CN*.p* vib_theta*.p* ../fig-vib




