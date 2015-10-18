#!/bin/sh
schemes="FE BE theta"

for scheme in $schemes; do

# Make demo: plug
rm -f tmp_*
python -c "from diffu1D_u0 import plug; plug('${scheme}', Fo=0.5, Nx=50);"
sh ../../make_moviefiles.sh tmp_frame%04d.png 20
name=../mov-diffu/diffu1D_u0_${scheme}_plug
if [ ! -d $name ]; then
mkdir $name
fi
mv -f tmp_frame*.png movie.* $name

# For the FE scheme, run the plug with Fo=0.25 to demonstrate a smooth
# solution
if [ ${scheme} == "FE" ]; then
rm -f tmp_*
python -c "from diffu1D_u0 import plug; plug('${scheme}', Fo=0.5, Nx=50);"
sh ../../make_moviefiles.sh tmp_frame%04d.png 40
name=../mov-diffu/diffu1D_u0_${scheme}_plug_Fo025
if [ ! -d $name ]; then
mkdir $name
fi
mv -f tmp_frame*.png movie.* $name
fi


# Make demo: peaked gaussian
rm -f tmp_*
python -c "from diffu1D_u0 import gaussian; gaussian('${scheme}', Fo=0.5, Nx=50, sigma=0.05);"
sh ../../make_moviefiles.sh tmp_frame%04d.png 20
name=../mov-diffu/diffu1D_u0_${scheme}_gaussian1
if [ ! -d $name ]; then
mkdir $name
fi
mv -f tmp_frame*.png movie.* $name

# Take a long time step for BE and theta=0.5
if [ ${scheme} != "FE" ]; then
rm -f tmp_*
python -c "from diffu1D_u0 import gaussian; gaussian('${scheme}', Fo=5, Nx=50, sigma=0.05);"
sh ../../make_moviefiles.sh tmp_frame%04d.png 2
name=../mov-diffu/diffu1D_u0_${scheme}_gaussian1_Fo5
if [ ! -d $name ]; then
mkdir $name
fi
mv -f tmp_frame*.png movie.* $name
fi


# Make demo: wide gaussian
rm -f tmp_*
python -c "from diffu1D_u0 import gaussian; gaussian('${scheme}', Fo=0.5, Nx=50, sigma=0.15);"
sh ../../make_moviefiles.sh tmp_frame%04d.png 20
name=../mov-diffu/diffu1D_u0_${scheme}_gaussian2
if [ ! -d $name ]; then
mkdir $name
fi
mv -f tmp_frame*.png movie.* $name

done
