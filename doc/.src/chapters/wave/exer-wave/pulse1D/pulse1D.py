import wave1D_dn_vc as wave
import os, sys, shutil, glob
import scitools.std as st

for pulse_tp in 'gaussian', 'cosinehat', 'half-cosinehat', 'plug':
    for Nx in 40, 80, 160:
        for sf in 1, 2, 4:
            if sf == 1 and Nx > 40:
                continue  # homogeneous medium with C=1: Nx=40 is sufficient
            print 'wave1D.pulse:', pulse_tp, Nx, sf
            wave.pulse(C=1, Nx=Nx, animate=False, # just hardcopies
                       version='vectorized',
                       T=2, loc='left', pulse_tp=pulse_tp,
                       slowness_factor=sf, medium=[7, 9],
                       every_frame=2*Nx/40)
            dir = 'pulse1D_%s_Nx%d_sf%g' % (pulse_tp, Nx, sf)
            if os.path.isdir(dir):
                shutil.rmtree(dir)
            os.mkdir(dir)
            for filename in glob.glob('frame*.png'):
                os.rename(filename, os.path.join(dir, filename))
            os.chdir(dir)
            st.movie('frame_*.png', encoder='html',
                     output_file='index.html', fps=4)
            os.chdir(os.pardir)

