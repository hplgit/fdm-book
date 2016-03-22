import wave1D_dn_vc as wave
import os, sys, shutil, glob

for pulse_tp in 'gaussian', 'cosinehat', 'half-cosinehat', 'plug':
    for Nx in 40, 80, 160:
        for sf in 2, 4:
            if sf == 1 and Nx > 40:
                continue  # homogeneous medium with C=1: Nx=40 enough
            print 'wave1D.pulse:', pulse_tp, Nx, sf

            wave.pulse(C=1, Nx=Nx, animate=False, # just hardcopies
                       version='vectorized',
                       T=2, loc='left', pulse_tp=pulse_tp,
                       slowness_factor=sf, medium=[0.7, 0.9],
                       skip_frame = 1,
                       sigma=0.05)
