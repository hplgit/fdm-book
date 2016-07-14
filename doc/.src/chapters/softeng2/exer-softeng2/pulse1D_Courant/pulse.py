import sys, os, glob
path = os.path.join(os.pardir, os.pardir,
                    os.pardir, os.pardir, 'wave', 'src-wave', 'wave1D')
#print path, glob.glob(path + '/w*')
sys.path.insert(0, path)
from wave1D_dn_vc import pulse
pulse_tp = sys.argv[1]
C = float(sys.argv[2])
pulse(pulse_tp=pulse_tp, C=C,
      Nx=100, animate=False, slowness_factor=4)
