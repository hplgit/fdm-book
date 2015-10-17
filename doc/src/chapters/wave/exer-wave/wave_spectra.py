import numpy as np
import matplotlib.pyplot as plt

def spectrum(f, x):
    # Discrete Fourier transform
    A = np.fft.rfft(f(x))
    A_amplitude = np.abs(A)

    # Compute the corresponding frequencies
    dx = x[1] - x[0]
    freqs = np.linspace(0, np.pi/dx, A_amplitude.size)

    plt.plot(freqs[:len(freqs)/2], A_amplitude[:len(freqs)/2])

# Mesh
L = 10; Nx = 100
x = np.linspace(0, L, Nx+1)

spectrum(lambda x: np.where(x < 5, 1, 0), x)
spectrum(lambda x: np.sin(np.pi*x/float(L)) + np.sin(np.pi*20*x/float(L)), x)
s = 0.5
spectrum(lambda x: 1./(np.sqrt(2*np.pi)*s)*np.exp(-0.5*((x-L/2.)/s)**2), x)

def f(x):
    r = np.zeros_like(x)
    r[len(x)/2] = 1
    return r

spectrum(f, x)

figfile = 'tmp'
plt.legend(['step', '2sin', 'gauss', 'peak'])
plt.savefig(figfile + '.pdf')
plt.savefig(figfile + '.png')
plt.show()
