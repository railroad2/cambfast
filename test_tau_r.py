import time

import numpy as np
import pylab as plt

import cambfast
from cambfast import get_spectrum_camb
from gbpipe.spectrum import get_spectrum_camb

fn = 'tau_r_lmax50_planck2018.npz'
cf = cambfast.CAMBfast(fn)

nsamples = 20
taus = np.linspace(0.04, 0.061, nsamples)
rs = np.linspace(-0.1, 0.11, nsamples)

fig1 = plt.figure()
ax1 = fig1.add_subplot() 
plt.xlabel('$l$')
plt.ylabel('$C_l (\mu \mathrm{K}^2)$')

fig2 = plt.figure()
ax2 = fig2.add_subplot()
plt.xlabel('$l$')
plt.ylabel(r'$\frac{C_l^\mathrm{camb}-C_l^{cf}}{C_l^\mathrm{camb}}$')

for tau, r in zip(taus, rs):
    cl_camb = get_spectrum_camb(lmax=50, 
                                ini_file='../ini/planck_2018.ini', 
                                tau=tau, r=r,
                                CMB_unit='muK', isDl=False)

    cl_cf   = cf.get_spectrum(lmax=50, 
                              tau=tau, r=r,
                              isDl=False)

    ax1.loglog(np.abs(cl_camb[0]), 'r-', linewidth=0.5)
    ax1.loglog(np.abs(cl_camb[1]), 'g-', linewidth=0.5)
    ax1.loglog(np.abs(cl_camb[2]), 'b-', linewidth=0.5)
    ax1.loglog(np.abs(cl_camb[3]), 'y-', linewidth=0.5)

    ax1.loglog(np.abs(cl_cf[0]), 'r--')
    ax1.loglog(np.abs(cl_cf[1]), 'g--')
    ax1.loglog(np.abs(cl_cf[2]), 'b--')
    ax1.loglog(np.abs(cl_cf[3]), 'y--')

    ax2.loglog(np.abs((cl_camb[0] - cl_cf[0])/cl_camb[0]), 'r')
    ax2.loglog(np.abs((cl_camb[1] - cl_cf[1])/cl_camb[1]), 'g')
    ax2.loglog(np.abs((cl_camb[2] - cl_cf[2])/cl_camb[2]), 'b')
    ax2.loglog(np.abs((cl_camb[3] - cl_cf[3])/cl_camb[3]), 'y')

ax1.loglog([], 'r-', label='TT')
ax1.loglog([], 'g-', label='EE')
ax1.loglog([], 'b-', label='BB')
ax1.loglog([], 'y-', label='TE')
ax1.legend()
plt.tight_layout()

ax2.loglog([], 'r-', label='TT')
ax2.loglog([], 'g-', label='EE')
ax2.loglog([], 'b-', label='BB')
ax2.loglog([], 'y-', label='TE')
ax2.legend()
plt.tight_layout()

plt.show()

