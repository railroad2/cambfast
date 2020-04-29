import numpy as np
import pylab as plt
import time
import cambfast

cf = cambfast.CAMBfast()
cf.load_funcs('./tau_r_lmax50_npts100.npz')
cfmc = cambfast.CAMBfast()
cfmc.load_MCdata('./tau_r_lmax50_npts5000_MC.npz', method='cubic')

t0 = time.time()
cl = cf.get_spectrum(lmax=50, tau=0.024827, r=0.02987)
print (f'interpolation : {time.time()-t0} s')
t0 = time.time()
clmc = cfmc.get_spectrum(lmax=50, tau=0.024827, r=0.02987)#*(2.7255*1e6)**2
print (f'interpolation MC: {time.time()-t0} s')
t0 = time.time()
cl1 = cambfast.get_spectrum_camb(lmax=50, tau=0.024827, r=0.02987, CMB_unit='muK')
print (f'CAMB : {time.time()-t0} s')

par0 = cf.pars[0]
par1 = cf.pars[1]


## grid vs random points
par0, par1 = np.mgrid[par0.min:par0.max:1j*par0.nstep, par1.min:par1.max:1j*par1.nstep]
rpts = cfmc.MCpoints
plt.plot(par0, par1, 'b.', markersize=1)
plt.plot(rpts.T[0], rpts.T[1], 'r.', markersize=1)

## spectra
plt.figure()
lines = plt.loglog(cl.T)
lines += plt.loglog(clmc.T)
lines += plt.loglog(cl1.T)
labels = ['grid TT', 'EE', 'BB', 'TE']
labels += ['MC TT', 'EE', 'BB', 'TE']
labels += ['camb TT', 'EE', 'BB', 'TE']
plt.legend(lines, labels)

## difference in spectra
plt.figure()
lines = plt.plot((cl1.T - cl.T)/cl1.T, label=['camb - grid TT', 'EE', 'BB', 'TE'])
linesmc = plt.plot((cl1.T - clmc.T)/cl1.T, label=['camb - mc TT', 'EE', 'BB', 'TE'])
labels = ['camb - grid, TT', 'EE', 'BB', 'TE']
labels += ['camb - mc, TT', 'EE', 'BB', 'TE']
plt.legend(lines+linesmc, labels)

plt.show()

