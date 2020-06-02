import time

import numpy as np

import cambfast
from cambfast import get_spectrum_camb

lmax = 50
npts = 100

cf = cambfast.CAMBfast()
cf.ini_file = '../ini/planck_2018.ini'
cf.add_parameter('tau', 0.05, 0.02, 0.08, npts, fixed=False)
cf.add_parameter('r', 0.00, -0.3, 0.3, npts, fixed=False)
cf.generate_interp(lmax, CMB_unit='muK')
#cf.generate_interp_MC(lmax, CMB_unit='muK', nsamples=npts)
#cf.write_funcs(f'tau_r_lmax{lmax}_npts{npts}_MC.npz')

