import time

import numpy as np

import cambfast
from cambfast import get_spectrum_camb

lmax = 50
npts = 50 

cf = cambfast.CAMBfast()
cf.inifile = '../ini/planck_2018_kmlee.ini'
cf.add_parameter('tau', 0.05, 0, 0.20, npts, fixed=False)
cf.add_parameter('r', 0.00, -0.3, 0.3, npts, fixed=True)
cf.generate_interp(lmax, CMB_unit='muK')
#cf.write_funcs(f'As_tau_r_lmax{lmax}_npts{npts}.npz')
#cf.write_funcs(f'tau_lmax{lmax}_npts{npts}_mintau0017.npz')
cf.write_funcs(f'tau_lmax{lmax}_npts{npts}_mintau0.npz')

