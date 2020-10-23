import time

import numpy as np

import cambfast
from cambfast import get_spectrum_camb

lmax = 50
npts = 200

cf = cambfast.CAMBfast()
cf.inifile = '../ini/planck_2018.ini'
cf.add_parameter('tau', 0.05, 0.003, 0.12, npts, fixed=False)
cf.add_parameter('r', 0.05, -0.3, 0.3, npts, fixed=True)
cf.generate_interp(lmax, CMB_unit='muK')
#cf.write_funcs(f'As_tau_r_lmax{lmax}_npts{npts}.npz')
cf.write_funcs(f'tau_lmax{lmax}_npts{npts}.npz')

